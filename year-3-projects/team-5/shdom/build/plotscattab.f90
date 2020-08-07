PROGRAM PLOTSCATTAB
! PLOTSCATTAB make a plotting file of the phase functions in a scattering 
! table file produced by make_mie_table.f90 or make_ice_table.f90.  
! The table contains the Legendre polynomial expansion of phase functions
! for several effective radii, which are converted to phase functions as
! a function of angle.  The program can also plot the phase functions in 
! the beginning of a tabulated phase function property file.
!
!  compile: pgf90 -fast -o plotscattab  plotscattab.f90
! 
!    Frank Evans    University of Colorado     May 2003

  IMPLICIT NONE
  INTEGER :: NANGLE, NOUT, NPHASE, MAXLEG
  INTEGER :: I, J, K, L
  LOGICAL :: FOUND
  REAL    :: PI, PL, PL1, PL2, SCAT, MU
  INTEGER, ALLOCATABLE :: NLEG(:)
  REAL,    ALLOCATABLE :: ANGLE(:), OUTREFF(:), REFF(:)
  REAL,    ALLOCATABLE :: LEGCOEF(:,:), PHASE(:,:)
  CHARACTER(LEN=1)  :: FILETYPE
  CHARACTER(LEN=72) :: INFILE, PLOTFILE

  
  WRITE (*,*) 'Scattering table or SHDOM property file input (S or P)'
  READ (*,'(A1)') FILETYPE
    WRITE (*,'(1X,A1)') FILETYPE

  WRITE (*,*) 'Input file name'
  READ (*,'(A)') INFILE
    WRITE (*,'(A)') INFILE

  WRITE (*,'(1X,A)') 'Number of angles (<0 to input |n| angles)'
  READ (*,*) NANGLE
    WRITE (*,*) NANGLE
  ALLOCATE (ANGLE(ABS(NANGLE)))
  IF (NANGLE <= 0) THEN
    NANGLE = ABS(NANGLE)
    WRITE (*,*) 'Input the angles (degrees)'
    DO J = 1, NANGLE
      WRITE (*,'(1X,I3,A)') J, ' : '
      READ (*,*) ANGLE(J)
    ENDDO
  ELSE
    ANGLE(:) = (/ (180.0*FLOAT(J-1)/(NANGLE-1), J=1,NANGLE) /)
  ENDIF


  IF (FILETYPE == 'P') THEN
    WRITE (*,'(1X,A)') 'Number of phase functions to output'
    READ (*,*) NOUT
      WRITE (*,*) NOUT
    ALLOCATE (OUTREFF(NOUT))
    WRITE (*,*) 'Input the tabulated phase function indices'
    READ (*,*) OUTREFF(:)
      WRITE (*,'(20(1X,F5.1))') OUTREFF(:)
  ELSE
    WRITE (*,'(1X,A)') 'Number of effective radii to output'
    READ (*,*) NOUT
      WRITE (*,*) NOUT
    ALLOCATE (OUTREFF(NOUT))
    WRITE (*,*) 'Input the effective radii (micron)'
    READ (*,*) OUTREFF(:)
      WRITE (*,'(20(1X,F6.2))') OUTREFF(:)
  ENDIF

  WRITE (*,*) 'Plotting output file name'
  READ (*,'(A)') PLOTFILE
    WRITE (*,'(A)') PLOTFILE


  IF (FILETYPE == 'P') THEN
    CALL READ_PROPERTY_SIZE (INFILE, NPHASE, MAXLEG)
    ALLOCATE (REFF(NPHASE), NLEG(NPHASE), LEGCOEF(0:MAXLEG,NPHASE))
    CALL READ_PROPERTY_PHASEFUNC (INFILE, NPHASE, MAXLEG, NLEG, LEGCOEF)
    REFF(:) = (/ (FLOAT(J), J=1,NPHASE) /)
  ELSE
    CALL READ_SCAT_TABLE_SIZE (INFILE, NPHASE, MAXLEG)
    ALLOCATE (REFF(NPHASE), NLEG(NPHASE), LEGCOEF(0:MAXLEG,NPHASE))
    CALL READ_SCAT_TABLE (INFILE, NPHASE, MAXLEG, REFF, NLEG, LEGCOEF)
  ENDIF

   ! Loop over each effective radius in Mie table (I), pulling out
   !   ones we want to output (K)
  ALLOCATE (PHASE(NANGLE,NOUT))
  PI = ACOS(-1.0)
  DO K = 1, NOUT
    FOUND = .FALSE.
    DO I = 1, NPHASE
      IF (REFF(I) == OUTREFF(K)) THEN
        FOUND = .TRUE.
         ! Sum the Legendre series for each angle in plot
        DO J = 1, NANGLE
          MU = COS(ANGLE(J)*PI/180.0)
          SCAT = 0.0
          PL1 = 1.0
          PL = 1.0
          DO L = 0, NLEG(I)
            IF (L .GT. 0)  PL = (2*L-1)*MU*PL1/L - (L-1)*PL2/L
            SCAT = SCAT + LEGCOEF(L,I)*PL
            PL2 = PL1
            PL1 = PL
          ENDDO
          PHASE(J,K) = SCAT
        ENDDO
      ENDIF
    ENDDO
    IF (.NOT. FOUND) THEN
      PRINT *, 'Phase function not found in input file:',OUTREFF(K)
      PHASE(:,K) = 1.0
    ENDIF
  ENDDO


   ! Output the phase functions
  OPEN (UNIT=1, FILE=PLOTFILE, STATUS='UNKNOWN')
  IF (FILETYPE == 'P') THEN
    WRITE (1,'(A,A40)') '!  Property file phase functions: ',INFILE
    WRITE (1,'(A)')'! Angle  cos(angle)  Phase functions with index numbers'
    WRITE (1,'(A,20(6X,I3,3X))') '!               ', (NINT(OUTREFF(K)), K=1,NOUT)
  ELSE
    WRITE (1,'(A,A40)') '!  Scattering table phase functions: ',INFILE
    WRITE (1,'(A)')'! Angle  cos(angle)  Phase functions for effective radii (um)'
    WRITE (1,'(A,20(6X,F6.2))') '!               ', (OUTREFF(K), K=1,NOUT)
  ENDIF
  DO J = 1, NANGLE
      WRITE (1,'(1X,F7.2,1X,F9.6,20(E12.4))') &
          ANGLE(J), COS(ANGLE(J)*PI/180.0), (PHASE(J,K), K=1,NOUT)
  ENDDO
  CLOSE (1)

  DEALLOCATE (PHASE, REFF, NLEG, LEGCOEF, OUTREFF, ANGLE)
END





SUBROUTINE READ_SCAT_TABLE_SIZE (SCATTABFILE, NRETAB, MAXLEG)
 ! Reads Legendre phase function coefficients from a table of Mie scattering 
 ! properties as a function of effective radius.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  INTEGER, INTENT(OUT) :: NRETAB, MAXLEG
  INTEGER :: I, L, NLEG
  REAL    :: REFF, EXT, SSALB, CHI

  OPEN (UNIT=3, FILE=SCATTABFILE, STATUS='OLD')
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*) NRETAB
  MAXLEG = 1
  DO I = 1, NRETAB
    READ (3,*) REFF, EXT, SSALB, NLEG
    READ (3,*) (CHI, L=0,NLEG) 
    MAXLEG = MAX(MAXLEG,NLEG)
  ENDDO
  CLOSE (3)
END SUBROUTINE READ_SCAT_TABLE_SIZE



SUBROUTINE READ_SCAT_TABLE (SCATTABFILE, NRETAB, MAXLEG, REFF, NLEG, LEGCOEF)
 ! Reads Legendre phase function coefficients from a table of Mie scattering 
 ! properties as a function of effective radius.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  INTEGER, INTENT(IN)  :: NRETAB, MAXLEG
  INTEGER, INTENT(OUT) :: NLEG(NRETAB)
  REAL,    INTENT(OUT) :: REFF(NRETAB), LEGCOEF(0:MAXLEG,NRETAB)
  INTEGER :: I, L
  REAL    :: EXT, SSALB

  OPEN (UNIT=3, FILE=SCATTABFILE, STATUS='OLD')
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  DO I = 1, NRETAB
    READ (3,*) REFF(I), EXT, SSALB, NLEG(I)
    READ (3,*) (LEGCOEF(L,I), L=0,NLEG(I)) 
  ENDDO
  CLOSE (3)
END SUBROUTINE READ_SCAT_TABLE



 
SUBROUTINE READ_PROPERTY_SIZE (PROPFILE, NUMPHASE, MAXLEG)
 ! Reads the header of the tabulated phase function property file to 
 ! get the maximum array sizes needed for allocatable arrays.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER, INTENT(OUT) :: NUMPHASE, MAXLEG
  INTEGER :: NPX, NPY, NPZ, NUML, I, K, L
  REAL    :: DELX, DELY, ZLEVELS, CHI
  CHARACTER(LEN=1) :: PROPTYPE
 
  ! Open the file, figure out the type, and get the grid size
  OPEN (UNIT=2, FILE=PROPFILE, STATUS='OLD')
  READ (2,'(A1)') PROPTYPE
  IF (PROPTYPE /= 'T') THEN
    PRINT *, 'Must be a tabulated phase function property file.'
    STOP
  ENDIF
  READ (2,*) NPX, NPY, NPZ
  READ (2,*) DELX, DELY, (ZLEVELS, K=1,NPZ)
  READ (2,*) NUMPHASE
  DO I = 1, NUMPHASE
    READ (2,*) NUML, (CHI, L=1,NUML)
    MAXLEG = MAX(NUML,MAXLEG)
  ENDDO
  CLOSE (2)
END SUBROUTINE READ_PROPERTY_SIZE



SUBROUTINE READ_PROPERTY_PHASEFUNC (PROPFILE, NUMPHASE, MAXLEG, NLEG, LEGCOEF)
 ! Reads the Legendre phase function coefficients from the header of a
 ! tabulated phase function file.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER, INTENT(IN)  :: NUMPHASE, MAXLEG
  INTEGER, INTENT(OUT) :: NLEG(NUMPHASE)
  REAL,    INTENT(OUT) :: LEGCOEF(0:MAXLEG,NUMPHASE)
  INTEGER :: NPX, NPY, NPZ, I, K, L
  REAL    :: DELX, DELY, ZLEVELS
 
  OPEN (UNIT=2, FILE=PROPFILE, STATUS='OLD')
  READ (2,*)
  READ (2,*) NPX, NPY, NPZ
  READ (2,*) DELX, DELY, (ZLEVELS, K=1,NPZ)
  READ (2,*)
  DO I = 1, NUMPHASE
    READ (2,*) NLEG(I), (LEGCOEF(L,I), L=1,NLEG(I))
    LEGCOEF(0,I) = 1.0
  ENDDO
  CLOSE (2)
END SUBROUTINE READ_PROPERTY_PHASEFUNC
 

