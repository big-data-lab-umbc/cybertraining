PROGRAM MAKE_ICE_TABLE

!   Reads the 0.25 degree gaussian forward peak smoothed version of 
! Ping Yang's shortwave ice crystal scattering database, and creates a
! scattering table for gamma size distributions of a sequence of effective 
! radii.  Only one particle shape is selected from the input database.  
! The scattering properties are averaged over the desired spectral range 
! with solar Planck function weighting.  The phase functions in the output
! scattering table are represented with Legendre series.
!
!  compile: pgf90 -fast -o make_ice_table  make_ice_table.f90
!
!    Frank Evans    University of Colorado     May 2003 

  IMPLICIT NONE
  INTEGER :: SHAPEIND, NRETAB
  REAL    :: WAVELEN1, WAVELEN2, SRETAB, ERETAB, ALPHA
  CHARACTER(LEN=80) :: ICESCATDB, ICETABFILE
  INTEGER, PARAMETER :: NSIZE=24   ! number of sizes in ice scattering db
  INTEGER, PARAMETER :: NANG=288   ! number of phase function angles in db
  INTEGER, PARAMETER :: NQUAD=2500 ! number of quad angles for Legendre conversion
  INTEGER, PARAMETER :: MAXLEG=NQUAD ! number of Legendre coeffs
  INTEGER :: I, J
  REAL    :: REFF, SCATTER, AREA, PI
  INTEGER, ALLOCATABLE :: NLEG(:)
  REAL, ALLOCATABLE :: DIAMAREA(:), DIAMVOL(:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:), QSCA(:)
  REAL, ALLOCATABLE :: ANGLES(:), PHASEFUNC(:,:), PHASEDIST(:)
  REAL, ALLOCATABLE :: EXTINCT(:), SSALB(:), LEGCOEF(:,:)
  REAL(8), ALLOCATABLE :: QUADMU(:), QUADWTS(:)


  CALL USER_INPUT (WAVELEN1, WAVELEN2, ICESCATDB, SHAPEIND, &
                   NRETAB, SRETAB, ERETAB, ALPHA, ICETABFILE)


  ! Allocate all the arrays here
  ALLOCATE (QUADMU(NQUAD), QUADWTS(NQUAD))
  ALLOCATE (DIAMAREA(NSIZE), DIAMVOL(NSIZE))
  ALLOCATE (QEXT(NSIZE), QSCA(NSIZE), ANGLES(NANG), PHASEFUNC(NANG,NSIZE))
  ALLOCATE (ND(NSIZE), PHASEDIST(NANG))
  ALLOCATE (EXTINCT(NRETAB), SSALB(NRETAB))
  ALLOCATE (NLEG(NRETAB), LEGCOEF(0:MAXLEG,NRETAB))

  ! Make the Gauss-Legendre quadrature abscissas and weights
  CALL GAUSQUAD (NQUAD, QUADMU, QUADWTS)


  ! Read in the ice scattering properties for the desired shape and 
  !   solar weighted average over the wavelength range.  The diameters
  !   are in microns.
  CALL READ_ICE_SCAT_DB (ICESCATDB, SHAPEIND, WAVELEN1, WAVELEN2, &
                         NSIZE, NANG, DIAMAREA, DIAMVOL, QEXT, QSCA, &
                         ANGLES, PHASEFUNC)

  IF (SRETAB < 2.0*0.5*DIAMVOL(1) .OR. ERETAB > 0.5*0.5*DIAMVOL(NSIZE)) THEN
    PRINT *, 'MAKE_ICE_TABLE: effective radius range outside possible range'
    STOP
  ENDIF

  PI = ACOS(-1.0)  
  ! Loop over the number of output tabulated effective radii
  DO I = 1, NRETAB
    ! Set tabulated effective radius
    IF (NRETAB <= 1) THEN
      REFF = SRETAB
    ELSE
      REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
    ENDIF

    ! Calculate the discrete size number concentrations (ND), which vary
    !   according to a gamma distribution in volume equivalent diameter,
    !   that give the desired effective radius (REFF) and IWC (1 g/m^3).
    CALL MAKE_SIZE_DIST (NSIZE, DIAMAREA, DIAMVOL, REFF, ALPHA, ND)

    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    PHASEDIST(:) = 0.0
    DO J = 1, NSIZE
      AREA = PI*(0.5*DIAMAREA(J))**2 
      EXTINCT(I) = EXTINCT(I) + QEXT(J)*AREA*ND(J)
      SCATTER = SCATTER + QSCA(J)*AREA*ND(J)
      PHASEDIST(:) = PHASEDIST(:) + QSCA(J)*AREA*ND(J)*PHASEFUNC(:,J)
    ENDDO
    PHASEDIST(:) = PHASEDIST(:)/SCATTER
    SSALB(I) = SCATTER/EXTINCT(I)
    EXTINCT(I) = 0.001*EXTINCT(I)

    ! Convert phase function from angle to Legendre coefficients
    CALL CONVERT_LEGENDRE (NANG, ANGLES, PHASEDIST, NQUAD, QUADMU, QUADWTS, &
                           MAXLEG, NLEG(I), LEGCOEF(:,I))

  ENDDO  ! end of effective radius loop


  CALL WRITE_SCAT_TABLE (ICETABFILE, WAVELEN1, WAVELEN2, SHAPEIND, &
                         NRETAB, SRETAB, ERETAB, ALPHA, &
                         EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)

  DEALLOCATE (DIAMAREA, DIAMVOL, QEXT, QSCA, ANGLES, PHASEFUNC)
  DEALLOCATE (ND, PHASEDIST, EXTINCT, SSALB, NLEG, LEGCOEF)
END





SUBROUTINE USER_INPUT (WAVELEN1, WAVELEN2, ICESCATDB, SHAPEIND, &
                       NRETAB, SRETAB, ERETAB, ALPHA, ICETABFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: SHAPEIND, NRETAB
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2, SRETAB, ERETAB, ALPHA
  CHARACTER(LEN=*), INTENT(OUT) :: ICESCATDB, ICETABFILE

  WRITE (*,*) 'Input ice scattering database name'
  READ (*,'(A)') ICESCATDB
    WRITE(*,'(1X,A70)') ICESCATDB

  WRITE (*,*) 'Wavelength range (micron)'
  READ (*,*) WAVELEN1, WAVELEN2
    WRITE (*,'(2(1X,F6.3))') WAVELEN1, WAVELEN2
  IF (WAVELEN1 > WAVELEN2) STOP 'USER_INPUT: wavelength1 must be <= wavelength2'

  WRITE(*,*) 'Index of ice crystal shape (1=hollow column, 2=solid column, 3=plate,'
  WRITE(*,*) '4=dendrite, 5=rough aggregate, 6=smooth agg., 7=4-bullet rosette, 8=6-bullet)'
  READ(*,*) SHAPEIND
    WRITE(*,'(1X,I1)') SHAPEIND

  WRITE(*,*) 'Number, starting, and ending tabulated effective radius'
  READ(*,*) NRETAB, SRETAB, ERETAB
    WRITE (*,'(1X,I3,2(1X,F7.2))') NRETAB, SRETAB, ERETAB

  WRITE(*,*) 'Gamma size distribution shape parameter (alpha)'
  READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA

  WRITE (*,*) 'Output scattering table name'
  READ (*,'(A)') ICETABFILE
    WRITE(*,'(1X,A70)') ICETABFILE
END SUBROUTINE USER_INPUT





SUBROUTINE READ_ICE_SCAT_DB (ICESCATDB, SHAPEIND, WAVELEN1, WAVELEN2, &
                         NSIZE, NANG, DIAMAREA, DIAMVOL, QEXT, QSCA, &
                         ANGLES, PHASEFUNC)
 ! Reads the 0.25 degree gaussian forward peak smoothed version of 
 ! Ping Yang's shortwave ice crystal scattering database.  The scattering
 ! properties are averaged over the desired spectral range with solar 
 ! Planck function weighting. 
 ! Inputs:
 !   ICESCATDB  database filename
 !   SHAPEIND   shape number (1-8)
 !   WAVELEN1, WAVELEN2  wavelength range (micron)
 !   NSIZE      number of particle sizes (24)
 !   NANG       number of phase function angles (288)
 ! Outputs: 
 !   DIAMAREA   equivalent area spherical diameter (microns)
 !   DIAMVOL    equivalent volume spherical diameter (microns)
 !   QEXT(:)    extinction efficiency (2 in geometric optics limit)
 !   QSCA(:)    scattering efficiency
 !   ANGLES(:)  angles of phase functions (degrees)
 !   PHASEFUNC(:,:) phase functions
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: ICESCATDB
  INTEGER, INTENT(IN)  :: NSIZE, NANG, SHAPEIND
  REAL,    INTENT(IN)  :: WAVELEN1, WAVELEN2
  REAL,    INTENT(OUT) :: DIAMAREA(NSIZE), DIAMVOL(NSIZE)
  REAL,    INTENT(OUT) :: QEXT(NSIZE), QSCA(NSIZE)
  REAL,    INTENT(OUT) :: ANGLES(NANG), PHASEFUNC(NANG,NSIZE)
  INTEGER, PARAMETER :: NWAVE=56
  REAL,    PARAMETER :: BBTEMP=5800., DELWAVE=0.002 ! for wavelen integration
  INTEGER :: I, J, K, N, ISH, D1, D2, D3
  REAL    :: WAVE, PLANCK, SUM, EXT, SSALB, PHASE
  REAL, ALLOCATABLE :: WAVEDB(:), WTWAVE(:)
  CHARACTER(LEN=288*3) :: PHASESTRING

  QEXT(:) = 0.0
  QSCA(:) = 0.0
  PHASEFUNC(:,:) = 0.0

  OPEN (UNIT=1, FILE=ICESCATDB, STATUS='OLD')  
  ! Read useful parts of header: wavelengths and phase function angles
  DO I = 1, 10
    READ (1,*)
  ENDDO
  READ (1,*) N, (ANGLES(I), I=1,N)
  IF (N /= NANG) STOP 'READ_ICE_SCAT_DB: wrong ice scattering database'
  READ (1,*)
  ALLOCATE (WAVEDB(NWAVE+1), WTWAVE(NWAVE))
  READ (1,*) (WAVEDB(I), I=1,NWAVE+1)
  READ (1,*)

  IF (WAVELEN1 < WAVEDB(1) .OR. WAVELEN2 > WAVEDB(NWAVE+1)) &
    STOP 'READ_ICE_SCAT_DB: Wavelengths out of range of ice scattering database'

  ! Make the solar Planck weighted weights for each db band in wavelength range
  SUM = 0.0
  WTWAVE(:) = 0.0
  WAVE = WAVELEN1
  K = 1
  DO WHILE (WAVE <= WAVELEN2)
    DO WHILE (WAVEDB(K+1) <= WAVE .AND. K < NWAVE)
      K = K + 1
    ENDDO
    PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
    SUM = SUM + PLANCK
    WTWAVE(K) = WTWAVE(K) + PLANCK
    WAVE = WAVE + DELWAVE
  ENDDO
  WTWAVE(:) = WTWAVE(:)/SUM
  print *
  print *,'Wavelen1 Wavelen2  Weight  <- summing database bands'
  do k = 1, nwave
    if (wtwave(k) > 0) then
      print '(2(2X,F5.3,2X),2X,F6.4)', wavedb(k:k+1), wtwave(k)
    endif
  enddo

  ! Skip over wrong particle shapes
  DO I = 1, SHAPEIND-1
    DO K = 1, NWAVE
      DO J = 1, NSIZE
        READ (1,*)
      ENDDO
    ENDDO
  ENDDO
  ! If this is a wavelength band we need then read the data in and do the
  !   weighted average of single scattering properties.
  DO K = 1, NWAVE
    IF (WTWAVE(K) > 0.0) THEN
      DO J = 1, NSIZE
        READ (1,'(I2,1X,F6.3,50X,2(1X,F8.3),2(1X,F8.6),10X,A)') &
            ISH, WAVE, DIAMAREA(J), DIAMVOL(J), EXT, SSALB, PHASESTRING
        IF (ISH /= SHAPEIND .OR. WAVE /= WAVEDB(K)) &
          STOP 'READ_ICE_SCAT_DB: Error reading scattering database'
        QEXT(J) = QEXT(J) + WTWAVE(K)*EXT
        QSCA(J) = QSCA(J) + WTWAVE(K)*EXT*SSALB
        DO I = 1, NANG
          D1 = ICHAR(PHASESTRING(3*I-2:3*I-2))-32
          D2 = ICHAR(PHASESTRING(3*I-1:3*I-1))-32
          D3 = ICHAR(PHASESTRING(3*I:3*I))-32
          PHASE = 10**(10.D0*(D1/95.D0+D2/95.D0**2+D3/95.D0**3) - 4.D0)
          PHASEFUNC(I,J) = PHASEFUNC(I,J) + WTWAVE(K)*EXT*SSALB*PHASE
        ENDDO
      ENDDO
    ELSE  ! Otherwise skip this wavelength
      DO J = 1, NSIZE
        READ (1,*)
      ENDDO
    ENDIF
  ENDDO

  DO J = 1, NSIZE
    PHASEFUNC(:,J) = PHASEFUNC(:,J)/QSCA(J)
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_ICE_SCAT_DB




SUBROUTINE MAKE_SIZE_DIST (NSIZE, DIAMAREA, DIAMVOL, REFF, ALPHA, ND)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle sizes of a gamma size distribution (in volume
 ! spherical equivalent diameter DIAMVOL) with a distribution 
 ! effective radius (0.75*volume/projected area) of REFF, gamma
 ! shape parameter ALPHA, and ice water content of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: DIAMAREA(NSIZE), DIAMVOL(NSIZE), REFF, ALPHA
  REAL,    INTENT(OUT) :: ND(NSIZE)
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
  INTEGER :: I
  REAL    :: DEHI, DELO, DEMID, RE


  ! Find Deff that gives Reff above desired value
  DEHI = 2.0*REFF
  I = 0
  RE = 0.5*REFF
  DO WHILE (RE <= REFF .AND. I < 4)
    DEHI = 2.0*DEHI
    I = I + 1
    CALL DO_GAMMA_DIST (NSIZE, DEHI, ALPHA, DIAMAREA, DIAMVOL, ND, RE)
  ENDDO
  IF (RE <= REFF) THEN
    PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,RE
    STOP
  ENDIF

  ! Find Deff that gives Reff below desired value
  DELO = 2.0*REFF
  I = 0
  RE = 2*REFF
  DO WHILE (RE >= REFF .AND. I < 4)
    DELO = 0.5*DELO
    I = I + 1
    CALL DO_GAMMA_DIST (NSIZE, DELO, ALPHA, DIAMAREA, DIAMVOL, ND, RE)
  ENDDO
  IF (RE >= REFF) THEN
    PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,RE
    STOP
  ENDIF

  ! Do bisection to get correct effective radius
  DO WHILE (ABS(RE-REFF) > TOL*REFF)
    DEMID = 0.5*(DELO+DEHI)
    CALL DO_GAMMA_DIST (NSIZE, DEMID, ALPHA, DIAMAREA, DIAMVOL, ND, RE)
    IF (RE < REFF) THEN
      DELO = DEMID
    ELSE
      DEHI = DEMID
    ENDIF
  ENDDO  
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_GAMMA_DIST (NSIZE, DEFF, ALPHA, DIAMAREA, DIAMVOL, ND, RE)
 ! For the input effective diameter (DEFF) and ALPHA, returns the 
 ! number concentrations ND [cm^-3] and the calculated effective radius 
 ! RE [um] for a gamma size distribution in DIAMVOL with total 
 ! IWC of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: DEFF, ALPHA, DIAMAREA(NSIZE), DIAMVOL(NSIZE)
  REAL,    INTENT(OUT) :: ND(NSIZE), RE
  REAL, PARAMETER :: DENSICE = 0.916
  INTEGER :: J
  REAL    :: PI, B, D, DELD, IWC, SUM2, SUM3

  PI = ACOS(-1.0)
  B = (ALPHA+3)/DEFF
  IWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    D = DIAMVOL(J)
    DELD = SQRT(DIAMVOL(J)*DIAMVOL(MIN(NSIZE,J+1))) &
         - SQRT(DIAMVOL(J)*DIAMVOL(MAX(1,J-1)))
    ND(J) = D**ALPHA *EXP(-B*D) * DELD
    IWC = IWC + 1.0E-6*DENSICE*ND(J)*(PI/6)*DIAMVOL(J)**3
    SUM2 = SUM2 + ND(J)*DIAMAREA(J)**2
    SUM3 = SUM3 + ND(J)*DIAMVOL(J)**3
  ENDDO
!  ND(:) = (1.0/IWC)*ND(:)
  DO J = 1, NSIZE
     ND(J) = (1.0/IWC)*ND(J)
   ENDDO
  RE = 0.5*SUM3/SUM2
END SUBROUTINE DO_GAMMA_DIST




SUBROUTINE CONVERT_LEGENDRE (NANG, ANGLES, PHASE, NQUAD, QUADMU, QUADWTS, &
                             MAXLEG, NLEG, LEGCOEF)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NANG, NQUAD, MAXLEG
  REAL,    INTENT(IN)  :: ANGLES(NANG), PHASE(NANG)
  REAL(8), INTENT(IN)  :: QUADMU(NQUAD), QUADWTS(NQUAD)
  INTEGER, INTENT(OUT) :: NLEG
  REAL,    INTENT(OUT) :: LEGCOEF(0:MAXLEG)
  INTEGER :: N, L
  REAL(8) :: INTEG, PL, PL1, PL2
  REAL(8), ALLOCATABLE :: PHASEQUAD(:), LEGEN(:)

  ALLOCATE (PHASEQUAD(NQUAD), LEGEN(0:MAXLEG))

  ! Interpolate phase function to the gaussian quadrature points
  CALL SPLINE_INTERP_PHASE (NANG, ANGLES, PHASE, &
                            NQUAD, QUADMU, PHASEQUAD, QUADWTS, INTEG)

  ! Do final small normalization by multiplication
  PHASEQUAD(:) = (2.D0/INTEG)*PHASEQUAD(:)

  ! Compute the Legendre coefficients for the smoothed phase func
  LEGEN(:) = 0.0
  DO N = 1, NQUAD
    ! Use upward recurrence to find Legendre polynomials
    PL1 = 1.0
    PL = 1.0 
    DO L = 0, MAXLEG
      IF (L .GT. 0) PL = (2*L-1)*QUADMU(N)*PL1/L-(L-1)*PL2/L
      LEGEN(L)=LEGEN(L) + PL*QUADWTS(N)*PHASEQUAD(N)
      PL2 = PL1
      PL1 = PL 
    ENDDO
  ENDDO
  ! Find last significant Legendre coefficient
  DO L = 0, MAXLEG
    LEGCOEF(L) = 0.5*(2*L+1)*LEGEN(L)
    IF (LEGCOEF(L) .GT. 1.0E-5) THEN   
      NLEG = L
    ENDIF
  ENDDO  

  DEALLOCATE (PHASEQUAD, LEGEN)
END SUBROUTINE CONVERT_LEGENDRE



SUBROUTINE SPLINE_INTERP_PHASE (NANG, ANGLES, PHASE, &
                                NQUAD, QUADMU, PHASEQUAD, QUADWTS, INTEG)
 ! Interpolates the phase function from ANGLES (degrees) to 
 ! quadrature points (in cos theta).  Also computes integral.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NANG, NQUAD
  REAL,    INTENT(IN) :: ANGLES(NANG), PHASE(NANG)
  REAL(8), INTENT(IN) :: QUADMU(NQUAD), QUADWTS(NQUAD)
  REAL(8), INTENT(OUT) :: PHASEQUAD(NQUAD), INTEG
  INTEGER :: I, J
  REAL(8) :: RD, SEVAL
  REAL(8), ALLOCATABLE :: X(:), Y(:), B(:), C(:), D(:)

  ALLOCATE (X(NANG), Y(NANG), B(NANG), C(NANG), D(NANG))

  RD = ACOS(-1.0D0)/180.D0
  DO I = 1, NANG
     X(I) = COS(RD*ANGLES(NANG+1-I))
     Y(I) = PHASE(NANG+1-I)
  ENDDO
  CALL SPLINE (NANG, X, Y, B, C, D)

  INTEG = 0.0D0
  DO J = 1, NQUAD
    PHASEQUAD(J) = SEVAL (NANG, QUADMU(J), X, Y, B, C, D)
    INTEG = INTEG + QUADWTS(J)*PHASEQUAD(J)
  ENDDO
END SUBROUTINE SPLINE_INTERP_PHASE



SUBROUTINE WRITE_SCAT_TABLE (ICETABFILE, WAVELEN1, WAVELEN2, SHAPEIND, &
                             NRETAB, SRETAB, ERETAB, ALPHA, &
                             EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)
 ! Writes the table of ice scattering properties as a function of 
 ! effective radius.  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: SHAPEIND, NRETAB, MAXLEG, NLEG(NRETAB)
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, SRETAB, ERETAB, ALPHA
  REAL,    INTENT(IN) :: EXTINCT(NRETAB), SSALB(NRETAB)
  REAL,    INTENT(IN) :: LEGCOEF(0:MAXLEG,NRETAB)
  CHARACTER(LEN=*), INTENT(IN) :: ICETABFILE
  INTEGER :: I, J, L, NL
  REAL    REFF
  CHARACTER(LEN=16), PARAMETER :: SHAPENAMES(8) = &
    (/ 'hollow column   ', 'solid column    ', 'plate           ', &
       'dendrite        ', 'rough aggregate ', 'smooth aggregate', &
       '4-bullet rosette', '6-bullet rosette' /)

  OPEN (UNIT=3, FILE=ICETABFILE, STATUS='REPLACE')
  WRITE (3,'(A)') '! Ice scattering table vs. effective radius (IWC=1 g/m^3)'
  WRITE (3,'(2(1X,F5.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  WRITE (3, '(A)') '! Ice crystal scattering from Ping Yang scattering data.'
  WRITE (3,'(I2,1X,A16,A)') SHAPEIND, SHAPENAMES(SHAPEIND), &
        ', particle shape number and name'
  WRITE (3,'(1X,F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NRETAB, SRETAB, ERETAB, &
        '  number, starting, ending effective radius'

  DO I = 1, NRETAB
    IF (NRETAB <= 1) THEN
      REFF = SRETAB
    ELSE
      REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
    ENDIF
    WRITE (3,'(1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)') &
        REFF, EXTINCT(I), SSALB(I), NLEG(I), '  Reff  Ext  Alb  Nleg'
    WRITE (3,'(2X,201(1X,F10.5))') (LEGCOEF(L,I), L=0,MIN(NLEG(I),200))
    DO J = 200, NLEG(I)-1, 200
      WRITE (3,'(2X,200(1X,F10.5))') (LEGCOEF(J+L,I),L=1,MIN(200,NLEG(I)-J))
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_SCAT_TABLE



SUBROUTINE GAUSQUAD (N, XA, WT)
 ! Generates the abscissas (X) and weights (W) for an N point
 ! Gauss-Legendre quadrature.  
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: XA(N), WT(N)
  INTEGER :: K, I, J, L
  REAL(8) :: X, XP, PL, PL1, PL2, DPL

  K = (N+1)/2
  DO J = 1, K
    X = COS(3.141592654*(J-.25)/(N+.5))
    I = 0
    DO WHILE (I < 10)
      PL1 = 1
      PL = X
      DO L = 2, N
         PL2 = PL1
         PL1 = PL
         PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
       ENDDO
       DPL = N*(X*PL-PL1)/(X*X-1)
       XP = X
       X = XP - PL/DPL
       I = I+1
       IF (ABS(X-XP) < 2*EPSILON(X)) EXIT
    ENDDO
    XA(J)     = -X
    XA(N-J+1) = X
    WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
    WT(N-J+1) = WT(J)
  ENDDO
END SUBROUTINE GAUSQUAD



subroutine spline (n, x, y, b, c, d)
  implicit none
  integer n
  real(8) x(n), y(n), b(n), c(n), d(n)

!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      integer nm1, ib, i
      double precision t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
end subroutine spline



double precision function seval(n, u, x, y, b, c, d)
  implicit none
  integer n
  real(8)  u, x(n), y(n), b(n), c(n), d(n)
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      real*8 dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
end function seval


