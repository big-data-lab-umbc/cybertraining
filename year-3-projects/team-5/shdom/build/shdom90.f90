!  SHDOM: Spherical harmonic discrete ordinate radiative transfer method.
!      See shdom.txt for documentation.
!      Fortran 90 version of the main program for using allocatable arrays.

      MODULE SHDOM_PROPERTY_ARRAYS
        INTEGER, SAVE              :: NPX, NPY, NPZ, NUMPHASE
        REAL,    SAVE              :: DELX, DELY, XSTART, YSTART
        REAL,    SAVE, ALLOCATABLE :: ZLEVELS(:)
        REAL,    SAVE, ALLOCATABLE :: TEMPP(:), EXTINCTP(:), ALBEDOP(:)
        REAL,    SAVE, ALLOCATABLE :: LEGENP(:), EXTDIRP(:)
        INTEGER*2,SAVE,ALLOCATABLE :: IPHASEP(:)
        INTEGER, SAVE              :: NZCKD
        REAL,    SAVE, ALLOCATABLE :: ZCKD(:), GASABS(:)
      END MODULE


      PROGRAM SHDOM
      USE SHDOM_PROPERTY_ARRAYS
      USE shdom_netcdf
      IMPLICIT NONE
!                         Number of output files and inputs parameter per file
      INTEGER, PARAMETER :: MAXOUT=30, MAXPAR=400
      REAL                 :: OUTPARMS(MAXPAR,MAXOUT)
      CHARACTER            :: OUTTYPES(MAXOUT)*1, OUTFILES(MAXOUT)*64 

      REAL,    ALLOCATABLE :: DELG(:), KABS(:,:)
      REAL,    ALLOCATABLE :: XGRID(:), YGRID(:), ZGRID(:)
      REAL,    ALLOCATABLE :: TEMP(:), PLANCK(:), EXTINCT(:), ALBEDO(:)
      REAL,    ALLOCATABLE :: LEGEN(:)
      INTEGER*2,ALLOCATABLE:: IPHASE(:)
      INTEGER, ALLOCATABLE :: GRIDPTR(:,:), NEIGHPTR(:,:), TREEPTR(:,:)
      INTEGER*2,ALLOCATABLE:: CELLFLAGS(:)
      REAL,    ALLOCATABLE :: GRIDPOS(:,:)
      INTEGER, ALLOCATABLE :: RSHPTR(:), SHPTR(:), OSHPTR(:)
      REAL,    ALLOCATABLE :: RADIANCE(:), SOURCE(:), DELSOURCE(:)
      INTEGER, ALLOCATABLE :: WORK1(:)
      REAL,    ALLOCATABLE :: WORK(:), WORK2(:)
      REAL,    ALLOCATABLE :: FLUXES(:), DIRFLUX(:)
      INTEGER, ALLOCATABLE :: NPHI0(:)
      REAL,    ALLOCATABLE :: MU(:), WTDO(:), PHI(:)
      REAL,    ALLOCATABLE :: CMU1(:), CMU2(:), CPHI1(:), CPHI2(:)
      INTEGER, ALLOCATABLE :: BCPTR(:)
      REAL,    ALLOCATABLE :: BCRAD(:)
      REAL,    ALLOCATABLE :: SFCPARMS(:)
      REAL,    ALLOCATABLE :: SFCGRIDPARMS(:)
      REAL,    ALLOCATABLE :: SUMFLUXES(:), SUMDIRFLUX(:), SUMFLUXDIV(:)
      REAL,    ALLOCATABLE :: SUMSHTERMS(:), SOURCE1OUT(:)
      REAL,    ALLOCATABLE :: SUMRADOUT(:), SUMVISOUT(:)
      REAL,    ALLOCATABLE :: TEMPPT(:), EXTINCTPT(:), ALBEDOPT(:)
      REAL,    ALLOCATABLE :: LEGENPT(:)
      INTEGER*2, ALLOCATABLE :: IPHASEPT(:)
      REAL,    ALLOCATABLE :: ALLFLUXES(:,:,:,:)
      REAL,    ALLOCATABLE :: ALLFLUXDIV(:,:,:), ALLSHTERMS(:,:,:,:)
      REAL,    ALLOCATABLE :: XALLGRID(:), YALLGRID(:)

      INTEGER :: MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO, MAXNZ, MAXNG
      INTEGER :: MAXLEG, MAXPGL, MAXIGL
      INTEGER :: MAXNBC, MAXBCRAD
      INTEGER :: MAXSFCPTS, MAXSFCPARS,  BIG_ARRAYS
      REAL    :: MEMWORD, WANTMEM
      REAL    :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
      REAL    :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
      REAL    :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

      INTEGER NX, NY, NZ, NX1, NY1
      INTEGER NXT, NYT, NPXT, NPYT, MAXPGT, MAXPGLT
      INTEGER NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
      INTEGER NG, SIDE
      INTEGER ML, MM, NLM, NLEG, NANG, NMU, NPHI, NPHI0MAX, NCS
      INTEGER NBPTS, NPTS, NCELLS, NSH, NBCELLS, OLDNPTS
      INTEGER MAXITER, ITER, TOTITER, IG, I, NUMOUT, BCFLAG, IPFLAG
      INTEGER NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
      INTEGER NCELLSTOT, NPTSTOT, NSHTOT
      LOGICAL KDIST, INRADFLAG, NEWGRIDFLAG, DELTAM, ACCELFLAG
      LOGICAL BASEOUT, HIGHORDERRAD, MASTERPROC, VALIDBEAM, BTEST
      LOGICAL RADIANCE_OUT, VISUAL_OUT
      LOGICAL FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
      REAL  SOLARFLUX, SOLARMU, SOLARAZ, SOLFLUX
      REAL  GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL  SOLACC, SOLCRIT, SPLITACC, SHACC
      REAL  DELXSFC, DELYSFC, MAXASYM, REDUCE
      REAL  MUOUT, PHIOUT
      REAL  UNIFZLEV, XO, YO, ZO, DIRPATH
      REAL  CPUTIME1, CPUTIME2, CPUTIMETOTAL, SUM_CPU_TIME
      CHARACTER SRCTYPE*1, UNITS*1
      CHARACTER GRIDTYPE*1, PROPTYPE*1, SFCTYPE*2
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
      CHARACTER*64 OutFileNC, RUNNAME
      EXTERNAL SUM_CPU_TIME


!          Choose the option to compute higher order radiance SH terms
      HIGHORDERRAD = .FALSE.
      NSHOUT = 4
      IF (HIGHORDERRAD) NSHOUT=5
!          Choose type of initialization/output for k-distributions:
!            True means initiatize with and output only the base grid.
      BASEOUT = .FALSE.

      CALL CPU_TIME(CPUTIME1) 

       ! Start up MPI for multiple processors if available, and
       !  return flag for master process (procnum=0)
      CALL START_MPI (MASTERPROC)


      IF (MASTERPROC) THEN 
        CALL USER_INPUT (RUNNAME, PROPFILE, SFCFILE, CKDFILE,        &
                       INSAVEFILE, OUTSAVEFILE,                      &
                       NX, NY, NZ, NMU, NPHI, BCFLAG,                &
                       IPFLAG, KDIST, DELTAM, GRIDTYPE,              &
                       SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                       GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,   &
                       ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,  &
                       MAXOUT,MAXPAR, NUMOUT,OUTTYPES,OUTPARMS,OUTFILES, &
                       OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,    &
                       NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
!        CALL NAMELIST_INPUT (RUNNAME, PROPFILE, SFCFILE, CKDFILE,    &
!                       INSAVEFILE, OUTSAVEFILE,                      &
!                       NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG,        &
!                       KDIST, DELTAM, GRIDTYPE,                      &
!                       SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
!                       GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,   &
!                       ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,  &
!                       MAXOUT,MAXPAR,NUMOUT,OUTTYPES,OUTPARMS,OUTFILES, &
!                       OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,    &
!                       NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
      ENDIF
      CALL BROADCAST_USER_INPUT (RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                       INSAVEFILE, OUTSAVEFILE,                      &
                       NX, NY, NZ, NMU, NPHI, BCFLAG,                &
                       IPFLAG, KDIST, DELTAM, GRIDTYPE,              &
                       SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                       GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,   &
                       ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,  &
                       MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,   &
                       MAX_TOTAL_MB, ADAPT_GRID_FACTOR,              &
                       NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

      NX = MAX(1,NX) ; NY = MAX(1,NY) ; NZ = MAX(2,NZ)
!         If single plane or column then force IP mode
      IF (NX .EQ. 1)  IPFLAG=IOR(IPFLAG,1)
      IF (NY .EQ. 1)  IPFLAG=IOR(IPFLAG,2)

!          Make ML and MM from NMU and NPHI
      NMU = MAX(2, 2*INT((NMU+1)/2) )
      NPHI = MAX(1,NPHI)
      ML = NMU-1
      MM = MAX(0,INT(NPHI/2)-1)
      NLEG = ML
      IF (DELTAM) NLEG = ML+1

!          Make ML and MM from NMU and NPHI
!          Compute NLM, NPHI0MAX.
      IF (SRCTYPE .EQ. 'T')  SOLARAZ = 0.0
      IF (BTEST(IPFLAG,1) .AND. ABS(SIN(SOLARAZ)) .LT. 1.0E-4) THEN
        NCS = 1
      ELSE
        NCS = 2
      ENDIF
      IF (NCS .EQ. 1) THEN
        NPHI0MAX = INT((NPHI+2)/2)
        NLM = (MM+1)*(ML+1) - (MM*(MM+1))/2
      ELSE
        NPHI0MAX = NPHI
        NLM = (2*MM+1)*(ML+1) - MM*(MM+1)
      ENDIF
      MEMWORD = NMU*(2+2*NPHI+2*NLM+2*33*32)

      ALLOCATE (MU(NMU), WTDO(NMU*NPHI), PHI(NMU*NPHI), NPHI0(NMU))
      ALLOCATE (CMU1(NLM*NMU), CMU2(NLM*NMU))
      ALLOCATE (CPHI1(NMU*33*32),CPHI2(NMU*33*32))


!         Find out the size of the input property file grid 
      IF (MASTERPROC) THEN
        CALL READ_PROPERTY_SIZE (PROPFILE, NLEG, NPX, NPY, NPZ, &
                                 NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
!        CALL read_property_size_netcdf(PROPFILE, NLEG, NPX, NPY, NPZ, &
!                                       NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
      ENDIF
      CALL BROADCAST_PROPERTY_SIZE (NPX, NPY, NPZ, DELX, DELY, &
                                    NUMPHASE, MAXLEG, MAXPGL)

       ! Set the full domain grid sizes (variables end in T)
      NXT=NX ; NYT=NY ; NPXT=NPX ; NPYT=NPY

       ! Map the full SHDOM domain onto the multiple processors:
       !   change NPX,NPY and NX,NY to the subdomain size, and output 
       !   the X and Y starting position of this subdomain.  Set bits
       !   2 and 3 in BCFLAG to indicate multiple processors in X and Y.
       !   Also open the slave processor log files.
      CALL MAP_SHDOM_MPI (BCFLAG, NPX,NPY, NX,NY,NZ, DELX, DELY, PROPFILE, &
                          XSTART, YSTART, RUNNAME)

!         Set up base grid point actual size (NX1xNY1xNZ)
      NX1 = NX+1
      IF (BTEST(IPFLAG,0) .OR. BTEST(BCFLAG,0) .OR. BTEST(BCFLAG,2)) NX1 = NX
      NY1 = NY+1
      IF (BTEST(IPFLAG,1) .OR. BTEST(BCFLAG,1) .OR. BTEST(BCFLAG,3)) NY1 = NY
      IF (NZ .LE. 1) THEN
        WRITE (6,*) 'NZ < 2'
        STOP
      ENDIF
      NBPTS = NX1*NY1*NZ
!        Calculate the number of base grid cells depending on the BCFLAG
      NBCELLS = (NZ-1)*(NX+IBITS(BCFLAG,0,1)-IBITS(BCFLAG,2,1)) &
                      *(NY+IBITS(BCFLAG,1,1)-IBITS(BCFLAG,3,1))

!         Determine the size of the property grid arrays so can allocate
      MAXPG = NPX*NPY*NPZ
      IF (MAXPGL .GT. MAXLEG*NUMPHASE) MAXPGL=MAXPG*MAXLEG
      ALLOCATE (ZLEVELS(NPZ))
      ALLOCATE (TEMPP(MAXPG), EXTINCTP(MAXPG), ALBEDOP(MAXPG))
      ALLOCATE (EXTDIRP(MAXPG), IPHASEP(MAXPG), LEGENP(MAXPGL))
      MEMWORD = MEMWORD + 4.5*MAXPG + MAXPGL + NUMPHASE*(MAXLEG+1)

      MAXNZ = NPZ
      IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
         ! Read the properties of the medium to the master processor and 
         !  distribute to all the processors
        IF (MASTERPROC) THEN
          MAXPGT = NPXT*NPYT*NPZ
          MAXPGLT = MAXPGL
          IF (MAXPGLT .GT. MAXLEG*NUMPHASE) MAXPGLT=MAXPGT*MAXLEG
          ALLOCATE (TEMPPT(MAXPGT), EXTINCTPT(MAXPGT), ALBEDOPT(MAXPGT))
          ALLOCATE (IPHASEPT(MAXPGT), LEGENPT(MAXPGLT))
          CALL READ_PROPERTIES (PROPFILE, NPXT, NPYT, NPZ, &
                 MAXLEG, NLEG, MAXNZ, MAXPGT, MAXPGLT, DELTAM, &
                 PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM, &
                 TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE, IPHASEPT)
!          CALL read_properties_netcdf(PROPFILE, NPXT, NPYT, NPZ, NUMPHASE,   &
!                                    MAXLEG, MAXPGT, MAXPGLT, DELTAM, NLEG, &
!                                    PROPTYPE, ZLEVELS, MAXASYM,          &
!                                    TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT) 
          IF (NUMPHASE == 0) THEN
            WRITE (6,*) 'Standard property file (with phase function Legendre series at every point)'
            WRITE (6,*) '  is not allowed with multiple processors.'
            STOP
          ENDIF
           ! Check values of medium properties
          CALL CHECK_PROPERTY_INPUT (NPXT, NPYT, NPZ, NLEG,  &
                      DELX, DELY, ZLEVELS, &
                      TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE, IPHASEPT)
        ENDIF
        CALL SCATTER_PROPERTIES (NPXT, NPYT, NPX, NPY, NPZ, NLEG, NUMPHASE, &
                       ZLEVELS, MAXASYM, &
                       TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, &
                       TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
        IF (MASTERPROC) DEALLOCATE (TEMPPT, EXTINCTPT, ALBEDOPT, IPHASEPT)
      ELSE
         ! Read the properties of the medium to the single processor
        CALL READ_PROPERTIES (PROPFILE, NPX, NPY, NPZ, &
                 MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM, &
                 PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM, &
                 TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
!          CALL read_properties_netcdf(PROPFILE, NPX, NPY, NPZ, NUMPHASE,   &
!                                    MAXLEG, MAXPG, MAXPGL, DELTAM, NLEG, &
!                                    PROPTYPE, ZLEVELS, MAXASYM,          &
!                                    TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP) 
         ! Check values of medium properties
        CALL CHECK_PROPERTY_INPUT (NPX, NPY, NPZ, NLEG,  &
                       DELX, DELY, ZLEVELS, &
                       TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
      ENDIF

      IF (MASTERPROC) THEN
        ! Check input parameters
        CALL CHECK_INPUT_PARMETERS (NMU, NPHI, DELTAM, MAXASYM, &
                        GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU, &
                        WAVENO, WAVELEN, GNDTEMP, GNDALBEDO, &
                        SPLITACC, SHACC)
      ENDIF


      IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
         ! Read the memory parameters from a file, if it exists, and send
         !  the parameters to all processors. Otherwise use the input ones.
        CALL READ_BROADCAST_MEM_PARAMS (MAX_TOTAL_MB, ADAPT_GRID_FACTOR, &
                            NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)
      ENDIF


!         Guess maximum number of grid points, cells, SH vector size needed
!           but don't let MAX_TOTAL_MB be exceeded
      IF (MAX_TOTAL_MB*1024.0**2 > 1.75*HUGE(MAXIV)) THEN
        MAX_TOTAL_MB = 1.75*HUGE(MAXIV)/1024.0**2
        WRITE (6,*) 'MAX_TOTAL_MB reduced to fit memory model: ', MAX_TOTAL_MB
      ENDIF
      IF (SPLITACC .LE. 0.0 .AND. INSAVEFILE(1:4) .EQ. 'NONE') &
        ADAPT_GRID_FACTOR = 1.0
      BIG_ARRAYS = 3
      IF (.NOT. ACCELFLAG) BIG_ARRAYS = 2
      WANTMEM = ADAPT_GRID_FACTOR*NBPTS *(28 + 16.5*CELL_TO_POINT_RATIO  &
                        + NPHI0MAX + NUM_SH_TERM_FACTOR*NLM*BIG_ARRAYS)
      REDUCE = MIN(1.0, ((MAX_TOTAL_MB*1024**2)/4 - MEMWORD) /WANTMEM)
      ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*REDUCE
      IF (ADAPT_GRID_FACTOR < 1.0) THEN
        WRITE (6,*) 'MAX_TOTAL_MB memory limit exceeded with just base grid.'
        STOP
      ENDIF
      IF (REDUCE < 1.0) THEN
         WRITE(6,*) 'ADAPT_GRID_FACTOR reduced to ',ADAPT_GRID_FACTOR
      ENDIF
      WANTMEM = REDUCE*WANTMEM
      IF (WANTMEM > HUGE(MAXIV)) THEN
        WRITE(6,*) 'Number of words of memory exceeds max integer size:', WANTMEM
        STOP
      ENDIF
      MAXIG = NINT(ADAPT_GRID_FACTOR*NBPTS)
      MAXIC = NINT(CELL_TO_POINT_RATIO*MAXIG)
      MAXIV = NINT(NUM_SH_TERM_FACTOR*NLM*MAXIG)
      MAXIDO = MAXIG*NPHI0MAX
      IF (NUMPHASE .GT. 0) THEN
        MAXIGL = NUMPHASE*(MAXLEG+1)
      ELSE
        MAXIGL = MAXIG*(MAXLEG+1)
      ENDIF
      IF (4.0*(MAXIV+MAXIG) > HUGE(MAXIV)) THEN
        WRITE(6,*) 'Size of big SH arrays (MAXIV) probably exceeds max integer number of bytes:',MAXIV
        STOP
      ENDIF
      IF (4.0*8.0*MAXIC > HUGE(MAXIV)) THEN
        WRITE(6,*) 'Size of GRIDPTR array (8*MAXIC) probably exceeds max integer number of bytes:',8*MAXIC
        STOP
      ENDIF

!         Allocate the internal adaptable grid arrays, including the BIG 3
      ALLOCATE (XGRID(NX+1), YGRID(NY+1), ZGRID(NZ))
      ALLOCATE (TEMP(MAXIG), PLANCK(MAXIG))
      ALLOCATE (EXTINCT(MAXIG), ALBEDO(MAXIG), IPHASE(MAXIG))
      ALLOCATE (LEGEN(MAXIGL))
      ALLOCATE (GRIDPTR(8,MAXIC), NEIGHPTR(6,MAXIC), TREEPTR(2,MAXIC))
      ALLOCATE (CELLFLAGS(MAXIC))
      ALLOCATE (GRIDPOS(3,MAXIG))
      ALLOCATE (RSHPTR(MAXIG+2), SHPTR(MAXIG+1), OSHPTR(MAXIG+1))
      ALLOCATE (RADIANCE(MAXIV+MAXIG), SOURCE(MAXIV))
      IF (ACCELFLAG)  ALLOCATE (DELSOURCE(MAXIV))
      ALLOCATE (WORK(MAXIDO), WORK1(8*MAXIG))
      ALLOCATE  (FLUXES(2*MAXIG), DIRFLUX(MAXIG))


      MAXNBC = MAXIG*3/NZ
!         Read in the surface properties      
      IF (SFCFILE(1:2) .EQ. 'NO') THEN
        SFCTYPE = 'FL'
      ELSE
        IF (MASTERPROC) THEN
          CALL READ_SURFACE_SIZE (SFCFILE, MAXSFCPTS, MAXSFCPARS)
        ENDIF
         ! Broadcast the surface array sizes
        CALL BROADCAST_SURFACE_SIZE (MAXSFCPTS, MAXSFCPARS)
        ALLOCATE (SFCPARMS(MAXSFCPARS*MAXSFCPTS))
        IF (MASTERPROC) THEN
          CALL READ_SURFACE (SFCFILE, MAXSFCPTS, MAXSFCPARS, &
                             SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC, &
                             NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
        ENDIF
         ! Broadcast the surface array data
        CALL BROADCAST_SURFACE_PARMS (SFCTYPE, NXSFC, NYSFC, NSFCPAR, &
                                 DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
        ALLOCATE (SFCGRIDPARMS(MAXSFCPARS*MAXNBC))
      ENDIF
      IF (SFCTYPE(2:2) .EQ. 'L') THEN
        MAXBCRAD=2*MAXNBC
      ELSE
        MAXBCRAD=(2+NMU*NPHI0MAX/2)*MAXNBC
      ENDIF
      ALLOCATE (BCPTR(2*MAXNBC), BCRAD(MAXBCRAD))




!          If doing a k-distribution then get the band info from the CKD file
      IF (KDIST) THEN
        IF (MASTERPROC) THEN
          CALL READ_CKD_SIZE (CKDFILE, WAVENO, NG, NZCKD)
        ENDIF
         ! Broadcast the k-distribution array sizes
        CALL BROADCAST_KDIST_SIZE (NG, NZCKD)
        ALLOCATE (ZCKD(NZCKD), GASABS(NZCKD), DELG(NG), KABS(NZCKD,NG))
        IF (MASTERPROC) THEN
          MAXNG = NG ; MAXNZ = NZCKD
          CALL READ_CKD (MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX, &
                         NG, DELG, NZCKD, ZCKD, KABS)
        ENDIF
        CALL BROADCAST_KDIST_PARMS (SOLFLUX, NG, DELG, NZCKD, ZCKD, KABS)
        SOLARFLUX = SOLFLUX*SOLARFLUX*ABS(SOLARMU)
      ELSE
        NG = 1
        ALLOCATE (DELG(NG))
        DELG(1) = 1.0
        NZCKD = 0
        BASEOUT = .FALSE.
      ENDIF


!         Figure out what type of output we will be generating
      NRAD = 0
      NVIS = 0
      RADIANCE_OUT = .FALSE.
      VISUAL_OUT = .FALSE.
      SOURCE_OUT = .FALSE.
      FLUX_OUT = .FALSE.
      FLUXDIV_OUT = .FALSE.
      SH_OUT = .FALSE.
      DO I = 1, NUMOUT
        IF (OUTTYPES(I) .EQ. 'R') THEN
          RADIANCE_OUT = .TRUE.
          IF (OUTPARMS(2,I) .LE. 0.0) THEN
            NXOUT = 1
          ELSE
            IF (BTEST(BCFLAG,2)) THEN
              NXOUT = MAX(1,NINT(NPXT*DELX/OUTPARMS(2,I)))
            ELSE
              NXOUT=MAX(1,NINT((NPX*DELX+2*OUTPARMS(4,I))/OUTPARMS(2,I)))
            ENDIF
          ENDIF
          IF (OUTPARMS(3,I) .LE. 0.0) THEN
            NYOUT = 1
          ELSE
            IF (BTEST(BCFLAG,3)) THEN
              NYOUT = MAX(1,NINT(NPYT*DELY/OUTPARMS(3,I)))
            ELSE
              NYOUT=MAX(1,NINT((NPY*DELY+2*OUTPARMS(5,I))/OUTPARMS(3,I)))
            ENDIF
          ENDIF
          NRAD = NRAD + NXOUT*NYOUT*NINT(OUTPARMS(6,I))
        ELSE IF (OUTTYPES(I) .EQ. 'V') THEN
          VISUAL_OUT = .TRUE.
          IF (NINT(OUTPARMS(1,I)) .EQ. 1) THEN
            NVIS = NVIS + NINT(OUTPARMS(10,I))*NINT(OUTPARMS(11,I))
          ELSE
            NSCAN = 1 + SQRT((OUTPARMS(4,I)-OUTPARMS(7,I))**2 &
                      +(OUTPARMS(5,I)-OUTPARMS(8,I))**2 &
                      +(OUTPARMS(6,I)-OUTPARMS(9,I))**2) /OUTPARMS(10,I)
            NPIX = 1 + ABS(OUTPARMS(12,I)-OUTPARMS(11,I))/OUTPARMS(13,I)
            NVIS = NVIS + NSCAN*NPIX
          ENDIF
        ELSE IF (OUTTYPES(I) .EQ. 'F') THEN
          IF (.NOT. FLUX_OUT) THEN 
            ALLOCATE (SUMFLUXES(2*MAXIG), SUMDIRFLUX(MAXIG))
          ENDIF
          FLUX_OUT = .TRUE.
        ELSE IF (OUTTYPES(I) .EQ. 'H') THEN
          IF (.NOT. FLUXDIV_OUT) THEN
            ALLOCATE (SUMFLUXDIV(MAXIG))
          ENDIF
          FLUXDIV_OUT = .TRUE.
        ELSE IF (OUTTYPES(I) .EQ. 'S') THEN
          IF (.NOT. SH_OUT) THEN
            ALLOCATE (SUMSHTERMS(NSHOUT*MAXIG))
          ENDIF
          SH_OUT = .TRUE.
        ELSE IF (OUTTYPES(I) .EQ. 'J') THEN
          IF (.NOT. SOURCE_OUT) THEN
            ALLOCATE (SOURCE1OUT(MAXIG))
          ENDIF
          SOURCE_OUT = .TRUE.
        ENDIF
      ENDDO 

      ALLOCATE (SUMRADOUT(NRAD), SUMVISOUT(NVIS), WORK2(MAX(NVIS,NRAD,MAXIG)))

       ! Notify user that certain outputs are not supported with multiple processors
      IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3) .AND. MASTERPROC) THEN
        DO I = 1, NUMOUT
          IF (OUTTYPES(I) .EQ. 'V') THEN
            VISUAL_OUT = .FALSE.
            WRITE (6,'(A,A)') 'V output format not supported for multiple processors: ',TRIM(OUTFILES(I))
          ELSE IF (OUTTYPES(I) .EQ. 'F') THEN
            IF (NINT(OUTPARMS(1,I)) .EQ. 2 .OR. NINT(OUTPARMS(1,I)) .EQ. 5) THEN
              WRITE (6,'(A,I1,A,A)') 'F output format ',NINT(OUTPARMS(1,I)), &
                  ' not supported for multiple processors: ',TRIM(OUTFILES(I))
            ENDIF
          ELSE IF (OUTTYPES(I) .EQ. 'H') THEN
            IF (NINT(OUTPARMS(1,I)) .EQ. 3) THEN
              WRITE (6,'(A,I1,A,A)') 'H output format ',NINT(OUTPARMS(1,I)), &
                  ' not supported for multiple processors: ',TRIM(OUTFILES(I))
            ENDIF
          ELSE IF (OUTTYPES(I) .EQ. 'S') THEN
            IF (NINT(OUTPARMS(1,I)) .EQ. 2) THEN
              WRITE (6,'(A,I1,A,A)') 'S output format ',NINT(OUTPARMS(1,I)), &
                  ' not supported for multiple processors: ',TRIM(OUTFILES(I))
            ENDIF
          ELSE IF (OUTTYPES(I) .EQ. 'J') THEN
            WRITE (6,'(A,A)') 'J output format not supported for multiple processors: ',TRIM(OUTFILES(I))
          ELSE IF (OUTTYPES(I) .EQ. 'M') THEN
            WRITE (6,'(A,A)') 'M output format not supported for multiple processors: ',TRIM(OUTFILES(I))
          ENDIF
        ENDDO 
      ENDIF



!         Read in previous run from binary file if desired
      CALL RESTORE_STATE (INSAVEFILE, NX, NY, NZ, ML, MM, NCS, NLM, &
            INRADFLAG, NEWGRIDFLAG,  XGRID, YGRID, ZGRID, &
            NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, &
            CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

      IF (NEWGRIDFLAG) THEN
!           Make the internal base grid lines
        CALL NEW_GRIDS (BCFLAG, GRIDTYPE, NPX, NPY, NPZ, NX, NY, NZ, &
                        XSTART, YSTART, DELX, DELY, ZLEVELS, &
                        XGRID, YGRID, ZGRID)
      ENDIF



      WRITE (6,*) 'Starting solution iterations for SHDOM run: ',TRIM(RUNNAME)
      TOTITER = 0
      OLDNPTS = 0

!         Loop over the k-distribution g's from high to low absorption
      DO IG = NG, 1, -1
!           If doing k-distribution then get the gas absorption profile.
        IF (KDIST) THEN
          DO I = 1, NZCKD
            GASABS(I) = KABS(I,IG)
          ENDDO
        ENDIF

!           If this is the first time through or we always want a base grid,
!           then make the base grid point positions and grid cell structure.
        IF (BASEOUT .OR. (IG .EQ. NG .AND. NEWGRIDFLAG)) THEN
          CALL INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, &
                   NX, NY, NZ, NX1, NY1, NPTS, NCELLS, &
                   XGRID, YGRID, ZGRID, GRIDPOS, &
                   GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
          NBCELLS = NCELLS
        ENDIF

!           Find the monochromatic SHDOM solution: adaptive grid output
!             in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        CALL SOLVE_RTE (NX, NY, NZ, NX1, NY1, NANG, &
                    ML, MM, NCS, NLM, NMU, NPHI, NLEG, NUMPHASE, &
                    NPHI0, MU, PHI, WTDO, &
                    MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG, &
                    BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD, &
                    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                    SFCTYPE, GNDTEMP, GNDALBEDO, &
                    NXSFC, NYSFC, DELXSFC, DELYSFC, &
                    NSFCPAR, SFCPARMS, SFCGRIDPARMS, &
                    UNITS, WAVENO, WAVELEN, &
                    ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER, &
                    SPLITACC, SHACC,  XGRID,YGRID,ZGRID, &
                    TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE, &
                    MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
                    CMU1, CMU2, CPHI1, CPHI2,  NPTS, GRIDPOS, &
                    NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                    RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2, &
                    SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        INRADFLAG = .FALSE.
        TOTITER = TOTITER + ITER

!         Compute the output and put in arrays.  First the radiance.
        IF (RADIANCE_OUT) THEN
          IRAD = 0
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'R') THEN
              IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
                CALL COMPUTE_RADIANCE_PAR (NX, NY, NZ, NPTS, NCELLS, &
                    ML, MM, NCS, NLEG, NUMPHASE, &
                    NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                    BCFLAG, NPXT*DELX, NPYT*DELY, IPFLAG, &
                    SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
                    SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
                    MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
                    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
                    XGRID, YGRID, ZGRID, GRIDPOS, &
                    GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                    EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
                    SHPTR, SOURCE, WORK, WORK(1+NPTS), &
                    OUTPARMS(1,I),  IRAD, WORK2)
              ELSE
                CALL COMPUTE_RADIANCE (NX, NY, NZ, NPTS, NCELLS, &
                    ML, MM, NCS, NLEG, NUMPHASE, &
                    NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                    BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
                    SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
                    MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
                    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
                    XGRID, YGRID, ZGRID, GRIDPOS, &
                    GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                    EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
                    SHPTR, SOURCE, WORK(1), WORK(1+NPTS), &
                    OUTPARMS(1,I),  IRAD, WORK2)
              ENDIF
            ENDIF
          ENDDO
          CALL SUM_OUTPUT (NG-IG, DELG(IG), NRAD, WORK2, SUMRADOUT)
        ENDIF

        IF (VISUAL_OUT) THEN
          IVIS = 0
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'V') THEN
              CALL VISUALIZE_RADIANCE (NX, NY, NZ, NPTS, NCELLS, &
                    ML, MM, NCS, NLM, NLEG, NUMPHASE, &
                    NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                    BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
                    SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
                    MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
                    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
                    XGRID, YGRID, ZGRID, GRIDPOS, &
                    GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                    EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
                    SHPTR, SOURCE, OUTPARMS(1,I),  IVIS, WORK2)
            ENDIF
          ENDDO
          CALL SUM_OUTPUT (NG-IG, DELG(IG), NVIS, WORK2, SUMVISOUT)
        ENDIF

        IF (SOURCE_OUT) THEN
!             Compute source function throughout grid for this one angle
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'J') THEN
              MUOUT = OUTPARMS(2,I)
              PHIOUT = OUTPARMS(3,I)*ACOS(-1.0)/180.0
              CALL COMPUTE_ONE_SOURCE (ML,MM,NCS, NLEG, NUMPHASE, &
                 NPTS, DELTAM, MUOUT, PHIOUT, &
                 SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE, &
                 DIRFLUX, SHPTR, SOURCE,  SOURCE1OUT)
            ENDIF
          ENDDO
        ENDIF



        IF (BASEOUT) THEN
!           If outputting only the base grid then throw out the other points
          CALL INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, &
                       NX, NY, NZ, NX1, NY1, &
                       NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS, &
                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
        ELSE
!           For the output from the grid point arrays (Flux, Net flux, SH)
!             we may need to interpolate the running k-distribution sum
!             to the newly created grid points.
          CALL INTERP_OUTPUT (OLDNPTS,NPTS,FLUX_OUT,FLUXDIV_OUT,SH_OUT, &
                  SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS, &
                  NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ENDIF

        IF (FLUX_OUT) THEN
          CALL SUM_OUTPUT (NG-IG,DELG(IG), 2*NPTS, FLUXES, SUMFLUXES)
          CALL SUM_OUTPUT (NG-IG,DELG(IG), NPTS, DIRFLUX, SUMDIRFLUX)
        ENDIF
        IF (FLUXDIV_OUT) THEN
!               Compute the net flux divergence at every point
          CALL COMPUTE_NETFLUXDIV (NPTS, RSHPTR, SRCTYPE, SOLARMU, &
                EXTINCT,ALBEDO, PLANCK,DIRFLUX, RADIANCE, WORK)
          CALL SUM_OUTPUT (NG-IG,DELG(IG),NPTS,WORK,SUMFLUXDIV)
        ENDIF

        IF (SH_OUT) THEN
          CALL COMPUTE_SH (NSHOUT, NPTS, SRCTYPE, SOLARMU, SOLARAZ, &
                        DIRFLUX, RSHPTR, ML, MM, NCS, RADIANCE, WORK)
          CALL SUM_OUTPUT (NG-IG,DELG(IG),NSHOUT*NPTS,WORK,SUMSHTERMS)
        ENDIF

        OLDNPTS = NPTS
      ENDDO


      CALL SAVE_STATE (OUTSAVEFILE, NX, NY, NZ, ML, MM, NCS, NLM, &
             WORK, XGRID, YGRID, ZGRID,  &
             NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, &
             CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)


      CALL CPU_TIME(CPUTIME2) 
      CPUTIMETOTAL = SUM_CPU_TIME(CPUTIME2 - CPUTIME1) 
      IF (MASTERPROC) THEN
        WRITE (6,'(A,F10.2)') 'Approximate total CPU time used (sec)', CPUTIMETOTAL 
      ENDIF

      NSH = SHPTR(NPTS+1)
      IRAD = 0
      IVIS = 0
      IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
        IF (MASTERPROC) THEN
          IF (FLUX_OUT)    ALLOCATE (ALLFLUXES(3,NZ,NYT,NXT))
          IF (FLUXDIV_OUT) ALLOCATE (ALLFLUXDIV(NZ,NYT,NXT))
          IF (SH_OUT)      ALLOCATE (ALLSHTERMS(NSHOUT,NZ,NYT,NXT))
          ALLOCATE (XALLGRID(NXT), YALLGRID(NYT))
        ENDIF
        CALL GATHER_OUTPUT (FLUX_OUT, FLUXDIV_OUT, SH_OUT, &
                            IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, &
                            SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS, &
                            NPXT, NPYT, DELX, DELY, XALLGRID, YALLGRID, &
                            ALLFLUXES, ALLFLUXDIV, ALLSHTERMS, &
                            NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

        IF (MASTERPROC) THEN
          DO I = 1, NUMOUT
            IF (OUTFILES(I) /= 'NONE') THEN
              CALL OUTPUT_RESULTS_PAR (NXT, NYT, NZ, NPTSTOT, NCELLSTOT, &
                  NSHTOT, ML,MM,NLM, NMU,NPHI,NANG, NG, &
                  PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE, &
                  BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                  SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                  SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                  SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, &
                  NPXT*DELX, NPYT*DELY, XALLGRID, YALLGRID, ZGRID, &
                  ALLFLUXES, ALLFLUXDIV, &
                  NSHOUT, ALLSHTERMS, IRAD, SUMRADOUT, &
                  OUTTYPES(I), OUTPARMS(1,I), OUTFILES(I))
            ENDIF
          ENDDO
          IF (TRIM(OutFileNC) /= 'NONE' .AND. TRIM(OutFileNC) /= '') THEN
            CALL output_results_netcdf_par (NXT, NYT, NZ, NPTSTOT, NCELLSTOT, &
                    NSHTOT, ML,MM,NLM, NMU,NPHI,NANG, NG, &
                    PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, &
                    BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                    SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                    SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,  CPUTIMETOTAL, &
                    NPXT*DELX, NPYT*DELY, XALLGRID, YALLGRID, ZGRID, &
                    ALLFLUXES, ALLFLUXDIV, NSHOUT, ALLSHTERMS, IRAD, SUMRADOUT, &
                    NUMOUT, OUTTYPES, OUTPARMS, OutFileNC)
          ENDIF
          IF (FLUX_OUT)    DEALLOCATE (ALLFLUXES)
          IF (FLUXDIV_OUT) DEALLOCATE (ALLFLUXDIV)
          IF (SH_OUT)      DEALLOCATE (ALLSHTERMS)
          DEALLOCATE (XALLGRID, YALLGRID)
        ENDIF

      ELSE

        DO I = 1, NUMOUT
          IF (OUTTYPES(I) .EQ. 'V') THEN
            CALL OUTPUT_IMAGE (NX, NY, NZ, NPTS, NCELLS, &
                  NSH, ML, MM, NLM, NMU, NPHI, NANG, NG, &
                  PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE, &
                  BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                  SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                  SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                  SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, &
                  IVIS, SUMVISOUT, OUTPARMS(1,I), OUTFILES(I))
          ELSE
            IF (OUTFILES(I) /= 'NONE') THEN
              CALL OUTPUT_RESULTS (NX, NY, NZ, NBPTS, NPTS, NCELLS, &
                NSH, ML, MM, NLM, NLEG, NUMPHASE, NMU, NPHI, NANG, NG, &
                PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, &
                BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                SPLITACC, SHACC, SOLACC, MAXITER,TOTITER, &
                XGRID, YGRID, ZGRID, GRIDPOS, &
                TREEPTR, GRIDPTR, CELLFLAGS, &
                EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, &
                SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, &
                IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT, &
                OUTTYPES(I), OUTPARMS(1,I), OUTFILES(I))
            ENDIF
          ENDIF
        ENDDO
        IF (TRIM(OutFileNC) /= 'NONE' .AND. TRIM(OutFileNC) /= '') THEN
          ALLOCATE (ALLFLUXES(3,NZ,NY1,NX1), ALLFLUXDIV(NZ,NY1,NX1))
          ALLOCATE (ALLSHTERMS(NSHOUT,NZ,NY1,NX1))
          IF (FLUX_OUT) THEN
            ALLFLUXES(1:2,:,:,:) = RESHAPE(SUMFLUXES(1:2*NBPTS), (/ 2, NZ, NY1, NX1 /) )
            ALLFLUXES(3,:,:,:) = RESHAPE(SUMDIRFLUX(1:NBPTS), (/ NZ, NY1, NX1 /) )
          ENDIF
          IF (FLUXDIV_OUT) THEN
            ALLFLUXDIV(:,:,:) = RESHAPE(SUMFLUXDIV(1:NBPTS), (/ NZ, NY1, NX1 /) )
          ENDIF
          IF (SH_OUT) THEN
            ALLSHTERMS(:,:,:,:) = RESHAPE(SUMSHTERMS(1:NSHOUT*NBPTS), (/ NSHOUT, NZ, NY1, NX1 /) )
          ENDIF
          CALL output_results_netcdf_par (NX, NY, NZ, NPTS, NCELLS, &
                    NSH, ML,MM,NLM, NMU,NPHI,NANG, NG, &
                    PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, &
                    BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                    SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                    SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, CPUTIMETOTAL, &
                    NPX*DELX, NPY*DELY, XGRID, YGRID, ZGRID, &
                    ALLFLUXES(:,:,1:NY,1:NX), ALLFLUXDIV(:,1:NY,1:NX), &
                    NSHOUT, ALLSHTERMS(:,:,1:NY,1:NX), IRAD, SUMRADOUT, &
                    NUMOUT, OUTTYPES, OUTPARMS, OutFileNC)
          DEALLOCATE (ALLFLUXES, ALLFLUXDIV, ALLSHTERMS)
        ENDIF
      ENDIF


       ! Calculate the actual memory parameters and output them to the log file
      MAXMB_OUT = 4*( NMU*(2+2*NPHI+2*NLM+2*33*32) &
                  + 4.5*MAXPG + MAXPGL + NUMPHASE*(MAXLEG+1) &
                  + 16.5*NCELLS + NPTS*(28+NPHI0MAX) + NSH*BIG_ARRAYS)/1024**2
      ADAPTGRIDFAC_OUT = FLOAT(NPTS)/NBPTS
      SHTERM_FAC_OUT = FLOAT(NSH)/(NLM*NPTS)
      CELL_POINT_OUT = FLOAT(NCELLS)/NPTS
      WRITE (6,'(F7.2,1X,F7.3,2(1X,F5.3),A)') &
          MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT, &
          '  Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'

      CALL END_SHDOM_MPI (NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, DELX, DELY,&
                          NPXT, NPYT, PROPFILE, RUNNAME)

      IF (SRCTYPE .NE. 'T') THEN
        CALL DIRECT_BEAM_PROP (9, 0.0, 0.0, 0.0, BCFLAG, IPFLAG, &
              DELTAM,ML,NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1), &
              UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
      ENDIF
      DEALLOCATE (SUMRADOUT, SUMVISOUT, WORK2)
      IF (FLUX_OUT) DEALLOCATE (SUMFLUXES, SUMDIRFLUX)
      IF (FLUXDIV_OUT) DEALLOCATE (SUMFLUXDIV)
      IF (SH_OUT) DEALLOCATE (SUMSHTERMS)
      IF (SOURCE_OUT) DEALLOCATE (SOURCE1OUT)
      IF (KDIST) DEALLOCATE (ZCKD, GASABS, KABS)
      DEALLOCATE (BCPTR, BCRAD, DELG)
      IF (ALLOCATED(SFCPARMS)) DEALLOCATE (SFCPARMS, SFCGRIDPARMS)
      DEALLOCATE (RADIANCE, SOURCE)
      IF (ACCELFLAG)  DEALLOCATE (DELSOURCE)
      DEALLOCATE (RSHPTR, SHPTR, OSHPTR, WORK, WORK1, FLUXES, DIRFLUX)
      DEALLOCATE (GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)
      DEALLOCATE (TEMP, PLANCK, EXTINCT, ALBEDO, IPHASE, LEGEN)
      DEALLOCATE (XGRID, YGRID, ZGRID)
      DEALLOCATE (ZLEVELS, TEMPP, EXTINCTP, ALBEDOP, EXTDIRP, IPHASEP, LEGENP)
      DEALLOCATE (MU, WTDO, PHI, NPHI0, CMU1, CMU2, CPHI1, CPHI2)

      END
 




 
      SUBROUTINE TRILIN_INTERP_PROP (X, Y, Z, INIT, NLEG, &
                     TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
!      Trilinearly interpolates the quantities on the input property
!     grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
!     ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the 
!     volume coefficients.  Also adds in the separate gaseous absorption.
!     Divides the phase function Legendre coefficients LEGEN by 2*l+1.
!     The phase function pointer IPHASE (for tabulated phase functions) 
!     is that of the maximum weighted scattering property grid point.
!     If INIT=.TRUE. then transfers the tabulated phase functions.

      USE SHDOM_PROPERTY_ARRAYS
      IMPLICIT NONE
      INTEGER NLEG
      INTEGER*2 IPHASE
      LOGICAL INIT
      REAL    X, Y, Z,  TEMP, EXTINCT, ALBEDO, LEGEN(0:NLEG,*)
      INTEGER IX, IXP, IY, IYP, IZ, L, IL, IM, IU
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I
      DOUBLE PRECISION U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
      DOUBLE PRECISION SCAT1,SCAT2,SCAT3,SCAT4,SCAT5,SCAT6,SCAT7,SCAT8
      DOUBLE PRECISION SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
      SAVE EXTMIN, SCATMIN

      IF (INIT) THEN
!         If there are tabulated phase functions, then transfer them
        DO I = 1, NUMPHASE
          LEGEN(0,I) = 1.0
          DO L = 1, NLEG
            LEGEN(L,I) = LEGENP(L+NLEG*(I-1))/(2*L+1)
          ENDDO
        ENDDO
        EXTMIN = 1.0E-5 / ( (ZLEVELS(NPZ)-ZLEVELS(1))/NPZ)
        SCATMIN = 0.1*EXTMIN
        RETURN
      ENDIF

!         Find the grid location and compute the interpolation factors
      IL=0
      IU=NPZ
      DO WHILE (IU-IL .GT. 1)
        IM = (IU+IL)/2
        IF (Z .GE. ZLEVELS(IM)) THEN
          IL = IM
        ELSE
          IU=IM
        ENDIF
      ENDDO
      IZ = MAX(IL,1)
      W = DBLE(Z - ZLEVELS(IZ))/(ZLEVELS(IZ+1) - ZLEVELS(IZ))
      W = MAX( MIN( W, 1.0D0), 0.0D0)
      IX = INT((X-XSTART)/DELX) + 1
      IF (ABS(X-XSTART-NPX*DELX) .LT. 0.01*DELX) IX = NPX
      IF (IX .LT. 1 .OR. IX .GT. NPX) THEN
        WRITE (6,*) 'TRILIN: Beyond X domain'
        STOP
      ENDIF
      IXP = MOD(IX,NPX) + 1
      U = DBLE(X-XSTART-DELX*(IX-1))/DELX
      U = MAX( MIN( U, 1.0D0), 0.0D0)
      IF (U .LT. 1.0D-5) U = 0.0D0
      IF (U .GT. 1.0D0-1.0D-5) U = 1.0D0
      IY = INT((Y-YSTART)/DELY) + 1
      IF (ABS(Y-YSTART-NPY*DELY) .LT. 0.01*DELY) IY = NPY
      IF (IY .LT. 1 .OR. IY .GT. NPY) THEN
        WRITE (6,*) 'TRILIN: Beyond Y domain'
        STOP
      ENDIF
      IYP = MOD(IY,NPY) + 1
      V = DBLE(Y-YSTART-DELY*(IY-1))/DELY
      V = MAX( MIN( V, 1.0D0), 0.0D0)
      IF (V .LT. 1.0D-5) V = 0.0D0
      IF (V .GT. 1.0D0-1.0D-5) V = 1.0D0

      F1 = (1-U)*(1-V)*(1-W)
      F2 =    U *(1-V)*(1-W)
      F3 = (1-U)*   V *(1-W)
      F4 =    U *   V *(1-W)
      F5 = (1-U)*(1-V)*   W
      F6 =    U *(1-V)*   W
      F7 = (1-U)*   V *   W
      F8 =    U *   V *   W
      I1 = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
      I2 = IZ + NPZ*(IY-1) + NPZ*NPY*(IXP-1)
      I3 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IX-1)
      I4 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IXP-1)
      I5 = I1+1
      I6 = I2+1
      I7 = I3+1
      I8 = I4+1

!         Trilinearly interpolate the temperature, extinction, scattering
      TEMP = F1*TEMPP(I1) + F2*TEMPP(I2) + F3*TEMPP(I3) + F4*TEMPP(I4) &
          + F5*TEMPP(I5) + F6*TEMPP(I6) + F7*TEMPP(I7) + F8*TEMPP(I8)
      EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3) &
              + F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6) &
              + F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
      SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
      SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
      SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
      SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
      SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
      SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
      SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
      SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
      SCATTER = SCAT1+SCAT2+SCAT3+SCAT4+SCAT5+SCAT6+SCAT7+SCAT8
      IF (EXTINCT .GT. EXTMIN) THEN
        ALBEDO = SCATTER/EXTINCT
      ELSE
        ALBEDO = SCATTER/EXTMIN
      ENDIF
!         For tabulated phase functions pick the one we are on top of
!         or the one with the most scattering weight.
      IF (NUMPHASE .GT. 0) THEN
        MAXSCAT = -1.0
        IF (SCAT1 .GT. MAXSCAT .OR. ABS(F1-1) .LT. 0.001) THEN
          MAXSCAT = SCAT1
          IPHASE = IPHASEP(I1)
        ENDIF
        IF (SCAT2 .GT. MAXSCAT .OR. ABS(F2-1) .LT. 0.001) THEN
          MAXSCAT = SCAT2
          IPHASE = IPHASEP(I2)
        ENDIF
        IF (SCAT3 .GT. MAXSCAT .OR. ABS(F3-1) .LT. 0.001) THEN
          MAXSCAT = SCAT3
          IPHASE = IPHASEP(I3)
        ENDIF
        IF (SCAT4 .GT. MAXSCAT .OR. ABS(F4-1) .LT. 0.001) THEN
          MAXSCAT = SCAT4
          IPHASE = IPHASEP(I4)
        ENDIF
        IF (SCAT5 .GT. MAXSCAT .OR. ABS(F5-1) .LT. 0.001) THEN
          MAXSCAT = SCAT5
          IPHASE = IPHASEP(I5)
        ENDIF
        IF (SCAT6 .GT. MAXSCAT .OR. ABS(F6-1) .LT. 0.001) THEN
          MAXSCAT = SCAT6
          IPHASE = IPHASEP(I6)
        ENDIF
        IF (SCAT7 .GT. MAXSCAT .OR. ABS(F7-1) .LT. 0.001) THEN
          MAXSCAT = SCAT7
          IPHASE = IPHASEP(I7)
        ENDIF
        IF (SCAT8 .GT. MAXSCAT .OR. ABS(F8-1) .LT. 0.001) THEN
          MAXSCAT = SCAT8
          IPHASE = IPHASEP(I8)
        ENDIF
      ELSE
        LEGEN(0,1) = 1.0
        DO L = 1, NLEG
          LEGEN(L,1) = SCAT1*LEGENP(L+NLEG*(I1-1)) &
                    + SCAT2*LEGENP(L+NLEG*(I2-1)) &
                    + SCAT3*LEGENP(L+NLEG*(I3-1)) &
                    + SCAT4*LEGENP(L+NLEG*(I4-1)) &
                    + SCAT5*LEGENP(L+NLEG*(I5-1)) &
                    + SCAT6*LEGENP(L+NLEG*(I6-1)) &
                    + SCAT7*LEGENP(L+NLEG*(I7-1)) &
                    + SCAT8*LEGENP(L+NLEG*(I8-1))
          IF (SCATTER .GT. SCATMIN) THEN
            LEGEN(L,1) = LEGEN(L,1)/SCATTER
          ELSE
            LEGEN(L,1) = LEGEN(L,1)/SCATMIN
          ENDIF
          LEGEN(L,1) = LEGEN(L,1)/(2*L+1)
        ENDDO
      ENDIF

!         Add in the gaseous absorption to extinction and albedo
      IF (NZCKD .GT. 0) THEN
        IL = 1
        IU = NZCKD
        DO WHILE (IU-IL .GT. 1)
          IM=(IU+IL)/2
          IF (Z .LE. ZCKD(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        I = MIN(MAX(IL,1),NZCKD-1)
        F = (Z-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
        F = MIN( MAX(F,0.0D0), 1.0D0)
        KG = (1.0-F)*GASABS(I) + F*GASABS(I+1)
        IF (EXTINCT+KG .GT. 0.0) THEN
          ALBEDO = ALBEDO*EXTINCT /(EXTINCT + KG)
        ELSE
          ALBEDO = 0.0
        ENDIF
        EXTINCT = EXTINCT + KG
      ENDIF
      RETURN
      END
 
 


 
      SUBROUTINE DIRECT_BEAM_PROP (INIT, XI, YI, ZI, BCFLAG, IPFLAG, &
                     DELTAM, ML, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, & 
                     DIRFLUX,  UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
!       Computes the direct beam flux at point (XI,YI,ZI) by integrating
!     the extinction through the property grid.  If called with 
!     INIT=1 then the property grid extinction array, solar direction
!     terms, and lowest uniform level are computed and saved.  If called
!     with INIT=2 then the input lowest uniform level (UNIFZLEV) is stored.
!     Call with INIT=0 to do the path integration and return the direct
!     beam flux on the horizontal (DIRFLUX) at XI,YI,ZI.
!     The DELTAM flag is used to determine whether to delta-M scale the 
!     extinction for the direct beam computation.  The extinction includes 
!     the gaseous absorption.  If the IPFLAG has bit 2 set then the direct
!     beam tracing is done in 3D, otherwise the lower two bits determine
!     the type of tracing: 0 for 3D tracing, 1 for XZ only tracing, 
!     2 for YZ only tracing, and 3 for Z only tracing.  If BCFLAG bits 0 
!     or 1 are set then have open boundary conditions in X and/or Y.  In 
!     this case when the ray tracing back to the sun reaches the boundary 
!     then independent pixel mode is entered so that only Z grid 
!     intersections occur.
!       For use with multiple processors (bit 2 or 3 set in BCFLAG), 
!     the VALIDBEAM flag is returned true if the ray to the sun made it 
!     to the top of the domain before hitting the subdomain side.  If the
!     flag is false SIDE is returned with the boundary hit (1=-X, 2=+X,
!     3=-Y, 4=+Y).  XE,YE,ZE returns the location of the exitting ray,
!     and path is the optical path from XI,YI,ZI to the sun.
      USE SHDOM_PROPERTY_ARRAYS
      IMPLICIT NONE
      INTEGER INIT, BCFLAG, IPFLAG, ML, NLEG, SIDE
      LOGICAL DELTAM, VALIDBEAM
      REAL    XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
      REAL    DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

      INTEGER IX, IY, IZ, JZ, IL, IM, IU
      INTEGER I, J, K, IP, JP, I1, I2, I3, I4
      INTEGER IPDIRECT, DI, DJ, DK
      LOGICAL CONSTX, CONSTY, HITBOUNDARY, BTEST
      DOUBLE PRECISION EXTINCT, ALBEDO, F
      DOUBLE PRECISION SUNMU, SUNAZ, PATH
      DOUBLE PRECISION EXTBEAMCUT, UNIFORMZLEV
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION XOFFS, YOFFS, DELXD, DELYD
      DOUBLE PRECISION X, Y, Z, XE, YE, ZE, XP, YP, ZP
      DOUBLE PRECISION X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
      DOUBLE PRECISION U0, V0, W0, U1, V1, W1, AX, AY, AZ
      DOUBLE PRECISION U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
      DOUBLE PRECISION E1,E2,E3,E4,E5,E6,E7,E8,  A, B, C, D
      DOUBLE PRECISION B1,B2,B3,B4,B5,B6,B7,B8,C1,C2,C3,C4,C5,C6,C7,C8 
      DOUBLE PRECISION UV,UMV,UVM,UMVM,UW,UMW,UWM,UMWM,VW,VMW,VWM,VMWM
      DOUBLE PRECISION VWU,VWUM,UWV,UWVM,UVW,UVWM
      REAL, SAVE, ALLOCATABLE :: GASEXT(:), EXTMIN(:), EXTMAX(:)

      SAVE  CX, CY, CZ, CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT
      SAVE  DELXD, DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV


      IF (INIT .EQ. 9) THEN
        DEALLOCATE (GASEXT, EXTMIN, EXTMAX)
        RETURN
      ENDIF

      IF (INIT .EQ. 1) THEN
        IF (.NOT. ALLOCATED(GASEXT)) ALLOCATE (GASEXT(NPZ), EXTMIN(NPZ), EXTMAX(NPZ))
!           Get the gaseous extinction at the property grid levels
        DO IZ = 1, NPZ
          IF (NZCKD .GT. 0) THEN
            IL = 1
            IU = NZCKD
            DO WHILE (IU-IL .GT. 1)
              IM=(IU+IL)/2
              IF (ZLEVELS(IZ) .LE. ZCKD(IM)) THEN
                IL = IM
              ELSE
                IU = IM
              ENDIF
            ENDDO
            I = MIN(MAX(IL,1),NZCKD-1)
            W0 = (ZLEVELS(IZ)-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
            W0 = MIN( MAX(W0,0.0D0), 1.0D0)
            GASEXT(IZ) = (1.0-W0)*GASABS(I) + W0*GASABS(I+1)
          ELSE
            GASEXT(IZ) = 0.0
          ENDIF
        ENDDO
!           First make the property grid extinction field, on which all
!           the direct beam paths will be computed.
        IP = 0
        DO IX = 1, NPX
          DO IY = 1, NPY
            DO IZ = 1, NPZ
              IP = IP + 1
              EXTINCT = EXTINCTP(IP)
              ALBEDO = ALBEDOP(IP)
!                 Add in the gaseous absorption to extinction
              IF (GASEXT(IZ) .GT. 0.0) THEN
                ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                EXTINCT = EXTINCT + GASEXT(IZ)
              ENDIF
!                 Do the Delta-M scaling if needed
              IF (DELTAM) THEN
                IF (NUMPHASE .GT. 0) THEN
                  F = LEGENP(ML+1+NLEG*(IPHASEP(IP)-1))/(2*ML+3)
                ELSE
                  F = LEGENP(ML+1+NLEG*(IP-1))/(2*ML+3)
                ENDIF
                EXTINCT = (1.0-ALBEDO*F)*EXTINCT
              ENDIF
              EXTDIRP(IP) = EXTINCT
            ENDDO
          ENDDO
        ENDDO

!           Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        IF (BTEST(IPFLAG,2)) IPDIRECT = 0

!           Make the ray direction (opposite to the solar beam)
!           SOLAR?? is direction of beam, SUN?? is direction to sun (opposite)
        SUNMU = -SOLARMU
        SUNAZ = SOLARAZ + ACOS(-1.0)
        CX = SQRT(1.0-SUNMU**2)*COS(SUNAZ)
        CY = SQRT(1.0-SUNMU**2)*SIN(SUNAZ)
        CZ = ABS(SUNMU)
        IF (ABS(CX) .GT. 1.0E-6) THEN
          CXINV = 1.0D0/CX
        ELSE
          CX = 0.0
          CXINV = 1.0E20
        ENDIF
        IF (ABS(CY) .GT. 1.0E-6) THEN
          CYINV = 1.0D0/CY
        ELSE
          CY = 0.0
          CYINV = 1.0E20
        ENDIF
        IF (ABS(CZ) .GT. 1.0E-6) THEN
          CZINV = 1.0D0/CZ
        ELSE
          CZ = 0.0
          CZINV = 1.0E20
        ENDIF
        DI = NINT(SIGN(1.0D0,CX))
        DJ = NINT(SIGN(1.0D0,CY))
        DK = NINT(SIGN(1.0D0,CZ))
        EPSZ = 1.0E-6*(ZLEVELS(NPZ)-ZLEVELS(1))
        EPSS = 1.0E-3*(ZLEVELS(NPZ)-ZLEVELS(1))/NPZ
        IF (.NOT. BTEST(IPDIRECT,0))  EPSS = MAX(EPSS,1.0D-4*DELX)
        IF (.NOT. BTEST(IPDIRECT,1))  EPSS = MAX(EPSS,1.0D-4*DELY)
        DELXD = DBLE(DELX)
        DELYD = DBLE(DELY)
        EPSS = MAX(0.001*DELXD,0.001*DELYD,EPSS)
        XDOMAIN = DELXD*NPX
        IF (BTEST(BCFLAG,2)) XDOMAIN = DELXD*(NPX-1)
        YDOMAIN = DELYD*NPY
        IF (BTEST(BCFLAG,3)) YDOMAIN = DELYD*(NPY-1)

!           Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0E-4
        DO IZ = 1, NPZ
          EXTMIN(IZ) = 1.0E20
          EXTMAX(IZ) = 0.0
        ENDDO
        IP = 0
        DO IX = 1, NPX
          DO IY = 1, NPY
            DO IZ = 1, NPZ
              IP = IP + 1
              EXTMIN(IZ) = MIN(EXTINCTP(IP), EXTMIN(IZ))
              EXTMAX(IZ) = MAX(EXTINCTP(IP), EXTMAX(IZ))
            ENDDO
          ENDDO
        ENDDO
        JZ = 0
        DO IZ = 1, NPZ
          IF (EXTMAX(IZ)-EXTMIN(IZ) .GT. EXTBEAMCUT)  JZ = IZ
        ENDDO
        JZ = MIN(NPZ,JZ+1)
        UNIFORMZLEV = ZLEVELS(JZ)
        UNIFZLEV = UNIFORMZLEV
        RETURN
!           Done with initialization
      ENDIF

      IF (INIT .EQ. 2) THEN
        UNIFORMZLEV = UNIFZLEV
        RETURN
      ENDIF



!         Here for computing the direct beam path for one starting point.
      Z = ZI
      X = XI - XSTART
      Y = YI - YSTART

!         Find the grid location 
      IL=0
      IU=NPZ
      DO WHILE (IU-IL .GT. 1)
        IM = (IU+IL)/2
        IF (Z .GE. ZLEVELS(IM)) THEN
          IL = IM
        ELSE
          IU=IM
        ENDIF
      ENDDO
      K = MAX(IL,1)
      I = INT(X/DELXD) + 1
      IF (I .GT. NPX .AND. ABS(X-XDOMAIN) .LT. 0.001*DELXD) I = NPX
      IF (I .LT. 1 .OR. I .GT. NPX) THEN
        WRITE (6,*) 'DIRECT_BEAM_PROP: Beyond X domain',I,XI,YI,ZI
        STOP
      ENDIF
      J = INT(Y/DELYD) + 1
      IF (J .GT. NPY .AND. ABS(Y-YDOMAIN) .LT. 0.001*DELYD) J = NPY
      IF (J .LT. 1 .OR. J .GT. NPY) THEN
        WRITE (6,*) 'DIRECT_BEAM_PROP: Beyond Y domain',J,XI,YI,ZI
        STOP
      ENDIF
      XE = X
      YE = Y
      ZE = Z
      XP = XE
      YP = YE
      ZP = ZE
      CONSTX = BTEST(IPDIRECT,0)
      CONSTY = BTEST(IPDIRECT,1)
      IF (CX .EQ. 0.0)  CONSTX = .TRUE.
      IF (CY .EQ. 0.0)  CONSTY = .TRUE.
      IF (BTEST(BCFLAG,0) .AND. (ABS(X) .LT. 0.01*DELXD  &
         .OR. ABS(X-(NPX-1)*DELXD) .LT. 0.01*DELXD))  CONSTX = .TRUE.
      IF (BTEST(BCFLAG,1) .AND. (ABS(Y) .LT. 0.01*DELYD  &
         .OR. ABS(Y-(NPY-1)*DELYD) .LT. 0.01*DELYD))  CONSTY = .TRUE.

       ! If have multiple subdomains (processors) and ray is going outwards 
       !   from a boundary then set the HITBOUNDARY flag
      HITBOUNDARY = .FALSE.
      IF (BTEST(BCFLAG,2)) THEN
        IF (CX .GT. 0.0 .AND. ABS(X-XDOMAIN) .LT. 0.001*DELXD) THEN
          SIDE = 2
          HITBOUNDARY = .TRUE.
        ELSE IF (CX .LT. 0.0 .AND. ABS(X) .LT.  0.001*DELXD) THEN
          SIDE = 1
          HITBOUNDARY = .TRUE.
        ENDIF
      ENDIF 
      IF (BTEST(BCFLAG,3)) THEN
        IF (CY .GT. 0.0 .AND. ABS(Y-YDOMAIN) .LT. 0.001*DELYD) THEN
          SIDE = 4
          HITBOUNDARY = .TRUE.
        ENDIF
        IF (CY .LT. 0.0 .AND. ABS(Y) .LT.  0.001*DELYD) THEN
          SIDE = 3
          HITBOUNDARY = .TRUE.
        ENDIF
      ENDIF 

!           Grid cell loop begin
      PATH = DIRPATH
      DO WHILE (.NOT. HITBOUNDARY .AND. ABS(ZE-ZLEVELS(NPZ)) .GT. EPSZ)
        IP = I + 1
        IF (I .EQ. NPX) THEN
          IF (BTEST(BCFLAG,0) .OR. BTEST(BCFLAG,2)) THEN
            IP=NPX
          ELSE
            IP = 1
          ENDIF   
        ENDIF     
        JP = J + 1
        IF (J .EQ. NPY) THEN
          IF (BTEST(BCFLAG,1) .OR. BTEST(BCFLAG,3)) THEN
            JP=NPY
          ELSE
            JP = 1
          ENDIF   
        ENDIF     
        X0 = DELXD*(I-1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J-1)
        Y1 = Y0 + DELYD
        Z0 = ZLEVELS(K)
        Z1 = ZLEVELS(K+1)
        IF (I .LT. 1 .OR. I .GT. NPX .OR. &
            J .LT. 1 .OR. J .GT. NPY .OR. &
            K .LT. 1 .OR. K .GE. NPZ) THEN
          WRITE (6,'(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
          WRITE (6,'(1(2X,3F9.5))') X, Y, Z
          STOP
        ENDIF
!           Get the eight corner extinction values
        I1 = K + NPZ*(J-1)  + NPZ*NPY*(I-1)
        I2 = K + NPZ*(J-1)  + NPZ*NPY*(IP-1)
        I3 = K + NPZ*(JP-1) + NPZ*NPY*(I-1)
        I4 = K + NPZ*(JP-1) + NPZ*NPY*(IP-1)
        E1 = EXTDIRP(I1)
        E2 = EXTDIRP(I2)
        E3 = EXTDIRP(I3)
        E4 = EXTDIRP(I4)
        E5 = EXTDIRP(I1+1)
        E6 = EXTDIRP(I2+1)
        E7 = EXTDIRP(I3+1)
        E8 = EXTDIRP(I4+1)

!           Compute the distance to the next grid plane in  X, Y, and Z
!             If in horizontal uniform region or doing IP then fix X and/or Y.
        IF (ZE .GE. UNIFORMZLEV) THEN
          CONSTX = .TRUE.
          CONSTY = .TRUE.
        ENDIF
        IF (CONSTX) THEN
          SOX = 1.0E30
        ELSE IF (CX .GT. 0.0) THEN
          SOX = (X1-XE)*CXINV
          XP = X1
        ELSE
          SOX = (X0-XE)*CXINV
          XP = X0
        ENDIF
        IF (CONSTY) THEN
          SOY = 1.0E30
        ELSE IF (CY .GT. 0.0) THEN
          SOY = (Y1-YE)*CYINV
          YP = Y1
        ELSE
          SOY = (Y0-YE)*CYINV
          YP = Y0
        ENDIF
        IF (CZ .GT. 0.0) THEN
          SOZ = (Z1-ZE)*CZINV
          ZP = Z1
        ELSE IF (CZ .LT. 0.0) THEN
          SOZ = (Z0-ZE)*CZINV
          ZP = Z0
        ELSE
          SOZ = 1.0E30
        ENDIF
        
!           The shortest distance is the plane we stop at:
!             get the exitting location and increment the cell
        XOFFS = 0.0
        YOFFS = 0.0
        IF (SOZ .LE. SOX .AND. SOZ .LE. SOY) THEN
          SO = SOZ
          IF (.NOT. CONSTX) XP = XE + SO*CX
          IF (.NOT. CONSTY) YP = YE + SO*CY
          K = K + DK
        ELSE IF (SOX .LE. SOY) THEN
          SO = SOX
          IF (.NOT. CONSTY) YP = YE + SO*CY
          ZP = ZE + SO*CZ
          I = I + DI
!             If have reached a horizontal boundary then either wrap around
!               (periodic) or go into IP mode (open boundaries).
          IF (I .EQ. 0) THEN
            IF (BTEST(BCFLAG,0)) THEN
              I = 1
              CONSTX = .TRUE.
            ELSE IF (BTEST(BCFLAG,2)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 1
            ELSE
              I = NPX
              XOFFS = XDOMAIN
            ENDIF
          ELSE IF (I .GE. NPX .AND. BTEST(BCFLAG,2)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 2
          ELSE IF (I .EQ. NPX+1) THEN
            IF (BTEST(BCFLAG,0)) THEN
              I = NPX
              CONSTX = .TRUE.
            ELSE
              I = 1
              XOFFS = -XDOMAIN
            ENDIF
          ENDIF
        ELSE
          SO = SOY
          IF (.NOT. CONSTX) XP = XE + SO*CX
          ZP = ZE + SO*CZ
          J = J + DJ
          IF (J .EQ. 0) THEN
            IF (BTEST(BCFLAG,1)) THEN
              J = 1
              CONSTY = .TRUE.
            ELSE IF (BTEST(BCFLAG,3)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 3
            ELSE
              J = NPY
              YOFFS = YDOMAIN
            ENDIF
          ELSE IF (J .GE. NPY .AND. BTEST(BCFLAG,3)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 4
          ELSE IF (J .EQ. NPY+1) THEN
            IF (BTEST(BCFLAG,1)) THEN
              J = NPY
              CONSTY = .TRUE.
            ELSE
              J = 1
              YOFFS = -YDOMAIN
            ENDIF
          ENDIF
        ENDIF
        IF (SO .LT. -EPSS) THEN
          WRITE (6,*) 'DIRECT_BEAM_PROP: SO<0', X,Y,Z, &
              XE,YE,ZE, XP,YP,ZP, CX,CY,CZ, SOX,SOY,SOZ
          STOP
        ENDIF
        SO = MAX(SO,0.0D0)
!           Make the starting and ending interpolation factors
        AX = 1.0D0/(X1-X0)
        AY = 1.0D0/(Y1-Y0)
        AZ = 1.0D0/(Z1-Z0)
        U0 = (XE-X0)*AX
        V0 = (YE-Y0)*AY
        W0 = (ZE-Z0)*AZ
        U1 = (XP-X0)*AX
        V1 = (YP-Y0)*AY
        W1 = (ZP-Z0)*AZ
!           Compute the cubic polynomial extinction coefficients
        U0M = 1.0-U0
        V0M = 1.0-V0
        W0M = 1.0-W0
        U1M = 1.0-U1
        V1M = 1.0-V1
        W1M = 1.0-W1
        DU = U1-U0
        DV = V1-V0
        DW = W1-W0
        UV   = U0 *V0
        UMV  = U0M*V0
        UVM  = U0 *V0M
        UMVM = U0M*V0M
        UW   = U0 *W0
        UMW  = U0M*W0
        UWM  = U0 *W0M
        UMWM = U0M*W0M
        VW   = V0 *W0
        VMW  = V0M*W0
        VWM  = V0 *W0M
        VMWM = V0M*W0M
        A =  (E1*U0M + E2*U0)*VMWM &
           + (E3*U0M + E4*U0)*VWM &
           + (E5*U0M + E6*U0)*VMW &
           + (E7*U0M + E8*U0)*VW
        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 =  DU*VMWM - DV*UWM  - DW*UVM
        B3 = -DU*VWM  + DV*UMWM - DW*UMV
        B4 =  DU*VWM  + DV*UWM  - DW*UV
        B5 = -DU*VMW  - DV*UMW  + DW*UMVM
        B6 =  DU*VMW  - DV*UW   + DW*UVM
        B7 = -DU*VW   + DV*UMW  + DW*UMV
        B8 =  DU*VW   + DV*UW   + DW*UV
        B = B1*E1 +B2*E2 +B3*E3 +B4*E4 +B5*E5 +B6*E6 +B7*E7 +B8*E8 
        VW = DV*DW
        VWU  = VW*U0
        VWUM = VW*U0M
        UW = DU*DW
        UWV  = UW*V0
        UWVM = UW*V0M
        UV = DU*DV
        UVW  = UV*W0
        UVWM = UV*W0M
        C1 = + VWUM + UWVM + UVWM
        C2 = + VWU  - UWVM - UVWM
        C3 = - VWUM + UWV  - UVWM
        C4 = - VWU  - UWV  + UVWM
        C5 = - VWUM - UWVM + UVW
        C6 = - VWU  + UWVM - UVW
        C7 = + VWUM - UWV  - UVW
        C8 = + VWU  + UWV  + UVW
        C = C1*E1 +C2*E2 +C3*E3 +C4*E4 +C5*E5 +C6*E6 +C7*E7 +C8*E8 
        D = DU*DV*DW*(E2+E3+E5+E8-E1-E4-E6-E7)
!            Compute the path through the cell: integration of extinction
        PATH = PATH + SO*(A +0.5D0*B +0.3333333333333333D0*C +0.25D0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
      ENDDO
      DIRFLUX = SOLARFLUX*EXP(-PATH)
      DIRPATH = PATH
      XO=XE+XSTART ; YO=YE+YSTART ; ZO=ZE
      VALIDBEAM = ABS(ZE-ZLEVELS(NPZ)) .LT. EPSZ
      IF (VALIDBEAM) SIDE=6
      RETURN
      END

