      PROGRAM SHDOM 
C      Spherical harmonic discrete ordinate radiative transfer method.
C      See shdom.doc for documentation.

C       All the large arrays are allocated here and are in two common blocks.
C       Don't forget the property array parameter statements
C       in TRILIN_INTERP_PROP and DIRECT_BEAM_PROP below.
      IMPLICIT NONE
      INTEGER    MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO
      PARAMETER (MAXPG=200000)
      PARAMETER (MAXIG=500000, MAXIC=1.5*MAXIG)
      PARAMETER (MAXIV=32000000, MAXIDO=32*MAXIG)
      INTEGER    MAXLEG, MAXPGL, MAXIGL
      PARAMETER (MAXLEG=3000)
      PARAMETER (MAXPGL=MAXLEG*300, MAXIGL=(MAXLEG+1)*300)
      INTEGER   MAXNBC, MAXBCRAD, MAXSFCPTS, MAXSFCPARS
      PARAMETER (MAXNBC=MAXIG/10, MAXBCRAD=(2+64)*MAXNBC)
      PARAMETER (MAXSFCPTS=MAXPG/10, MAXSFCPARS=4)
      INTEGER   MAXNCMU, MAXNCPHI
      PARAMETER (MAXNCMU=64*4096, MAXNCPHI=64*33*32)
      INTEGER   MAXNMU, MAXNPHI
      PARAMETER (MAXNMU=64, MAXNPHI=128)
      INTEGER   MAXNXY, MAXNZ, MAXRADOUT, MAXOUT, MAXPAR, MAXNG
      PARAMETER (MAXNXY=1000, MAXNZ=501)
      PARAMETER (MAXRADOUT=100000, MAXOUT=20, MAXPAR=400, MAXNG=16)
      INTEGER BCPTR(MAXNBC*2), NPHI0(MAXNMU)
      INTEGER RSHPTR(MAXIG+2), SHPTR(MAXIG+1), OSHPTR(MAXIG+1)
      INTEGER GRIDPTR(8*MAXIC)
      INTEGER NEIGHPTR(6*MAXIC), TREEPTR(2*MAXIC)
      INTEGER WORK1(8*MAXIG)
      INTEGER*2 CELLFLAGS(MAXIC)
      INTEGER*2 IPHASEP(MAXPG), IPHASE(MAXIG)
      REAL  OUTPARMS(MAXPAR,MAXOUT)
      REAL  ZCKD(MAXNZ), DELG(MAXNG), KABS(MAXNZ*MAXNG), GASABS(MAXNZ)
      REAL  MU(MAXNMU), WTDO(MAXNMU*MAXNPHI), PHI(MAXNMU*MAXNPHI)
      REAL  CMU1(MAXNCMU),CMU2(MAXNCMU),CPHI1(MAXNCPHI),CPHI2(MAXNCPHI)
      REAL  SFCPARMS(MAXSFCPARS*MAXSFCPTS), BCRAD(MAXBCRAD)
      REAL  SFCGRIDPARMS(MAXSFCPARS*MAXNBC)
      REAL  ZLEVELS(MAXNZ)
      REAL  XGRID(MAXNXY+1), YGRID(MAXNXY+1), ZGRID(MAXNZ)
      REAL  TEMPP(MAXPG),EXTINCTP(MAXPG),ALBEDOP(MAXPG), EXTDIRP(MAXPG)
      REAL  LEGENP(MAXPGL)
      REAL  GRIDPOS(3*MAXIG)
      REAL  TEMP(MAXIG), PLANCK(MAXIG)
      REAL  EXTINCT(MAXIG), ALBEDO(MAXIG), LEGEN(MAXIGL)
      REAL  FLUXES(2*MAXIG),    DIRFLUX(MAXIG)
      REAL  SUMFLUXES(2*MAXIG), SUMDIRFLUX(MAXIG)
      REAL  SUMFLUXDIV(MAXIG),  SUMSHTERMS(4*MAXIG), SOURCE1OUT(MAXIG)
      REAL  SUMRADOUT(MAXRADOUT), SUMVISOUT(MAXRADOUT)
      REAL  WORK2(MAXRADOUT+MAXIG), WORK(MAXIDO)
      REAL  RADIANCE(MAXIV+MAXIG), SOURCE(MAXIV), DELSOURCE(MAXIV)
      CHARACTER  OUTTYPES(MAXOUT)*1, OUTFILES(MAXOUT)*64

      INTEGER NPX, NPY, NPZ, NX, NY, NZ, NX1, NY1
      INTEGER NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
      INTEGER NZCKD, NG, NUMPHASE
      INTEGER ML, MM, NLM, NLEG, NANG, NMU, NPHI, NPHI0MAX, NCS
      INTEGER NBPTS, NPTS, NCELLS, NSH, NBCELLS, OLDNPTS
      INTEGER MAXITER, ITER, TOTITER, IG, I, NUMOUT, BCFLAG, IPFLAG
      INTEGER NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
      LOGICAL KDIST, INRADFLAG, NEWGRIDFLAG, DELTAM, ACCELFLAG
      LOGICAL BASEOUT, HIGHORDERRAD, BTEST
      LOGICAL RADIANCE_OUT, VISUAL_OUT
      LOGICAL FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
      REAL  SOLARFLUX, SOLARMU, SOLARAZ, SOLFLUX
      REAL  GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL  SOLACC, SOLCRIT, SPLITACC, SHACC
      REAL  DELX, DELY, DELXSFC, DELYSFC, MAXASYM
      REAL  MUOUT, PHIOUT
      REAL  MAX_TOTAL_MB, ADAPT_GRID_FACTOR
      REAL  NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
      CHARACTER SRCTYPE*1, UNITS*1
      CHARACTER GRIDTYPE*1, PROPTYPE*1, SFCTYPE*2
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
      CHARACTER*64 OUTFILENC



C         Common block for the input medium properties
      COMMON /SHDOMPROP/ NPX, NPY, NPZ, DELX, DELY, ZLEVELS,
     .        TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP,
     .        EXTDIRP, NZCKD, ZCKD, GASABS



C          Choose the option to compute higher order radiance SH terms
C            Change to SUMSHTERMS(5*MAXIG) for this option.
      HIGHORDERRAD = .FALSE.
      NSHOUT = 4
      IF (HIGHORDERRAD) NSHOUT=5
C          Choose type of initialization/output for k-distributions:
C            True means initiatize with and output only the base grid.
      BASEOUT = .FALSE.
 
      CALL USER_INPUT (PROPFILE, SFCFILE, CKDFILE, 
     .                 INSAVEFILE, OUTSAVEFILE,
     .                 NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG, 
     .                 KDIST, DELTAM, GRIDTYPE, 
     .                 SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .                 GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, 
     .                 ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,
     .                 MAXOUT,MAXPAR,NUMOUT,OUTTYPES,OUTPARMS,OUTFILES,
     .                 OUTFILENC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,
     .                 NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
c      CALL NAMELIST_INPUT (PROPFILE, SFCFILE, CKDFILE, 
c     .                 INSAVEFILE, OUTSAVEFILE,
c     .                 NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG, 
c     .                 KDIST, DELTAM, GRIDTYPE, 
c     .                 SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
c     .                 GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, 
c     .                 ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,
c     .                 MAXOUT,MAXPAR,NUMOUT,OUTTYPES,OUTPARMS,OUTFILES,
c     .                 OUTFILENC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,
c     .                 NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

C          Make ML and MM from NMU and NPHI
      NMU = MAX(2, 2*INT((NMU+1)/2) )
      NPHI = MAX(1,NPHI)
      ML = NMU-1
      MM = MAX(0,INT(NPHI/2)-1)


      NLEG = ML
      IF (DELTAM) NLEG = ML+1

C         Read in the properties of the medium
      CALL READ_PROPERTIES (PROPFILE, NPX, NPY, NPZ, 
     .         MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM, 
     .         PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,
     .         TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)

      
      IF (SFCFILE(1:2) .EQ. 'NO') THEN
        SFCTYPE = 'FL'
      ELSE
        CALL READ_SURFACE (SFCFILE, MAXSFCPTS, MAXSFCPARS, 
     .                     SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .                     NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
      ENDIF

C          If doing a k-distribution then get the band info from the CKD file
      IF (KDIST) THEN
        CALL READ_CKD (MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX, 
     .                 NG, DELG, NZCKD, ZCKD, KABS)
        SOLARFLUX = SOLFLUX*SOLARFLUX*ABS(SOLARMU)
      ELSE
        NG = 1
        DELG(1) = 1.0
        NZCKD = 0
        BASEOUT = .FALSE.
      ENDIF


C         Check input parameters
      CALL CHECK_INPUT_PARMETERS (NMU, NPHI, DELTAM, MAXASYM,
     .           GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU, 
     .           WAVENO, WAVELEN, GNDTEMP, GNDALBEDO, 
     .           SPLITACC, SHACC)

C         Check values of medium properties
      CALL CHECK_PROPERTY_INPUT (NPX, NPY, NPZ, NLEG,  
     .         DELX, DELY, ZLEVELS,
     .         TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)


      NX = MAX(1,NX)
      NY = MAX(1,NY)
      NZ = MAX(2,NZ)
C         If single plane or column then force IP mode
      IF (NX .EQ. 1)  IPFLAG=IOR(IPFLAG,1)
      IF (NY .EQ. 1)  IPFLAG=IOR(IPFLAG,2)
C         IP mode takes precedence over open boundary conditions
      IF (BTEST(IPFLAG,0))  BCFLAG=IAND(BCFLAG,2)
      IF (BTEST(IPFLAG,1))  BCFLAG=IAND(BCFLAG,1)
C         Set up base grid point actual size (NX1xNY1xNZ);  IP mode NX,NY=1
      NX1 = NX+1
      IF (BTEST(IPFLAG,0) .OR. BTEST(BCFLAG,0)) NX1 = NX
      NY1 = NY+1
      IF (BTEST(IPFLAG,1) .OR. BTEST(BCFLAG,1)) NY1 = NY
      NBPTS = NX1*NY1*NZ
      NBCELLS = (NZ-1)*(NX+IBITS(BCFLAG,0,1))*(NY+IBITS(BCFLAG,1,1))

C          Make ML and MM from NMU and NPHI
C          Compute NLM, NPHI0MAX, and NLEG.
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


C           Check to see that arrays sizes are not exceeded
      IF (NLM*NMU .GT. MAXNCMU) STOP 'MAXNCMU exceeded'
      IF (NMU*33*32 .GT. MAXNCPHI) STOP 'MAXNCPHI exceeded'
      IF (NBPTS .GT. MAXIG) STOP 'MAXIG exceeded'
      IF (NBCELLS .GT. MAXIC) STOP 'MAXIC exceeded'
      IF (NBPTS*NPHI0MAX .GT. MAXIDO) STOP 'MAXIDO exceeded'
      IF (NBPTS*(NCS+2) .GT. MAXIV) STOP 'MAXIV may be exceeded'
      IF (PROPTYPE .EQ. 'E' .OR. PROPTYPE .EQ. 'T') THEN
        IF (NUMPHASE*(NLEG+1) .GT. MAXIGL) STOP 'MAXIGL exceeded'
      ELSE
        IF (NBPTS*(NLEG+1) .GT. MAXIGL) STOP 'MAXIGL exceeded'
      ENDIF
      IF (NX .GT. MAXNXY .OR. NY .GT. MAXNXY) STOP 'MAXNXY exceeded'
      IF (NZ .GT. MAXNZ) STOP 'MAXNZ exceeded'
      IF (NZ .LE. 1) STOP 'NZ < 2'
      IF (NUMOUT .GT. MAXOUT) STOP 'MAXOUT exceeded'


 
C         Figure out what type of output we will be generating
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
          IF (OUTPARMS(2,I) .EQ. 0.0) THEN
            NXOUT = 1
          ELSE
            NXOUT=MAX(1,NINT((NPX*DELX+2*OUTPARMS(4,I))/OUTPARMS(2,I)))
          ENDIF
          IF (OUTPARMS(3,I) .EQ. 0.0) THEN
            NYOUT = 1
          ELSE
            NYOUT=MAX(1,NINT((NPY*DELY+2*OUTPARMS(5,I))/OUTPARMS(3,I)))
          ENDIF
          NRAD = NRAD + NXOUT*NYOUT*NINT(OUTPARMS(6,I))
          IF (NRAD .GT. MAXRADOUT) 
     .      STOP 'Main: Radiance output array exceeded'
        ELSE IF (OUTTYPES(I) .EQ. 'V') THEN
          VISUAL_OUT = .TRUE.
          IF (NINT(OUTPARMS(1,I)) .EQ. 1) THEN
            NVIS = NVIS + NINT(OUTPARMS(10,I))*NINT(OUTPARMS(11,I))
          ELSE
            NSCAN = 1 + SQRT((OUTPARMS(4,I)-OUTPARMS(7,I))**2
     .               +(OUTPARMS(5,I)-OUTPARMS(8,I))**2 
     .               +(OUTPARMS(6,I)-OUTPARMS(9,I))**2) /OUTPARMS(10,I)
            NPIX = 1 + ABS(OUTPARMS(12,I)-OUTPARMS(11,I))/OUTPARMS(13,I)
            NVIS = NVIS + NSCAN*NPIX
          ENDIF
          IF (NVIS .GT. MAXRADOUT) 
     .      STOP 'Main: Visualization output array exceeded'
        ENDIF
        IF (OUTTYPES(I) .EQ. 'J')  SOURCE_OUT = .TRUE.
        IF (OUTTYPES(I) .EQ. 'F')  FLUX_OUT = .TRUE.
        IF (OUTTYPES(I) .EQ. 'H')  FLUXDIV_OUT = .TRUE.
        IF (OUTTYPES(I) .EQ. 'S')  SH_OUT = .TRUE.
      ENDDO 
 

C         Read in previous run from binary file if desired
      CALL RESTORE_STATE (INSAVEFILE, NX, NY, NZ, ML, MM, NCS, NLM,
     .       INRADFLAG, NEWGRIDFLAG,  XGRID, YGRID, ZGRID,
     .       NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, 
     .       CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

      IF (NEWGRIDFLAG) THEN
C           Make the internal base grid lines
        CALL NEW_GRIDS (BCFLAG, GRIDTYPE, NPX, NPY, NPZ, NX, NY, NZ,
     .                  0.0, 0.0, DELX, DELY, ZLEVELS, 
     .                  XGRID, YGRID, ZGRID)
      ENDIF



      TOTITER = 0
      OLDNPTS = 0

C         Loop over the k-distribution g's from high to low absorption
      DO IG = NG, 1, -1
C           If doing k-distribution then get the gas absorption profile.
        IF (KDIST) THEN
          DO I = 1, NZCKD
            GASABS(I) = KABS(I+NZCKD*(IG-1))
          ENDDO
        ENDIF

C           If this is the first time through or we always want a base grid,
C           then make the base grid point positions and grid cell structure.
        IF (BASEOUT .OR. (IG .EQ. NG .AND. NEWGRIDFLAG)) THEN
          CALL INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, 
     .              NX, NY, NZ, NX1, NY1, NPTS, NCELLS, 
     .              XGRID, YGRID, ZGRID, GRIDPOS,
     .              GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
          NBCELLS = NCELLS
        ENDIF

C           Find the monochromatic SHDOM solution: adaptive grid output
C             in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        CALL SOLVE_RTE (NX, NY, NZ, NX1, NY1, NANG, 
     .               ML, MM, NCS, NLM, NMU, NPHI, NLEG, NUMPHASE,
     .               NPHI0, MU, PHI, WTDO,
     .               MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG, 
     .               BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD,
     .               SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .               SFCTYPE, GNDTEMP, GNDALBEDO, 
     .               NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .               NSFCPAR, SFCPARMS, SFCGRIDPARMS, 
     .               UNITS, WAVENO, WAVELEN, 
     .               ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,
     .               SPLITACC, SHACC,  XGRID,YGRID,ZGRID,
     .               TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,
     .               MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .               CMU1, CMU2, CPHI1, CPHI2,  NPTS, GRIDPOS, 
     .               NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,
     .               SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        INRADFLAG = .FALSE.
        TOTITER = TOTITER + ITER

C         Compute the output and put in arrays
        IF (RADIANCE_OUT) THEN
          IRAD = 0
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'R') THEN
              CALL COMPUTE_RADIANCE (NX, NY, NZ, NPTS, NCELLS,
     .               ML, MM, NCS, NLEG, NUMPHASE,
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .               BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ, 
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .               MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .               GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .               XGRID, YGRID, ZGRID, GRIDPOS, 
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .               SHPTR, SOURCE, WORK(1), WORK(1+NPTS),
     .               OUTPARMS(1,I),  IRAD, WORK2)
            ENDIF
          ENDDO
          CALL SUM_OUTPUT (NG-IG, DELG(IG), NRAD, WORK2, SUMRADOUT)
        ENDIF

        IF (VISUAL_OUT) THEN
          IVIS = 0
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'V') THEN
              CALL VISUALIZE_RADIANCE (NX, NY, NZ, NPTS, NCELLS,
     .               ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .               BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .               MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .               GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .               XGRID, YGRID, ZGRID, GRIDPOS,
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .               SHPTR, SOURCE, OUTPARMS(1,I),  IVIS, WORK2)
            ENDIF
          ENDDO
          CALL SUM_OUTPUT (NG-IG, DELG(IG), NVIS, WORK2, SUMVISOUT)
        ENDIF

        IF (SOURCE_OUT) THEN
C             Compute source function throughout grid for this one angle
          DO I = 1, NUMOUT
            IF (OUTTYPES(I) .EQ. 'J') THEN
              MUOUT = OUTPARMS(2,I)
              PHIOUT = OUTPARMS(3,I)*ACOS(-1.0)/180.0
              CALL COMPUTE_ONE_SOURCE (ML,MM,NCS, NLEG, NUMPHASE, 
     .            NPTS, DELTAM, MUOUT, PHIOUT, 
     .            SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE, 
     .            DIRFLUX, SHPTR, SOURCE,  SOURCE1OUT)
            ENDIF
          ENDDO
        ENDIF



        IF (BASEOUT) THEN
C           If outputting only the base grid then throw out the other points
          CALL INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, 
     .                  NX, NY, NZ, NX1, NY1, 
     .                  NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,
     .                  GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
        ELSE
C           For the output from the grid point arrays (Flux, Net flux, SH)
C             we may need to interpolate the running k-distribution sum
C             to the newly created grid points.
          CALL INTERP_OUTPUT (OLDNPTS,NPTS,FLUX_OUT,FLUXDIV_OUT,SH_OUT,
     .             SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,
     .             NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ENDIF

        IF (FLUX_OUT) THEN
          CALL SUM_OUTPUT (NG-IG,DELG(IG), 2*NPTS, FLUXES, SUMFLUXES)
          CALL SUM_OUTPUT (NG-IG,DELG(IG), NPTS, DIRFLUX, SUMDIRFLUX)
        ENDIF
        IF (FLUXDIV_OUT) THEN
C               Compute the net flux divergence at every point
          CALL COMPUTE_NETFLUXDIV (NPTS, RSHPTR, SRCTYPE, SOLARMU,
     .           EXTINCT,ALBEDO, PLANCK,DIRFLUX, RADIANCE, WORK)
         CALL SUM_OUTPUT (NG-IG,DELG(IG),NPTS,WORK,SUMFLUXDIV)
        ENDIF

        IF (SH_OUT) THEN
          CALL COMPUTE_SH (NSHOUT, NPTS, SRCTYPE, SOLARMU, SOLARAZ,
     .                   DIRFLUX, RSHPTR, ML, MM, NCS, RADIANCE, WORK)
          CALL SUM_OUTPUT (NG-IG,DELG(IG),NSHOUT*NPTS,WORK,SUMSHTERMS)
        ENDIF

        OLDNPTS = NPTS
      ENDDO
 

      CALL SAVE_STATE (OUTSAVEFILE, NX, NY, NZ, ML, MM, NCS, NLM,
     .       WORK, XGRID, YGRID, ZGRID, 
     .       NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, 
     .       CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

      NSH = SHPTR(NPTS+1)
      IRAD = 0
      IVIS = 0
      DO I = 1, NUMOUT
        IF (OUTTYPES(I) .EQ. 'V') THEN
          CALL OUTPUT_IMAGE (NX, NY, NZ, NPTS, NCELLS,
     .             NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,
     .             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,
     .             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,
     .             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,
     .             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,
     .             IVIS, SUMVISOUT, OUTPARMS(1,I), OUTFILES(I))
        ELSE
          CALL OUTPUT_RESULTS (NX, NY, NZ, NBPTS, NPTS, NCELLS, 
     .         NSH, ML, MM, NLM, NLEG, NUMPHASE, NMU, NPHI, NANG, NG,
     .         PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,
     .         BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, 
     .         SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .         SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, 
     .         SPLITACC, SHACC, SOLACC, MAXITER,TOTITER, 
     .         XGRID, YGRID, ZGRID, GRIDPOS,
     .         TREEPTR, GRIDPTR, CELLFLAGS,
     .         EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, 
     .         SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,
     .         IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,
     .         OUTTYPES(I), OUTPARMS(1,I), OUTFILES(I))
        ENDIF
      ENDDO


      END
 







 
      SUBROUTINE TRILIN_INTERP_PROP (X, Y, Z, INIT, NLEG,
     .               TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
C      Trilinearly interpolates the quantities on the input property
C     grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
C     ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the 
C     volume coefficients.  Also adds in the separate gaseous absorption.
C     Divides the phase function Legendre coefficients LEGEN by 2*l+1.
C     The phase function pointer IPHASE (for tabulated phase functions) 
C     is that of the maximum weighted scattering property grid point.
C     If INIT=.TRUE. then transfers the tabulated phase functions.
      INTEGER NLEG
      INTEGER*2 IPHASE
      LOGICAL INIT
      REAL    X, Y, Z,  TEMP, EXTINCT, ALBEDO, LEGEN(0:NLEG,*)
      INTEGER IX, IXP, IY, IYP, IZ, L, IL, IM, IU
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I
      DOUBLE PRECISION U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
      DOUBLE PRECISION SCAT1,SCAT2,SCAT3,SCAT4,SCAT5,SCAT6,SCAT7,SCAT8
      DOUBLE PRECISION SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN

      INTEGER NPX, NPY, NPZ, NZCKD, NUMPHASE
      INTEGER    MAXPG, MAXLEG, MAXPGL, MAXNZ
      PARAMETER (MAXPG=200000, MAXLEG=3000, MAXPGL=MAXLEG*300)
      PARAMETER (MAXNZ=501)
      INTEGER*2 IPHASEP(MAXPG)
      REAL  DELX, DELY, ZLEVELS(MAXNZ)
      REAL  TEMPP(MAXPG),EXTINCTP(MAXPG),ALBEDOP(MAXPG),LEGENP(MAXPGL)
      REAL  EXTDIRP(MAXPG), ZCKD(MAXNZ), GASABS(MAXNZ)
      COMMON /SHDOMPROP/ NPX, NPY, NPZ, DELX, DELY, ZLEVELS,
     .        TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP,
     .        EXTDIRP, NZCKD, ZCKD, GASABS
      SAVE EXTMIN, SCATMIN

      IF (INIT) THEN
C         If there are tabulated phase functions, then transfer them
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

C         Find the grid location and compute the interpolation factors
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
      IX = INT(X/DELX) + 1
      IF (ABS(X-NPX*DELX) .LT. 0.01*DELX) IX = NPX
      IF (IX .LT. 1 .OR. IX .GT. NPX)  STOP 'TRILIN: Beyond X domain'
      IXP = MOD(IX,NPX) + 1
      U = DBLE(X-DELX*(IX-1))/DELX
      U = MAX( MIN( U, 1.0D0), 0.0D0)
      IF (U .LT. 1.0D-5) U = 0.0D0
      IF (U .GT. 1.0D0-1.0D-5) U = 1.0D0
      IY = INT(Y/DELY) + 1
      IF (ABS(Y-NPY*DELY) .LT. 0.01*DELY) IY = NPY
      IF (IY .LT. 1 .OR. IY .GT. NPY)  STOP 'TRILIN: Beyond Y domain'
      IYP = MOD(IY,NPY) + 1
      V = DBLE(Y-DELY*(IY-1))/DELY
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

C         Trilinearly interpolate the temperature, extinction, scattering
      TEMP = F1*TEMPP(I1) + F2*TEMPP(I2) + F3*TEMPP(I3) + F4*TEMPP(I4)
     .     + F5*TEMPP(I5) + F6*TEMPP(I6) + F7*TEMPP(I7) + F8*TEMPP(I8)
      EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3)
     .        + F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6)
     .        + F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
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
C         For tabulated phase functions pick the one we are on top of
C         or the one with the most scattering weight.
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
          LEGEN(L,1) = SCAT1*LEGENP(L+NLEG*(I1-1))
     .               + SCAT2*LEGENP(L+NLEG*(I2-1))
     .               + SCAT3*LEGENP(L+NLEG*(I3-1))
     .               + SCAT4*LEGENP(L+NLEG*(I4-1))
     .               + SCAT5*LEGENP(L+NLEG*(I5-1))
     .               + SCAT6*LEGENP(L+NLEG*(I6-1))
     .               + SCAT7*LEGENP(L+NLEG*(I7-1))
     .               + SCAT8*LEGENP(L+NLEG*(I8-1))
          IF (SCATTER .GT. SCATMIN) THEN
            LEGEN(L,1) = LEGEN(L,1)/SCATTER
          ELSE
            LEGEN(L,1) = LEGEN(L,1)/SCATMIN
          ENDIF
          LEGEN(L,1) = LEGEN(L,1)/(2*L+1)
        ENDDO
      ENDIF

C         Add in the gaseous absorption to extinction and albedo
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
 
 


 
      SUBROUTINE DIRECT_BEAM_PROP (X, Y, Z, INIT, BCFLAG, IPFLAG,
     .                DELTAM, ML, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, 
     .                DIRFLUX)
C       Computes the direct beam flux at point (X,Y,Z) by integrating
C     the extinction through the property grid.  If called with 
C     INIT=.TRUE. then the property grid extinction array, solar 
C     direction terms, and lowest uniform level are computed and saved.  
C     The DELTAM flag is used to determine whether to delta-M scale the 
C     extinction for the direct beam computation.  The extinction includes 
C     the gaseous absorption.  If the IPFLAG has bit 2 set then the direct
C     beam tracing is done in 3D, otherwise the lower two bits determine
C     the type of tracing: 0 for 3D tracing, 1 for XZ only tracing, 
C     2 for YZ only tracing, and 3 for Z only tracing.  If BCFLAG>0 then
C     have open boundary conditions in X and/or Y. In this case when the
C     ray tracing back to the sun reaches the boundary then independent
C     pixel mode is entered so that only Z grid intersections occur.

      INTEGER BCFLAG, IPFLAG, ML, NLEG
      LOGICAL INIT, DELTAM
      REAL    X, Y, Z, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX

      INTEGER IX, IY, IZ, JZ, IL, IM, IU
      INTEGER I, J, K, IP, JP, I1, I2, I3, I4
      INTEGER IPDIRECT, DI, DJ, DK
      LOGICAL CONSTX, CONSTY, BTEST
      DOUBLE PRECISION EXTINCT, ALBEDO, F
      DOUBLE PRECISION SUNMU, SUNAZ, PATH
      DOUBLE PRECISION EXTBEAMCUT, UNIFORMZLEV
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION XOFFS, YOFFS, DELXD, DELYD
      DOUBLE PRECISION XE, YE, ZE, XP, YP, ZP
      DOUBLE PRECISION X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
      DOUBLE PRECISION U0, V0, W0, U1, V1, W1, AX, AY, AZ
      DOUBLE PRECISION U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
      DOUBLE PRECISION E1,E2,E3,E4,E5,E6,E7,E8,  A, B, C, D
      DOUBLE PRECISION B1,B2,B3,B4,B5,B6,B7,B8,C1,C2,C3,C4,C5,C6,C7,C8 
      DOUBLE PRECISION UV,UMV,UVM,UMVM,UW,UMW,UWM,UMWM,VW,VMW,VWM,VMWM
      DOUBLE PRECISION VWU,VWUM,UWV,UWVM,UVW,UVWM

      INTEGER NPX, NPY, NPZ, NZCKD, NUMPHASE
      INTEGER    MAXPG, MAXLEG, MAXPGL, MAXNZ
      PARAMETER (MAXPG=200000, MAXLEG=3000, MAXPGL=MAXLEG*300)
      PARAMETER (MAXNZ=501)
      INTEGER*2 IPHASEP(MAXPG)
      REAL  DELX, DELY, ZLEVELS(MAXNZ)
      REAL  TEMPP(MAXPG),EXTINCTP(MAXPG),ALBEDOP(MAXPG),LEGENP(MAXPGL)
      REAL  EXTDIRP(MAXPG), ZCKD(MAXNZ), GASABS(MAXNZ)
      COMMON /SHDOMPROP/ NPX, NPY, NPZ, DELX, DELY, ZLEVELS,
     .        TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP,
     .        EXTDIRP, NZCKD, ZCKD, GASABS

      REAL  GASEXT(MAXNZ), EXTMIN(MAXNZ), EXTMAX(MAXNZ)

      SAVE  GASEXT, EXTMIN, EXTMAX
      SAVE  CX, CY, CZ, CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT
      SAVE  DELXD, DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV


      IF (INIT) THEN
        IF (NPZ .GT. MAXNZ)  STOP 'DIRECT_BEAM_PROP: MAXNZ exceeded'
C           Get the gaseous extinction at the property grid levels
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
C           First make the property grid extinction field, on which all
C           the direct beam paths will be computed.
        IP = 0
        DO IX = 1, NPX
          DO IY = 1, NPY
            DO IZ = 1, NPZ
              IP = IP + 1
              EXTINCT = EXTINCTP(IP)
              ALBEDO = ALBEDOP(IP)
C                 Add in the gaseous absorption to extinction
              IF (GASEXT(IZ) .GT. 0.0) THEN
                ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                EXTINCT = EXTINCT + GASEXT(IZ)
              ENDIF
C                 Do the Delta-M scaling if needed
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

C           Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        IF (BTEST(IPFLAG,2)) IPDIRECT = 0

C           Make the ray direction (opposite to the solar beam)
C           SOLAR?? is direction of beam, SUN?? is direction to sun (opposite)
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
        EPSS = 1.0E-4*(ZLEVELS(NPZ)-ZLEVELS(1))/NPZ
        IF (.NOT. BTEST(IPDIRECT,0))  EPSS = MAX(EPSS,1.0D-4*DELX)
        IF (.NOT. BTEST(IPDIRECT,1))  EPSS = MAX(EPSS,1.0D-4*DELY)
        DELXD = DBLE(DELX)
        DELYD = DBLE(DELY)
        XDOMAIN = DELXD*NPX
        YDOMAIN = DELYD*NPY

C           Find the Z level above which the medium is plane-parallel
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
        RETURN
C           Done with initialization
      ENDIF


C         Here for computing the direct beam path for one starting point.
C         Find the grid location 
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
      IF (ABS(X-XDOMAIN) .LT. 0.01*DELXD) I = NPX
      IF (I .LT. 1 .OR. I .GT. NPX)
     .  STOP 'DIRECT_BEAM_PROP: Beyond X domain'
      J = INT(Y/DELYD) + 1
      IF (ABS(Y-YDOMAIN) .LT. 0.01*DELYD) J = NPY
      IF (J .LT. 1 .OR. J .GT. NPY)
     .  STOP 'DIRECT_BEAM_PROP: Beyond Y domain'
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
      IF (BTEST(BCFLAG,0) .AND. (ABS(X) .LT. 0.01*DELXD 
     .    .OR. ABS(X-(NPX-1)*DELXD) .LT. 0.01*DELXD))  CONSTX = .TRUE.
      IF (BTEST(BCFLAG,1) .AND. (ABS(Y) .LT. 0.01*DELYD 
     .    .OR. ABS(Y-(NPY-1)*DELYD) .LT. 0.01*DELYD))  CONSTY = .TRUE.
 

C           Grid cell loop begin
      PATH = 0.0
      DO WHILE (ABS(ZE-ZLEVELS(NPZ)) .GT. EPSZ)
        IP = I + 1
        IF (I .EQ. NPX) THEN
          IF (BTEST(BCFLAG,0)) THEN
            IP=NPX
          ELSE
            IP = 1
          ENDIF   
        ENDIF     
        JP = J + 1
        IF (J .EQ. NPY) THEN
          IF (BTEST(BCFLAG,1)) THEN
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
        IF (I .LT. 1 .OR. I .GT. NPX .OR.
     .      J .LT. 1 .OR. J .GT. NPY .OR.
     .      K .LT. 1 .OR. K .GE. NPZ) THEN
          WRITE (*,'(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
          WRITE (*,'(1(2X,3F9.5))') X, Y, Z
          STOP
        ENDIF
C           Get the eight corner extinction values
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

C           Compute the distance to the next grid plane in  X, Y, and Z
C             If in horizontal uniform region or doing IP then fix X and/or Y.
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
        
C           The shortest distance is the plane we stop at:
C             get the exitting location and increment the cell
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
C             If have reached a horizontal boundary then either wrap around
C               (periodic) or go into IP mode (open boundaries).
          IF (I .EQ. 0) THEN
            IF (BTEST(BCFLAG,0)) THEN
              I = 1
              CONSTX = .TRUE.
            ELSE
              I = NPX
              XOFFS = XDOMAIN
            ENDIF
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
            ELSE
              J = NPY
              YOFFS = YDOMAIN
            ENDIF
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
          WRITE (*,*) 'DIRECT_BEAM_PROP: SO<0', X, Y, Z,
     .        XE, YE, ZE, XP, YP, ZP, SOX, SOY, SOZ
          STOP
        ENDIF
        SO = MAX(SO,0.0D0)
C           Make the starting and ending interpolation factors
        AX = 1.0D0/(X1-X0)
        AY = 1.0D0/(Y1-Y0)
        AZ = 1.0D0/(Z1-Z0)
        U0 = (XE-X0)*AX
        V0 = (YE-Y0)*AY
        W0 = (ZE-Z0)*AZ
        U1 = (XP-X0)*AX
        V1 = (YP-Y0)*AY
        W1 = (ZP-Z0)*AZ
C           Compute the cubic polynomial extinction coefficients
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
        A =  (E1*U0M + E2*U0)*VMWM 
     .     + (E3*U0M + E4*U0)*V WM
     .     + (E5*U0M + E6*U0)*VMW
     .     + (E7*U0M + E8*U0)*V W
        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 =  DU*VMWM - DV*U WM - DW*U VM
        B3 = -DU*V WM + DV*UMWM - DW*UMV
        B4 =  DU*V WM + DV*U WM - DW*U V
        B5 = -DU*VMW  - DV*UMW  + DW*UMVM
        B6 =  DU*VMW  - DV*U W  + DW*U VM
        B7 = -DU*V W  + DV*UMW  + DW*UMV
        B8 =  DU*V W  + DV*U W  + DW*U V
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
C            Compute the path through the cell: integration of extinction
        PATH = PATH + SO*(A +0.5D0*B +0.3333333333333333D0*C +0.25D0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
      ENDDO
      DIRFLUX = SOLARFLUX*EXP(-PATH)
      RETURN
      END


 
