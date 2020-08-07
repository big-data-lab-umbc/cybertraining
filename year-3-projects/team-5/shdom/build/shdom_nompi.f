 ! Dummy routines to allow parallel SHDOM to compile without MPI.



      SUBROUTINE START_MPI (MASTERPROC)
      LOGICAL MASTERPROC

      MASTERPROC = .TRUE.
      RETURN
      END



      SUBROUTINE MAP_SHDOM_MPI (BCFLAG, NPX,NPY, NX,NY,NZ, DELX,DELY, 
     .                          PROPFILE, XSTART, YSTART)
      INTEGER BCFLAG, NPX, NPY, NX, NY, NZ
      REAL    DELX, DELY
      CHARACTER PROPFILE*64
      REAL    XSTART, YSTART

      XSTART = 0.0
      YSTART = 0.0
      RETURN
      END




      SUBROUTINE BROADCAST_USER_INPUT (PROPFILE, SFCFILE, CKDFILE,
     .                 INSAVEFILE, OUTSAVEFILE,
     .                 NX, NY, NZ, NMU, NPHI, BCFLAG,
     .                 IPFLAG, KDIST, DELTAM, GRIDTYPE,
     .                 SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .                 GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,
     .                 ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,
     .                 MAXOUT,MAXPAR,NUMOUT,OUTTYPES,OUTPARMS,OUTFILES)
      INTEGER NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
      INTEGER MAXOUT, MAXPAR, MAXITER, NUMOUT
      LOGICAL KDIST, DELTAM, ACCELFLAG
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL    SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR,MAXOUT)
      CHARACTER SRCTYPE*1, GRIDTYPE*1, UNITS*1, OUTTYPES(*)*1
      CHARACTER PROPFILE*64, SFCFILE*64, CKDFILE*64
      CHARACTER OUTFILES(*)*64, INSAVEFILE*64, OUTSAVEFILE*64

      RETURN
      END 



      SUBROUTINE BROADCAST_PROPERTY_SIZE (NPX, NPY, NPZ, DELX, DELY,
     .                                    NUMPHASE, MAXLEG, MAXPGL)
      INTEGER NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
      REAL    DELX, DELY

      RETURN
      END


      SUBROUTINE SCATTER_PROPERTIES (NPXT,NPYT, NPX,NPY,NPZ, NLEG, 
     .                    NUMPHASE, ZLEVELS, MAXASYM, TEMPPT, 
     .                    EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,
     .                    TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
      INTEGER   NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE, NLEG
      REAL      ZLEVELS(NPZ), MAXASYM
      REAL      TEMPPT(NPZ,NPYT,NPXT), EXTINCTPT(NPZ,NPYT,NPXT)
      REAL      ALBEDOPT(NPZ,NPYT,NPXT), LEGENPT(*)
      INTEGER*2 IPHASEPT(NPZ,NPYT,NPXT)
      REAL     TEMPP(NPZ,NPY,NPX), EXTINCTP(NPZ,NPY,NPX) 
      REAL     ALBEDOP(NPZ,NPY,NPX), LEGENP(*)
      INTEGER*2 IPHASEP(NPZ,NPY,NPX)

      RETURN
      END


      SUBROUTINE READ_BROADCAST_MEM_PARAMS (MAX_TOTAL_MB, 
     .                ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, 
     .                CELL_TO_POINT_RATIO, RUNNAME)
      REAL MAX_TOTAL_MB, ADAPT_GRID_FACTOR
      REAL NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
      CHARACTER(*) RUNNAME

      RETURN
      END


      SUBROUTINE BROADCAST_SURFACE_SIZE (MAXSFCPTS, MAXSFCPARS)
      INTEGER MAXSFCPTS, MAXSFCPARS

      RETURN
      END


      SUBROUTINE BROADCAST_SURFACE_PARMS (SFCTYPE, NXSFC,NYSFC,NSFCPAR,
     .                  DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
      CHARACTER SFCTYPE*2
      INTEGER   NXSFC, NYSFC, NSFCPAR
      REAL      DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO

      RETURN
      END




      SUBROUTINE BROADCAST_KDIST_SIZE (NG, NZCKD)
      INTEGER NG, NZCKD

      RETURN
      END


      SUBROUTINE BROADCAST_KDIST_PARMS (SOLFLUX, NG, DELG, 
     .                                  NZCKD, ZCKD, KABS)
      INTEGER NG, NZCKD
      REAL    SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG)

      RETURN
      END




      SUBROUTINE GATHER_OUTPUT (FLUX_OUT, FLUXDIV_OUT, SH_OUT,
     .                   IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,
     .                   SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,
     .                   NPX, NPY, DELX, DELY, XALLGRID, YALLGRID,
     .                   ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,
     .                   NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)
      LOGICAL FLUX_OUT, FLUXDIV_OUT, SH_OUT
      INTEGER IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPX, NPY
      REAL    DELX, DELY
      REAL    SUMFLUXES(2,NBPTS), SUMDIRFLUX(NBPTS)
      REAL    SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT,NBPTS)
      REAL    XALLGRID(NXT), YALLGRID(NYT)
      REAL    ALLFLUXES(3,NZ,NYT,NXT), ALLFLUXDIV(NZ,NYT,NXT)
      REAL    ALLSHTERMS(NSHOUT,NZ,NYT,NXT)
      INTEGER NCELLS, NPTS, NSH
      INTEGER NCELLSTOT, NPTSTOT, NSHTOT

      RETURN
      END


      REAL FUNCTION SUM_CPU_TIME (cpuTimes)
        real cpuTimes
        SUM_CPU_TIME = cpuTimes
        RETURN
      END 


      SUBROUTINE TOTAL_ALBEDO_MAX (ALBMAX)
      REAL ALBMAX
      RETURN
      END


      SUBROUTINE UNIFY_SPLITTING (DOSPLIT, STARTSPLITACC)
      LOGICAL DOSPLIT
      REAL    STARTSPLITACC
      RETURN
      END

      SUBROUTINE TOTAL_SPLITCRIT_MAX (SPLITCRIT)
      REAL SPLITCRIT
      RETURN
      END


      SUBROUTINE MAKE_DIRECT_PAR (SPT, NPTS, BCFLAG, IPFLAG, DELTAM, 
     .                ML, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,
     .                NX, XGRID, NY, YGRID,  DIRFLUX)
      INTEGER SPT, NPTS, BCFLAG, IPFLAG, ML, NLEG, NX, NY
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GRIDPOS(3,NPTS), XGRID(NX+1), YGRID(NY+1)
      REAL    DIRFLUX(NPTS)

      RETURN
      END




      SUBROUTINE FIND_BOUNDARY_POINTS (BCFLAG, IPFLAG, NPTS, SWEEPORD,
     .               GRIDPTR, GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)
      INTEGER BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS,*), GRIDPTR(8,*)
      INTEGER NX, NY, NZ
      REAL GRIDPOS(3,NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)

      RETURN
      END



      SUBROUTINE CALC_BOUNDARY_RADIANCES (BCFLAG, IPFLAG, JOCT, IZ,
     .                           NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                           NA, NPTS, NCELLS, GRIDPTR,
     .                           NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,
     .                           MU, PHI, EXTINCT, SOURCE,
     .                           KANG, GRIDRAD)
      INTEGER BCFLAG, IPFLAG, JOCT, IZ
      INTEGER NX, NY, NZ, NA, NPTS, NCELLS, KANG
      REAL    XGRID(NX), YGRID(NY), ZGRID(NZ)
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(NCELLS)
      REAL    GRIDPOS(3,NPTS), MU, PHI
      REAL    EXTINCT(NPTS), SOURCE(NA,NPTS), GRIDRAD(NPTS)

      RETURN
      END



      SUBROUTINE COMPUTE_RADIANCE_PAR (NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NCS, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, XDOMAIN, YDOMAIN, IPFLAG,
     .             SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, SOURCE1, GRIDRAD,
     .             OUTPARMS,  NRAD, RADOUT)
      INTEGER NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NRAD
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), BCPTR(MAXNBC,2)
      INTEGER*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL    MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
      REAL    XDOMAIN, YDOMAIN, XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
      REAL    GRIDPOS(3,NPTS)
      REAL    SFCGRIDPARMS(*), BCRAD(*)
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), FLUXES(2,NPTS), SOURCE(*)
      REAL    SOURCE1(NPTS), GRIDRAD(NPTS), OUTPARMS(*), RADOUT(*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1

      RETURN
      END



      SUBROUTINE CALC_ACCEL_SOLCRIT (DOACCEL, DELJDOT, DELJOLD, DELJNEW, 
     .                               JNORM, ACCELPAR, SOLCRIT)
C     Calculates the acceleration parameter and solution criterion from
C     the delta source function vector dot  products.
      LOGICAL DOACCEL
      REAL    DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR, SOLCRIT
      REAL    R, THETA, A
      SAVE   A
      DATA   A/0.0/

C       Accelerate if desired, didn't last time, and things are converging. 
      IF (DOACCEL .AND. A .EQ. 0.0 .AND. DELJNEW .LT. DELJOLD) THEN
C       Compute the acceleration extrapolation factor and apply it.
        R = SQRT(DELJNEW/DELJOLD)
        THETA = ACOS(DELJDOT/SQRT(DELJOLD*DELJNEW))
        A = (1 - R*COS(THETA) + R**(1+0.5*3.14159/THETA))
     .         /(1 + R**2  - 2*R*COS(THETA))  - 1.0
        A = MIN(10.0,MAX(0.0,A))
C         WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
      ELSE
        A = 0.0
      ENDIF
      ACCELPAR = A

      IF (JNORM .GT. 0.0) THEN
        SOLCRIT = SQRT(DELJNEW/JNORM)
      ELSE IF (DELJNEW .EQ. 0.0) THEN
        SOLCRIT = 0.0
      ENDIF
      RETURN
      END




      SUBROUTINE END_SHDOM_MPI (NPTS, GRIDPOS, NPX,NPY, XSTART,YSTART,
     .                          DELX, DELY, NPXT, NPYT, PROPFILE)
      INTEGER NPTS, NPX, NPY, NPXT, NPYT
      REAL    GRIDPOS(3,NPTS), XSTART, YSTART, DELX, DELY
      CHARACTER PROPFILE*64

      RETURN
      END


      SUBROUTINE ABORT_SHDOM_MPI (ERRSTR)
      CHARACTER(*) ERRSTR

      WRITE (6,*) ERRSTR
      stop
      END

