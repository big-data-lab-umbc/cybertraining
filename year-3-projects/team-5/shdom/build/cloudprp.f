       PROGRAM CLOUDPRP
C       Makes an SHDOM tabulated phase function format property file
C     from an input file containing liquid water content and (optionally)
C     effective radius.  The extinction, single scattering albedo, and
C     phase functions are made by assuming a gamma distribution of
C     droplets and using Lorenz-Mie theory.  If the effective radius
C     is not specified then the droplet concentration is assumed fixed
C     at some specified value.  The index of refraction is obtained
C     by averaging over the specified wavelength range for water,
C     while the central wavelength is used in the Mie computations.
C     The input cloud LWC field is for a regularly spaced grid.  This
C     grid is used as the output property file grid, except that 
C     extra Z levels may be added.  A vertical profile of horizontally
C     uniform aerosol scattering may be specified.  The aerosol properties
C     (0.55 um extinction, effective radius, distribution width, and 
C     index of refraction) are input from a file as a function of height.  
C     A Mie calculation gives the aerosol optical properties for each height.
C     The aerosol and molecular Rayleigh scattering phase functions at the 
C     aerosol input levels are put at the beginning of the tabulated phase 
C     function output.
C
C     MODIFICATIONS: 
C       SEP 1998  Faster Mie table generation. Added aerosols.  KFE
C       JAN 1999  Choice of gamma or lognormal particle size 
C                   distributions  -  LHC 1/8/1999
C       JUL 1999  Changed aerosol input file.  KFE
C     
C     compile with f77 -O2 -o cloudprp cloudprp.f mieindsub.f

      IMPLICIT NONE
      INTEGER MAXNPTS, MAXNTAB, MAXLEG, MAXNZ, MAXLEGA, MAXNAERO
      PARAMETER (MAXNPTS=2000000, MAXNTAB=300, MAXLEG=2000, MAXNZ=512)
      PARAMETER (MAXLEGA=200, MAXNAERO=50)
      INTEGER NX, NY, NZC, NZO, NZT, I, MAXTABLEG, MAXOUTLEG
      INTEGER NRETAB, NLEG(MAXNTAB), ICLDLEV(MAXNZ)
      INTEGER NZA, NTABAERO, ITABAERO(MAXNZ), NLEGAERO(MAXNAERO)
      LOGICAL REFFLAG
      REAL    DROPCONC, ALPHA, WAVELEN1, WAVELEN2, WAVELENCEN
      REAL    SRETAB, ERETAB,  DELX, DELY, RAYLCOEF
      REAL    ZLEVELS(MAXNZ), ZCLD(MAXNZ), ZOTHER(MAXNZ)
      REAL    TEMP(MAXNZ), TEMPCLD(MAXNZ), TEMPOTHER(MAXNZ)
      REAL    EXTRAYL(MAXNZ) 
      REAL    LWC(MAXNPTS), REFF(MAXNPTS)
      REAL    LEGEN(0:MAXLEG,MAXNTAB)
      REAL    EXTINCT(MAXNTAB), ALBEDO(MAXNTAB)
      REAL    ZAERO(MAXNZ), EXT55AERO(MAXNZ)
      REAL    REFFAERO(MAXNZ), WIDTHAERO(MAXNZ)
      REAL    EXTAERO(MAXNZ), SSALBAERO(MAXNZ)
      REAL    LEGENAERO(0:MAXLEGA,MAXNAERO)
      COMPLEX AERO55INDEX(MAXNZ), AEROINDEX(MAXNZ)
      COMPLEX RINDEX
      CHARACTER  CLOUDFILE*72, MIEFILE*72, PROPFILE*72, AEROFILE*72
      CHARACTER  MIEFLAG*1, PARTTYPE*1, distflag*1, aerodist*1


      WRITE (*,*) 'Mie table file flag (Input, Output)'
      READ (*,'(A1)') MIEFLAG
        WRITE (*,*) MIEFLAG
      WRITE (*,*) 'Mie table file name'
      READ (*,'(A)') MIEFILE
        WRITE (*,*) MIEFILE

      IF (MIEFLAG .EQ. 'O') THEN
        WRITE (*,*) 'Water, ice droplets, or aerosols (W,I,A)'
        READ (*,'(A1)') PARTTYPE
          WRITE (*,*) PARTTYPE
        IF (PARTTYPE .EQ. 'A') THEN
          READ (*,*) RINDEX
            WRITE (*,*) RINDEX
        ENDIF
        WRITE (*,*) 'Wavelength range (micron)'
        READ (*,*) WAVELEN1, WAVELEN2
          WRITE (*,*) WAVELEN1, WAVELEN2
        WRITE (*,*) 'Cloud particle distribution type: ',
     &              'G = Gamma, L = Lognormal'
        READ (*,'(A1)') distflag
          WRITE (*,*) distflag
        WRITE (*,*) 'Size distribution shape parameter - alpha',
     &     'for lognormal, this is log standard deviation'  
        READ (*,*) ALPHA
          WRITE (*,*) ALPHA
        WRITE (*,*) 
     .       'Number, starting, and ending tabulated effective radius'
        READ (*,*) NRETAB, SRETAB, ERETAB
          WRITE (*,*) NRETAB, SRETAB, ERETAB
        WRITE (*,*) 'Maximum Legendre order in Mie table'
        READ (*,*) MAXTABLEG
          WRITE (*,*) MAXTABLEG


C           Get the Planck weighted central wavelength for Mie calculations
        CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)

C           Get the Planck weighted index of refraction
        IF (PARTTYPE .NE. 'A') THEN
          CALL GET_REFRACT_INDEX (PARTTYPE, WAVELEN1, WAVELEN2, RINDEX)
          WRITE (*,'(A,F6.3,A,2E12.4)') ' Wavelength: ', WAVELENCEN,
     .       '   Average index: ', RINDEX
        ENDIF

C           Make the table of scattering properties as a function of
C             effective radius for LWC=1 g/m^3 for droplet distributions
        CALL MAKE_MIE_TABLE (NRETAB, SRETAB, ERETAB, 
     .                       WAVELENCEN, PARTTYPE, RINDEX, ALPHA, 
     .                  EXTINCT, ALBEDO, MAXLEG, NLEG, LEGEN, distflag)
        CALL WRITE_MIE_TABLE (MIEFILE, NRETAB, SRETAB, ERETAB, 
     .                        WAVELENCEN, RINDEX, ALPHA, 
     .                        EXTINCT, ALBEDO, MAXLEG, MAXTABLEG,
     .                        NLEG, LEGEN, distflag)

      ELSE


        WRITE (*,*) 'Output property file name'
        READ (*,'(A)') PROPFILE
          WRITE (*,*) PROPFILE
        WRITE (*,*) 'Input cloud file name (or NONE)'
        READ (*,'(A)') CLOUDFILE
          WRITE (*,*) CLOUDFILE
        WRITE (*,*) 'Droplet concentration (/cm^3) '
        READ (*,*) DROPCONC
          WRITE (*,*) DROPCONC
        WRITE (*,*) 'Molecular scattering coefficient (K/(km mb))'
        READ (*,*) RAYLCOEF
          WRITE (*,*) RAYLCOEF
        WRITE (*,*) 'Input aerosol property file (or NONE)'
        READ (*,'(A)') AEROFILE
          WRITE (*,*) AEROFILE
        WRITE (*,*) 'Aerosol distribution type: ',
     &              'G = Gamma, L = Lognormal'
        READ (*,'(a)') aerodist
          WRITE (*,*) aerodist
        WRITE (*,*) 'Maximum Legendre order in output'
        READ (*,*) MAXOUTLEG
          WRITE (*,*) MAXOUTLEG
        WRITE (*,*) 'Number of extra Z levels'
        READ (*,*) NZO
          WRITE (*,*) NZO
        DO I = 1, NZO
          WRITE (*,*) 'Height (km) and temperature (K) :', I
          READ (*,*) ZOTHER(I), TEMPOTHER(I)
            WRITE (*,*) ZOTHER(I), TEMPOTHER(I)
        ENDDO


C           Read in the scattering table
        CALL READ_MIE_TABLE (MIEFILE, MAXNTAB, NRETAB, SRETAB, ERETAB, 
     .                       WAVELENCEN, RINDEX, ALPHA, distflag,
     .                       EXTINCT, ALBEDO, MAXLEG, NLEG, LEGEN)

C           Read in the cloud LWC file
        CALL READ_CLOUD_FILE (CLOUDFILE, MAXNPTS, MAXNZ, REFFLAG, 
     .           NX, NY, NZC, DELX, DELY, ZCLD, TEMPCLD, LWC, REFF)

C           If the file had no effective radius then make them from LWC
        IF (.NOT. REFFLAG .AND. DROPCONC .GT. 0.0) THEN
          CALL MAKE_REFF (NX*NY*NZC, DROPCONC, ALPHA, LWC, REFF, 
     &                    distflag)
        ENDIF

C           Combine the cloud and extra levels 
        CALL ORGANIZE_LEVELS (NZC,ZCLD,TEMPCLD, NZO,ZOTHER,TEMPOTHER,
     .                        NZT, ZLEVELS, TEMP, ICLDLEV)

C           Compute the molecular Rayleigh scattering from density
        CALL RAYLEIGH_EXTINCT (NZT, ZLEVELS, TEMP, RAYLCOEF, EXTRAYL)

        IF (AEROFILE .NE. 'NONE') THEN
C             Read in the aerosol property file if there is one
          CALL READ_AERO_FILE (AEROFILE, MAXNZ, NZA, ZAERO,
     .                         EXT55AERO, REFFAERO, WIDTHAERO, 
     .                         AERO55INDEX, AEROINDEX)
        ELSE
          NZA = 0
        ENDIF

C             Do Mie calculations for the aerosols, make the phase function
C               table including the Rayleigh phase function, and 
C               interpolate properties to the output vertical grid
        CALL CALC_MIE_AEROSOLS (WAVELENCEN, NZA, ZAERO, 
     .                EXT55AERO, REFFAERO, WIDTHAERO, 
     .                AERO55INDEX, AEROINDEX,
     .                NZT, ZLEVELS, EXTRAYL, 
     .                EXTAERO, SSALBAERO, ITABAERO,
     .                NTABAERO, MAXLEGA, NLEGAERO, LEGENAERO, aerodist)

C           Interpolate in the scattering property table and output
C            the atmosphere property file
        CALL OUTPUT_PROPERTIES (PROPFILE, NX,NY,NZT, 
     .           DELX,DELY,ZLEVELS, TEMP, EXTAERO, SSALBAERO, 
     .           ITABAERO, NTABAERO, NLEGAERO, MAXLEGA, LEGENAERO,
     .           NZC, ICLDLEV, LWC, REFF,
     .           NRETAB, SRETAB, ERETAB, EXTINCT, ALBEDO, 
     .           MAXLEG, MAXOUTLEG, NLEG, LEGEN)

      ENDIF

      END





      SUBROUTINE READ_CLOUD_FILE (CLOUDFILE, MAXNPTS, MAXNZ, REFFLAG,
     .                            NX, NY, NZC, DELX, DELY, 
     .                            ZCLD, TEMPCLD, LWC, REFF)
C       Reads the liquid water content file.  If it has effective radius
C     (2 in the first place in file) then the effective radius is also
C     read in and the REFFLAG is returned true.  The header containing 
C     the array size (NX,NY,NZC) and heights (ZCLD) and temperature
C     profile (TEMPCLD) are also read in.
      IMPLICIT NONE
      INTEGER MAXNPTS, MAXNZ, NX, NY, NZC
      LOGICAL REFFLAG
      REAL    DELX, DELY, ZCLD(*), TEMPCLD(*), LWC(*), REFF(*)
      CHARACTER CLOUDFILE*(*)
      INTEGER NP, IX, IY, IZ, K, NPTS
      REAL    LWCT, REFFT

C         Read in header of LWC file
      OPEN (UNIT=1, FILE=CLOUDFILE, STATUS='OLD')
      READ (1,*) NP
      REFFLAG = (NP .EQ. 2)
      READ (1,*) NX, NY, NZC
      NPTS = NX*NY*NZC
      IF (NPTS .GT. MAXNPTS)  STOP 'READ_CLOUD_FILE: MAXNPTS exceeded'
      IF (NZC .GT. MAXNZ)  STOP 'READ_CLOUD_FILE: MAXNZ exceeded'
      READ (1,*) DELX, DELY
      READ (1,*) (ZCLD(IZ), IZ=1, NZC)
      READ (1,*) (TEMPCLD(IZ), IZ=1, NZC)
C         Initialize the output arrays (in case there are missing pixels)
      DO K = 1, NPTS
        LWC(K) = 0.0
        REFF(K) = 0.0
      ENDDO
C         Read the data in
      REFFT = 10.0
      DO WHILE (.TRUE.)
        IF (REFFLAG) THEN
          READ (1,*,END=190) IX, IY, IZ, LWCT, REFFT
        ELSE
          READ (1,*,END=190) IX, IY, IZ, LWCT
        ENDIF
        IF (IX .GE. 1 .AND. IX .LE. NX .AND.
     .      IY .GE. 1 .AND. IY .LE. NY .AND.
     .      IZ .GE. 1 .AND. IZ .LE. NZC) THEN
          K = IZ + NZC*(IY-1) + NZC*NY*(IX-1)
          LWC(K) = LWCT
          REFF(K) = REFFT
        ENDIF
      ENDDO
190   CONTINUE
      CLOSE (1)
      RETURN
      END




      SUBROUTINE READ_AERO_FILE (AEROFILE, MAXNZ, NZA, ZAERO,
     .                           EXT55AERO, REFFAERO, WIDTHAERO, 
     .                           AERO55INDEX, AEROINDEX)
C       Reads in the aerosol property file.  The 0.55 micron extinction,
C     the effective radius, and distribution width (alpha or sigma), and 
C     index of refraction at 0.55 microns and the current wavelength are 
C     defined at particular heights.
      IMPLICIT NONE
      INTEGER MAXNZ, NZA
      REAL    ZAERO(*), EXT55AERO(*), REFFAERO(*), WIDTHAERO(*)
      COMPLEX AERO55INDEX(*), AEROINDEX(*)
      CHARACTER*(*) AEROFILE
      INTEGER I
      REAL    MR, MI, MR0, MI0

      OPEN (UNIT=1, FILE=AEROFILE, STATUS='OLD')
      I = 1
      DO WHILE (.TRUE.)
        READ (1,*,END=190) ZAERO(I), EXT55AERO(I), 
     .       REFFAERO(I), WIDTHAERO(I), MR0, MI0, MR, MI
        AERO55INDEX(I) = CMPLX(MR0,MI0)
        AEROINDEX(I) = CMPLX(MR,MI)
        I = I + 1
        IF (I .GT. MAXNZ) STOP 'READ_AERO_FILE: MAXNZ exceeded'
      ENDDO
190   CONTINUE
      CLOSE (1)
      NZA = I - 1
      IF (NZA .LT. 2) STOP 'READ_AERO_FILE: must be at least two levels'
      RETURN
      END



      SUBROUTINE OUTPUT_PROPERTIES (PROPFILE, NX, NY, NZT, 
     .             DELX,DELY,ZLEVELS, TEMP, EXTAERO, SSALBAERO, 
     .             ITABAERO, NTABAERO, NLEGAERO, MAXLEGA, LEGENAERO,
     .             NZC, ICLDLEV, LWC, REFF, NRETAB, SRETAB, ERETAB, 
     .             EXTINCT, ALBEDO, MAXLEG, MAXOUTLEG, NLEG, LEGEN)
C       Interpolates the effective radius the scattering property table 
C     and output the atmosphere property file.  NRETAB+NTABAERO+1 tabulated
C     phase functions are output: 1 is isotropic (for no scattering),
C     2 to NTABAERO+1 is aerosol+Rayleigh scattering (for above the 
C     extinction cutoff), and the rest are for Mie cloud scattering 
C     (for above the cutoff).
C     For in cloud scattering, the extinction and albedo are linearly
C     interpolated in effective radius, while the closest Legendre phase 
C     function expansion is used.
      IMPLICIT NONE
      INTEGER NRETAB, NX, NY, NZT, NZC, MAXLEG, MAXOUTLEG
      INTEGER ICLDLEV(NZT), NLEG(NRETAB)
      INTEGER NTABAERO, MAXLEGA, NLEGAERO(NTABAERO), ITABAERO(NZT)
      REAL    DELX, DELY, ZLEVELS(NZT), TEMP(NZT)
      REAL    LWC(NZC,NY,NX), REFF(NZC,NY,NX)
      REAL    SRETAB, ERETAB
      REAL    EXTINCT(NRETAB), ALBEDO(NRETAB), LEGEN(0:MAXLEG,NRETAB)
      REAL    EXTAERO(NZT), SSALBAERO(NZT)
      REAL    LEGENAERO(0:MAXLEGA,NTABAERO)
      CHARACTER PROPFILE*(*)
      INTEGER IX, IY, IZ, IZC, I, J, L, IPH, NL, NL1
      REAL    EXTCUT, EXT, ALB, EXTC, ALBC, RI, F, REDELINV

      EXTCUT = 1.0E-4

C         Open the property file and output the header including the
C         tabulated phase functions
      OPEN (UNIT=2, FILE=PROPFILE, STATUS='UNKNOWN')
      WRITE (2,'(A)') 'Tabulated phase function property file'
      WRITE (2,'(3(1X,I5))') NX, NY, NZT
      WRITE (2,'(2(1X,F7.4),202(1X,F8.4))') DELX, DELY, 
     .           (ZLEVELS(I),I=1,NZT)
      WRITE (2,'(1X,I4)') NRETAB+NTABAERO+1
      WRITE (2,'(1X,I4)') 0
      DO I = 1, NTABAERO
        NL = NLEGAERO(I)
        WRITE (2,'(1X,I4,5X,200(1X,F9.5))') 
     .       NL,(LEGENAERO(L,I), L=1,MIN(NL,200))
        DO J = 200, NL-1, 200
          WRITE (2,'(200(1X,F9.5))') 
     .          (LEGENAERO(J+L,I), L=1,MIN(NL-J,200))
        ENDDO

      ENDDO
      DO I = 1, NRETAB
        NL = MIN(MAXOUTLEG,NLEG(I))
        WRITE (2,'(1X,I4,5X,200(1X,F9.5))') 
     .       NL,(LEGEN(L,I), L=1,MIN(NL,200))
        DO J = 200, NL-1, 200
          WRITE (2,'(200(1X,F9.5))') (LEGEN(J+L,I), L=1,MIN(NL-J,200))
        ENDDO
      ENDDO

      IF (NRETAB .LE. 1) THEN
        REDELINV = 0.0
      ELSE
        REDELINV = (NRETAB-1)/(ERETAB-SRETAB)
      ENDIF
      DO IX = 1, NX
        DO IY = 1, NY
          DO IZ = 1, NZT
            IZC = ICLDLEV(IZ)
            IF (IZC .GT. 0) THEN
              RI = (REFF(IZC,IY,IX)-SRETAB)*REDELINV + 1
              IF (RI .GT. NRETAB) THEN
                WRITE (*,'(A,3I4,F12.6,F12.6)')
     .              'Warning: effective radius beyond table :', 
     .              IX, IY, IZ, LWC(IZC,IY,IX), REFF(IZC,IY,IX)
              ENDIF
              RI = MAX(1.0,MIN(FLOAT(NRETAB),RI))
              I = INT(RI)
              F = RI - I
              IPH = MAX(1,MIN(NRETAB,NINT(RI))) + 1 + NTABAERO
              EXTC = (1-F)*EXTINCT(I) + F*EXTINCT(I+1)
              EXTC = LWC(IZC,IY,IX)*EXTC
              EXT = EXTAERO(IZ) + EXTC
              ALBC = (1-F)*ALBEDO(I) + F*ALBEDO(I+1)
              IF (EXT .GT. 0.0) THEN
                ALB = (ALBC*EXTC+SSALBAERO(IZ)*EXTAERO(IZ))/EXT
              ELSE
                ALB = 1.0
              ENDIF
            ELSE
              EXTC = 0.0
            ENDIF
            IF (EXTC .GT. EXTCUT .AND. EXTC .GT. EXTAERO(IZ)) THEN
              WRITE (2,'(3(1X,I4),1X,F5.1,1X,F9.4,1X,F8.6,1X,I3)') 
     .            IX, IY, IZ, TEMP(IZ), EXT, ALB, IPH
            ELSE IF (EXTAERO(IZ) .GT. EXTCUT) THEN
              WRITE (2,'(3(1X,I4),1X,F5.1,1X,F8.4,1X,F8.6,1X,I3)') 
     .            IX, IY, IZ, TEMP(IZ), 
     .            EXTAERO(IZ), SSALBAERO(IZ), ITABAERO(IZ)+1
            ELSE
              WRITE (2,'(3(1X,I4),1X,F5.1,1X,F8.4,1X,F8.6,1X,I3)') 
     .          IX, IY, IZ, TEMP(IZ), 0.0, 0.0, 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CLOSE (2)
      RETURN
      END






      SUBROUTINE MAKE_REFF (NPTS, DROPCONC, ALPHA, LWC, REFF, distflag)
C       Makes the effective radius REFF [micron] from the liquid water 
C     content LWC [g/m^3] array and the constant droplet concentration
C     DROPCONC [/cm^3] for a gamma droplet size distribution with
C     shape parameter ALPHA or for a lognormal distribution, per distflag.
C     Uses a density of 1 g/cm^3 appropriate for water droplets.
      IMPLICIT NONE
      INTEGER NPTS
      CHARACTER distflag*1
      REAL    DROPCONC, ALPHA, LWC(NPTS), REFF(NPTS)
      INTEGER I

      DO I = 1, NPTS
        IF (distflag .EQ. 'G') THEN
          REFF(I) = 100* ( LWC(I) *(0.75*(ALPHA+3)*(ALPHA+3))
     .            /(3.14159*(ALPHA+1)*(ALPHA+2)*DROPCONC) )**(1.0/3)
        ELSEIF (distflag .EQ. 'L') THEN
          REFF(I) = 100.*( 3.*LWC(I)*EXP(3.*ALPHA**2)/
     .                     (4.*3.14159*DROPCONC) )**(1.0/3)
        ELSE
          WRITE (*,*) 'MAKE_REFF: Unrecognized distflag'
          STOP
        ENDIF
      ENDDO
      RETURN
      END





      SUBROUTINE ORGANIZE_LEVELS (NZC, ZCLD, TEMPCLD, 
     .             NZO, ZOTHER, TEMPOTHER, NZT, ZLEVELS, TEMP, ICLDLEV)
C       Combines the cloud and extra levels, returning the resulting
C     heights (ZLEVELS) and temperature profile (TEMP), and the 
C     pointers to the cloud levels (ICLDLEV).
      IMPLICIT NONE
      INTEGER NZC, NZO, NZT, ICLDLEV(*), nzou
      REAL    ZCLD(NZC), TEMPCLD(NZC), ZOTHER(NZO), TEMPOTHER(NZO)
      REAL    ZLEVELS(*), TEMP(*)
      INTEGER J, K

C          Put the cloud and other Z levels (and temperatures) in one array
      K = 0
      DO J = 1, NZC
        K = K + 1
        ZLEVELS(K) = ZCLD(J)
        TEMP(K) = TEMPCLD(J)
      ENDDO
      nzou=0
      DO J = 1, NZO
C         Do not use extra levels within cloud layer - 
C           cloud properties do not get applied
        IF (ZOTHER(J).LT.ZCLD(1) .OR. ZOTHER(J).GT.ZCLD(NZC)) THEN
          K = K + 1
          ZLEVELS(K) = ZOTHER(J)
          TEMP(K) = TEMPOTHER(J)
          nzou=nzou+1
        ENDIF
      ENDDO
      NZT = NZC+nzou
C          Sort the levels
      CALL SSORT (ZLEVELS, TEMP, NZT, 2)
C          Check to see if two Z levels are same
      DO K = 1, NZT-1
        IF (ZLEVELS(K) .EQ. ZLEVELS(K+1)) THEN
          WRITE (*,*) 'Same Z levels: ', ZLEVELS(K)
          STOP
        ENDIF
      ENDDO
C          Make the cloud level pointers
      J = 1
      DO K = 1, NZT
        ICLDLEV(K) = 0
        IF (ZLEVELS(K) .EQ. ZCLD(J)) THEN
          ICLDLEV(K) = J
          J = MIN(NZC,J+1)
        ENDIF
      ENDDO

      RETURN
      END






      SUBROUTINE RAYLEIGH_EXTINCT (NZT, ZLEVELS,TEMP, RAYLCOEF, EXTRAYL)
C       Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
C     from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
C     a linear lapse rate between levels to compute the pressure at
C     each level.  The Rayleigh extinction is proportional to air
C     density, with the coefficient RAYLCOEF in [K/(mb km)].
      IMPLICIT NONE
      INTEGER NZT
      REAL    ZLEVELS(NZT), TEMP(NZT), RAYLCOEF, EXTRAYL(NZT)
      INTEGER I
      REAL    PRES, LAPSE, TS, DZ

C           Find surface pressure by integrating hydrostatic relation
C           for a dry atmosphere up to surface height.
      PRES = 1013.
      TS = TEMP(1)
      LAPSE = 6.5*0.001
      PRES = PRES*(TS/(TS+LAPSE*ZLEVELS(1)*1000.))**(9.8/(287.*LAPSE))

C         Use layer mean temperature to compute fractional pressure change.
      DO I = 1, NZT-1
        EXTRAYL(I) = RAYLCOEF*PRES/TEMP(I)
        DZ = 1000.*(ZLEVELS(I+1)-ZLEVELS(I))
        LAPSE = (TEMP(I)-TEMP(I+1))/DZ
        IF (ABS(LAPSE) .GT. 0.00001) THEN
          PRES = PRES*(TEMP(I+1)/TEMP(I))**(9.8/(287.*LAPSE))
        ELSE
          PRES = PRES*EXP(-9.8*DZ/(287.*TEMP(I)))
        ENDIF
      ENDDO  
      EXTRAYL(NZT) = RAYLCOEF*PRES/TEMP(NZT)
      RETURN
      END




      SUBROUTINE CALC_MIE_AEROSOLS (WAVELEN, NZA, ZAERO, 
     .                EXT55AERO, REFFAERO, WIDTHAERO, 
     .                AERO55INDEX, AEROINDEX,
     .                NZT, ZLEVELS, EXTRAYL, 
     .                EXTAERO, SSALBAERO, ITABAERO, 
     .                NTABAERO, MAXLEGA, NLEGAERO, LEGENAERO, aerodist)
C      Does the Mie calculations for each input aerosol level.
C      The Mie calculation is done for the desired wavelength and
C      then normalized by the Mie extinction at 0.55 micron.  The
C      aerosol optical properties are interpolated to the output grid.
C      The molecular Rayleigh scattering and aerosol scattering
C      properties are combined for each output level.  The tabulated
C      phase functions are combined Rayleigh and aerosol for the
C      input aerosol levels.
      IMPLICIT NONE
      INTEGER NZA, NZT, ITABAERO(NZT)
      INTEGER NTABAERO, MAXLEGA
      INTEGER NLEGAERO(*)
      REAL    ZAERO(NZA), EXT55AERO(NZA)
      REAL    REFFAERO(NZA), WIDTHAERO(NZA)
      REAL    ZLEVELS(NZT), EXTRAYL(NZT), WAVELEN
      REAL    EXTAERO(NZT), SSALBAERO(NZT)
      REAL    LEGENAERO(0:MAXLEGA,*)
      COMPLEX AERO55INDEX(NZA), AEROINDEX(NZA)
      INTEGER MAXNZA
      PARAMETER (MAXNZA=100)
      REAL    EXT(MAXNZA), ALB(MAXNZA)
      INTEGER I, J, L, NUMRAD
      REAL    ALPHA, SIGMA, A, B, RMIN,RMAX, GAMMLN,TWOPI, XMAX, EXT55
      REAL    EXTR, EXTA, ALBA, SCATA, F
      CHARACTER aerodist*1

 
      IF (NZA .EQ. 0) THEN
C           Do the Rayleigh only case
        NTABAERO = 1
        NLEGAERO(1) = 2
        LEGENAERO(0,1) = 1.0
        LEGENAERO(1,1) = 0.0
        LEGENAERO(2,1) = 0.5
        DO I = 1, NZT
          EXTAERO(I) = EXTRAYL(I)
          SSALBAERO(I) = 1.0
          ITABAERO(I) = 1
        ENDDO
      ELSE
C           Do the Mie scattering for the input aerosol layers.
C             Run Mie code without phase function for the 0.55 um normalization.
        IF (NZA .GT. MAXNZA) STOP 'CALC_MIE_AEROSOLS: MAXNZA exceeded'
        NTABAERO = NZA
        TWOPI = 2.0*ACOS(-1.0)
        DO J = 1, NZA
          ALPHA = WIDTHAERO(J)
          IF (aerodist .EQ. 'L') THEN
            SIGMA = ALPHA
            B = REFFAERO(J)*EXP(-2.5*SIGMA**2)
            A = 1000.0/( (2*TWOPI/3.)* SQRT(TWOPI)*SIGMA 
     .                 * B**3 *EXP(4.5*SIGMA**2) )
          ELSE
            B = (ALPHA+3)/REFFAERO(J)
            A = 1000.*(0.75/3.1416)*B**(ALPHA+4)/EXP(GAMMLN(ALPHA+4))
          ENDIF
          RMIN = REFFAERO(J)*0.01
          RMAX = REFFAERO(J)*10.0
          XMAX = TWOPI*RMAX/WAVELEN
          NUMRAD = MIN(1000,MAX(100,NINT(XMAX*3.0)))
          CALL MIE (0.55, AERO55INDEX(J), RMIN, RMAX, NUMRAD, MAXLEGA,
     .              A, B, ALPHA, 1.0, .FALSE.,
     .             EXT55, ALB(J), NLEGAERO(J), LEGENAERO(0,J) ,aerodist)
          CALL MIE (WAVELEN, AEROINDEX(J), RMIN, RMAX, NUMRAD, MAXLEGA,
     .              A, B, ALPHA, 1.0, .TRUE.,
     .            EXT(J), ALB(J), NLEGAERO(J), LEGENAERO(0,J) ,aerodist)
          EXT(J) = EXT(J)*EXT55AERO(J)/EXT55
        ENDDO
C           Interpolate the Rayleigh extinction to the aerosol input grid
C             so the phase functions can be averaged for the table.
        I = 1
        DO J = 1, NZA
          DO WHILE (ZAERO(J) .GE. ZLEVELS(I+1) .AND. I .LT. NZT-1)
            I = I + 1
          ENDDO
          F = (ZAERO(J)-ZLEVELS(I))/(ZLEVELS(I+1)-ZLEVELS(I))
          F = MIN(MAX(F,0.0),1.0)
          EXTR = (1-F)*EXTRAYL(I) + F*EXTRAYL(I+1)
          SCATA = EXT(J)*ALB(J)
          NLEGAERO(J) = MAX(2,NLEGAERO(J))
          DO L = 1, NLEGAERO(J)
            LEGENAERO(L,J) = LEGENAERO(L,J)*SCATA/(SCATA+EXTR)
          ENDDO
          LEGENAERO(2,J) = LEGENAERO(2,J) + 0.5*EXTR/(SCATA+EXTR)
        ENDDO        
C           Interpolate the aerosol properties to the output grid and
C             add in the Rayleigh scattering at the same time.
        J = 1
        DO I = 1, NZT
          DO WHILE (ZLEVELS(I) .GE. ZAERO(J+1) .AND. J .LT. NZA-1)
            J = J + 1
          ENDDO
          F = (ZLEVELS(I)-ZAERO(J))/(ZAERO(J+1)-ZAERO(J))
          F = MIN(MAX(F,0.0),1.0)
          EXTA = (1-F)*EXT(J) + F*EXT(J+1)
          ALBA = (1-F)*ALB(J) + F*ALB(J+1)
          EXTAERO(I) =  EXTA + EXTRAYL(I)
          SSALBAERO(I) = (EXTA*ALBA+EXTRAYL(I))/EXTAERO(I)
          ITABAERO(I) = NINT(J+F)
        ENDDO
      ENDIF

      RETURN
      END



      SUBROUTINE READ_MIE_TABLE (MIEFILE, MAXNTAB,NRETAB,SRETAB,ERETAB,
     .                           WAVELEN, RINDEX, ALPHA, distflag,
     .                           EXTINCT, ALBEDO, MAXLEG, NLEG, LEGEN)
C       Reads a table of Mie scattering properties as a function of 
C     effective radius.  This includes reading in the parameters in
C     the header, such as wavelength and index of refraction.
      IMPLICIT NONE
      INTEGER MAXNTAB, NRETAB, MAXLEG, NLEG(MAXNTAB)
      REAL    SRETAB, ERETAB, WAVELEN, ALPHA
      REAL    EXTINCT(MAXNTAB),ALBEDO(MAXNTAB), LEGEN(0:MAXLEG,MAXNTAB)
      COMPLEX RINDEX
      CHARACTER*1 distflag
      CHARACTER MIEFILE*(*)
      INTEGER I, L
      REAL    REFF, MR, MI

      OPEN (UNIT=3, FILE=MIEFILE, STATUS='OLD')
      READ (3,*)
      READ (3,*) WAVELEN
      READ (3,*) MR, MI
      RINDEX = CMPLX(MR,MI)
      READ (3,'(A1)') distflag
      READ (3,*) ALPHA
      READ (3,*) NRETAB, SRETAB, ERETAB
      IF (NRETAB .GT. MAXNTAB) 
     .    STOP 'READ_MIE_TABLE: MAXNTAB exceeded'
      DO I = 1, NRETAB
        READ (3,*) REFF, EXTINCT(I), ALBEDO(I), NLEG(I)
        IF (NLEG(I) .GT. MAXLEG) STOP 'READ_MIE_TABLE: MAXLEG exceeded'
        READ (3,*) (LEGEN(L,I), L=0,NLEG(I)) 
      ENDDO
      CLOSE (3)
      RETURN
      END





      SUBROUTINE WRITE_MIE_TABLE (MIEFILE, NRETAB, SRETAB, ERETAB, 
     .                            WAVELEN, RINDEX, ALPHA, 
     .                            EXTINCT, ALBEDO, MAXLEG, MAXTABLEG, 
     .                            NLEG, LEGEN, distflag)
C       Writes the table of Mie scattering properties as a function of 
C     effective radius.  Includes a header for verifying the Mie table
C     parameters.
      IMPLICIT NONE
      INTEGER NRETAB, MAXLEG, MAXTABLEG, NLEG(NRETAB)
      REAL    SRETAB, ERETAB, WAVELEN, ALPHA
      REAL    EXTINCT(NRETAB), ALBEDO(NRETAB), LEGEN(0:MAXLEG,NRETAB)
      COMPLEX RINDEX
      CHARACTER MIEFILE*(*), distflag*1
      INTEGER I, J, L, NL
      REAL    REFF

      OPEN (UNIT=3, FILE=MIEFILE, STATUS='UNKNOWN')
      WRITE (3,'(A)') '! Mie table vs. effective radius (LWC=1 g/m^3)'
      WRITE (3,'(E13.6,A)') WAVELEN, '  wavelength (micron)'
      WRITE (3,'(2(1X,E13.6),A)') RINDEX, '  index of refraction'
      WRITE (3,'(A1,1X,A)') distflag, ' distribution type'
      WRITE (3,'(F7.5,A)') ALPHA, '  distribution shape parameter'
      WRITE (3,'(1X,I3,2(1X,F8.4),A)') NRETAB, SRETAB, ERETAB, 
     .  '  number, starting, ending effective radius'
      DO I = 1, NRETAB
        IF (NRETAB .LE. 1) THEN
          REFF = SRETAB
        ELSE
          REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
        ENDIF
        NL = MIN(MAXTABLEG,NLEG(I))
        WRITE (3,'(1X,F8.4,1X,E12.5,1X,F8.6,1X,I4,A)') REFF, 
     .    EXTINCT(I), ALBEDO(I), NL, '  Reff  Ext  Alb  Nleg'
        WRITE (3,'(2X,201(1X,F9.5))') (LEGEN(L,I), L=0,MIN(NL,200))
        DO J = 200, NL-1, 200
          WRITE (3,'(2X,200(1X,F9.5))') (LEGEN(J+L,I),L=1,MIN(200,NL-J))
        ENDDO
      ENDDO
      CLOSE (3)
      RETURN
      END




      SUBROUTINE MAKE_MIE_TABLE (NRETAB, SRETAB, ERETAB, 
     .                           WAVELEN, PARTTYPE, RINDEX, ALPHA, 
     .                           EXTINCT, ALBEDO, MAXLEG, 
     .                           NLEG, LEGEN, distflag)
C       Makes the table of scattering properties as a function of effective
C     radius.  There are NRETAB entries in the table from SRETAB to ERETAB.
C     The table is made for a gamma or lognormal droplet size distribution
C     with LWC=1 g/m^3, the tabulated effective radius [micron], and shape
C     parameter ALPHA.  The Mie computations are done for a wavelength
C     of WAVELEN [micron] and refractive index RINDEX.
      IMPLICIT NONE
      INTEGER NRETAB, MAXLEG, NLEG(NRETAB)
      REAL    SRETAB, ERETAB, WAVELEN, ALPHA
      REAL    EXTINCT(NRETAB), ALBEDO(NRETAB), LEGEN(0:MAXLEG,NRETAB)
      COMPLEX RINDEX
      CHARACTER PARTTYPE*1, distflag*1
      INTEGER MAXNRAD, MAXLEGR
      PARAMETER (MAXNRAD=5000, MAXLEGR=2500)
      INTEGER I, K, L, NRAD, NLEG0, NL, NLEGR(MAXNRAD)
      REAL    TWOPI, DENS, X, DELX, RAD, DELRAD, RADMIN, RADMAX
      REAL    REFF, LWC, A, B, GAMMLN, CONC, SCATTER
      REAL    RADIUS(MAXNRAD), WTS(MAXNRAD)
      REAL    EXTINCTR(MAXNRAD), SCATTERR(MAXNRAD)
      REAL    LEGENR(0:MAXLEGR,MAXNRAD)

      TWOPI = 2.0*ACOS(-1.0)
       
C         Call the Mie code for all the radii. The spacing of the radii
C           is determined from a power law in the size parameter x
      WRITE (*,*) 'Mie computations:'
      RADMIN = 0.02*SRETAB
      RADMAX = 5.0*ERETAB
      RAD = RADMIN
      K = 1
      DO WHILE (RAD .LE. RADMAX)
        X = TWOPI*RAD/WAVELEN
        DELX = MAX(0.03,0.03*X**0.5)
        DELRAD = DELX*WAVELEN/TWOPI
        RAD = RAD + DELRAD
        K = K + 1
        IF (K .GT. MAXNRAD) STOP 'MAKE_MIE_TABLE: MAXNRAD exceeded'
      ENDDO

      WRITE (*,*) 'Index  Radius  SizePar'
      RAD = RADMIN
      K = 1
      DO WHILE (RAD .LE. RADMAX)
        RADIUS(K) = RAD
        CALL MIE_ONE (WAVELEN, RINDEX, RAD, MAXLEGR,
     .                EXTINCTR(K), SCATTERR(K), NLEGR(K), LEGENR(0,K) )
        X = TWOPI*RAD/WAVELEN
        IF (MOD(K,50) .EQ. 0) THEN
          WRITE (*,'(1X,I4,1X,F7.2,1X,F8.3)') K, RAD, X
        ENDIF
        DELX = MAX(0.03,0.03*X**0.5)
c        DELX = 0.1
        DELRAD = DELX*WAVELEN/TWOPI
        RAD = RAD + DELRAD
        K = K + 1
      ENDDO
      NRAD = K-1
      DO K = 2, NRAD-1
        WTS(K) = (RADIUS(K+1)-RADIUS(K-1))/2.0
      ENDDO
      WTS(1) = (RADIUS(1)+RADIUS(2))/2.0 - RADMIN
      WTS(NRAD) = RADMAX - (RADIUS(NRAD-1)+RADIUS(NRAD))/2.0
      WRITE (*,*) NRAD, ' radii'


      LWC = 1.0
      DENS=1.0
      IF (PARTTYPE .EQ. 'I') DENS=0.916

C         Integrate over the particle size distributions for each Reff
      DO I = 1, NRETAB
        EXTINCT(I) = 0.0
        SCATTER = 0.0   
        DO L = 0, MAXLEG
          LEGEN(L,I) = 0.0
        ENDDO
        NL = 0
C         Make the effective radius and the gamma distribution parameters
C         Set up for Mie code gamma distribution  n(r)=a * r^alpha * exp(-b*r)
C         Reff =  (alpha+3)/b   M = 4/3 pi rho a b^-(alpha+4) * Gamma(alpha+4)
        IF (NRETAB .LE. 1) THEN
          REFF = SRETAB
        ELSE
          REFF = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
        ENDIF
        IF (distflag .EQ. 'G') THEN
          B = (ALPHA+3)/REFF
          A = (0.75/3.14159)*B**(ALPHA+4)*1000./EXP(GAMMLN(ALPHA+4.))
        ELSEIF (distflag .EQ. 'L') THEN
          B=REFF*EXP(-2.5*ALPHA**2)
          A = 1000./( (2*TWOPI/3.)* SQRT(TWOPI)*ALPHA 
     .                 * B**3 *EXP(4.5*ALPHA**2) )
        ELSE
          WRITE (*,*) 'MAKE_MIE_TABLE: unrecognized distflag'
          STOP
        ENDIF
        A = A*LWC/DENS
        DO K = 1, NRAD
          IF (distflag .EQ. 'G') THEN
C             Evaluate gamma distribution at each radius
            CONC = A* RADIUS(K)**ALPHA * EXP(-B*RADIUS(K))
          ELSEIF (distflag .EQ. 'L') THEN
            CONC = A/RADIUS(K)*EXP(-.5*(ALOG(RADIUS(K)/B))**2/ALPHA**2)
          ELSE
            WRITE (*,*) 'MAKE_MIE_TABLE: unrecognized distflag'
            STOP
          ENDIF
          EXTINCT(I) = EXTINCT(I) + CONC*WTS(K)*EXTINCTR(K)
          SCATTER = SCATTER + CONC*WTS(K)*SCATTERR(K)
          NLEG0 = MIN(MAXLEG,NLEGR(K))
          DO L = 0, NLEG0
            LEGEN(L,I) = LEGEN(L,I) + CONC*WTS(K)*LEGENR(L,K)
          ENDDO
          NL = MAX(NL, NLEG0)
        ENDDO
        ALBEDO(I) = SCATTER/EXTINCT(I)
        DO L = 0, NL
          LEGEN(L,I) = LEGEN(L,I)/SCATTER
          IF (LEGEN(L,I) .GT. 0.5E-5) NLEG(I) = L
        ENDDO
        IF (LEGEN(0,I) .LT. 0.9999) THEN
          WRITE (*,*) 'Phase function not normalized for Reff=',
     .                REFF,LEGEN(0,I)
          STOP
        ENDIF
      ENDDO  

      RETURN
      END





      SUBROUTINE GET_REFRACT_INDEX (PARTTYPE,WAVELEN1,WAVELEN2, RINDEX)
C       Returns the index of refraction for water or ice averaged over
C     the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
C     averaging is done at 0.05 micron intervals and is weighted by
C     a Planck function, using a solar temperature for central wavelengths
C     less than 4 microns and otherwise a cloud temperature.
      IMPLICIT NONE
      REAL    WAVELEN1, WAVELEN2
      COMPLEX RINDEX
      CHARACTER PARTTYPE*1
      REAL    WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
      REAL    MRE, MIM, SUMP, SUMMR, SUMMI, A

      DELWAVE = 0.05
      WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
      WAVECUT = 4.0
      IF (WAVECEN .LT. WAVECUT) THEN
        BBTEMP = 5800
      ELSE
        BBTEMP = 270
      ENDIF

      SUMP = 0.0
      SUMMR = 0.0
      SUMMI = 0.0
      WAVE = WAVELEN1
      DO WHILE (WAVE .LE. WAVELEN2)
        PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
        SUMP = SUMP + PLANCK
        IF (PARTTYPE .EQ. 'I') THEN
          CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
        ELSE
          CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
        ENDIF
        SUMMR = SUMMR + PLANCK*MRE
        SUMMI = SUMMI + PLANCK*MIM
        WAVE = WAVE + DELWAVE
      ENDDO
      MRE = SUMMR/SUMP
      MIM = SUMMI/SUMP
      RINDEX = CMPLX(MRE,-MIM)
      RETURN
      END




      SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
C       Returns the Planck weighted center wavelength for water or ice 
C     averaged over the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).
C     Averaging is done using a solar temperature for wavelengths less 
C     than 4 microns and otherwise a cloud temperature.
      IMPLICIT NONE
      REAL    WAVELEN1, WAVELEN2, WAVELENCEN
      REAL    WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
      REAL    SUMP, SUMW

      IF (WAVELEN1 .EQ. WAVELEN2) THEN
        WAVELENCEN = WAVELEN1

      ELSE
        WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
        DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
        WAVECUT = 4.0
        IF (WAVECEN .LT. WAVECUT) THEN
          BBTEMP = 5800
        ELSE
          BBTEMP = 270
        ENDIF
        SUMP = 0.0
        SUMW = 0.0
        WAVE = WAVELEN1
        DO WHILE (WAVE .LE. WAVELEN2)
          PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
          SUMP = SUMP + PLANCK
          SUMW = SUMW + PLANCK*WAVE
          WAVE = WAVE + DELWAVE
        ENDDO
        WAVELENCEN = SUMW/SUMP
        IF (WAVELEN2-WAVELEN1 .GT. 0.2*WAVELENCEN) THEN
          IF (WAVECEN .LT. WAVECUT) THEN
            WAVELENCEN = NINT(100*WAVELENCEN)/100.0
          ELSE
            WAVELENCEN = NINT(10*WAVELENCEN)/10.0
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END




      SUBROUTINE MIE (WAVELENGTH, MINDEX, RAD1, RAD2, NUMRAD, MAXLEG,
     .                AD, BD, ALPHA, GAMMA, PHASEFLAG,
     .                EXTINCTION, ALBEDO, NLEGEN, LEGEN, aerodist)
C       Computes the Mie scattering properties for a gamma or lognormal
C     distribution of spheres.
      IMPLICIT NONE
      INTEGER     MAXLEG, NLEGEN, NUMRAD
      LOGICAL     PHASEFLAG
      REAL        WAVELENGTH, RAD1, RAD2
      REAL        AD, BD, ALPHA, GAMMA
      COMPLEX     MINDEX
      REAL        EXTINCTION, ALBEDO, LEGEN(*)
      INTEGER     MAXN
      PARAMETER   (MAXN=5000)
      REAL*8      PI
      PARAMETER   (PI = 3.14159265358979D0)
      INTEGER     NTERMS, NQUAD, NMIE, NLEG
      INTEGER     I, L, M, IR, MIN0
      REAL*8      X, DELRAD, RADIUS, NDENS, TMP
      REAL*8      QEXT, QSCAT, SCATTER
      REAL*8      DISTRIBUTION
      REAL*8      MU(MAXN), WTS(MAXN)
      REAL*8      P1, PL, PL1, PL2
      REAL*8      SUMQE, SUMQS, ONE
      REAL*8      SUMP1(MAXN), COEF1(MAXN)
      COMPLEX*16  A(MAXN), B(MAXN), MSPHERE
      CHARACTER aerodist*1
 
 
C           Find the maximum number of terms required in the Mie series,
      MSPHERE = MINDEX
      X = 2.0D0*PI*RAD2/WAVELENGTH
      NTERMS = 0
      CALL MIECALC (NTERMS, X, MSPHERE, A, B)
      NLEGEN = 2*NTERMS
      NLEGEN = MIN0(MAXLEG, NLEGEN)
      NQUAD  = (NLEGEN + 2*NTERMS + 2)/2
      IF (NQUAD .GT. MAXN)  STOP 'MIE: MAXN exceeded' 
 
C           Get the Gauss-Legendre quadrature abscissas and weights
      CALL GAUSQUAD (NQUAD, MU, WTS)
 
      SUMQE = 0.0
      SUMQS = 0.0
      DO I = 1, NQUAD
        SUMP1(I) = 0.0
      ENDDO
 
C               Integration loop over radius of spheres
      IF (NUMRAD .GT. 0)  DELRAD = (RAD2-RAD1)/NUMRAD
      DO IR = 1, NUMRAD+1
          RADIUS = RAD1 + (IR-1)*DELRAD
          NDENS = DISTRIBUTION (DBLE(AD), DBLE(BD), 
     .                     DBLE(ALPHA), DBLE(GAMMA), RADIUS, aerodist)
          IF ((IR .EQ. 1 .OR. IR .EQ. NUMRAD+1)
     .                      .AND. NUMRAD .GT. 0) THEN
              NDENS = 0.5*NDENS
          ENDIF
          X = 2.0D0*PI*RADIUS/WAVELENGTH
          NMIE = 0
          CALL MIECALC (NMIE, X, MSPHERE, A, B)
          CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT)
          SUMQE = SUMQE + QEXT*NDENS*RADIUS**2
          SUMQS = SUMQS + QSCAT*NDENS*RADIUS**2
          IF (PHASEFLAG) THEN
            NMIE = MIN0(NMIE, NTERMS)
            DO I = 1, NQUAD
              CALL MIEANGLE (NMIE, A, B, MU(I), P1)
              SUMP1(I) = SUMP1(I) + P1*NDENS
            ENDDO
          ENDIF
        ENDDO
 
 
C           Multiply the sums by the integration delta and other constants
C             Put quadrature weights in angular array for later
      IF (NUMRAD .EQ. 0) DELRAD = 1.0
 
      EXTINCTION = PI*SUMQE*DELRAD
      SCATTER = PI*SUMQS*DELRAD
      ALBEDO = SCATTER/EXTINCTION
 
C         If the phase function is not desired then leave now
      IF (.NOT. PHASEFLAG) RETURN
 
      TMP = (WAVELENGTH**2/(PI*SCATTER)) * DELRAD
      DO I = 1, NQUAD
        SUMP1(I) = TMP*SUMP1(I) *WTS(I)
      ENDDO
 
C           Integrate the angular scattering functions times Legendre
C             polynomials to find the Legendre coefficients
      DO M = 1, NLEGEN+1
        COEF1(M) = 0.0
      ENDDO
C           Use upward recurrence to find Legendre polynomials
      DO I = 1, NQUAD
          PL1 = 1.0
          PL = 1.0
          DO L = 0, NLEGEN
              M = L + 1
              IF (L .GT. 0)  PL = (2*L-1)*MU(I)*PL1/L - (L-1)*PL2/L
              COEF1(M) = COEF1(M) + SUMP1(I)*PL
              PL2 = PL1
              PL1 = PL
          ENDDO
      ENDDO
      NLEG = NLEGEN
      DO L = 0, NLEG
          M = L + 1
          LEGEN(M) = (2*L+1)/2.0 *COEF1(M)
          IF (LEGEN(M) .GT. 1.0E-5)  NLEGEN = L
      ENDDO
       
      RETURN
      END
 
 


 
      REAL*8 FUNCTION DISTRIBUTION (A, B, ALPHA, GAMMA, R, distflag)
C        DISTRIBUTION returns the particle density for a given radius R
C      for a modified gamma distribution specified by A, B, ALPHA, GAMMA:
C           N(r) = a * r^alpha * exp(-b * r^gamma)     .
C      or a log-normal distribution:
C           N(r) = a/r * exp(- ln(r/b)^2 / (2*alpha^2) )     .
C      depending in distflag.
      IMPLICIT NONE
      REAL*8   A, B, ALPHA, GAMMA, R
      CHARACTER distflag*1
 
      IF (DISTFLAG .EQ. 'G') THEN
C           Modified gamma distibution
        DISTRIBUTION = A* R**ALPHA * DEXP(-B*R**GAMMA)
      ELSEIF (DISTFLAG .EQ. 'L') THEN
C           Log-normal distibution
        DISTRIBUTION = A/R *DEXP(-.5*(DLOG(R/B))**2/ALPHA**2)
      ELSE
        WRITE (*,*) 'Unrecognized distflag in DISTRIBUTION'
      ENDIF
 
      RETURN
      END





      SUBROUTINE SSORT (X, Y, N, KFLAG)
C***BEGIN PROLOGUE  SSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SSORT sorts array X and optionally makes the same interchanges in
C   array Y.  The array X may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      Y - array to be (optionally) carried along
C      N - number of values in array X to be sorted
C      KFLAG - control parameter
C            =  2  means sort X in increasing order and carry Y along.
C            =  1  means sort X in increasing order (ignoring Y)
C            = -1  means sort X in decreasing order (ignoring Y)
C            = -2  means sort X in decreasing order and carry Y along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***END PROLOGUE  SSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
c      REAL X(*), Y(*)
      REAL X(*)
      INTEGER Y(*)
C     .. Local Scalars ..
c      REAL R, T, TT, TTY, TY
      REAL R, T, TT
      INTEGER TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(51), IU(51)
C     .. External Subroutines ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***First executable statement  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         STOP 'The number of values to be sorted is not positive.'
      ENDIF

      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        STOP 'The sort control parameter, K, is not 2, 1, -1, or -2.'
      ENDIF

C     Alter array X to get decreasing order if needed
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

C     Sort X only
      M = 1
      I = 1
      J = NN
      R = 0.375E0

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   30 K = I

C     Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)

C     If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J

C     If last element of array is less than than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

C        If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF

C     Find an element in the second half of the array which is smaller
C     than T
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40

C     Find an element in the first half of the array which is greater
C     than T
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50

C     Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF

C     Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

C     Begin again on another portion of the unsorted array
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I

   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80

C     Sort X and carry Y along
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I

C     Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)

C     If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J

C     If last element of array is less than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

C        If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF

C     Find an element in the second half of the array which is smaller
C     than T
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130

C     Find an element in the first half of the array which is greater
C     than T
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140

C     Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF

C     Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

C     Begin again on another portion of the unsorted array
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I

  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170

C     Clean up
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END


