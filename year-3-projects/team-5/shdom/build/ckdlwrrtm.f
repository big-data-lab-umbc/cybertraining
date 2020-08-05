C ************************  PROGRAM CKDLWRRTM  *******************************
C
C  This program and its subroutines create a CKD file for use with SHDOM,
C    using the accurate RRTM correlated k-distribution parameterization
C    for 16 bands in the longwave.
C    
C  They take the following steps:
C
C    a)  read in the user's atmospheric profile
C    b)  set the concentrations of extra molecules
C    c)  compute the appropriate absorption coefficients for all bands
C    d)  find the subinterval weights for all bands
C    e)  write the profile, weight, and absorption information to a CKD file
C
C  Many thanks to the authors of RRTM for developing this parameterization:
C    Eli J. Mlawer, Steven J. Taubman, and Shepard A. Clough
C
C  The RRTM reference is:
C      Mlawer, E. J., S. J. Taubman, P. D. Brown, M. J. Iacono, 
C    S. A. Clough, 1997: Radiative transfer for inhomogeneous atmospheres: 
C    RRTM, a validated correlated-k model of the longwave.  
C    J. Geophys. Res., 102, 16663-16682.
C
C  Much of this code is taken directly from the portions of RRTM that
C    compute absorption.
C
C  There are, however, some appropriate simplifications, such as fewer
C    inputs and only one set of weights per band.  This code does not
C    operate exactly the same as RRTM, but the upwelling and downwelling
C    fluxes from SHDOM agree with RRTM to within 2 W/m^2 at all levels
C    in standard atmospheres.
C
C  Written by Timothy C. Benner with the guidance of K. Franklin Evans
C    at the University of Colorado, July 1999.
C
C ****************************************************************************

      PROGRAM CKDRRTM

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS)
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PLANKG/    FRACS(NBANDS,MAXLEV,MG)
      COMMON /PRECISE/   ONEMINUS
      COMMON /PROFILE/   NLEVELS,HEIGHT(MAXLEV),
     &                   PRES(MAXLEV),TEMP(MAXLEV),
     &                   RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /TAUGCOM/   TAUG(MAXLEV,MG)

      DATA WAVENUM1(1) /10./,    WAVENUM2(1) /250./
      DATA WAVENUM1(2) /250./,   WAVENUM2(2) /500./
      DATA WAVENUM1(3) /500./,   WAVENUM2(3) /630./
      DATA WAVENUM1(4) /630./,   WAVENUM2(4) /700./
      DATA WAVENUM1(5) /700./,   WAVENUM2(5) /820./
      DATA WAVENUM1(6) /820./,   WAVENUM2(6) /980./
      DATA WAVENUM1(7) /980./,   WAVENUM2(7) /1080./
      DATA WAVENUM1(8) /1080./,  WAVENUM2(8) /1180./
      DATA WAVENUM1(9) /1180./,  WAVENUM2(9) /1390./
      DATA WAVENUM1(10) /1390./, WAVENUM2(10) /1480./
      DATA WAVENUM1(11) /1480./, WAVENUM2(11) /1800./
      DATA WAVENUM1(12) /1800./, WAVENUM2(12) /2080./
      DATA WAVENUM1(13) /2080./, WAVENUM2(13) /2250./
      DATA WAVENUM1(14) /2250./, WAVENUM2(14) /2380./
      DATA WAVENUM1(15) /2380./, WAVENUM2(15) /2600./
      DATA WAVENUM1(16) /2600./, WAVENUM2(16) /3000./

      DATA NG  /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /1, 1,10, 9, 9, 1, 9, 1,11, 1, 1, 9, 9, 1, 9, 9/
      DATA NSPB /1, 1, 5, 6, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0/

      CHARACTER  ATMFILE*72, OUTFILE*72
      REAL       UMCO2, UMCH4, UMN2O
      INTEGER    IB, IG, L
      REAL       SOLARFLUX(NBANDS)
      REAL       ABSORB(NBANDS,MAXLEV,MG), DELG(NBANDS,MG)
      REAL       LAYABS(MAXLEV), TOTABS, DG

      DATA       SOLARFLUX /NBANDS*0.0/

      ONEMINUS = 1. - 1.E-6

C     Accept user inputs - same as ckdfu.

      UMCO2 = 360
      UMCH4 = 1.7
      UMN2O = 0.3

      WRITE (*,*) 'Atmospheric profile file '
      READ (*,'(A)') ATMFILE
      WRITE (*,*) '   ', ATMFILE

      WRITE (*,*) 'Output k-distribution file '
      READ (*,'(A)') OUTFILE
      WRITE (*,*) '   ', OUTFILE

      WRITE (*,'(A,A)') ' Surface concentration of CO2, CH4, and N2O',
     .                  ' (ppmv) (try: 360, 1.7, 0.3)'
      READ (*,*) UMCO2, UMCH4, UMN2O
      WRITE (*,*) UMCO2, UMCH4, UMN2O

C     Read atmospheric profile from ATMFILE.

      WRITE (*,*) 'Reading atmospheric profile'

      CALL READ_ATMOS(ATMFILE,UMCO2,UMCH4,UMN2O)

C     Set the concentrations of "cross-section" gas molecules.

      WRITE (*,*) 'Setting cross-section concentrations'

      CALL SET_XSECT

C     Calculate information needed by the radiative transfer routine
C     that is specific to this atmosphere, especially some of the 
C     coefficients and indices needed to compute the optical depths
C     by interpolating data from stored reference atmospheres. 

      WRITE (*,*) 'Setting coefficients'

      CALL SETCOEF

C     Get the absorption coefficients.

      WRITE (*,*) 'Getting absorption coefficients'

      CALL TAUGB1
      DO L = 1, NLEVELS
          DO IG = 1, NG(1)
              ABSORB(1,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB2
      DO L = 1, NLEVELS
          DO IG = 1, NG(2)
              ABSORB(2,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB3
      DO L = 1, NLEVELS
          DO IG = 1, NG(3)
              ABSORB(3,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB4
      DO L = 1, NLEVELS
          DO IG = 1, NG(4)
              ABSORB(4,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB5
      DO L = 1, NLEVELS
          DO IG = 1, NG(5)
              ABSORB(5,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB6
      DO L = 1, NLEVELS
          DO IG = 1, NG(6)
              ABSORB(6,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB7
      DO L = 1, NLEVELS
          DO IG = 1, NG(7)
              ABSORB(7,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB8
      DO L = 1, NLEVELS
          DO IG = 1, NG(8)
              ABSORB(8,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB9
      DO L = 1, NLEVELS
          DO IG = 1, NG(9)
              ABSORB(9,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB10
      DO L = 1, NLEVELS
          DO IG = 1, NG(10)
              ABSORB(10,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB11
      DO L = 1, NLEVELS
          DO IG = 1, NG(11)
              ABSORB(11,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB12
      DO L = 1, NLEVELS
          DO IG = 1, NG(12)
              ABSORB(12,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB13
      DO L = 1, NLEVELS
          DO IG = 1, NG(13)
              ABSORB(13,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB14
      DO L = 1, NLEVELS
          DO IG = 1, NG(14)
              ABSORB(14,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB15
      DO L = 1, NLEVELS
          DO IG = 1, NG(15)
              ABSORB(15,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO
      CALL TAUGB16
      DO L = 1, NLEVELS
          DO IG = 1, NG(16)
              ABSORB(16,L,IG) = TAUG(L,IG)
          ENDDO
      ENDDO

C     Compute the subinterval weights for each band by using an
C     absorption-weighted average of the weights at each level.

      DO IB = 1, 16
         TOTABS = 0.0
         DO L = 1, NLEVELS
            LAYABS(L) = 0.0
            DO IG = 1, NG(IB)
               LAYABS(L) = LAYABS(L)+FRACS(IB,L,IG)*ABSORB(IB,L,IG)
            ENDDO
            TOTABS = TOTABS+LAYABS(L)
         ENDDO
         DO IG = 1, NG(IB)
            DG = 0.0
            DO L=1, NLEVELS
               DG = DG+FRACS(IB,L,IG)*LAYABS(L)
            ENDDO
            DELG(IB,IG) = DG/TOTABS
         ENDDO
      ENDDO

C     Write CKD information to OUTFILE.
C     Transform back to SHDOM format of height decreasing with level.

      WRITE (*,*) 'Writing CKD information'

      OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')

      WRITE (2,'(A)') '! Correlated k-distribution file: method of RRTM'
      WRITE (2,'(I3,A)') NBANDS, '  ! number of bands'
      WRITE (2,'(A,A)') '!Band  Wavenums    SolF  Ng  Delta g'

      DO IB = 1, NBANDS
        WRITE (2,'(1X,I2,2(1X,F6.0),1X,F6.2,1X,I2,16(1X,F7.5))')
     .         IB, WAVENUM1(IB), WAVENUM2(IB), SOLARFLUX(IB),
     .         NG(IB), (DELG(IB,IG), IG=1, NG(IB))
      ENDDO

      WRITE (2,'(I3,A)') NLEVELS, '  ! number of levels'
      WRITE (2,'(1X,F5.1,2(1X,F5.3),A)')  UMCO2, UMCH4, UMN2O,  
     .       '  ! concentration of CO2, CH4, N2O in ppmv'
      WRITE (2,'(A)') '!  Z     Pres  Temp      qH2O        qO3'
      DO L = NLEVELS, 1, -1
        WRITE (2,'(1X,F7.3,1X,F6.1,1X,F5.1,2(1X,E11.4))')
     .        HEIGHT(L), PRES(L), TEMP(L), 
     .        RHOH2O(L)/RHOAIR(L), RHOO3(L)/RHOAIR(L)
      ENDDO

      WRITE (2,'(A)') '!IB IZ    Kabs (km^-1)'
      BAD = 0
      DO IB = 1, NBANDS
        DO L = NLEVELS, 1, -1
          WRITE (2,'(1X,I2,1X,I3,16(1X,E12.5))')
     .        IB, NLEVELS-L+1, (ABSORB(IB,L,IG), IG=1, NG(IB))
          DO IG = 1, NG(IB)
            IF ((ABSORB(IB,L,IG) .LT. 0) .AND. (BAD .LT. 10)) THEN
              WRITE (*,*) '***** Negative absorption(s):  ', IB, IG
              BAD = BAD+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CLOSE (2)

      END


C ************************  SUBROUTINE READ_ATMOS  ***************************
C
C  Read in the atmospheric profile from a user-provided file.
C
C      1  -  H2O
C      2  -  CO2
C      3  -  O3
C      4  -  N2O
C      5  -  CO
C      6  -  CH4
C      7  -  O2
C
C ****************************************************************************

      SUBROUTINE READ_ATMOS(ATMFILE,UMCO2,UMCH4,UMN2O)

      IMPLICIT NONE

C     Maximum number of atmospheric levels
      
      INTEGER    MAXLEV
      PARAMETER (MAXLEV = 201)

C     Molecular weight of dry air, water, ozone [g/mole]
      
      REAL       MDA, MDW, MDO
      PARAMETER (MDA = 28.964, MDW = 18.015, MDO = 47.997)
      
C     Avagadro's # [# mol/mole]  

      REAL       NA  
      PARAMETER (NA = 6.022E+23)

C     Common block

      COMMON /PROFILE/   NLEVELS,HEIGHT(MAXLEV),
     &                   PRES(MAXLEV),TEMP(MAXLEV),
     &                   RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /SPECIES/   COLDRY(MAXLEV),WKL(35,MAXLEV)

      INTEGER    NLEVELS
      REAL       HEIGHT, PRES, TEMP, RHOAIR, RHOH2O, RHOO3
      REAL       COLDRY, WKL

C     Other variables

      CHARACTER*(*) ATMFILE
      REAL          UMCO2, UMCH4, UMN2O, BUFFER(6,MAXLEV)
      INTEGER       L

C     Open and read the ATMFILE file.

      OPEN (UNIT=1, FILE=ATMFILE, STATUS='OLD')
      L = 1
100   CONTINUE
      READ (1,*,END=120) BUFFER(1,L), BUFFER(2,L), BUFFER(3,L),
     .                   BUFFER(4,L), BUFFER(5,L), BUFFER(6,L)
      L = L + 1
      IF (L .LT. MAXLEV) GOTO 100 
120   CONTINUE
      NLEVELS = L - 1
      CLOSE (1)

C     Transform to RRTM format of height increasing with level.

      DO L = 1, NLEVELS
         HEIGHT(L) = BUFFER(1,NLEVELS-L+1)
         PRES(L) = BUFFER(2,NLEVELS-L+1)
         TEMP(L) = BUFFER(3,NLEVELS-L+1)
         RHOAIR(L) = BUFFER(4,NLEVELS-L+1)  
         RHOH2O(L) = BUFFER(5,NLEVELS-L+1)
         RHOO3(L) = BUFFER(6,NLEVELS-L+1)
      ENDDO

C     Convert ATMFILE gas densities from [g/m^3] to [# mol/cm^2/km].

      DO L = 1, NLEVELS
         COLDRY(L) = RHOAIR(L)*NA/(MDA*10.00)
         WKL(1,L) = RHOH2O(L)*NA/(MDW*10.00)
         WKL(2,L) = COLDRY(L)*UMCO2*1.0E-6
         WKL(3,L) = RHOO3(L)*NA/(MDO*10.00)
         WKL(4,L) = COLDRY(L)*UMN2O*1.0E-6
         WKL(5,L) = 0.0
         WKL(6,L) = COLDRY(L)*UMCH4*1.0E-6
         WKL(7,L) = 0.0
      ENDDO

      WRITE (*,*) '  Number of levels = ', NLEVELS

      RETURN
      END 


C ************************  SUBROUTINE SET_XSECT  ****************************
C
C  Set the concentrations of the "cross-section" gas molecules.
C
C      1  -  CCL4
C      2  -  CFC11
C      3  -  CFC12
C      4  -  CFC22
C
C  Mixing ratios come from RRTM's RRTATM.F data blocks.
C  Set the standard altitude and the mixing ratios [ppmv] for each gas, 
C    then interpolate them to the user's altitude levels,
C    while multiplying by COLDRY to get concentrations [mol/cm^2/km].
C
C ****************************************************************************

      SUBROUTINE SET_XSECT

      IMPLICIT NONE

C     Parameters

      INTEGER    MAXLEV, MAXXSEC, MAXPROD, NUMDATA, NUMEXTRA
      PARAMETER (MAXLEV = 201)
      PARAMETER (MAXXSEC = 4)
      PARAMETER (MAXPROD = MAXLEV*MAXXSEC)
      PARAMETER (NUMDATA = 50)
      PARAMETER (NUMEXTRA = MAXLEV-NUMDATA)

C     Common block

      COMMON /PROFILE/   NLEVELS,HEIGHT(MAXLEV),
     &                   PRES(MAXLEV),TEMP(MAXLEV),
     &                   RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /SPECIES/   COLDRY(MAXLEV),WKL(35,MAXLEV)
      COMMON /XSEC/      WX(MAXXSEC,MAXLEV)

      INTEGER    NLEVELS
      REAL       HEIGHT, PRES, TEMP, RHOAIR, RHOH2O, RHOO3
      REAL       COLDRY, WKL, WX

C     Data

      DATA WX /MAXPROD*0.0/

      REAL ALTZ(MAXLEV)
      REAL XMOL1(MAXLEV),XMOL2(MAXLEV),XMOL3(MAXLEV),XMOL4(MAXLEV)

C     Set the standard altitude array [km].

      DATA ALTZ /
     *       0.0,       1.0,       2.0,       3.0,       4.0,
     *       5.0,       6.0,       7.0,       8.0,       9.0,
     *      10.0,      11.0,      12.0,      13.0,      14.0,
     *      15.0,      16.0,      17.0,      18.0,      19.0,
     *      20.0,      21.0,      22.0,      23.0,      24.0,
     *      25.0,      27.5,      30.0,      32.5,      35.0,
     *      37.5,      40.0,      42.5,      45.0,      47.5,
     *      50.0,      55.0,      60.0,      65.0,      70.0,
     *      75.0,      80.0,      85.0,      90.0,      95.0,
     *     100.0,     105.0,     110.0,     115.0,     120.0,
     *     NUMEXTRA*0.0 /

C     Set CCL4 mixing ratios [ppmv].

      DATA XMOL1 /
     *  1.300E-04, 1.300E-04, 1.299E-04, 1.299E-04, 1.298E-04,
     *  1.297E-04, 1.296E-04, 1.295E-04, 1.294E-04, 1.293E-04,
     *  1.292E-04, 1.289E-04, 1.285E-04, 1.266E-04, 1.247E-04,
     *  1.187E-04, 1.127E-04, 1.026E-04, 9.256E-05, 8.037E-05,
     *  6.817E-05, 5.611E-05, 4.405E-05, 3.395E-05, 2.385E-05,
     *  1.701E-05, 5.027E-06, 8.202E-07, 1.204E-07, 1.304E-08,
     *  1.050E-09, 4.864E-11, 5.081E-12, 5.372E-13, 5.548E-14,
     *  5.688E-15, 2.281E-16, 5.092E-18, 1.699E-19, 3.184E-21,
     *  9.600E-23, 1.638E-24, 4.605E-26, 6.985E-28, 1.743E-29,
     *  2.224E-31, 4.283E-33, 0.       , 0.       , 0.       ,
     *  NUMEXTRA*0.0 /

C     Set CFC11 mixing ratios [ppmv].

      DATA XMOL2 /
     *  1.400E-04, 1.400E-04, 1.399E-04, 1.399E-04, 1.398E-04,
     *  1.397E-04, 1.396E-04, 1.396E-04, 1.395E-04, 1.394E-04,
     *  1.392E-04, 1.389E-04, 1.386E-04, 1.368E-04, 1.349E-04,
     *  1.292E-04, 1.234E-04, 1.138E-04, 1.041E-04, 9.216E-05,
     *  8.021E-05, 6.799E-05, 5.576E-05, 4.480E-05, 3.384E-05,
     *  2.550E-05, 9.634E-06, 2.441E-06, 5.553E-07, 1.024E-07,
     *  1.581E-08, 1.939E-09, 3.811E-10, 7.716E-11, 1.585E-11,
     *  3.658E-12, 4.173E-13, 3.465E-14, 3.353E-15, 2.383E-16,
     *  2.084E-17, 1.346E-18, 1.080E-19, 6.099E-21, 4.246E-22,
     *  1.923E-23, 1.110E-24, 5.158E-26, 3.393E-27, 3.738E-28,
     *  NUMEXTRA*0.0 /

C     Set CFC12 mixing ratios [ppmv].

      DATA XMOL3 /
     *  2.400E-04, 2.400E-04, 2.399E-04, 2.399E-04, 2.398E-04,
     *  2.398E-04, 2.397E-04, 2.396E-04, 2.395E-04, 2.394E-04,
     *  2.393E-04, 2.390E-04, 2.387E-04, 2.370E-04, 2.353E-04,
     *  2.300E-04, 2.247E-04, 2.157E-04, 2.066E-04, 1.952E-04,
     *  1.838E-04, 1.712E-04, 1.585E-04, 1.452E-04, 1.319E-04,
     *  1.183E-04, 8.552E-05, 5.683E-05, 3.498E-05, 2.013E-05,
     *  1.111E-05, 6.014E-06, 3.446E-06, 1.998E-06, 1.181E-06,
     *  7.687E-07, 3.876E-07, 1.818E-07, 8.265E-08, 3.432E-08,
     *  1.380E-08, 4.984E-09, 1.704E-09, 4.917E-10, 1.272E-10,
     *  2.351E-11, 3.640E-12, 4.251E-13, 4.981E-14, 8.792E-15,
     *  NUMEXTRA*0.0 /

C     Set CFC22 mixing ratios [ppmv].

      DATA XMOL4 /
     *  6.000E-05, 5.994E-05, 5.987E-05, 5.982E-05, 5.977E-05,
     *  5.974E-05, 5.970E-05, 5.968E-05, 5.966E-05, 5.963E-05,
     *  5.960E-05, 5.955E-05, 5.949E-05, 5.921E-05, 5.893E-05,
     *  5.808E-05, 5.723E-05, 5.582E-05, 5.441E-05, 5.265E-05,
     *  5.089E-05, 4.897E-05, 4.705E-05, 4.502E-05, 4.298E-05,
     *  4.084E-05, 3.548E-05, 3.021E-05, 2.514E-05, 2.062E-05,
     *  1.686E-05, 1.392E-05, 1.184E-05, 1.036E-05, 9.356E-06,
     *  8.784E-06, 8.163E-06, 7.741E-06, 7.449E-06, 7.201E-06,
     *  6.919E-06, 6.524E-06, 5.872E-06, 4.867E-06, 3.396E-06,
     *  1.808E-06, 6.935E-07, 2.066E-07, 5.485E-08, 1.930E-08,
     *  NUMEXTRA*0.0 /

C     Other variables

      INTEGER L, L1, L2
      REAL    IFRAC, C0

C     Interpolate each gas at each user level.

      DO L = 1, NLEVELS
C         Get lower and upper indices.
          L1 = 1
          DO WHILE (ALTZ(L1) .LT. HEIGHT(L) .AND. L1 .LT. NUMDATA)
              L1 = L1 + 1
          ENDDO
          L1 = MAX(L1 - 1, 1)
          L2 = L1 + 1
C         Compute how far between indices to go (i.e. what fraction).
          IFRAC = (HEIGHT(L)-ALTZ(L1))/(ALTZ(L2)-ALTZ(L1))
C         Interpolate gases, multiply to get concentrations.
C         Multiply by 1.0E-06 to convert ppmv -> true mixing ratio.
C         Multiply by 1.0E-20 to get column amounts, as SETCOEF does for WKL
          C0 = COLDRY(L)*1.0E-06*1.0E-20
          WX(1,L) = C0*(XMOL1(L1)+(XMOL1(L2)-XMOL1(L1))*IFRAC)
          WX(2,L) = C0*(XMOL2(L1)+(XMOL2(L2)-XMOL2(L1))*IFRAC)
          WX(3,L) = C0*(XMOL3(L1)+(XMOL3(L2)-XMOL3(L1))*IFRAC)
          WX(4,L) = C0*(XMOL4(L1)+(XMOL4(L2)-XMOL4(L1))*IFRAC)
      ENDDO

      RETURN
      END


C ************************  SUBROUTINE SETCOEF  ******************************
C
C  For a given atmosphere, calculate the indices and fractions
C    related to the pressure and temperature interpolations.
C
C ****************************************************************************

      SUBROUTINE SETCOEF

      INTEGER    MAXLEV
      PARAMETER (MAXLEV=201)

C  Input      
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /SPECIES/  COLDRY(MAXLEV),WKL(35,MAXLEV)

C  Output
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC, SELFFRAC, INDSELF
      COMMON /FOREIGN/  FORFAC(MAXLEV)

      DIMENSION SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      DIMENSION PREF(59),PREFLOG(59),TREF(59)

C     These pressures are chosen such that the ln of the first pressure
C     has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
C     each subsequent ln(pressure) differs from the previous one by 0.2.
      DATA PREF /
     &    1.05363E+03,8.62642E+02,7.06272E+02,5.78246E+02,4.73428E+02,
     &    3.87610E+02,3.17348E+02,2.59823E+02,2.12725E+02,1.74164E+02,
     &    1.42594E+02,1.16746E+02,9.55835E+01,7.82571E+01,6.40715E+01,
     &    5.24573E+01,4.29484E+01,3.51632E+01,2.87892E+01,2.35706E+01,
     &    1.92980E+01,1.57998E+01,1.29358E+01,1.05910E+01,8.67114E+00,
     &    7.09933E+00,5.81244E+00,4.75882E+00,3.89619E+00,3.18993E+00,
     &    2.61170E+00,2.13828E+00,1.75067E+00,1.43333E+00,1.17351E+00,
     &    9.60789E-01,7.86628E-01,6.44036E-01,5.27292E-01,4.31710E-01,
     &    3.53455E-01,2.89384E-01,2.36928E-01,1.93980E-01,1.58817E-01,
     &    1.30029E-01,1.06458E-01,8.71608E-02,7.13612E-02,5.84256E-02,
     &    4.78349E-02,3.91639E-02,3.20647E-02,2.62523E-02,2.14936E-02,
     &    1.75975E-02,1.44076E-02,1.17959E-02,9.65769E-03/
      DATA PREFLOG /
     &     6.9600E+00, 6.7600E+00, 6.5600E+00, 6.3600E+00, 6.1600E+00,
     &     5.9600E+00, 5.7600E+00, 5.5600E+00, 5.3600E+00, 5.1600E+00,
     &     4.9600E+00, 4.7600E+00, 4.5600E+00, 4.3600E+00, 4.1600E+00,
     &     3.9600E+00, 3.7600E+00, 3.5600E+00, 3.3600E+00, 3.1600E+00,
     &     2.9600E+00, 2.7600E+00, 2.5600E+00, 2.3600E+00, 2.1600E+00,
     &     1.9600E+00, 1.7600E+00, 1.5600E+00, 1.3600E+00, 1.1600E+00,
     &     9.6000E-01, 7.6000E-01, 5.6000E-01, 3.6000E-01, 1.6000E-01,
     &    -4.0000E-02,-2.4000E-01,-4.4000E-01,-6.4000E-01,-8.4000E-01,
     &    -1.0400E+00,-1.2400E+00,-1.4400E+00,-1.6400E+00,-1.8400E+00,
     &    -2.0400E+00,-2.2400E+00,-2.4400E+00,-2.6400E+00,-2.8400E+00,
     &    -3.0400E+00,-3.2400E+00,-3.4400E+00,-3.6400E+00,-3.8400E+00,
     &    -4.0400E+00,-4.2400E+00,-4.4400E+00,-4.6400E+00/
C     These are the temperatures associated with the respective 
C     pressures for the MLS standard atmosphere. 
      DATA TREF /
     &     2.9420E+02, 2.8799E+02, 2.7894E+02, 2.6925E+02, 2.5983E+02,
     &     2.5017E+02, 2.4077E+02, 2.3179E+02, 2.2306E+02, 2.1578E+02,
     &     2.1570E+02, 2.1570E+02, 2.1570E+02, 2.1706E+02, 2.1858E+02,
     &     2.2018E+02, 2.2174E+02, 2.2328E+02, 2.2479E+02, 2.2655E+02,
     &     2.2834E+02, 2.3113E+02, 2.3401E+02, 2.3703E+02, 2.4022E+02,
     &     2.4371E+02, 2.4726E+02, 2.5085E+02, 2.5457E+02, 2.5832E+02,
     &     2.6216E+02, 2.6606E+02, 2.6999E+02, 2.7340E+02, 2.7536E+02,
     &     2.7568E+02, 2.7372E+02, 2.7163E+02, 2.6955E+02, 2.6593E+02,
     &     2.6211E+02, 2.5828E+02, 2.5360E+02, 2.4854E+02, 2.4348E+02,
     &     2.3809E+02, 2.3206E+02, 2.2603E+02, 2.2000E+02, 2.1435E+02,
     &     2.0887E+02, 2.0340E+02, 1.9792E+02, 1.9290E+02, 1.8809E+02,
     &     1.8329E+02, 1.7849E+02, 1.7394E+02, 1.7212E+02/

      STPFAC = 296./1013.

      LAYTROP = 0
      LAYSWTCH = 0
      LAYLOW = 0
      DO 7000 LAY = 1, NLEVELS

C        Find the two reference pressures on either side of the
C        layer pressure.  Store them in JP and JP1.  Store in FP the
C        fraction of the difference (in ln(pressure)) between these
C        two values that the layer pressure lies.
         PLOG = ALOG(PRES(LAY))
         JP(LAY) = INT(36. - 5*(PLOG+0.04))
         IF (JP(LAY) .LT. 1) THEN
            JP(LAY) = 1
         ELSEIF (JP(LAY) .GT. 58) THEN
            JP(LAY) = 58
         ENDIF
         JP1 = JP(LAY) + 1
         FP = 5. * (PREFLOG(JP(LAY)) - PLOG)

C        Determine, for each reference pressure (JP and JP1), which
C        reference temperature (these are different for each  
C        reference pressure) is nearest the layer temperature but does
C        not exceed it.  Store these indices in JT and JT1, resp.
C        Store in FT (resp. FT1) the fraction of the way between JT
C        (JT1) and the next highest reference temperature that the 
C        layer temperature falls.
         JT(LAY) = INT(3. + (TEMP(LAY)-TREF(JP(LAY)))/15.)
         IF (JT(LAY) .LT. 1) THEN
            JT(LAY) = 1
         ELSEIF (JT(LAY) .GT. 4) THEN
            JT(LAY) = 4
         ENDIF
         FT = ((TEMP(LAY)-TREF(JP(LAY)))/15.) - FLOAT(JT(LAY)-3)
         JT1(LAY) = INT(3. + (TEMP(LAY)-TREF(JP1))/15.)
         IF (JT1(LAY) .LT. 1) THEN
            JT1(LAY) = 1
         ELSEIF (JT1(LAY) .GT. 4) THEN
            JT1(LAY) = 4
         ENDIF
         FT1 = ((TEMP(LAY)-TREF(JP1))/15.) - FLOAT(JT1(LAY)-3)

         WATER = WKL(1,LAY)/COLDRY(LAY)
         SCALEFAC = PRES(LAY) * STPFAC / TEMP(LAY)

C        If the pressure is less than ~100mb, perform a different
C        set of species interpolations.
         IF (PLOG .LE. 4.56) GO TO 5300
         LAYTROP =  LAYTROP + 1
C        For one band, the "switch" occurs at ~300 mb. 
         IF (PLOG .GE. 5.76) LAYSWTCH = LAYSWTCH + 1
         IF (PLOG .GE. 6.62) LAYLOW = LAYLOW + 1

         FORFAC(LAY) = SCALEFAC / (1.+WATER)
C        Set up factors needed to separately include the water vapor
C        self-continuum in the calculation of absorption coefficient.
         SELFFAC(LAY) = WATER * FORFAC(LAY)
         FACTOR = (TEMP(LAY)-188.0)/7.2
         INDSELF(LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
         SELFFRAC(LAY) = FACTOR - FLOAT(INDSELF(LAY) + 7)

C        Calculate needed column amounts.
         COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
         COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
         COLO3(LAY) = 1.E-20 * WKL(3,LAY)
         COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
         COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
         IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
C        Using E = 1334.2 cm-1.
         CO2REG = 3.55E-24 * COLDRY(LAY)
         CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *
     &        272.63*EXP(-1919.4/TEMP(LAY))/(8.7604E-4*TEMP(LAY))
         GO TO 5400

C        Above LAYTROP.
 5300    CONTINUE

         FORFAC(LAY) = SCALEFAC / (1.+WATER)
C        Calculate needed column amounts.
         COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
         COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
         COLO3(LAY) = 1.E-20 * WKL(3,LAY)
         COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
         COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
         IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
         CO2REG = 3.55E-24 * COLDRY(LAY)
         CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *
     &        272.63*EXP(-1919.4/TEMP(LAY))/(8.7604E-4*TEMP(LAY))
 5400    CONTINUE

C        We have now isolated the layer ln pressure and temperature,
C        between two reference pressures and two reference temperatures 
C        (for each reference pressure).  We multiply the pressure 
C        fraction FP with the appropriate temperature fractions to get 
C        the factors that will be needed for the interpolation that yields
C        the optical depths (performed in routines TAUGBn for band n).

         COMPFP = 1. - FP
         FAC10(LAY) = COMPFP * FT
         FAC00(LAY) = COMPFP * (1. - FT)
         FAC11(LAY) = FP * FT1
         FAC01(LAY) = FP * (1. - FT1)

 7000 CONTINUE
      RETURN
      END


C ************************  SUBROUTINES TAUGBN  ******************************
C
C     TAUGBn calculates the optical depths and Planck fractions
C     per g-value and layer for band n.
C
C  Output:  absorption coefficients (km^-1) and subinterval weights
C
C     COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
C     COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)
C
C  Input
C
C     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
C     COMMON /PRECISE/  ONEMINUS
C     COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
C    &                  PRES(MAXLEV),TEMP(MAXLEV),
C    &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
C     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
C    &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
C    &                  COLO3(MAXLEV),COLN2O(MAXLEV),
C    &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
C     COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
C    &                  FAC10(MAXLEV),FAC11(MAXLEV)
C     COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
C     COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
C
C  Description
C
C     NG(IBAND) - number of g-values in band IBAND
C     NSPA(IBAND) - for the lower atmosphere, the number of reference
C                   atmospheres that are stored for band IBAND per
C                   pressure level and temperature.  Each of these
C                   atmospheres has different relative amounts of the
C                   key species for the band (i.e. different binary
C                   species parameters).
C     NSPB(IBAND) - same for upper atmosphere
C     ONEMINUS - since problems are caused in some cases by interpolation
C                parameters equal to or greater than 1, for these cases
C                these parameters are set to this value, slightly < 1.
C     PRES - level pressures (mb)
C     TEMP - level temperatures (degrees K)
C     LAYTROP - layer at which switch is made from one combination of
C               key species to another
C     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - concentrations of water
C               vapor,carbon dioxide, ozone, nitrous ozide, methane,
C               respectively (molecules/cm^2/km)
C     CO2MULT - for bands in which carbon dioxide is implemented as a
C               trace species, this is the factor used to multiply the
C               band's average CO2 absorption coefficient to get the added
C               contribution to the optical depth relative to 355 ppm.
C     FACij(LAY) - for layer LAY, these are factors that are needed to
C                  compute the interpolation factors that multiply the
C                  appropriate reference k-values.  A value of 0 (1) for
C                  i,j indicates that the corresponding factor multiplies
C                  reference k-value for the lower (higher) of the two
C                  appropriate temperatures, and altitudes, respectively.
C     JP - the index of the lower (in altitude) of the two appropriate
C          reference pressure levels needed for interpolation
C     JT, JT1 - the indices of the lower of the two appropriate reference
C               temperatures needed for interpolation (for pressure
C               levels JP and JP+1, respectively)
C     SELFFAC - scale factor needed to water vapor self-continuum, equals
C               (water vapor density)/(atmospheric density at 296K and
C               1013 mb)
C     SELFFRAC - factor needed for temperature interpolation of reference
C                water vapor self-continuum data
C     INDSELF - index of the lower of the two appropriate reference
C               temperatures needed for the self-continuum interpolation
C
C  Data input
C
C     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG)
C        (note:  n is the band number)
C
C  Description
C
C     KA - k-values for low reference atmospheres (no water vapor
C          self-continuum) (units: cm^2/molecule)
C     KB - k-values for high reference atmospheres (all sources)
C          (units: cm^2/molecule)
C     SELFREF - k-values for water vapor self-continuum for reference
C               atmospheres (used below LAYTROP)
C               (units: cm^2/molecule)
C
C     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)
C     EQUIVALENCE (KA,ABSA),(KB,ABSB)
C
C ****************************************************************************

      SUBROUTINE TAUGB1

C     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /FOREIGN/  FORFAC(MAXLEV)
      COMMON /K1/       KA(5,13,MG), KB(5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(65,MG),ABSB(235,MG),FORREF(MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

      DATA FRACREFA/
     &    0.08452097,0.17952873,0.16214369,0.13602182,
     &    0.12760490,0.10302561,0.08392423,0.06337652,
     &    0.04206551,0.00487497,0.00410743,0.00344421,
     &    0.00285731,0.00157327,0.00080648,0.00012406/
      DATA FRACREFB/
     &    0.15492001,0.17384727,0.15165100,0.12675308,
     &    0.10986247,0.09006091,0.07584465,0.05990077,
     &    0.04113461,0.00438638,0.00374754,0.00313924,
     &    0.00234381,0.00167167,0.00062744,0.00010889/
      DATA FORREF/
     &     -4.50470E-02,-1.18908E-01,-7.21730E-02,-2.83862E-02,
     &     -3.01961E-02,-1.56877E-02,-1.53684E-02,-1.29135E-02,
     &     -1.27963E-02,-1.81742E-03, 4.40008E-05, 1.05260E-02,
     &      2.17290E-02, 1.65571E-02, 7.60751E-02, 1.47405E-01/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

      WRITE (*,*) '  Computing absorption for band 1'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(1) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(1) + 1
         INDS = INDSELF(LAY)
         DO IG = 1, NG(1)
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG) +
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) +
     &           FORFAC(LAY) * FORREF(IG))
            FRACS(1,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(1) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(1) + 1
         DO IG = 1, NG(1)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG) + 
     &           FORFAC(LAY) * FORREF(IG))
            FRACS(1,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB2

C     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /SPECIES/  COLDRY(MAXLEV),WKL(35,MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /FOREIGN/  FORFAC(MAXLEV)
      COMMON /K2/       KA(5,13,MG), KB(5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(65,MG),ABSB(235,MG),FORREF(MG)
      DIMENSION FC00(MAXLEV),FC01(MAXLEV),FC10(MAXLEV),FC11(MAXLEV)
      DIMENSION FRACREFA(MG,13),FRACREFB(MG),REFPARAM(13)

C     These are the mixing ratios for H2O for a MLS atmosphere at the 
C     13 RRTM reference pressure levels:  1.8759999E-02, 1.2223309E-02, 
C     5.8908667E-03, 2.7675382E-03, 1.4065107E-03, 7.5969833E-04, 
C     3.8875898E-04, 1.6542293E-04, 3.7189537E-05, 7.4764857E-06, 
C     4.3081886E-06, 3.3319423E-06, 3.2039343E-06/        

C     The following are parameters related to the reference water vapor
C     mixing ratios by REFPARAM(I) = REFH2O(I) / (.002+REFH2O(I)).
C     These parameters are used for the Planck function interpolation.
      DATA REFPARAM/
     &  0.903661, 0.859386, 0.746542, 0.580496, 0.412889, 0.275283, 
     &  0.162745, 7.63929E-02, 1.82553E-02, 3.72432E-03, 
     &  2.14946E-03, 1.66320E-03, 1.59940E-03/    

C     The ith set of reference fractions are from the ith reference
C     pressure level.
      DATA FRACREFA/
     &    0.18068060,0.16803175,0.15140158,0.12221480,
     &    0.10240850,0.09330297,0.07518960,0.05611294,
     &    0.03781487,0.00387192,0.00321285,0.00244440,
     &    0.00179546,0.00107704,0.00038798,0.00005060,
     &    0.17927621,0.16731168,0.15129538,0.12328085,
     &    0.10243484,0.09354796,0.07538418,0.05633071,
     &    0.03810832,0.00398347,0.00320262,0.00250029,
     &    0.00178666,0.00111127,0.00039438,0.00005169,
     &    0.17762886,0.16638555,0.15115446,0.12470623,
     &    0.10253213,0.09383459,0.07560240,0.05646568,
     &    0.03844077,0.00409142,0.00322521,0.00254918,
     &    0.00179296,0.00113652,0.00040169,0.00005259,
     &    0.17566043,0.16539773,0.15092199,0.12571971,
     &    0.10340609,0.09426189,0.07559051,0.05678188,
     &    0.03881499,0.00414102,0.00328551,0.00258795,
     &    0.00181648,0.00115145,0.00040969,0.00005357,
     &    0.17335825,0.16442548,0.15070701,0.12667464,
     &    0.10452303,0.09450833,0.07599410,0.05706393,
     &    0.03910370,0.00417880,0.00335256,0.00261708,
     &    0.00185491,0.00116627,0.00041759,0.00005464,
     &    0.17082544,0.16321516,0.15044247,0.12797612,
     &    0.10574646,0.09470057,0.07647423,0.05738756,
     &    0.03935621,0.00423789,0.00342651,0.00264549,
     &    0.00190188,0.00118281,0.00042592,0.00005583,
     &    0.16809277,0.16193336,0.15013184,0.12937409,
     &    0.10720784,0.09485368,0.07692636,0.05771774,
     &    0.03966988,0.00427754,0.00349696,0.00268946,
     &    0.00193536,0.00120222,0.00043462,0.00005712,
     &    0.16517997,0.16059248,0.14984852,0.13079269,
     &    0.10865030,0.09492947,0.07759736,0.05812201,
     &    0.03997169,0.00432356,0.00355308,0.00274031,
     &    0.00197243,0.00122401,0.00044359,0.00005849,
     &    0.16209179,0.15912023,0.14938223,0.13198245,
     &    0.11077233,0.09487948,0.07831636,0.05863440,
     &    0.04028239,0.00436804,0.00360407,0.00279885,
     &    0.00200364,0.00124861,0.00045521,0.00005996,
     &    0.15962425,0.15789343,0.14898103,0.13275230,
     &    0.11253940,0.09503502,0.07884382,0.05908009,
     &    0.04053524,0.00439971,0.00364269,0.00284965,
     &    0.00202758,0.00127076,0.00046408,0.00006114,
     &    0.15926200,0.15770932,0.14891729,0.13283882,
     &    0.11276010,0.09507311,0.07892222,0.05919230,
     &    0.04054824,0.00440833,0.00365575,0.00286459,
     &    0.00203786,0.00128405,0.00046504,0.00006146,
     &    0.15926351,0.15770483,0.14891177,0.13279966,
     &    0.11268171,0.09515216,0.07890341,0.05924807,
     &    0.04052851,0.00440870,0.00365425,0.00286878,
     &    0.00205747,0.00128916,0.00046589,0.00006221,
     &    0.15937765,0.15775780,0.14892603,0.13273248,
     &    0.11252731,0.09521657,0.07885858,0.05927679,
     &    0.04050184,0.00440285,0.00365748,0.00286791,
     &    0.00207507,0.00129193,0.00046679,0.00006308/
C     From P = 0.432 mb.
      DATA FRACREFB/
     &    0.17444289,0.16467269,0.15021490,0.12460902,
     &    0.10400643,0.09481928,0.07590704,0.05752856,
     &    0.03931715,0.00428572,0.00349352,0.00278938,
     &    0.00203448,0.00130037,0.00051560,0.00006255/
      DATA FORREF/
     &     -2.34550E-03,-8.42698E-03,-2.01816E-02,-5.66701E-02,
     &     -8.93189E-02,-6.37487E-02,-4.56455E-02,-4.41417E-02,
     &     -4.48605E-02,-4.74696E-02,-5.16648E-02,-5.63099E-02,
     &     -4.74781E-02,-3.84704E-02,-2.49905E-02, 2.02114E-03/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

      WRITE (*,*) '  Computing absorption for band 2'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum is 
C     interpolated (in temperature) separately.
      DO LAY = 1, LAYTROP
         WATER = 1.E20 * COLH2O(LAY) / COLDRY(LAY)
         H2OPARAM = WATER/(WATER +.002)
         DO IFRAC = 2, 12
            IF (H2OPARAM .GE. REFPARAM(IFRAC)) GO TO 1900
         ENDDO
 1900    CONTINUE

         FP = FAC11(LAY) + FAC01(LAY)
         IF (FP .EQ. 1. .OR. FP .LE. 0. ) THEN
            CORR1 = 1
            CORR2 = 1
         ELSE
            RTFP = SQRT(MAX(0.0,FP))
            CORR1 = RTFP/FP
            CORR2 = (1.-RTFP)/(1.-FP)
         ENDIF
         FC00(LAY) = FAC00(LAY) * CORR2 
         FC10(LAY) = FAC10(LAY) * CORR2 
         FC01(LAY) = FAC01(LAY) * CORR1 
         FC11(LAY) = FAC11(LAY) * CORR1 
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(2) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(2) + 1
         INDS = INDSELF(LAY)
         DO IG = 1, NG(2)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FC00(LAY) * ABSA(IND0,IG) +
     &           FC10(LAY) * ABSA(IND0+1,IG) +
     &           FC01(LAY) * ABSA(IND1,IG) + 
     &           FC11(LAY) * ABSA(IND1+1,IG) + 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) +
     &           FORFAC(LAY) * FORREF(IG))
            FRACS(2,LAY,IG) = FRACREFA(IG,IFRAC) + FRACINT *
     &           (FRACREFA(IG,IFRAC-1)-FRACREFA(IG,IFRAC))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         FP = FAC11(LAY) + FAC01(LAY)
         RTFP = SQRT(MAX(FP,0.0))
         CORR1 = RTFP/FP
         CORR2 = (1.-RTFP)/(1.-FP)
         FC00(LAY) = FAC00(LAY) * CORR2
         FC10(LAY) = FAC10(LAY) * CORR2
         FC01(LAY) = FAC01(LAY) * CORR1
         FC11(LAY) = FAC11(LAY) * CORR1
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(2) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(2) + 1
         DO IG = 1, NG(2)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FC00(LAY) * ABSB(IND0,IG) +
     &           FC10(LAY) * ABSB(IND0+1,IG) +
     &           FC01(LAY) * ABSB(IND1,IG) + 
     &           FC11(LAY) * ABSB(IND1+1,IG) +
     &           FORFAC(LAY) * FORREF(IG))
            FRACS(2,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB3

C     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /FOREIGN/  FORFAC(MAXLEV)
      COMMON /K3/       KA(10,5,13,MG), KB(5,5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(650,MG),ABSB(1175,MG),FORREF(MG)
      DIMENSION ABSN2OA(MG),ABSN2OB(MG)
      DIMENSION FRACREFA(MG,10),FRACREFB(MG,5)
      DIMENSION N2OREF(59),H2OREF(59),CO2REF(59), ETAREF(10)

C     From P = 1053.6 mb.
      DATA FRACREFA/
     &    0.15116400,0.14875700,0.14232300,0.13234501,
     &    0.11881600,0.10224100,0.08345580,0.06267490,
     &    0.04250650,0.00462650,0.00382259,0.00302600,
     &    0.00222004,0.00141397,0.00053379,0.00007421,
     &    0.15266000,0.14888400,0.14195900,0.13179500,
     &    0.11842700,0.10209000,0.08336130,0.06264370,
     &    0.04247660,0.00461946,0.00381536,0.00302601,
     &    0.00222004,0.00141397,0.00053302,0.00007498,
     &    0.15282799,0.14903000,0.14192399,0.13174300,
     &    0.11835300,0.10202700,0.08329830,0.06264830,
     &    0.04246910,0.00460242,0.00381904,0.00301573,
     &    0.00222004,0.00141397,0.00053379,0.00007421,
     &    0.15298399,0.14902800,0.14193401,0.13173500,
     &    0.11833300,0.10195800,0.08324730,0.06264770,
     &    0.04246490,0.00460489,0.00381123,0.00301893,
     &    0.00221093,0.00141397,0.00053379,0.00007421,
     &    0.15307599,0.14907201,0.14198899,0.13169800,
     &    0.11827300,0.10192300,0.08321600,0.06263490,
     &    0.04245600,0.00460846,0.00380836,0.00301663,
     &    0.00221402,0.00141167,0.00052807,0.00007376,
     &    0.15311401,0.14915401,0.14207301,0.13167299,
     &    0.11819300,0.10188900,0.08318760,0.06261960,
     &    0.04243890,0.00461584,0.00380929,0.00300815,
     &    0.00221736,0.00140588,0.00052776,0.00007376,
     &    0.15316001,0.14925499,0.14213000,0.13170999,
     &    0.11807700,0.10181400,0.08317400,0.06260300,
     &    0.04242720,0.00461520,0.00381381,0.00301285,
     &    0.00220275,0.00140371,0.00052776,0.00007376,
     &    0.15321200,0.14940999,0.14222500,0.13164200,
     &    0.11798200,0.10174500,0.08317500,0.06253640,
     &    0.04243130,0.00461724,0.00381534,0.00300320,
     &    0.00220091,0.00140364,0.00052852,0.00007300,
     &    0.15312800,0.14973100,0.14234400,0.13168900,
     &    0.11795200,0.10156100,0.08302990,0.06252240,
     &    0.04240980,0.00461035,0.00381381,0.00300176,
     &    0.00220160,0.00140284,0.00052774,0.00007376,
     &    0.15292500,0.14978001,0.14242400,0.13172600,
     &    0.11798800,0.10156400,0.08303050,0.06251670,
     &    0.04240970,0.00461302,0.00381452,0.00300250,
     &    0.00220126,0.00140324,0.00052850,0.00007300/
C     From P = 64.1 mb.
      DATA FRACREFB/
     &    0.16340201,0.15607700,0.14601400,0.13182700,
     &    0.11524700,0.09666570,0.07825360,0.05849780,
     &    0.03949650,0.00427980,0.00353719,0.00279303,
     &    0.00204788,0.00130139,0.00049055,0.00006904,
     &    0.15762900,0.15494700,0.14659800,0.13267800,
     &    0.11562700,0.09838360,0.07930420,0.05962700,
     &    0.04036360,0.00438053,0.00361463,0.00285723,
     &    0.00208345,0.00132135,0.00050528,0.00008003,
     &    0.15641500,0.15394500,0.14633600,0.13180400,
     &    0.11617100,0.09924170,0.08000510,0.06021420,
     &    0.04082730,0.00441694,0.00365364,0.00287723,
     &    0.00210914,0.00135784,0.00054651,0.00008003,
     &    0.15482700,0.15286300,0.14392500,0.13244100,
     &    0.11712000,0.09994920,0.08119200,0.06104360,
     &    0.04135600,0.00446685,0.00368377,0.00290767,
     &    0.00215445,0.00142865,0.00056142,0.00008003,
     &    0.15975100,0.15653500,0.14214399,0.12892200,
     &    0.11508400,0.09906020,0.08087940,0.06078190,
     &    0.04140530,0.00452724,0.00374558,0.00295328,
     &    0.00218509,0.00138644,0.00056018,0.00008003/
      DATA ABSN2OA/
     &     1.50387E-01,2.91407E-01,6.28803E-01,9.65619E-01,
     &     1.15054E-00,2.23424E-00,1.83392E-00,1.39033E-00,
     &     4.28457E-01,2.73502E-01,1.84307E-01,1.61325E-01,
     &     7.66314E-02,1.33862E-01,6.71196E-07,1.59293E-06/
      DATA ABSN2OB/
     &     9.37044E-05,1.23318E-03,7.91720E-03,5.33005E-02,
     &     1.72343E-01,4.29571E-01,1.01288E+00,3.83863E+00,
     &     1.15312E+01,1.08383E+00,2.24847E+00,1.51268E+00,
     &     3.33177E-01,7.82102E-01,3.44631E-01,1.61039E-03/
      DATA ETAREF/
     &     0.,0.125,0.25,0.375,0.5,0.625,0.75,0.875,0.9875,1.0/
      DATA H2OREF/
     &     1.87599E-02,1.22233E-02,5.89086E-03,2.76753E-03,1.40651E-03, 
     &     7.59698E-04,3.88758E-04,1.65422E-04,3.71895E-05,7.47648E-06, 
     &     4.30818E-06,3.33194E-06,3.20393E-06,3.16186E-06,3.25235E-06, 
     &     3.42258E-06,3.62884E-06,3.91482E-06,4.14875E-06,4.30810E-06,
     &     4.44204E-06,4.57783E-06,4.70865E-06,4.79432E-06,4.86971E-06, 
     &     4.92603E-06,4.96688E-06,4.99628E-06,5.05266E-06,5.12658E-06, 
     &     5.25028E-06,5.35708E-06,5.45085E-06,5.48304E-06,5.50000E-06, 
     &     5.50000E-06,5.45359E-06,5.40468E-06,5.35576E-06,5.25327E-06,
     &     5.14362E-06,5.03396E-06,4.87662E-06,4.69787E-06,4.51911E-06, 
     &     4.33600E-06,4.14416E-06,3.95232E-06,3.76048E-06,3.57217E-06, 
     &     3.38549E-06,3.19881E-06,3.01212E-06,2.82621E-06,2.64068E-06, 
     &     2.45515E-06,2.26962E-06,2.08659E-06,1.93029E-06/
      DATA N2OREF/
     &     3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,
     &     3.19652E-07,3.15324E-07,3.03830E-07,2.94221E-07,2.84953E-07,
     &     2.76714E-07,2.64709E-07,2.42847E-07,2.09547E-07,1.71945E-07,
     &     1.37491E-07,1.13319E-07,1.00354E-07,9.12812E-08,8.54633E-08,
     &     8.03631E-08,7.33718E-08,6.59754E-08,5.60386E-08,4.70901E-08,
     &     3.99774E-08,3.29786E-08,2.60642E-08,2.10663E-08,1.65918E-08,
     &     1.30167E-08,1.00900E-08,7.62490E-09,6.11592E-09,4.66725E-09,
     &     3.28574E-09,2.84838E-09,2.46198E-09,2.07557E-09,1.85507E-09,
     &     1.65675E-09,1.45843E-09,1.31948E-09,1.20716E-09,1.09485E-09,
     &     9.97803E-10,9.31260E-10,8.64721E-10,7.98181E-10,7.51380E-10,
     &     7.13670E-10,6.75960E-10,6.38250E-10,6.09811E-10,5.85998E-10,
     &     5.62185E-10,5.38371E-10,5.15183E-10,4.98660E-10/
      DATA CO2REF/
     &     53*3.55E-04, 3.5470873E-04, 3.5427220E-04, 3.5383567E-04,
     &     3.5339911E-04, 3.5282588E-04, 3.5079606E-04/  
      DATA FORREF/
     &      1.76842E-04, 1.77913E-04, 1.25186E-04, 1.07912E-04,
     &      1.05217E-04, 7.48726E-05, 1.11701E-04, 7.68921E-05,
     &      9.87242E-05, 9.85711E-05, 6.16557E-05,-1.61291E-05,
     &     -1.26794E-04,-1.19011E-04,-2.67814E-04, 6.95005E-05/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB,N2OMULT
      STRRAT = 1.19268

      WRITE (*,*) '  Computing absorption for band 3'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         IF (JS .EQ. 8) THEN
            IF (FS .GE. 0.9) THEN
               JS = 9
               FS = 10. * (FS - 0.9)
            ELSE
               FS = FS/0.9
            ENDIF
         ENDIF
         NS = JS + INT(FS + 0.5)
         FP = FAC01(LAY) + FAC11(LAY)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(3) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(3) + JS
         INDS = INDSELF(LAY)
         COLREF1 = N2OREF(JP(LAY))
         COLREF2 = N2OREF(JP(LAY)+1)
         IF (NS .EQ. 10) THEN
            WCOMB1 = H2OREF(JP(LAY))
            WCOMB2 = H2OREF(JP(LAY)+1)
         ELSE
            WCOMB1 = STRRAT * CO2REF(JP(LAY))/(1.-ETAREF(NS))
            WCOMB2 = STRRAT * CO2REF(JP(LAY)+1)/(1.-ETAREF(NS))
         ENDIF
         RATIO = (COLREF1/WCOMB1)+FP*((COLREF2/WCOMB2)-(COLREF1/WCOMB1))
         CURRN2O = SPECCOMB * RATIO
         N2OMULT = COLN2O(LAY) - CURRN2O 
         DO IG = 1, NG(3)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+10,IG) +
     &           FAC110 * ABSA(IND0+11,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+10,IG) +
     &           FAC111 * ABSA(IND1+11,IG)) +
     &           COLH2O(LAY) * 
     &           (SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))) +
     &           FORFAC(LAY) * FORREF(IG))
     &           + N2OMULT * ABSN2OA(IG)
            FRACS(3,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         SPECCOMB = COLH2O(LAY) + STRRAT*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         NS = JS + INT(FS + 0.5)
         FP = FAC01(LAY) + FAC11(LAY)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(3) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(3) + JS
         COLREF1 = N2OREF(JP(LAY)) 
         COLREF2 = N2OREF(JP(LAY)+1) 
         IF (NS .EQ. 5) THEN
            WCOMB1 = H2OREF(JP(LAY))
            WCOMB2 = H2OREF(JP(LAY)+1)
         ELSE
            WCOMB1 = STRRAT * CO2REF(JP(LAY))/(1.-ETAREF(NS))
            WCOMB2 = STRRAT * CO2REF(JP(LAY)+1)/(1.-ETAREF(NS))
         ENDIF
         RATIO = (COLREF1/WCOMB1)+FP*((COLREF2/WCOMB2)-(COLREF1/WCOMB1))
         CURRN2O = SPECCOMB * RATIO
         N2OMULT = COLN2O(LAY) - CURRN2O 
         DO IG = 1, NG(3)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &           FAC100 * ABSB(IND0+1,IG) +
     &           FAC010 * ABSB(IND0+5,IG) +
     &           FAC110 * ABSB(IND0+6,IG) +
     &           FAC001 * ABSB(IND1,IG) + 
     &           FAC101 * ABSB(IND1+1,IG) +
     &           FAC011 * ABSB(IND1+5,IG) +
     &           FAC111 * ABSB(IND1+6,IG)) +
     &           COLH2O(LAY) * FORFAC(LAY) * FORREF(IG) 
     &           + N2OMULT * ABSN2OB(IG)
            FRACS(3,LAY,IG) = FRACREFB(IG,JS) + FS *
     &           (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB4

C     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K4/       KA(9,5,13,MG), KB(6,5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(585,MG),ABSB(1410,MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG,6)

      DATA FRACREFA/
     &    0.15579100,0.14918099,0.14113800,0.13127001,
     &    0.11796300,0.10174300,0.08282370,0.06238150,
     &    0.04213440,0.00458968,0.00377949,0.00298736,
     &    0.00220743,0.00140644,0.00053024,0.00007459,
     &    0.15292799,0.15004000,0.14211500,0.13176700,
     &    0.11821100,0.10186300,0.08288040,0.06241390,
     &    0.04220720,0.00459006,0.00377919,0.00298743,
     &    0.00220743,0.00140644,0.00053024,0.00007459,
     &    0.14386199,0.15125300,0.14650001,0.13377000,
     &    0.11895900,0.10229400,0.08312110,0.06239520,
     &    0.04225560,0.00459428,0.00378865,0.00298860,
     &    0.00220743,0.00140644,0.00053024,0.00007459,
     &    0.14359100,0.14561599,0.14479300,0.13740200,
     &    0.12150100,0.10315400,0.08355480,0.06247240,
     &    0.04230980,0.00459916,0.00378373,0.00300063,
     &    0.00221111,0.00140644,0.00053024,0.00007459,
     &    0.14337599,0.14451601,0.14238000,0.13520500,
     &    0.12354200,0.10581200,0.08451810,0.06262440,
     &    0.04239590,0.00460297,0.00378701,0.00300466,
     &    0.00221899,0.00141020,0.00053024,0.00007459,
     &    0.14322001,0.14397401,0.14117201,0.13401900,
     &    0.12255500,0.10774100,0.08617650,0.06296420,
     &    0.04249590,0.00463406,0.00378241,0.00302037,
     &    0.00221583,0.00141103,0.00053814,0.00007991,
     &    0.14309500,0.14364301,0.14043900,0.13348100,
     &    0.12211600,0.10684700,0.08820590,0.06374610,
     &    0.04264730,0.00464231,0.00384022,0.00303427,
     &    0.00221825,0.00140943,0.00055564,0.00007991,
     &    0.15579100,0.14918099,0.14113800,0.13127001,
     &    0.11796300,0.10174300,0.08282370,0.06238150,
     &    0.04213440,0.00458968,0.00377949,0.00298736,
     &    0.00220743,0.00140644,0.00053024,0.00007459,
     &    0.15937001,0.15159500,0.14242800,0.13078900,
     &    0.11671300,0.10035700,0.08143450,0.06093850,
     &    0.04105320,0.00446233,0.00369844,0.00293784,
     &    0.00216425,0.00143403,0.00054571,0.00007991/
C     From P = 1.17 mb.
      DATA FRACREFB/
     &    0.15558299,0.14930600,0.14104301,0.13124099,
     &    0.11792900,0.10159200,0.08314130,0.06240450,
     &    0.04217020,0.00459313,0.00379798,0.00299835,
     &    0.00218950,0.00140615,0.00053010,0.00007457,
     &    0.15592700,0.14918999,0.14095700,0.13115700,
     &    0.11788900,0.10158000,0.08313780,0.06240240,
     &    0.04217000,0.00459313,0.00379798,0.00299835,
     &    0.00218950,0.00140615,0.00053010,0.00007457,
     &    0.15949000,0.15014900,0.14162201,0.13080800,
     &    0.11713500,0.10057100,0.08170080,0.06128110,
     &    0.04165600,0.00459202,0.00379835,0.00299717,
     &    0.00218958,0.00140616,0.00053010,0.00007457,
     &    0.15967900,0.15038200,0.14196999,0.13074800,
     &    0.11701700,0.10053000,0.08160790,0.06122690,
     &    0.04128310,0.00456598,0.00379486,0.00299457,
     &    0.00219016,0.00140619,0.00053011,0.00007456,
     &    0.15989800,0.15057300,0.14207700,0.13068600,
     &    0.11682900,0.10053900,0.08163610,0.06121870,
     &    0.04121690,0.00449061,0.00371235,0.00294207,
     &    0.00217778,0.00139877,0.00053011,0.00007455,
     &    0.15950100,0.15112500,0.14199100,0.13071300,
     &    0.11680800,0.10054600,0.08179050,0.06120910,
     &    0.04126050,0.00444324,0.00366843,0.00289369,
     &    0.00211550,0.00134746,0.00050874,0.00007863/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      STRRAT1 = 850.577
      STRRAT2 = 35.7416

      WRITE (*,*) '  Computing absorption for band 4'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(4) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(4) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(4)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            FRACS(4,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         SPECCOMB = COLO3(LAY) + STRRAT2*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         IF (JS .GT. 1) THEN
            JS = JS + 1
         ELSEIF (FS .GE. 0.0024) THEN
            JS = 2
            FS = (FS - 0.0024)/0.9976
         ELSE
            JS = 1
            FS = FS/0.0024
         ENDIF
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(4) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(4) + JS
         DO IG = 1, NG(4)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &           FAC100 * ABSB(IND0+1,IG) +
     &           FAC010 * ABSB(IND0+6,IG) +
     &           FAC110 * ABSB(IND0+7,IG) +
     &           FAC001 * ABSB(IND1,IG) + 
     &           FAC101 * ABSB(IND1+1,IG) +
     &           FAC011 * ABSB(IND1+6,IG) +
     &           FAC111 * ABSB(IND1+7,IG))
            FRACS(4,LAY,IG) = FRACREFB(IG,JS) + FS *
     &           (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB5

C     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)

      INTEGER MAXLEV, MG, NBANDS, MAXXSEC
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16, MAXXSEC=4)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /XSEC/     WX(MAXXSEC,MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K5/       KA(9,5,13,MG), KB(5,5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(585,MG),ABSB(1175,MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG,5),CCL4(MG)

C     From P = 387.6 mb.
      DATA FRACREFA/
     &    0.13966499,0.14138900,0.13763399,0.13076700,
     &    0.12299100,0.10747700,0.08942000,0.06769200,   
     &    0.04587610,0.00501173,0.00415809,0.00328398,
     &    0.00240015,0.00156222,0.00059104,0.00008323,
     &    0.13958199,0.14332899,0.13785399,0.13205400,
     &    0.12199700,0.10679600,0.08861080,0.06712320,
     &    0.04556030,0.00500863,0.00416315,0.00328629,
     &    0.00240023,0.00156220,0.00059104,0.00008323,
     &    0.13907100,0.14250501,0.13889600,0.13297300,
     &    0.12218700,0.10683800,0.08839260,0.06677310,
     &    0.04538570,0.00495402,0.00409863,0.00328219,
     &    0.00240805,0.00156266,0.00059104,0.00008323,
     &    0.13867700,0.14190100,0.13932300,0.13327099,
     &    0.12280800,0.10692500,0.08844510,0.06658510,
     &    0.04519340,0.00492276,0.00408832,0.00323856,  
     &    0.00239289,0.00155698,0.00059104,0.00008323,
     &    0.13845000,0.14158800,0.13929300,0.13295600,
     &    0.12348300,0.10736700,0.08859480,0.06650610,
     &    0.04498230,0.00491335,0.00406968,0.00322901,
     &    0.00234666,0.00155235,0.00058813,0.00008323,
     &    0.13837101,0.14113200,0.13930500,0.13283101,
     &    0.12349200,0.10796400,0.08890490,0.06646480,
     &    0.04485990,0.00489554,0.00405264,0.00320313,
     &    0.00234742,0.00151159,0.00058438,0.00008253,
     &    0.13834500,0.14093500,0.13896500,0.13262001,
     &    0.12326900,0.10828900,0.08950050,0.06674610,
     &    0.04476560,0.00489624,0.00400962,0.00317423,
     &    0.00233479,0.00148249,0.00058590,0.00008253,
     &    0.13831300,0.14069000,0.13871400,0.13247600,
     &    0.12251400,0.10831300,0.08977090,0.06776920,
     &    0.04498390,0.00484111,0.00398948,0.00316069,
     &    0.00229741,0.00150104,0.00058608,0.00008253,
     &    0.14027201,0.14420401,0.14215700,0.13446601,
     &    0.12303700,0.10596100,0.08650370,0.06409570,
     &    0.04312310,0.00471110,0.00393954,0.00310850,
     &    0.00229588,0.00146366,0.00058194,0.00008253/
C     From P = 1.17 mb.
      DATA FRACREFB/
     &    0.14339100,0.14358699,0.13935301,0.13306700,
     &    0.12135700,0.10590600,0.08688240,0.06553220,
     &    0.04446740,0.00483580,0.00399413,0.00316225,
     &    0.00233007,0.00149135,0.00056246,0.00008059,
     &    0.14330500,0.14430299,0.14053699,0.13355300,
     &    0.12151200,0.10529100,0.08627630,0.06505230,
     &    0.04385850,0.00476555,0.00395010,0.00313878,
     &    0.00232273,0.00149354,0.00056246,0.00008059,
     &    0.14328399,0.14442700,0.14078601,0.13390100,
     &    0.12132600,0.10510600,0.08613660,0.06494630,
     &    0.04381310,0.00475378,0.00394166,0.00313076,
     &    0.00231235,0.00149159,0.00056301,0.00008059,
     &    0.14326900,0.14453100,0.14114200,0.13397101,
     &    0.12127200,0.10493400,0.08601380,0.06483360,
     &    0.04378900,0.00474655,0.00393549,0.00312583,
     &    0.00230686,0.00148433,0.00056502,0.00008059,  
     &    0.14328900,0.14532700,0.14179000,0.13384600,
     &    0.12093700,0.10461500,0.08573010,0.06461340,
     &    0.04366570,0.00473087,0.00392539,0.00311238,
     &    0.00229865,0.00147572,0.00056517,0.00007939/
      DATA CCL4/
     &     26.1407,  53.9776,  63.8085,  36.1701,
     &     15.4099, 10.23116,  4.82948,  5.03836,
     &     1.75558,0.,0.,0.,
     &     0.,0.,0.,0./

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      STRRAT1 = 90.4894
      STRRAT2 = 0.900502

      WRITE (*,*) '  Computing absorption for band 5'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(5) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(5) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(5)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
     &           + WX(1,LAY) * CCL4(IG)
            FRACS(5,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         SPECCOMB = COLO3(LAY) + STRRAT2*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(5) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(5) + JS
         DO IG = 1, NG(5)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &           FAC100 * ABSB(IND0+1,IG) +
     &           FAC010 * ABSB(IND0+5,IG) +
     &           FAC110 * ABSB(IND0+6,IG) +
     &           FAC001 * ABSB(IND1,IG) + 
     &           FAC101 * ABSB(IND1+1,IG) +
     &           FAC011 * ABSB(IND1+5,IG) +
     &           FAC111 * ABSB(IND1+6,IG))
     &           + WX(1,LAY) * CCL4(IG)
            FRACS(5,LAY,IG) = FRACREFB(IG,JS) + FS *
     &           (FRACREFB(IG,JS+1) - FRACREFB(IG,JS))
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB6

C     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

      INTEGER MAXLEV, MG, NBANDS, MAXXSEC
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16, MAXXSEC=4)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /XSEC/     WX(MAXXSEC,MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K6/       KA(5,13,MG), SELFREF(10,MG)

      DIMENSION ABSA(65,MG)
      DIMENSION FRACREFA(MG),ABSCO2(MG), CFC11ADJ(MG), CFC12(MG)

C     From P = 706 mb.
      DATA FRACREFA/
     &    0.13739009,0.14259538,0.14033118,0.13547136,
     &    0.12569460,0.11028396,0.08626066,0.06245148,
     &    0.04309394,0.00473551,0.00403920,0.00321695,
     &    0.00232470,0.00147662,0.00056095,0.00007373/
C      DATA CFC11/
C     &     0., 0., 26.5435, 108.850,
C     &     58.7804, 54.0875, 41.1065, 35.6120,
C     &     41.2328, 47.7402, 79.1026, 64.3005,
C     &     108.206, 141.617, 186.565, 58.4782/
C     CFC11 is multiplied by 1.385 to account for the 1060-1107 cm-1 band.
      DATA CFC11ADJ/
     &     0.,  0., 36.7627,    150.757,    
     &     81.4109, 74.9112, 56.9325, 49.3226,  
     &     57.1074, 66.1202, 109.557, 89.0562,  
     &     149.865, 196.140, 258.393, 80.9923/   
      DATA CFC12/
     &     62.8368, 43.2626, 26.7549, 22.2487,
     &     23.5029, 34.8323, 26.2335, 23.2306,
     &     18.4062, 13.9534, 22.6268, 24.2604,
     &     30.0088, 26.3634, 15.8237, 57.5050/

      DATA ABSCO2/
     &     7.44852E-05, 6.29208E-05, 7.34031E-05, 6.65218E-05,
     &     7.87511E-05, 1.22489E-04, 3.39785E-04, 9.33040E-04,
     &     1.54323E-03, 4.07220E-04, 4.34332E-04, 8.76418E-05,
     &     9.80381E-05, 3.51680E-05, 5.31766E-05, 1.01542E-05/

      EQUIVALENCE (KA,ABSA)
      REAL KA

      WRITE (*,*) '  Computing absorption for band 6'

C     Compute the optical depth by interpolating in ln(pressure) and
C     temperature. The water vapor self-continuum is interpolated
C     (in temperature) separately.  
      DO LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(6) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(6) + 1
         INDS = INDSELF(LAY)
         DO IG = 1, NG(6)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG) + 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY)*
     &           (SELFREF(INDS+1,IG)-SELFREF(INDS,IG))))
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
     &           + CO2MULT(LAY) * ABSCO2(IG)
            FRACS(6,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

C     Nothing important goes on above LAYTROP in this band.
      DO LAY = LAYTROP+1, NLEVELS
         DO IG = 1, NG(6)
            TAUG(LAY,IG) = 0.0 
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
            FRACS(6,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB7

C     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K7/       KA(9,5,13,MG),KB(5,13:59,MG),SELFREF(10,MG)

      DIMENSION ABSA(585,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG),ABSCO2(MG)

      DATA FRACREFA/
     &  0.16461779, 0.14889984, 0.14233345, 0.13156526,
     &  0.11679733, 0.09988949, 0.08078653, 0.06006384,
     &  0.04028391, 0.00435899, 0.00359173, 0.00281707,
     &  0.00206767, 0.00135012, 0.00050720, 0.00007146,
     &  0.16442357, 0.14944240, 0.14245804, 0.13111183,
     &  0.11688625, 0.09983791, 0.08085148, 0.05993948,
     &  0.04028057, 0.00435939, 0.00358708, 0.00284036,
     &  0.00208869, 0.00133256, 0.00049260, 0.00006931,
     &  0.16368519, 0.15018989, 0.14262174, 0.13084342,
     &  0.11682195, 0.09996257, 0.08074036, 0.05985692,
     &  0.04045362, 0.00436208, 0.00358257, 0.00287122,
     &  0.00211004, 0.00133804, 0.00049260, 0.00006931,
     &  0.16274056, 0.15133780, 0.14228874, 0.13081114,
     &  0.11688486, 0.09979610, 0.08073687, 0.05996741,
     &  0.04040616, 0.00439869, 0.00368910, 0.00293041,
     &  0.00211604, 0.00133536, 0.00049260, 0.00006931,
     &  0.16176532, 0.15207882, 0.14226955, 0.13079646,
     &  0.11688191, 0.09966998, 0.08066384, 0.06020275,
     &  0.04047901, 0.00446696, 0.00377456, 0.00294410,
     &  0.00211082, 0.00133536, 0.00049260, 0.00006931,
     &  0.15993737, 0.15305527, 0.14259829, 0.13078023,
     &  0.11686983, 0.09980131, 0.08058286, 0.06031430,
     &  0.04082833, 0.00450509, 0.00377574, 0.00294823,
     &  0.00210977, 0.00133302, 0.00049260, 0.00006931,
     &  0.15371189, 0.15592396, 0.14430280, 0.13076764,
     &  0.11720382, 0.10023471, 0.08066396, 0.06073554,
     &  0.04121581, 0.00451202, 0.00377832, 0.00294609,
     &  0.00210943, 0.00133336, 0.00049260, 0.00006931,
     &  0.14262275, 0.14572631, 0.14560597, 0.13736825,
     &  0.12271351, 0.10419556, 0.08294533, 0.06199794,
     &  0.04157615, 0.00452842, 0.00377704, 0.00293852,
     &  0.00211034, 0.00133278, 0.00049259, 0.00006931,
     &  0.14500433, 0.14590444, 0.14430299, 0.13770708,
     &  0.12288283, 0.10350952, 0.08269450, 0.06130579,
     &  0.04144571, 0.00452096, 0.00377382, 0.00294532,
     &  0.00210943, 0.00133228, 0.00049260, 0.00006931/
      DATA FRACREFB/
     &    0.15355594,0.15310939,0.14274909,0.13129812,
     &    0.11736792,0.10118213,0.08215259,0.06165591,
     &    0.04164486,0.00451141,0.00372837,0.00294095,
     &    0.00215259,0.00136792,0.00051233,0.00007075/
      DATA ABSCO2/
     &     9.30038E-05, 1.74061E-04, 2.09293E-04, 2.52360E-04,
     &     3.13404E-04, 4.16619E-04, 6.27394E-04, 1.29386E-03,
     &     4.05192E-03, 3.97050E-03, 7.00634E-04, 6.06617E-04,
     &     7.66978E-04, 6.70661E-04, 7.89971E-04, 7.55709E-04/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      STRRAT1 = 8.21104E4

      WRITE (*,*) '  Computing absorption for band 7'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLO3(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(7) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(7) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(7)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
     &           + CO2MULT(LAY) * ABSCO2(IG)
            FRACS(7,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(7) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(7) + 1
         DO IG = 1, NG(7)
            TAUG(LAY,IG) = COLO3(LAY) *
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + CO2MULT(LAY) * ABSCO2(IG)
            FRACS(7,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB8

C     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)

      INTEGER MAXLEV, MG, NBANDS, MAXXSEC
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16, MAXXSEC=4)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /XSEC/     WX(MAXXSEC,MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K8/       KA(5,7,MG), KB(5,7:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(35,MG),ABSB(265,MG),CFC12(MG),CFC22ADJ(MG)
      DIMENSION ABSN2OA(MG),ABSN2OB(MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG),ABSCO2A(MG),ABSCO2B(MG)
      DIMENSION N2OREF(59),H2OREF(59),O3REF(59)

C     From P = 1053.6 mb.
      DATA FRACREFA/
     &    0.15309700,0.15450300,0.14458799,0.13098200,
     &    0.11817900,0.09953490,0.08132080,0.06139960,
     &    0.04132010,0.00446788,0.00372533,0.00294053,
     &    0.00211371,0.00128122,0.00048050,0.00006759/
C     From P = 28.9 mb.
      DATA FRACREFB/
     &    0.14105400,0.14728899,0.14264800,0.13331699,
     &    0.12034100,0.10467000,0.08574980,0.06469390,
     &    0.04394640,0.00481284,0.00397375,0.00315006,
     &    0.00228636,0.00144606,0.00054604,0.00007697/
      DATA CFC12/
     &     85.4027, 89.4696, 74.0959, 67.7480,
     &     61.2444, 59.9073, 60.8296, 63.0998,
     &     59.6110, 64.0735, 57.2622, 58.9721,
     &     43.5505, 26.1192, 32.7023, 32.8667/
C     Original CFC22 is multiplied by 1.485 to account for the 780-850 cm-1 
C     and 1290-1335 cm-1 bands.
      DATA CFC22ADJ/
     &     135.335, 89.6642, 76.2375, 65.9748,
     &     63.1164, 60.2935, 64.0299, 75.4264,
     &     51.3018, 7.07911, 5.86928, 0.398693,
     &     2.82885, 9.12751, 6.28271, 0./

      DATA ABSCO2A/
     &     1.11233E-05, 3.92400E-05, 6.62059E-05, 8.51687E-05,
     &     7.79035E-05, 1.34058E-04, 2.82553E-04, 5.41741E-04,
     &     1.47029E-05, 2.34982E-05, 6.91094E-08, 8.48917E-08,
     &     6.58783E-08, 4.64849E-08, 3.62742E-08, 3.62742E-08/
      DATA ABSCO2B/
     &     4.10977E-09, 5.65200E-08, 1.70800E-07, 4.16840E-07,
     &     9.53684E-07, 2.36468E-06, 7.29502E-06, 4.93883E-05, 
     &     5.10440E-04, 9.75248E-04, 1.36495E-03, 2.40451E-03,
     &     4.50277E-03, 2.24486E-02, 4.06756E-02, 2.17447E-10/
      DATA ABSN2OA/
     &     1.28527E-02,5.28651E-02,1.01668E-01,1.57224E-01,
     &     2.76947E-01,4.93048E-01,6.71387E-01,3.48809E-01,
     &     4.19840E-01,3.13558E-01,2.44432E-01,2.05108E-01,
     &     1.21423E-01,1.22158E-01,1.49702E-01,1.47799E-01/
      DATA ABSN2OB/
     &     3.15864E-03,4.87347E-03,8.63235E-03,2.16053E-02,
     &     3.63699E-02,7.89149E-02,3.53807E-01,1.27140E-00,
     &     2.31464E-00,7.75834E-02,5.15063E-02,4.07059E-02,
     &     5.91947E-02,5.83546E-02,3.12716E-01,1.47456E-01/

      DATA H2OREF/
     &     1.87599E-02,1.22233E-02,5.89086E-03,2.76753E-03,1.40651E-03, 
     &     7.59698E-04,3.88758E-04,1.65422E-04,3.71895E-05,7.47648E-06, 
     &     4.30818E-06,3.33194E-06,3.20393E-06,3.16186E-06,3.25235E-06, 
     &     3.42258E-06,3.62884E-06,3.91482E-06,4.14875E-06,4.30810E-06,
     &     4.44204E-06,4.57783E-06,4.70865E-06,4.79432E-06,4.86971E-06, 
     &     4.92603E-06,4.96688E-06,4.99628E-06,5.05266E-06,5.12658E-06, 
     &     5.25028E-06,5.35708E-06,5.45085E-06,5.48304E-06,5.50000E-06, 
     &     5.50000E-06,5.45359E-06,5.40468E-06,5.35576E-06,5.25327E-06,
     &     5.14362E-06,5.03396E-06,4.87662E-06,4.69787E-06,4.51911E-06, 
     &     4.33600E-06,4.14416E-06,3.95232E-06,3.76048E-06,3.57217E-06, 
     &     3.38549E-06,3.19881E-06,3.01212E-06,2.82621E-06,2.64068E-06, 
     &     2.45515E-06,2.26962E-06,2.08659E-06,1.93029E-06/
      DATA N2OREF/
     &     3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,
     &     3.19652E-07,3.15324E-07,3.03830E-07,2.94221E-07,2.84953E-07,
     &     2.76714E-07,2.64709E-07,2.42847E-07,2.09547E-07,1.71945E-07,
     &     1.37491E-07,1.13319E-07,1.00354E-07,9.12812E-08,8.54633E-08,
     &     8.03631E-08,7.33718E-08,6.59754E-08,5.60386E-08,4.70901E-08,
     &     3.99774E-08,3.29786E-08,2.60642E-08,2.10663E-08,1.65918E-08,
     &     1.30167E-08,1.00900E-08,7.62490E-09,6.11592E-09,4.66725E-09,
     &     3.28574E-09,2.84838E-09,2.46198E-09,2.07557E-09,1.85507E-09,
     &     1.65675E-09,1.45843E-09,1.31948E-09,1.20716E-09,1.09485E-09,
     &     9.97803E-10,9.31260E-10,8.64721E-10,7.98181E-10,7.51380E-10,
     &     7.13670E-10,6.75960E-10,6.38250E-10,6.09811E-10,5.85998E-10,
     &     5.62185E-10,5.38371E-10,5.15183E-10,4.98660E-10/
      DATA O3REF/
     &     3.01700E-08,3.47254E-08,4.24769E-08,5.27592E-08,6.69439E-08,
     &     8.71295E-08,1.13911E-07,1.56771E-07,2.17878E-07,3.24430E-07,
     &     4.65942E-07,5.68057E-07,6.96065E-07,1.11863E-06,1.76175E-06,
     &     2.32689E-06,2.95769E-06,3.65930E-06,4.59503E-06,5.31891E-06,
     &     5.96179E-06,6.51133E-06,7.06350E-06,7.69169E-06,8.25771E-06,
     &     8.70824E-06,8.83245E-06,8.71486E-06,8.09434E-06,7.33071E-06,
     &     6.31014E-06,5.36717E-06,4.48289E-06,3.83913E-06,3.28270E-06,
     &     2.82351E-06,2.49061E-06,2.16453E-06,1.83845E-06,1.66182E-06,
     &     1.50517E-06,1.34852E-06,1.19718E-06,1.04822E-06,8.99264E-07,
     &     7.63432E-07,6.53806E-07,5.44186E-07,4.34564E-07,3.64210E-07,
     &     3.11938E-07,2.59667E-07,2.07395E-07,1.91456E-07,1.93639E-07,
     &     1.95821E-07,1.98004E-07,2.06442E-07,2.81546E-07/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB,N2OMULT

      WRITE (*,*) '  Computing absorption for band 8'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  
      DO LAY = 1, LAYSWTCH
         FP = FAC01(LAY) + FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(8) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(8) + 1
         INDS = INDSELF(LAY)
         COLREF1 = N2OREF(JP(LAY))
         COLREF2 = N2OREF(JP(LAY)+1)
         WCOMB1 = H2OREF(JP(LAY))
         WCOMB2 = H2OREF(JP(LAY)+1)
         RATIO = (COLREF1/WCOMB1)+FP*((COLREF2/WCOMB2)-(COLREF1/WCOMB1))
         CURRN2O = COLH2O(LAY) * RATIO
         N2OMULT = COLN2O(LAY) - CURRN2O 
         DO IG = 1, NG(8)
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG) + 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
     &           + CO2MULT(LAY) * ABSCO2A(IG)
     &           + N2OMULT * ABSN2OA(IG)
            FRACS(8,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      DO LAY = LAYSWTCH+1, NLEVELS
         FP = FAC01(LAY) + FAC11(LAY)
         IND0 = ((JP(LAY)-7)*5+(JT(LAY)-1))*NSPB(8) + 1
         IND1 = ((JP(LAY)-6)*5+(JT1(LAY)-1))*NSPB(8) + 1
         COLREF1 = N2OREF(JP(LAY))
         COLREF2 = N2OREF(JP(LAY)+1)
         WCOMB1 = O3REF(JP(LAY))
         WCOMB2 = O3REF(JP(LAY)+1)
         RATIO = (COLREF1/WCOMB1)+FP*((COLREF2/WCOMB2)-(COLREF1/WCOMB1))
         CURRN2O = COLO3(LAY) * RATIO
         N2OMULT = COLN2O(LAY) - CURRN2O 
         DO IG = 1, NG(8)
            TAUG(LAY,IG) = COLO3(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
     &           + CO2MULT(LAY) * ABSCO2B(IG)
     &           + N2OMULT * ABSN2OB(IG)
            FRACS(8,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB9

C     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K9/       KA(11,5,13,MG),KB(5,13:59,MG),SELFREF(10,MG)

      DIMENSION ABSA(715,MG),ABSB(235,MG),ABSN2O(3*MG)
      DIMENSION FRACREFA(MG,9),FRACREFB(MG)
      DIMENSION N2OREF(13),H2OREF(13),CH4REF(13),ETAREF(11)

C     From P = 1053.6 mb.
      DATA FRACREFA/
     &    0.16898900,0.15898301,0.13575301,0.12600900,
     &    0.11545800,0.09879170,0.08106830,0.06063440,
     &    0.03988780,0.00421760,0.00346635,0.00278779,
     &    0.00206225,0.00132324,0.00050033,0.00007038,
     &    0.18209399,0.15315101,0.13571000,0.12504999,
     &    0.11379100,0.09680810,0.08008570,0.05970280,
     &    0.03942860,0.00413383,0.00343186,0.00275558,
     &    0.00204657,0.00130219,0.00045454,0.00005664,
     &    0.18459500,0.15512000,0.13395500,0.12576801,
     &    0.11276800,0.09645190,0.07956650,0.05903340,
     &    0.03887050,0.00412226,0.00339453,0.00273518,
     &    0.00196922,0.00119411,0.00040263,0.00005664,
     &    0.18458800,0.15859900,0.13278100,0.12589300,   
     &    0.11272700,0.09599660,0.07903030,0.05843600,
     &    0.03843400,0.00405181,0.00337980,0.00263818,
     &    0.00186869,0.00111807,0.00040263,0.00005664,
     &    0.18459301,0.16176100,0.13235000,0.12528200,
     &    0.11237100,0.09618840,0.07833760,0.05800770,
     &    0.03787610,0.00408253,0.00330363,0.00250445,
     &    0.00176725,0.00111753,0.00040263,0.00005664,
     &    0.18454400,0.16505300,0.13221300,0.12476600,
     &    0.11158300,0.09618120,0.07797340,0.05740380,
     &    0.03742820,0.00392691,0.00312208,0.00246306,
     &    0.00176735,0.00111721,0.00040263,0.00005664,
     &    0.18452001,0.16697501,0.13445500,0.12391300,
     &    0.11059100,0.09596890,0.07761050,0.05643200,
     &    0.03686520,0.00377086,0.00309351,0.00246297,
     &    0.00176765,0.00111700,0.00040263,0.00005664,
     &    0.18460999,0.16854499,0.13922299,0.12266400,
     &    0.10962200,0.09452030,0.07653800,0.05551340,
     &    0.03609660,0.00377043,0.00309367,0.00246304,
     &    0.00176749,0.00111689,0.00040263,0.00005664,
     &    0.18312500,0.16787501,0.14720701,0.12766500,
     &    0.10890900,0.08935530,0.07310870,0.05443140,
     &    0.03566380,0.00376446,0.00309521,0.00246510,   
     &    0.00176139,0.00111543,0.00040263,0.00005664/
C     From P = 0.071 mb.
      DATA FRACREFB/
     &    0.20148601,0.15252700,0.13376500,0.12184600,
     &    0.10767800,0.09307410,0.07674570,0.05876940,
     &    0.04001480,0.00424612,0.00346896,0.00269954,
     &    0.00196864,0.00122562,0.00043628,0.00004892/

      DATA N2OREF/
     &     3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,3.20000E-07,
     &     3.19652E-07,3.15324E-07,3.03830E-07,2.94221E-07,2.84953E-07,
     &     2.76714E-07,2.64709E-07,2.42847E-07/
      DATA H2OREF/
     &     1.8759999E-02, 1.2223309E-02, 5.8908667E-03, 2.7675382E-03,  
     &     1.4065107E-03, 7.5969833E-04, 3.8875898E-04, 1.6542293E-04,  
     &     3.7189537E-05, 7.4764857E-06, 4.3081886E-06, 3.3319423E-06,  
     &     3.2039343E-06/  
      DATA CH4REF/
     &     1.7000001E-06, 1.7000001E-06, 1.6998713E-06, 1.6904165E-06,  
     &     1.6671424E-06, 1.6350652E-06, 1.6097551E-06, 1.5590465E-06,  
     &     1.5119849E-06, 1.4741138E-06, 1.4384609E-06, 1.4002215E-06,  
     &     1.3573376E-06/
      DATA ETAREF/
     &     0.,0.125,0.25,0.375,0.5,0.625,0.75,0.875,0.96,0.99,1.0/
      DATA ABSN2O/
C     From P = 952.
     &     3.26267E-01,2.42869E-00,1.15455E+01,7.39478E-00,
     &     5.16550E-00,2.54474E-00,3.53082E-00,3.82278E-00,
     &     1.81297E-00,6.65313E-01,1.23652E-01,1.83895E-03,
     &     1.70592E-03,2.68434E-09,0.,0.,
C     From P = 620.
     &     2.08632E-01,1.11865E+00,4.95975E+00,8.10907E+00,
     &     1.10408E+01,5.45460E+00,4.18611E+00,3.53422E+00,
     &     2.54164E+00,3.65093E-01,5.84480E-01,2.26918E-01,
     &     1.36230E-03,5.54400E-10,6.83703E-10,0.,
C     From P=313.
     &     6.20022E-02,2.69521E-01,9.81928E-01,1.65004E-00,
     &     3.08089E-00,5.38696E-00,1.14600E+01,2.41211E+01,
     &     1.69655E+01,1.37556E-00,5.43254E-01,3.52079E-01,
     &     4.31888E-01,4.82523E-06,5.74747E-11,0./

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      STRRAT = 21.6282
      IOFF = 0

      WRITE (*,*) '  Computing absorption for band 9'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         JFRAC = JS
         FS = AMOD(SPECMULT,1.0)
         FFRAC = FS
         IF (JS .EQ. 8) THEN
            IF (FS. LE. 0.68) THEN
               FS = FS/0.68
            ELSEIF (FS .LE. 0.92) THEN
               JS = JS + 1
               FS = (FS-0.68)/0.24
            ELSE
               JS = JS + 2
               FS = (FS-0.92)/0.08
            ENDIF
         ELSEIF (JS .EQ.9) THEN
            JS = 10
            FS = 1.
         ENDIF
         FP = FAC01(LAY) + FAC11(LAY)
         NS = JS + INT(FS + 0.5)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(9) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(9) + JS
         INDS = INDSELF(LAY)
         IF (LAY .EQ. LAYLOW) IOFF = 16
         IF (LAY .EQ. LAYSWTCH) IOFF = 32
         COLREF1 = N2OREF(JP(LAY))
         COLREF2 = N2OREF(JP(LAY)+1)
         IF (NS .EQ. 11) THEN
            WCOMB1 = H2OREF(JP(LAY))
            WCOMB2 = H2OREF(JP(LAY)+1)
         ELSE
            WCOMB1 = STRRAT * CH4REF(JP(LAY))/(1.-ETAREF(NS))
            WCOMB2 = STRRAT * CH4REF(JP(LAY)+1)/(1.-ETAREF(NS))
         ENDIF
         RATIO = (COLREF1/WCOMB1)+FP*((COLREF2/WCOMB2)-(COLREF1/WCOMB1))
         CURRN2O = SPECCOMB * RATIO
         N2OMULT = COLN2O(LAY) - CURRN2O 
         DO IG = 1, NG(9)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+11,IG) +
     &           FAC110 * ABSA(IND0+12,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+11,IG) +
     &           FAC111 * ABSA(IND1+12,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
     &           + N2OMULT * ABSN2O(IG+IOFF)
            FRACS(9,LAY,IG) = FRACREFA(IG,JFRAC) + FFRAC *
     &           (FRACREFA(IG,JFRAC+1) - FRACREFA(IG,JFRAC))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(9) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(9) + 1
         DO IG = 1, NG(9)
            TAUG(LAY,IG) = COLCH4(LAY) *
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
            FRACS(9,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB10

C     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /K10/      KA(5,13,MG), KB(5,13:59,MG)

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C     From P = 473 mb.
      DATA FRACREFA/
     &    0.16271301,0.15141940,0.14065412,0.12899506,
     &    0.11607002,0.10142808,0.08116794,0.06104711,
     &    0.04146209,0.00447386,0.00372902,0.00287258,
     &    0.00206028,0.00134634,0.00049232,0.00006927/
C     From P = 1.17 mb.
      DATA FRACREFB/
     &    0.16571465,0.15262246,0.14036226,0.12620729,
     &    0.11477834,0.09967982,0.08155201,0.06159503,
     &    0.04196607,0.00453940,0.00376881,0.00300437,
     &    0.00223034,0.00139432,0.00051516,0.00007095/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

      WRITE (*,*) '  Computing absorption for band 10'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  
      DO LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(10) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(10) + 1
         DO IG = 1, NG(10)
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
            FRACS(10,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(10) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(10) + 1
         DO IG = 1, NG(10)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
            FRACS(10,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB11

C     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K11/      KA(5,13,MG), KB(5,13:59,MG) , SELFREF(10,MG)

      DIMENSION ABSA(65,MG),ABSB(235,MG)

      DIMENSION FRACREFA(MG),FRACREFB(MG)

C     From P = 473 mb.
      DATA FRACREFA/
     &    0.14152819,0.13811260,0.14312185,0.13705885,
     &    0.11944738,0.10570189,0.08866373,0.06565409,
     &    0.04428961,0.00481540,0.00387058,0.00329187,
     &    0.00238294,0.00150971,0.00049287,0.00005980/
C     From P = 1.17 mb.
      DATA FRACREFB/
     &    0.10874039,0.15164889,0.15149839,0.14515044,
     &    0.12486220,0.10725017,0.08715712,0.06463144,
     &    0.04332319,0.00441193,0.00393819,0.00305960,
     &    0.00224221,0.00145100,0.00055586,0.00007934/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

      WRITE (*,*) '  Computing absorption for band 11'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(11) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(11) + 1
         INDS = INDSELF(LAY)
         DO IG = 1, NG(11)
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG) +
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))
            FRACS(11,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(11) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(11) + 1
         DO IG = 1, NG(11)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
            FRACS(11,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB12

C     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K12/      KA(9,5,13,MG),SELFREF(10,MG)

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C     From P = 706.3 mb.
      DATA FRACREFA/
     &    0.21245100,0.15164700,0.14486700,0.13075501,
     &    0.11629600,0.09266050,0.06579930,0.04524000,
     &    0.03072870,0.00284297,0.00234660,0.00185208,
     &    0.00133978,0.00082214,0.00031016,0.00004363,
     &    0.14703900,0.16937999,0.15605700,0.14159000,
     &    0.12088500,0.10058500,0.06809110,0.05131470,
     &    0.03487040,0.00327281,0.00250183,0.00190024,
     &    0.00133978,0.00082214,0.00031016,0.00004363,
     &    0.13689300,0.16610400,0.15723500,0.14299500,
     &    0.12399400,0.09907820,0.07169690,0.05367370,
     &    0.03671630,0.00378148,0.00290510,0.00221076,
     &    0.00142810,0.00093527,0.00031016,0.00004363,
     &    0.13054299,0.16273800,0.15874299,0.14279599,
     &    0.12674300,0.09664900,0.07462200,0.05620080,
     &    0.03789090,0.00411690,0.00322920,0.00245036,
     &    0.00178303,0.00098595,0.00040802,0.00010150,
     &    0.12828299,0.15824600,0.15688400,0.14449100,
     &    0.12787800,0.09517830,0.07679350,0.05890820,
     &    0.03883570,0.00442304,0.00346796,0.00255333,
     &    0.00212519,0.00116168,0.00067065,0.00010150,
     &    0.12649800,0.15195100,0.15646499,0.14569700,
     &    0.12669300,0.09653520,0.07887920,0.06106920,
     &    0.04043910,0.00430390,0.00364453,0.00314360,
     &    0.00203206,0.00187787,0.00067075,0.00010150,
     &    0.12500300,0.14460599,0.15672199,0.14724600,
     &    0.11978900,0.10190200,0.08196710,0.06315770,
     &    0.04240100,0.00433645,0.00404097,0.00329466,
     &    0.00288491,0.00187803,0.00067093,0.00010150,
     &    0.12317200,0.14118700,0.15242000,0.13794300,
     &    0.12119200,0.10655400,0.08808350,0.06521370,
     &    0.04505680,0.00485949,0.00477105,0.00401468,
     &    0.00288491,0.00187786,0.00067110,0.00010150,
     &    0.10193600,0.11693000,0.13236099,0.14053200,
     &    0.13749801,0.12193100,0.10221000,0.07448910,
     &    0.05205320,0.00572312,0.00476882,0.00403380,
     &    0.00288871,0.00187396,0.00067218,0.00010150/

      EQUIVALENCE (KA,ABSA)
      REAL KA
      STRRAT1 = 0.009736757

      WRITE (*,*) '  Computing absorption for band 12'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(12) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(12) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(12)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            FRACS(12,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         DO IG = 1, NG(12)
            TAUG(LAY,IG) = 0.0
            FRACS(12,LAY,IG) = 0.0
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB13

C     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K13/      KA(9,5,13,MG),SELFREF(10,MG)

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C     From P = 706.3 mb.
      DATA FRACREFA/
     &    0.17683899,0.17319500,0.15712699,0.13604601,
     &    0.10776200,0.08750010,0.06808820,0.04905150,
     &    0.03280360,0.00350836,0.00281864,0.00219862,
     &    0.00160943,0.00101885,0.00038147,0.00005348,
     &    0.17535400,0.16999300,0.15610200,0.13589200,
     &    0.10842100,0.08988550,0.06943920,0.04974900,
     &    0.03323400,0.00352752,0.00289402,0.00231003,
     &    0.00174659,0.00101884,0.00038147,0.00005348,
     &    0.17409500,0.16846400,0.15641899,0.13503000,
     &    0.10838600,0.08985800,0.07092720,0.05075710,
     &    0.03364180,0.00354241,0.00303507,0.00243391,
     &    0.00177502,0.00114638,0.00043585,0.00005348,
     &    0.17248300,0.16778600,0.15543500,0.13496999,
     &    0.10826300,0.09028740,0.07156720,0.05187120,
     &    0.03424890,0.00363933,0.00324715,0.00255030,
     &    0.00187380,0.00116978,0.00051229,0.00009768,
     &    0.17061099,0.16715799,0.15405200,0.13471501,
     &    0.10896400,0.09069460,0.07229760,0.05218280,
     &    0.03555340,0.00379576,0.00330240,0.00274693,
     &    0.00201587,0.00119598,0.00061885,0.00009768,
     &    0.16789700,0.16629100,0.15270300,0.13360199,
     &    0.11047200,0.09151080,0.07325000,0.05261450,
     &    0.03657990,0.00450092,0.00349537,0.00283321,
     &    0.00208396,0.00140354,0.00066587,0.00009768,
     &    0.16412200,0.16387400,0.15211500,0.13062200,
     &    0.11325100,0.09348130,0.07381380,0.05434740,
     &    0.03803160,0.00481346,0.00393592,0.00296633,
     &    0.00222532,0.00163762,0.00066648,0.00009768,
     &    0.15513401,0.15768200,0.14850400,0.13330200,
     &    0.11446500,0.09868230,0.07642050,0.05624170,
     &    0.04197810,0.00502288,0.00429452,0.00315347,
     &    0.00263559,0.00171772,0.00066860,0.00009768,
     &    0.15732600,0.15223300,0.14271900,0.13563600,
     &    0.11859600,0.10274200,0.07934560,0.05763410,
     &    0.03921740,0.00437741,0.00337921,0.00280212,
     &    0.00200156,0.00124812,0.00064664,0.00009768/

      EQUIVALENCE (KA,ABSA)
      REAL KA
      STRRAT1 = 16658.87

      WRITE (*,*) '  Computing absorption for band 13'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLN2O(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(13) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(13) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(13)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            FRACS(13,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         DO IG = 1, NG(13)
            TAUG(LAY,IG) = 0.0
            FRACS(13,LAY,IG) = 0.0
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB14

C     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K14/      KA(5,13,MG), KB(5,13:59,MG), SELFREF(10,MG)

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

C     From P = 1053.6 mb.
      DATA FRACREFA/
     &    0.18446200,0.16795200,0.14949700,0.12036000,
     &    0.10440100,0.09024280,0.07435880,0.05629380,
     &    0.03825420,0.00417276,0.00345278,0.00272949,
     &    0.00200378,0.00127404,0.00050721,0.00004141/
C     From P = 0.64 mb.
      DATA FRACREFB/
     &    0.19128500,0.16495700,0.14146100,0.11904500,
     &    0.10350200,0.09151190,0.07604270,0.05806020,
     &    0.03979950,0.00423959,0.00357439,0.00287559,
     &    0.00198860,0.00116529,0.00043616,0.00005987/

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

      WRITE (*,*) '  Computing absorption for band 14'

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(14) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(14) + 1
         INDS = INDSELF(LAY)
         DO IG = 1, NG(14)
            TAUG(LAY,IG) = COLCO2(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG) +
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))
            FRACS(14,LAY,IG) = FRACREFA(IG)
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(14) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(14) + 1
         DO IG = 1, NG(14)
            TAUG(LAY,IG) = COLCO2(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
            FRACS(14,LAY,IG) = FRACREFB(IG)
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB15

C     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K15/      KA(9,5,13,MG),SELFREF(10,MG)

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C     From P = 1053.6 mb.
      DATA FRACREFA/
     &    0.11287100,0.12070200,0.12729000,0.12858100,
     &    0.12743001,0.11961800,0.10290400,0.07888980,
     &    0.05900120,0.00667979,0.00552926,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.13918801,0.16353001,0.16155800,0.14090499,
     &    0.11322300,0.08757720,0.07225720,0.05173390,
     &    0.04731360,0.00667979,0.00552926,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.14687300,0.17853101,0.15664500,0.13351700,
     &    0.10791200,0.08684320,0.07158090,0.05198410,
     &    0.04340110,0.00667979,0.00552926,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.15760700,0.17759100,0.15158001,0.13193300,
     &    0.10742800,0.08693760,0.07159490,0.05196250,
     &    0.04065270,0.00667979,0.00552926,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.16646700,0.17299300,0.15018500,0.13138700,
     &    0.10735900,0.08713110,0.07130330,0.05279420,
     &    0.03766730,0.00667979,0.00552926,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.17546000,0.16666500,0.14969499,0.13105400,
     &    0.10782500,0.08718610,0.07156770,0.05308320,
     &    0.03753960,0.00432465,0.00509623,0.00436993,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.18378501,0.16064601,0.14940400,0.13146400,
     &    0.10810300,0.08775740,0.07115360,0.05400040,
     &    0.03689970,0.00388333,0.00323610,0.00353414,
     &    0.00320611,0.00204765,0.00077371,0.00010894,
     &    0.18966800,0.15744300,0.14993000,0.13152599,
     &    0.10899200,0.08858690,0.07142920,0.05399600,
     &    0.03433460,0.00374886,0.00302066,0.00240653,
     &    0.00199205,0.00204765,0.00077371,0.00010894,
     &    0.11887100,0.12479600,0.12569501,0.12839900,
     &    0.12473500,0.12012800,0.11086700,0.08493590,
     &    0.05063770,0.00328723,0.00266849,0.00210232,
     &    0.00152114,0.00095635,0.00035374,0.00004980/

      EQUIVALENCE (KA,ABSA)
      REAL KA
      STRRAT1 = 0.2883201

      WRITE (*,*) '  Computing absorption for band 15'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLN2O(LAY) + STRRAT1*COLCO2(LAY)
         SPECPARM = COLN2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(15) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(15) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(15)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            FRACS(15,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         DO IG = 1, NG(15)
            TAUG(LAY,IG) = 0.0
            FRACS(15,LAY,IG) = 0.0
         ENDDO
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB16

C     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)

      INTEGER MAXLEV, MG, NBANDS
      PARAMETER (MAXLEV=201, MG=16, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MAXLEV,MG)
      COMMON /PLANKG/   FRACS(NBANDS,MAXLEV,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLEVELS,HEIGHT(MAXLEV),
     &                  PRES(MAXLEV),TEMP(MAXLEV),
     &                  RHOAIR(MAXLEV),RHOH2O(MAXLEV),RHOO3(MAXLEV)
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,
     &                  COLH2O(MAXLEV),COLCO2(MAXLEV),
     &                  COLO3(MAXLEV),COLN2O(MAXLEV),
     &                  COLCH4(MAXLEV),CO2MULT(MAXLEV)
      COMMON /INTFAC/   FAC00(MAXLEV),FAC01(MAXLEV),
     &                  FAC10(MAXLEV),FAC11(MAXLEV)
      COMMON /INTIND/   JP(MAXLEV),JT(MAXLEV),JT1(MAXLEV)
      COMMON /SELF/     SELFFAC(MAXLEV),SELFFRAC(MAXLEV),INDSELF(MAXLEV)
      COMMON /K16/      KA(9,5,13,MG),SELFREF(10,MG)

      DIMENSION ABSA(585,MG)
      DIMENSION FRACREFA(MG,9)

C     From P = 862.6 mb.
      DATA FRACREFA/
     &    0.17356300,0.18880001,0.17704099,0.13661300,
     &    0.10691600,0.08222480,0.05939860,0.04230810,
     &    0.02526330,0.00244532,0.00193541,0.00150415,
     &    0.00103528,0.00067068,0.00024951,0.00003348,
     &    0.17779499,0.19837400,0.16557600,0.13470000,
     &    0.11013600,0.08342720,0.05987030,0.03938700,
     &    0.02293650,0.00238849,0.00192400,0.00149921,
     &    0.00103539,0.00067150,0.00024822,0.00003348,
     &    0.18535601,0.19407199,0.16053200,0.13300700,
     &    0.10779000,0.08408500,0.06480450,0.04070160,
     &    0.02203590,0.00227779,0.00189074,0.00146888,
     &    0.00103147,0.00066770,0.00024751,0.00003348,
     &    0.19139200,0.18917400,0.15748601,0.13240699,
     &    0.10557300,0.08383260,0.06724060,0.04364450,
     &    0.02175820,0.00225436,0.00184421,0.00143153,
     &    0.00103027,0.00066066,0.00024222,0.00003148,
     &    0.19547801,0.18539500,0.15442000,0.13114899,
     &    0.10515600,0.08350350,0.06909780,0.04671630,
     &    0.02168820,0.00224400,0.00182009,0.00139098,
     &    0.00102582,0.00065367,0.00023202,0.00003148,
     &    0.19757500,0.18266800,0.15208900,0.12897800,
     &    0.10637200,0.08391220,0.06989830,0.04964120,
     &    0.02155800,0.00224310,0.00177358,0.00138184,
     &    0.00101538,0.00063370,0.00023227,0.00003148,
     &    0.20145500,0.17692900,0.14940600,0.12690400,  
     &    0.10828800,0.08553720,0.07004940,0.05153430,
     &    0.02268740,0.00216943,0.00178603,0.00137754,
     &    0.00098344,0.00063165,0.00023218,0.00003148,
     &    0.20383500,0.17047501,0.14570600,0.12679300,
     &    0.11043100,0.08719150,0.07045440,0.05345420,
     &    0.02448340,0.00215839,0.00175893,0.00138296,
     &    0.00098318,0.00063188,0.00023199,0.00003148,
     &    0.18680701,0.15961801,0.15092900,0.13049100,
     &    0.11418400,0.09380540,0.07093450,0.05664280,
     &    0.02938410,0.00217751,0.00176766,0.00138275,
     &    0.00098377,0.00063181,0.00023193,0.00003148/

      EQUIVALENCE (KA,ABSA)
      REAL KA
      STRRAT1 = 830.411

      WRITE (*,*) '  Computing absorption for band 16'

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO LAY = 1, LAYTROP
         SPECCOMB = COLH2O(LAY) + STRRAT1*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB 
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*SPECPARM
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)
         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS) * FAC01(LAY)
         FAC011 = (1. - FS) * FAC11(LAY)
         FAC101 = FS * FAC01(LAY)
         FAC111 = FS * FAC11(LAY)
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(16) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(16) + JS
         INDS = INDSELF(LAY)
         DO IG = 1, NG(16)
            TAUG(LAY,IG) = SPECCOMB * 
     &          (FAC000 * ABSA(IND0,IG) +
     &           FAC100 * ABSA(IND0+1,IG) +
     &           FAC010 * ABSA(IND0+9,IG) +
     &           FAC110 * ABSA(IND0+10,IG) +
     &           FAC001 * ABSA(IND1,IG) + 
     &           FAC101 * ABSA(IND1+1,IG) +
     &           FAC011 * ABSA(IND1+9,IG) +
     &           FAC111 * ABSA(IND1+10,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            FRACS(16,LAY,IG) = FRACREFA(IG,JS) + FS *
     &           (FRACREFA(IG,JS+1) - FRACREFA(IG,JS))
         ENDDO
      ENDDO

      DO LAY = LAYTROP+1, NLEVELS
         DO IG = 1, NG(16)
            TAUG(LAY,IG) = 0.0
            FRACS(16,LAY,IG) = 0.0
         ENDDO
      ENDDO

      RETURN
      END
