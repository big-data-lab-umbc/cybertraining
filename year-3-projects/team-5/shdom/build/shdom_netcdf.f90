module shdom_netcdf
  ! 
  ! Netcdf interface for SHDOM input and output. 
  !
  ! Robert Pincus, University of Colordado, Dec 2008
  !
  use netcdf
  implicit none
  
  private :: IntToChar, radiance_Name
  
contains
  ! -------------------------------------------------------------------------------------------
  subroutine convert_Prp_to_netcdf(asciiFileName, netcdfFileName)
    character(len = *), intent(in) :: asciiFileName, netcdfFileName
    !
    ! Convert an ASCII shdom property file of type 'T' to a netcdf version. 
    !
    character(len=1) :: proptype
    integer          :: NPX, NPY, NPZ, numPhase, numL, MAXLEG, I, j, IX, IY, IZ, iph
    real             :: delx, dely, zlevs, chi, ext, alb, tmp

    real,    dimension(:, :, :), allocatable :: extinction, singleScatteringAlbedo, temperature
    integer, dimension(:, :, :), allocatable :: phaseFunctionIndex
    real,    dimension(:, :),    allocatable :: phaseFunctionCoefficients
    real,    dimension(:),       allocatable :: zLevels
    
      !
      ! Open the ASCII file, ensure it's of type "T", and get the grid sizes, 
      !   number of phase functions, and maximum number of terms in the phase function series
      !
      OPEN (UNIT=1, FILE=trim(asciiFileName), STATUS='OLD')
      READ (1,'(A1)') PROPTYPE
      ! Property file type T is for tabulated phase function format
      IF (PROPTYPE /= 'T') stop "Can only convert property file format 'T'" 
      READ (1,*) NPX, NPY, NPZ
      READ (1,*) DELX, DELY, (zlevs, j=1,NPZ)
      MAXLEG = 0
      READ (1,*) NUMPHASE
      DO I = 1, NUMPHASE
        READ (1,*) NUML, (CHI, j=1,NUML)
        MAXLEG = MAX(NUML, MAXLEG)
      ENDDO
      
      !
      !  Allocate 3D arrays for the properties. Initialize to zero since the ASCII format
      !     doesn't require every cell to be filled in
      !
      allocate(extinction(NPX, NPY, NPZ),  singleScatteringAlbedo(NPX, NPY, NPZ), &
               temperature(NPX, NPY, NPZ), phaseFunctionIndex(NPX, NPY, NPZ),     &
               zLevels(NPZ))
      allocate(phaseFunctionCoefficients(maxleg, numPhase))
      extinction(:, :, :) = 0.; singleScatteringAlbedo(:, :, :) = 0.
      temperature(:, :, :) = 0.; phaseFunctionIndex(:, :, :) = 0
      zLevels(:) = 0. 
      phaseFunctionCoefficients(:, :) = 0.    
      
      !
      ! Starting at the top, read in the phase functions... 
      !
      rewind(1)
      READ (1,'(A1)') PROPTYPE
      READ (1,*) NPX, NPY, NPZ
      READ (1,*) DELX, DELY, zLevels(:) 
      READ (1,*) NUMPHASE
      DO I = 1, NUMPHASE
        READ (1,*) NUML, phaseFunctionCoefficients(:numL, i) 
      ENDDO
      !
      ! ... and then the properties at each grid point until we hit the end of the file
      !
      DO WHILE (.TRUE.)
        READ (1,*,END=290) IX, IY, IZ, TMP, EXT, ALB, IPH
        IF ( all((/ IX >= 1, IX <= NPX, &
                    IY >= 1, IY <= NPY, &
                    IZ >= 1, IZ <= NPZ /)) ) THEN
          temperature(ix, iy, iz) = tmp
          extinction( ix, iy, iz) = ext
          singleScatteringAlbedo(ix, iy, iz) = alb
          phaseFunctionIndex(    ix, iy, iz) = IPH
        ENDIF
      END DO
 290  close(1)
 
    call write_properties_netcdf(extinction, singleScatteringAlbedo, phaseFunctionIndex, &
                                 temperature, phaseFunctionCoefficients,                 &
                                 delX, delY, zLevels, netcdfFileName)
                                 
    deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex, &
               phaseFunctionCoefficients, temperature, zLevels)
  end subroutine convert_Prp_to_netcdf
! -------------------------------------------------------------------------------------------
  subroutine write_properties_netcdf(extinction, singleScatteringAlbedo, phaseFunctionIndex, &
                                     temperature, phaseFunctionCoefficients,                 &
                                     delX, delY, zLevels, fileName)
    real,    dimension(:, :, :), intent(in) :: extinction, singleScatteringAlbedo, temperature
    integer, dimension(:, :, :), intent(in) :: phaseFunctionIndex
    real,    dimension(:, :),    intent(in) :: phaseFunctionCoefficients
    real,                        intent(in) :: delX, delY
    real,    dimension(:),       intent(in) :: zLevels
    character(len = *),          intent(in) :: fileName
    !
    ! Writes a netcdf property file for shdom_netcdf
    !   Dimensions are x, y, z; also legendreTerm (MAXLEG), phaseFunction (NUMPHASE)
    !   Variables are 
    !    real: extinction(x, y, z); singleScatteringAlbedo(x, y, z), temperature(x, y, z); 
    !    integer: phaseFunctionIndex(x, y, z)  
    !    real: phaseFunctionCoefficients(legendreTerm, phaseFunction)
    !    real: z(z) - dimension variable
    !    The file also needs to have global attibutes DelX, DelY OR contain 
    !    dimension variables for x, y; these must be evenly spaced
    !    
    integer                :: numX, numY, numZ, j
    integer                :: ncFileId, xDimId, yDimId, zDimId, termDimId, phaseFunctionDimId, ncVarID
    integer, dimension(24) :: status
    
    numX = size(extinction, 1)
    numY = size(extinction, 2)
    numZ = size(extinction, 3)
    !
    ! Ensure that all the 3D arrays are the same shape
    !
    if(any((/ size(singleScatteringAlbedo, 1) /= numX,      &
              size(singleScatteringAlbedo, 2) /= numY,      &
              size(singleScatteringAlbedo, 3) /= numZ /)) ) &
      stop "write_properties_netcdf: singleScattingAlbedo isn't the same shape as extinction" 
    if(any((/ size(phaseFunctionIndex, 1) /= numX,      &
              size(phaseFunctionIndex, 2) /= numY,      &
              size(phaseFunctionIndex, 3) /= numZ /)) ) &
      stop "write_properties_netcdf: phaseFunctionIndex isn't the same shape as extinction" 
    if(any((/ size(temperature, 1) /= numX,      &
              size(temperature, 2) /= numY,      &
              size(temperature, 3) /= numZ /)) ) &
      stop "write_properties_netcdf: temperature isn't the same shape as extinction" 
    !
    ! ... and that the zLevels is consistent with the 
    !
    if(size(zLevels) /= numZ) stop "write_properties_netcdf: zLevels isn't the same extent as extinction"

    !
    ! Create the netcdf file
    !
    status( :) = nf90_NoErr
    status( 1) = nf90_create(trim(fileName), nf90_clobber, ncFileId) 
    !
    ! Define the dimensions
    !
    status( 2) = nf90_def_dim(ncFileId, 'x', numX, xDimId) 
    status( 3) = nf90_def_dim(ncFileId, 'y', numY, yDimId) 
    status( 4) = nf90_def_dim(ncFileId, 'z', numZ, zDimId) 
    status( 5) = nf90_def_dim(ncFileId, 'legendreTerm',  size(phaseFunctionCoefficients, 1), &
                              termDimId) 
    status( 6) = nf90_def_dim(ncFileId, 'phaseFunction', size(phaseFunctionCoefficients, 2), &
                              phaseFunctionDimId) 
    !
    ! Define the variables - first the one dimension variable, then the fields themselves 
    !
    status( 7) = nf90_def_var(ncFileId, 'x', nf90_float, xDimId, ncVarId)
    status( 8) = nf90_def_var(ncFileId, 'y', nf90_float, yDimId, ncVarId)
    status( 9) = nf90_def_var(ncFileId, 'z', nf90_float, zDimId, ncVarId)
    status(10) = nf90_def_var(ncFileId, 'extinction',             nf90_float, &
                              (/ xDimId, yDimId, zDimId /), ncVarId)
    status(11) = nf90_def_var(ncFileId, 'singleScatteringAlbedo', nf90_float, &
                              (/ xDimId, yDimId, zDimId /), ncVarId)
    status(12) = nf90_def_var(ncFileId, 'temperature',            nf90_float, &
                              (/ xDimId, yDimId, zDimId /), ncVarId)
    status(13) = nf90_def_var(ncFileId, 'phaseFunctionIndex',     nf90_short, &
                              (/ xDimId, yDimId, zDimId /), ncVarId)
                              
    status(14) = nf90_def_var(ncFileId, 'phaseFunctionCoefficients', nf90_float, &
                              (/ termDimId, phaseFunctionDimId /), ncVarId)
    !
    ! DelX, DelY attributes
    !
    status(15) = nf90_put_att(ncFileId, nf90_global, "DelX", delx)
    status(16) = nf90_put_att(ncFileId, nf90_global, "DelY", dely)
    status(17) = nf90_endDef(ncFileId) 
    if(any(status(:) /= nf90_noErr)) print *, "Error defining netcdf property file " // trim(fileName)  
    !
    ! Now write the fields out and close the file 
    !
    status( 1) = nf90_inq_varid(ncFileId, 'x', ncVarId) 
    status( 2) = nf90_put_var(ncFileId, ncVarId, (/ (j * delX, j = 0, numX - 1) /)) 
    status( 3) = nf90_inq_varid(ncFileId, 'y', ncVarId) 
    status( 4) = nf90_put_var(ncFileId, ncVarId, (/ (j * delY, j = 0, numY - 1) /)) 
    status( 5) = nf90_inq_varid(ncFileId, 'z', ncVarId) 
    status( 6) = nf90_put_var(ncFileId, ncVarId, zLevels) 
    status( 7) = nf90_inq_varid(ncFileId, 'extinction', ncVarId) 
    status( 8) = nf90_put_var(ncFileId, ncVarId, extinction) 
    status( 9) = nf90_inq_varid(ncFileId, 'singleScatteringAlbedo', ncVarId) 
    status(10) = nf90_put_var(ncFileId, ncVarId, singleScatteringAlbedo) 
    status(11) = nf90_inq_varid(ncFileId, 'temperature', ncVarId) 
    status(12) = nf90_put_var(ncFileId, ncVarId, temperature) 
    status(13) = nf90_inq_varid(ncFileId, 'phaseFunctionIndex', ncVarId) 
    status(14) = nf90_put_var(ncFileId, ncVarId, phaseFunctionIndex) 
    status(15) = nf90_inq_varid(ncFileId, 'phaseFunctionCoefficients', ncVarId) 
    status(16) = nf90_put_var(ncFileId, ncVarId, phaseFunctionCoefficients) 
    status(17) = nf90_close(ncFileId) 
    if(any(status(:) /= nf90_noErr)) print *, "Error writing netcdf property file " // trim(fileName)  
    
  end subroutine write_properties_netcdf
  ! -------------------------------------------------------------------------------------------
  subroutine read_property_size_netcdf(PROPFILE, NLEG, NPX, NPY, NPZ, &
                                       NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
    character(len = *), intent(in)  :: propfile
    INTEGER,            intent(in)  :: NLEG
    INTEGER,            intent(out) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    REAL,               intent(out) :: DELX, DELY
    
    integer                            :: ncFileId, ncDimId, ncVarId
    integer                            :: i
    integer,           dimension(16)   :: status 
    real, allocatable, dimension(:)    :: dimVals
  
    status( :) = nf90_noErr
    status( 1) = nf90_open(trim(propfile), NF90_NOWRITE, ncFileId) 
    !
    ! Get the length of each dimension
    !
    status( 2) = nf90_inq_dimid(ncFileId, 'x', ncDimId)
    status( 3) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = NPX)
    status( 4) = nf90_inq_dimid(ncFileId, 'y', ncDimId)
    status( 5) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = NPY)
    status( 6) = nf90_inq_dimid(ncFileId, 'z', ncDimId)
    status( 7) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = NPZ)
    status( 8) = nf90_inq_dimid(ncFileId, 'legendreTerm', ncDimId)
    status( 9) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = MAXLEG)
    status(10) = nf90_inq_dimid(ncFileId, 'phaseFunction', ncDimId)
    status(11) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = NUMPHASE)
    if(any(status(:) /= nf90_NoErr)) stop "Error reading dimension sizes from netcdf file" 
    MAXLEG = max(MAXLEG, NLEG) 
    
    !
    ! Determine DELX,Y. These are determined by default from the dimension variable 
    !   but can alternatively be set as a global attribute
    !
    ! X 
    if(nf90_inq_varid(ncFileId, 'x', ncVarId) == nf90_NoErr) then 
      allocate(dimVals(NPX))
      if(nf90_get_var(ncFileID, ncVarId, dimVals) /= nf90_noErr) &
        stop "Error reading x dimension values from netcdf file" 
      delx = dimvals(2) - dimvals(1) 
      if(any(abs((dimvals(2:NPX) - dimvals(1:NPX-1)) - delx) > 2. * spacing(dimvals(2:NPX)))) &
        stop "X values are not evenly spaced in netcdf file"
      deallocate(dimvals) 
    else
      if(nf90_get_att(ncFileID, nf90_global, "DelX", DELX) /= nf90_NoErr) &
        stop "Error reading DELX attributes from netcdf file" 
    end if
    
    ! Y
    if(nf90_inq_varid(ncFileId, 'y', ncVarId) == nf90_NoErr) then 
      allocate(dimVals(NPY))
      if(nf90_get_var(ncFileID, ncVarId, dimVals) /= nf90_noErr) &
        stop "Error reading y dimension values from netcdf file" 
      dely = dimvals(2) - dimvals(1) 
      if(any(abs((dimvals(2:NPY) - dimvals(1:NPY-1)) - delx) > 2. * spacing(dimvals(2:NPY)))) &
        stop "Y values are not evenly spaced in netcdf file"
      deallocate(dimvals) 
    else
      if(nf90_get_att(ncFileID, nf90_global, "DelY", DELY) /= nf90_NoErr) &
        stop "Error reading DELY attributes from netcdf file" 
    end if 
    
    MAXPGL = MAXLEG * NUMPHASE
    status( 1) = nf90_close(ncFileId)
    
  end subroutine read_property_size_netcdf
  ! -------------------------------------------------------------------------------------------
  subroutine read_properties_netcdf(PROPFILE, NPX, NPY, NPZ, NUMPHASE,   &
                                    MAXLEG, MAXPG, MAXPGL, DELTAM, NLEG, &
                                    PROPTYPE, ZLEVELS, MAXASYM,          &
                                    TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP) 
    character(len = *), intent(in   )  :: propfile
    INTEGER,            intent(in   )  :: NPX, NPY, NPZ, NUMPHASE
    INTEGER,            intent(in   )  :: MAXPG, MAXPGL, MAXLEG
    LOGICAL,            intent(in   )  :: DELTAM
    integer,            intent(inout) :: NLEG
    character,          intent(  out) :: PROPTYPE
    real, dimension(:), intent(  out) :: zLevels
    real,               intent(  out) :: MAXASYM
    real, dimension(:), intent(  out) :: TEMPP, EXTINCTP, ALBEDOP, LEGENP
    INTEGER*2, &
          dimension(:), intent(  out) :: IPHASEP
! 
!   If doing delta-M then NLEG is the minimum number of Legendre terms on
!   input and the actual maximum number of terms on output, except that
!   it may not exceed MAXLEG; otherwise NLEG is the number of Legendre
!   terms to be used (regardless of what is in property file).

    integer                                    :: ncFileId, ncDimId, ncVarId, numTerms, i, j
    integer,                dimension(8   )    :: status
    real,      allocatable, dimension(:, :, :) :: tempArray
    integer*2, allocatable, dimension(:, :, :) :: tempIntArray

    IF (NPX*NPY*NPZ > MAXPG)  STOP 'READ_PROPERTIES: MAXPG exceeded'
    PROPTYPE = 'T' 
    status( :) = nf90_noErr
    status( 1) = nf90_open(trim(propfile), NF90_NOWRITE, ncFileId) 
    
    !
    ! 3D  arrays - real for singleScatteringAlbedo, extinction, and (optionally) temperature; 
    !  integer for phaseFunctionIndex
    !
    allocate(tempArray(NPZ, 1, 1))
    
    status(1) = nf90_inq_varid(ncFileId, "z", ncVarId)
    status(2) = nf90_get_var(ncFileId, ncVarId, tempArray(:, 1, 1))
    if(any(status(:) /= nf90_NoErr)) stop "Error reading z-Levels from netcdf file" 
    zLevels = tempArray(:, 1, 1)
    
    deallocate(tempArray)
    allocate(tempArray(NPX, NPY, NPZ)) 
    
    status(1) = nf90_inq_varid(ncFileId, "extinction", ncVarId)
    status(2) = nf90_get_var(ncFileId, ncVarId, tempArray)
    if(any(status(:) /= nf90_NoErr)) stop "Error reading extinction from netcdf file" 
    EXTINCTP(:NPX*NPY*NPZ) = pack(reshape(tempArray, shape = (/ NPZ, NPY, NPX /), order = (/ 3, 2, 1 /) ), &
                                  mask = .true.) 

    status(1) = nf90_inq_varid(ncFileId, "singleScatteringAlbedo", ncVarId)
    status(2) = nf90_get_var(ncFileId, ncVarId, tempArray)
    if(any(status(:) /= nf90_NoErr)) stop "Error reading singleScatteringAlbedo from netcdf file" 
    ALBEDOP(:NPX*NPY*NPZ) = pack(reshape(tempArray, shape = (/ NPZ, NPY, NPX /), order = (/ 3, 2, 1 /) ), &
                                 mask = .true.) 

    !
    ! Temperature may or may not be present - return 0s if not
    !
    if(nf90_inq_varid(ncFileId, "temperature", ncVarId) == nf90_NoErr) then 
      if(nf90_get_var(ncFileId, ncVarId, tempArray) /= nf90_noErr) &
        stop "Error reading temperature from netcdf file"
      TEMPP(:NPX*NPY*NPZ) = pack(reshape(tempArray, shape = (/ NPZ, NPY, NPX /), order = (/ 3, 2, 1 /) ), &
                                 mask = .true.) 
    else
      tempp(:NPX*NPY*NPZ) = 250. 
    end if 
    
    deallocate(tempArray)
    allocate(tempIntArray(NPX, NPY, NPZ)) 
    
    status(1) = nf90_inq_varid(ncFileId, "phaseFunctionIndex", ncVarId)
    status(2) = nf90_get_var(ncFileId, ncVarId, tempIntArray)
    if(any(status(:) /= nf90_NoErr)) stop "Error reading phaseFunctionIndex from netcdf file" 
    IPHASEP(:NPX*NPY*NPZ) = pack(reshape(tempIntArray, shape = (/ NPZ, NPY, NPX /), order = (/ 3, 2, 1 /) ), &
                                 mask = .true.) 
    
    deallocate(tempIntArray)
    
    !
    ! 2D array for the phase functions
    ! The number of phase function terms returned depends on whether we're doing delta-M or not, 
    !    and may be larger than the number of terms in the file
    !
    status(1) = nf90_inq_dimid(ncFileId, 'legendreTerm', ncDimId)
    status(2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = numTerms)
    allocate(tempArray(numTerms, NUMPHASE, 1))
    status(3) = nf90_inq_varid(ncFileId, "phaseFunctionCoefficients", ncVarId)
    status(4) = nf90_get_var(ncFileId, ncVarId, tempArray(:, :, 1))
    if(any(status(:) /= nf90_NoErr)) stop "Error reading phaseFunctionCoefficients from netcdf file" 
    If(deltaM) then 
      NLEG = numTerms
      IF (NUMPHASE*NLEG > MAXPGL) STOP 'READ_PROPERTIES: MAXPGL exceeded'
      LEGENP(:numTerms * NUMPHASE) = pack(tempArray(:numTerms, :, 1), mask = .true.)
    else if(numTerms >= NLEG) then
      IF (NUMPHASE*NLEG > MAXPGL) STOP 'READ_PROPERTIES: MAXPGL exceeded'
      LEGENP(:NLEG * NUMPHASE) = pack(tempArray(:NLEG, :, 1), mask = .true.)
    else 
      IF (NUMPHASE*NLEG > MAXPGL) STOP 'READ_PROPERTIES: MAXPGL exceeded'
      !
      ! Need to pad with 0s.
      !
      do i = 1, numPhase
        LEGENP(i*(NLEG-1)+1:i*NLEG) = (/ tempArray(:, i, 1), (0., j = numTerms + 1, NLEG) /) 
      end do 
    end if 
    MAXASYM = maxval(tempArray(1, :, 1))/3.0
    
    deallocate(tempArray)

  end subroutine read_properties_netcdf
  ! -------------------------------------------------------------------------------------------
  subroutine output_results_netcdf_par(NX, NY, NZ, NPTS, NCELLS,                            &
                                       NSH, ML,MM,NLM, NMU,NPHI,NANG, NG,                   &
                                       PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,  &
                                       BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,           &
                                       SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,                 &
                                       SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                                       SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, CPUTIME,  &
                                       XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID,               &
                                       FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,      &
                                       NUMOUT, OUTTYPES, OUTPARMS, OutFileNC)
    implicit none    
    integer, intent(in) :: NX, NY, NZ, NPTS, NCELLS, &
                           NSH, ML, MM, NLM, NMU, NPHI, NANG, NG
    integer, intent(in) :: BCFLAG, IPFLAG, MAXITER, TOTITER, IRAD, NSHOUT, NUMOUT
    logical, intent(in) :: deltaM
    real,    intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ, GNDTEMP, GNDALBEDO, SKYRAD, &
                           SOLACC, SPLITACC, SHACC, CPUTIME, XDOMAIN, YDOMAIN, WAVELEN
    character (len = *), intent(in) :: SRCTYPE, SFCTYPE, UNITS, GRIDTYPE, & 
                                       PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    real, dimension(:, :, :, :), &
                         intent(in) :: fluxes, & ! dim(3, nz, ny, nx)
                                       SHTERMS     ! dim(nshout, nz, ny, nx)
    real, dimension(:, :, :), &
                         intent(in) :: fluxdiv   ! dim(nz, ny, nx)
    real, dimension(:),  intent(in) :: XGRID, YGRID, ZGRID, WAVENO, radout
    character(len = 1), dimension(:), &          ! dim(numout)
                         intent(in) :: OUTTYPES
    real, dimension(:, :),  &                    ! dim(:, numout)
                         intent(in) :: OUTPARMS
    character (len = *), intent(in) :: OutFileNC                    
    !
    ! Write SHDOM output to a netcdf file. A single file holds all data sets. 
    !   Only the output type/format combinations supported by write_results_par
    !   are included. 
    ! Desired outputs are specifed in the same way as the ASCII outputs. It's possible
    !   to request the same output/type combination more than once. With the exception of 
    !   radiance ("R") output, requests after the first are redundent and ignored. 
    !   Radiance output dimensions and variables for the second and higher request have 
    !   ".N" added to their name
    ! When the calculation includes a solar component the downwelling fluxes are broken down into 
    !   total and diffuse components - this differs from the SHDOM convention of direct/diffuse
    
    ! Local variables
    CHARACTER(len = 32) :: GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME
    INTEGER :: NANGOUT, NXOUT, NYOUT
    REAL    :: STARTX, STARTY
    
    integer :: i, j, numRadianceOutputs, radianceStartPos
    real    :: c, z0
    logical :: outputDiffuse, &
               ! Catch duplicated output types
               F1done, F3done, F4done, H1done, H2Done, Sdone
    
    integer :: ncFileId, ncvarid
    integer :: xDimId, YDimId, zDimId, radxDimId, radyDimId, radDirDimId
    integer, dimension(64) :: status

    ! ---------------------------------------
    !
    ! Expand grid, surface, and source names
    !
    select case(gridtype(1:1)) 
      case('p', 'P')
        gridname = 'Property-File (Z)'
      case('f', 'F')
        gridname = 'Input-File (Z)'
      case default
        gridname = 'Even (Z)' 
    end select
    gridname = 'Even (X,Y) ' // trim(gridname) 

    select case(srctype(1:1)) 
      case('s', 'S')
        SOURCENAME = 'Solar'
        outputDiffuse = .true. 
      case('t', 'T')
        SOURCENAME = 'Thermal'
        outputDiffuse = .false. 
      case ('b', 'B')
        SOURCENAME = 'Solar/Thermal'
        outputDiffuse = .true. 
      case default
        SOURCENAME = "UNKNOWN"
        outputDiffuse = .false. 
    end select
    if(deltaM) sourcename = trim(sourcename) // "  Delta-M method"
    C = 1
    if (UNITS(1:1) == 'T') C = 1./acos(-1.)
    
    select case(SFCTYPE(1:2)) 
      case('fl', 'FL')
        SFCNAME = 'Fixed Lambertian'
      case('vl', 'VL')
        SFCNAME = 'Variable Lambertian'
      case('vf', 'VF')
        SFCNAME = 'Variable Fresnel'
      case('vr', 'VR')
        SFCNAME = 'Variable RPV'
      case default
        SFCNAME = "UNKNOWN"
    end select
    IF (SFCTYPE(1:1) .EQ. 'V') SFCNAME = 'Variable other'
    
    status( :) = nf90_noErr
    status( 1) = nf90_create(trim(outfileNC), nf90_clobber, ncFileId)
    ! ------------------------------------------------------
    !
    ! Metadata - this replicates what's in the header of each ASCII SHDOM 
    !   output file
    !
    status( 2) = nf90_put_att(ncFileId, NF90_Global, "description", &
      'Spherical Harmonic Discrete Ordinate Radiative Transfer Output') 

    !   Simultation parameters
    status( 3) = nf90_put_att(ncFileId, NF90_Global, 'L',    ML)
    status( 4) = nf90_put_att(ncFileId, NF90_Global, 'M',    MM) 
    status( 5) = nf90_put_att(ncFileId, NF90_Global, 'NLM',  NLM) 
    status( 6) = nf90_put_att(ncFileId, NF90_Global, 'NMu',  NMU) 
    status( 7) = nf90_put_att(ncFileId, NF90_Global, 'NPhi', NPHI) 
    status( 8) = nf90_put_att(ncFileId, NF90_Global, 'Nang', NANG) 
    status( 9) = nf90_put_att(ncFileId, NF90_Global, 'NSH',  NSH) 
    ! We're skipping nx, ny, nz beacuse they'll be in the output
    status(10) = nf90_put_att(ncFileId, NF90_Global, 'NPts',   NPTS) 
    status(11) = nf90_put_att(ncFileId, NF90_Global, 'NCells', NCELLS) 
    
    !   File names
    status(12) = nf90_put_att(ncFileId, NF90_Global, 'Property_File',    trim(PROPFILE))
    status(13) = nf90_put_att(ncFileId, NF90_Global, 'Input_Save_File',  trim(INSAVEFILE))
    status(14) = nf90_put_att(ncFileId, NF90_Global, 'Output_Save_File', trim(OUTSAVEFILE))
    status(15) = nf90_put_att(ncFileId, NF90_Global, 'Correlated_K_Dist_File', trim(CKDFILE))
    status(16) = nf90_put_att(ncFileId, NF90_Global, 'NUM_G',                  NG)
    
    status(17) = nf90_put_att(ncFileId, NF90_Global, "Source_Type",         trim(sourcename)) 
    status(18) = nf90_put_att(ncFileId, NF90_Global, "Grid_Type",           trim(GRIDNAME)) 
    status(19) = nf90_put_att(ncFileId, NF90_Global, "Independent_Pixel",   IPFLAG) 
    status(20) = nf90_put_att(ncFileId, NF90_Global, "Surface_Type",        trim(SFCNAME)) 
    status(21) = nf90_put_att(ncFileId, NF90_Global, "Horiz_Boundary_Condition", BCFLAG) 
    
    !   Solution crieria at end 
    status(22) = nf90_put_att(ncFileId, NF90_Global, 'Splitting_Accuracy',          SPLITACC) 
    status(23) = nf90_put_att(ncFileId, NF90_Global, 'Spherical_Harmonic_Accuracy', SHACC) 
    status(24) = nf90_put_att(ncFileId, NF90_Global, 'Solution_Accuracy',           SOLACC) 
    status(25) = nf90_put_att(ncFileId, NF90_Global, 'Maximum_Iterations',          MAXITER) 
    status(26) = nf90_put_att(ncFileId, NF90_Global, 'Number_Iterations',           TOTITER) 
    
    ! Boundary conditions 
    status(27) = nf90_put_att(ncFileId, NF90_Global, 'Sky_rad', SKYRAD)
    if(SOURCENAME(1:1) == 'T' .or. SOURCENAME(1:1) == 'B') & 
      status(28)  = nf90_put_att(ncFileId, NF90_Global, 'Ground_Temp', GNDTEMP)
    if(SOURCENAME(1:1) == 'S' .or. SOURCENAME(1:1) == 'B') then  
      status(29)  = nf90_put_att(ncFileId, NF90_Global, 'Solar_Flux', SOLARFLUX)
      status(30)  = nf90_put_att(ncFileId, NF90_Global, 'Solar_Mu',   SOLARMU)
      status(31)  = nf90_put_att(ncFileId, NF90_Global, 'Solar_Az',   SOLARAZ*180.0/ACOS(-1.0))
    end if 
    if(SFCNAME(1:1) == 'V') then 
      status(32)  = nf90_put_att(ncFileId, NF90_Global, 'Surface_File', trim(SFCFILE))
    else
      if(SOURCENAME(1:1) == 'T' .or. SOURCENAME(1:1) == 'B') & 
        status(33)  = nf90_put_att(ncFileId, NF90_Global, 'Ground_Emis', 1.0-GNDALBEDO)
      if(SOURCENAME(1:1) == 'S' .or. SOURCENAME(1:1) == 'B') & 
        status(34)  = nf90_put_att(ncFileId, NF90_Global, 'Ground_Albedo', GNDALBEDO)
    end if 
  
    status(35) = nf90_put_att(ncFileId, NF90_Global, "Cpu_time_total", CPUTIME)
    if(any(status(:) /= nf90_NoErr)) print *, "Error defining metadata for netcdf output" 
    
    ! ------------------------------------------------------
    !
    ! Define dimensions
    !
    ! x, y, and z on the base grid
    status( 1) = nf90_def_dim(ncFileId, "x", NX, xDimId)
    status( 2) = nf90_def_dim(ncFileId, "y", NY, yDimId)
    status( 3) = nf90_def_dim(ncFileId, "z", NZ, zDimId)
    !
    ! Dimension variables
    !
    status( 4) = nf90_def_var(ncFileId, "x", nf90_float, xDimId, ncVarId)
    status( 5) = nf90_def_var(ncFileId, "y", nf90_float, yDimId, ncVarId)
    status( 6) = nf90_def_var(ncFileId, "z", nf90_float, zDimId, ncVarId)
    if(any(status(:) /= nf90_NoErr)) print *, "Error defining default dimensions for netcdf output" 

    ! ------------------------------------------------------
    !
    ! Define variables (including units) for each output 
    !   Defer writing the data because it's substantially more efficient to define all variables 
    !   before writing any of them out
    !   If we're writing out radiance we have to define more dimensions and dim. variables
    !
    F1done = .false.; F3done = .false.; F4done = .false. 
    H1done = .false.; H2Done = .false.; Sdone = .false. 
    numRadianceOutputs = 0
    do i = 1, numOut
      select case (outtypes(i))
        case ('R') ! Radiance 
          numRadianceOutputs = numRadianceOutputs + 1 
          ! 
          ! Radiance dimensions and dimension variables 
          !
          STARTX = -outparms(4, i)
          IF (outparms(2, i) == 0.0) THEN
            NXOUT = 1
          ELSE
            NXOUT = MAX(1,NINT(XDOMAIN/outparms(2, i)))
          ENDIF
          STARTY = -outparms(5, i)
          IF (outparms(3, i) == 0.0) THEN
            NYOUT = 1
          ELSE
            NYOUT = MAX(1,NINT(YDOMAIN/outparms(3, i)))
          ENDIF
          Z0 = MIN( MAX(outparms(1, i),ZGRID(1)), ZGRID(NZ))
          NANGOUT = NINT(outparms(6, i))

          status( 1) = nf90_def_dim(ncFileId, trim(radiance_Name("x", numRadianceOutputs)), &
                                    NXOUT, radXDimID)
          status( 2) = nf90_def_dim(ncFileId, trim(radiance_Name("y", numRadianceOutputs)), &
                                    NYOUT, radYDimId)
          status( 3) = nf90_def_dim(ncFileId, trim(radiance_Name("direction", numRadianceOutputs)), &
                                    NANGOUT, radDirDimId)
          status( 4) = nf90_def_var(ncFileId, trim(radiance_Name("x", numRadianceOutputs)), &
                                    nf90_float, radXDimID, ncVarId)
          status( 5) = nf90_def_var(ncFileId, trim(radiance_Name("y", numRadianceOutputs)), &
                                    nf90_float, radYDimId, ncVarId)
          !
          ! Variables and dimension
          !
          status( 7) = nf90_def_var(ncFileId, trim(radiance_Name("mu",  numRadianceOutputs)), &
                                    nf90_float, radDirDimId, ncVarId)
          status( 8) = nf90_def_var(ncFileId, trim(radiance_Name("phi", numRadianceOutputs)), &
                                              nf90_float, radDirDimId, ncVarId)
          
          select case(units) 
            case('T')
              UNITSNAME = 'Kelvin'
            case('B') 
              UNITSNAME = 'Watts/(M^2 Ster)'
            case default
              UNITSNAME = 'Watts/(M^2 Micron Ster)'
          end select
          status( 9) = nf90_def_var(ncFileId, trim(radiance_Name("", numRadianceOutputs)), nf90_float, &
                                      (/radXDimID, radYDimID, radDirDimId/), ncVarId)
          status(10) = nf90_put_att(ncFileId, ncVarId, "Z_level", Z0)
          status(11) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          
        case ('F')  ! Fluxes 
          ! Possible formats: 
          !   1 - flux at top and bottom of medium at grid points
          !   2 - flux at a given level at regular locations (not supported)
          !   3 - domain averaged vertical profile
          !   4 - fluxes at every base grid point
          !   5 - fluxes at every grid point (not supported)
          select case(units) 
            case('T')
              UNITSNAME = 'Kelvin'
            case('B') 
              UNITSNAME = 'Watts/(M^2)'
            case default
              UNITSNAME = 'Watts/(M^2 Micron)'
          end select
          
          select case(NINT(outparms(1, i)))
            case (1) 
              if(F1done) cycle 
              status( 1) = nf90_def_var(ncFileId, "fluxUp_Top",   nf90_float, &
                                         (/ xDimId, yDimId /), ncVarId)
              status( 2) = nf90_put_att(ncFileId, ncVarId, "Z_level", ZGRID(NZ))
              status( 3) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))

              status( 4) = nf90_def_var(ncFileId, "fluxDown_Bottom", nf90_float, &
                                         (/ xDimId, yDimId /), ncVarId)
              status( 5) = nf90_put_att(ncFileId, ncVarId, "Z_level", ZGRID(1))
              status( 6) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              
              if(outputDiffuse) then 
                !
                ! Diffuse component for solar RT
                !
                status( 7) = nf90_def_var(ncFileId, "fluxDown_Diffuse_Bottom", nf90_float, &
                                           (/ xDimId, yDimId /), ncVarId)
                status( 8) = nf90_put_att(ncFileId, ncVarId, "Z_level", ZGRID(1))
                status( 9) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              end if 
              f1Done = .true.
            case (3) 
              if (F3Done) cycle
              status( 1) = nf90_def_var(ncFileId, "fluxUp_DomainMean",   nf90_float, zDimId, ncVarId)
              status( 2) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))

              status( 3) = nf90_def_var(ncFileId, "fluxDown_DomainMean", nf90_float, zDimId, ncVarId)
              status( 4) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              
              if(outputDiffuse) then 
                ! Diffuse component 
                status( 5) = nf90_def_var(ncFileId, "fluxDown_Diffuse_DomainMean", &
                                                                           nf90_float, zDimId, ncVarId)
                status( 6) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              end if 
              F3Done = .true. 
            case (4, 5) 
              if(F4done) cycle 
              status( 1) = nf90_def_var(ncFileId, "fluxUp",   nf90_float, &
                                         (/ xDimId, yDimId, zDimId /), ncVarId)
              status( 2) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))

              status( 3) = nf90_def_var(ncFileId, "fluxDown", nf90_float, &
                                         (/ xDimId, yDimId, zDimId /), ncVarId)
              status( 4) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              
              if(outputDiffuse) then 
                !
                ! Diffuse component for solar RT
                !
                status( 5) = nf90_def_var(ncFileId, "fluxDown_Diffuse", nf90_float, &
                                           (/ xDimId, yDimId, zDimId /), ncVarId)
                status( 6) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              end if 
              F4done = .true. 
            case default 
              ! Not supported or an error 
          end select
          
        case ('H') 
          !
          ! Heating output: net flux convergence 
          !
          select case(units) 
            case('T')
              UNITSNAME = 'Kelvin'
            case('B') 
              UNITSNAME = 'Watts/(M^2 Km)'
            case default
              UNITSNAME = 'Watts/(M^2 Micron Km)'
          end select
          ! Possible formats: 
          !   1 - domain averaged vertical profile
          !   2 - flux convergence for every base grid point
          !   3 - flux convergence for every grid point (not supported) 
          select case(NINT(outparms(1, i))) 
            case (1) 
              if(H1Done) cycle
              status( 1) = nf90_def_var(ncFileId, "fluxDivergence_DomainMean", nf90_float, zDimId, ncVarId)
              status( 2) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              H1Done = .true. 
            case (2, 3) 
              if(H2Done) cycle
              status( 1) = nf90_def_var(ncFileId, "fluxDivergence", nf90_float, &
                                        (/ xDimId, YDimId, zDimId/), ncVarId)
              status( 2) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
              H2Done = .true. 
            case default
          end select
        
        case ('S') 
          if(SDone) cycle 
          ! Spherical Harmonic output: 
          !   Output mean intensity and net flux (x,y,z components).
          !      and maybe normalized rms of higher order terms 
          select case(units) 
            case('T')
              UNITSNAME = 'Kelvin'
            case('B') 
              UNITSNAME = 'Watts/(M^2 Ster)'
            case default
              UNITSNAME = 'Watts/(M^2 Micron)'
          end select
          status( 1) = nf90_def_var(ncFileId, "meanIntensity",   nf90_float, &
                                     (/ xDimId, yDimId, zDimId /), ncVarId)
          status( 2) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          status( 3) = nf90_def_var(ncFileId, "netFlux_xDir",   nf90_float, &
                                     (/ xDimId, yDimId, zDimId /), ncVarId)
          status( 4) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          status( 5) = nf90_def_var(ncFileId, "netFlux_yDir",   nf90_float, &
                                     (/ xDimId, yDimId, zDimId /), ncVarId)
          status( 6) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          status( 7) = nf90_def_var(ncFileId, "netFlux_zDir",   nf90_float, &
                                     (/ xDimId, yDimId, zDimId /), ncVarId)
          status( 8) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          if(nshout == 5) then 
            status( 9) = nf90_def_var(ncFileId, "HigherOrder_SHTerms_RMS",   nf90_float, &
                                     (/ xDimId, yDimId, zDimId /), ncVarId)
            status(10) = nf90_put_att(ncFileId, ncVarId, "units", trim(UNITSNAME))
          end if 
          Sdone = .true. 

        case ('J', 'M')
          ! not supported
        case default 
          ! not supported - but shouldn't ever happen
      end select
      if(any(status(:) /= nf90_NoErr)) print *, "Error setting up output number ", i, status(1:11)
      status(:) = nf90_NoErr
    end do 
    
    if(nf90_enddef(ncFileId) /= nf90_NoErr) print *, "Error ending definition of netcdf file"
    
    ! ------------------------------------------------------
    !
    ! Write the various fields to the netcdf file 
    ! We could keep track of which outputs have been written and avoid re-writing, but we 
    !   won't bother for now. 
    !
     ! Write the base grid positions
    status( 1) = nf90_inq_varid(ncFileId, "x", ncVarId)
    status( 2) = nf90_put_var(ncFileId, ncVarId, XGRID(1:NX)) 
    status( 3) = nf90_inq_varid(ncFileId, "y", ncVarId)
    status( 4) = nf90_put_var(ncFileId, ncVarId, YGRID(1:NY)) 
    status( 5) = nf90_inq_varid(ncFileId, "z", ncVarId)
    status( 6) = nf90_put_var(ncFileId, ncVarId, ZGRID) 
    if(any(status(1:6) /= nf90_NoErr)) print *, "Error writing grid positions for netcdf output", status(1:6)

    status(:) = nf90_noerr
    numRadianceOutputs = 0; radianceStartPos = 1
    do i = 1, numOut
      select case (outtypes(i))
        case ('R') ! Radiance 
          numRadianceOutputs = numRadianceOutputs + 1 
          STARTX = -outparms(4, i)
          IF (outparms(2, i) == 0.0) THEN
            NXOUT = 1
          ELSE
            NXOUT = MAX(1,NINT(XDOMAIN/outparms(2, i)))
          ENDIF
          STARTY = -outparms(5, i)
          IF (outparms(3, i) == 0.0) THEN
            NYOUT = 1
          ELSE
            NYOUT = MAX(1,NINT(YDOMAIN/outparms(3, i)))
          ENDIF
          Z0 = MIN( MAX(outparms(1, i),ZGRID(1)), ZGRID(NZ))
          NANGOUT = NINT(outparms(6, i))
          
          status( 1) = nf90_inq_varid(ncFileId, trim(radiance_Name("x", numRadianceOutputs)), ncVarId)
          status( 2) = nf90_put_var(ncFileId, ncVarId, STARTX + outparms(2, i) * (/ (j, j = 0, NXOUT-1) /)) 

          status( 3) = nf90_inq_varid(ncFileId, trim(radiance_Name("y", numRadianceOutputs)), ncVarId)
          status( 4) = nf90_put_var(ncFileId, ncVarId, STARTY + outparms(3, i) * (/ (j, j = 0, NYOUT-1) /)) 

          status( 5) = nf90_inq_varid(ncFileId, trim(radiance_Name("mu", numRadianceOutputs)), ncVarId)
          status( 6) = nf90_put_var(ncFileId, ncVarId, outparms(5+2:5+NANGOUT*2:2, i))

          status( 7) = nf90_inq_varid(ncFileId, trim(radiance_Name("phi", numRadianceOutputs)), ncVarId)
          status( 8) = nf90_put_var(ncFileId, ncVarId, outparms(6+2:6+NANGOUT*2:2, i))

          status( 9) = nf90_inq_varid(ncFileId, trim(radiance_Name("", numRadianceOutputs)), ncVarId)
          status(10) = nf90_put_var(ncFileId, ncVarId, & 
                                    reshape(RADOUT(radianceStartPos:radianceStartPos+NXOUT*NYOUT*NANGOUT-1), &
                                            (/ NXOUT, NYOUT, NANGOUT /) ) ) 
          radianceStartPos = radianceStartPos + NXOUT*NYOUT*NANGOUT
          
        case ('F')  ! Fluxes 
          select case(NINT(outparms(1, i)))
            case (1) 
              status( 1) = nf90_inq_varid(ncFileId, "fluxUp_Top", ncVarId)
              status( 2) = nf90_put_var(ncFileId, ncVarId, C*transpose(FLUXES(2,NZ,:,:)))
              
              status( 3) = nf90_inq_varid(ncFileId, "fluxDown_Bottom", ncVarId)
              status( 4) = nf90_put_var(ncFileId, ncVarId, C*transpose(FLUXES(1,1,:,:) + FLUXES(3,1,:,:)))
              
              if(outputDiffuse) then 
                status( 5) = nf90_inq_varid(ncFileId, "fluxDown_Diffuse_Bottom", ncVarId)
                status( 6) = nf90_put_var(ncFileId, ncVarId, C*transpose(FLUXES(1,1,:,:)))
              end if 
              
            case (3) 
              status( 1) = nf90_inq_varid(ncFileId, "fluxUp_DomainMean",   ncVarId)
              status( 2) = nf90_put_var(ncFileId, ncVarId, &
                           C*sum(sum(FLUXES(2,:,:,:), dim = 3), dim = 2)/(nx*ny))

              status( 3) = nf90_inq_varid(ncFileId, "fluxDown_DomainMean", ncVarId)
              status( 4) = nf90_put_var(ncFileId, ncVarId, &
                           C*sum(sum(FLUXES(1,:,:,:) + FLUXES(3,:,:,:), dim = 3), dim = 2)/(nx*ny))
             
              if(outputDiffuse) then 
                ! Diffuse component 
                status( 5) = nf90_inq_varid(ncFileId, "fluxDown_Diffuse_DomainMean", ncVarId)
                status( 6) = nf90_put_var(ncFileId, ncVarId, &
                                            C*sum(sum(FLUXES(1,:,:,:), dim = 3), dim = 2)/(nx*ny))
              end if 
           case (4) 
              status( 1) = nf90_inq_varid(ncFileId, "fluxUp",   ncVarId)
              status( 2) = nf90_put_var(ncFileId, ncVarId,         &
                                        C*reshape(FLUXES(2,:,:,:), &
                                                  shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )

              status( 3) = nf90_inq_varid(ncFileId, "fluxDown", ncVarId)
              status( 4) = nf90_put_var(ncFileId, ncVarId,                           &
                                        C*reshape(FLUXES(1,:,:,:) + FLUXES(3,:,:,:), &
                                                  shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
              
              if(outputDiffuse) then 
                status( 5) = nf90_inq_varid(ncFileId, "fluxDown_Diffuse", ncVarId)
                status( 6) = nf90_put_var(ncFileId, ncVarId,         &
                                          C*reshape(FLUXES(1,:,:,:), &
                                                    shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
              end if 
            case default 
              ! Not supported or an error 
              print *, "Flux output format", NINT(outparms(1, i)), "is not supported in netcdf files." 
          end select
          
        case ('H') 
          select case(NINT(outparms(1, i))) 
            case (1) 
              status( 1) = nf90_inq_varid(ncFileId, "fluxDivergence_DomainMean", ncVarId)
              status( 2) = nf90_put_var(ncFileId, ncVarId, &
                                          -sum(sum(FLUXDIV(:,:,:), dim = 3), dim = 2)/(nx*ny))
            case (2) 
              status( 1) = nf90_inq_varid(ncFileId, "fluxDivergence", ncVarId)
              status( 2) = nf90_put_var(ncFileId, ncVarId, &
                                        reshape(-FLUXDIV(:,:,:), &
                                                shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
            case default
              print *, "Heating rate output format", NINT(outparms(1, i)), "is not supported." 
          end select
        
        case ('S') 
          ! Spherical Harmonic output: 
          !   Output mean intensity and net flux (x,y,z components).
          !      and maybe normalized rms of higher order terms 
          status( 1) = nf90_inq_varid(ncFileId, "meanIntensity", ncVarId)
          status( 2) = nf90_put_var(ncFileId, ncVarId,        & 
                                    reshape(SHTERMS(1,:,:,:), &
                                            shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
          
          status( 3) = nf90_inq_varid(ncFileId, "netFlux_xDir", ncVarId)
          status( 4) = nf90_put_var(ncFileId, ncVarId,        &
                                    reshape(SHTERMS(2,:,:,:), &
                                            shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
          
          status( 5) = nf90_inq_varid(ncFileId, "netFlux_yDir", ncVarId)
          status( 6) = nf90_put_var(ncFileId, ncVarId,        &
                                    reshape(SHTERMS(3,:,:,:), &
                                            shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )

          status( 7) = nf90_inq_varid(ncFileId, "netFlux_zDir", ncVarId)
          status( 8) = nf90_put_var(ncFileId, ncVarId,        &
                                    reshape(SHTERMS(4,:,:,:), &
                                            shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
          
          if(nshout == 5) then 
            status( 9) = nf90_inq_varid(ncFileId, "HigherOrder_SHTerms_RMS", ncVarId)
            status(10) = nf90_put_var(ncFileId, ncVarId,      &
                                      reshape(SHTERMS(5,:,:,:)/max(1.e-20,  SHTERMS(1,:,:,:)), &
                                              shape = (/ NX, NY, NZ /), order = (/ 3, 2, 1 /)) )
          end if 

        case ('J', 'M')
          ! not supported
          print *, "Output format " // outtypes(i) // " is not supported in netcdf files." 
        case default 
          ! not supported - but shouldn't ever happen
          print *, "Output format " // outtypes(i) // " is unknown." 
      end select
      if(any(status(:) /= nf90_NoErr)) print *, "Error writing output number", i,status(1:10)
      status(:) = nf90_noerr
    end do 
    
    if(nf90_close(ncFileId) /= nf90_NoErr) print *, "Error closing netcdf file"
    
  end subroutine output_results_netcdf_par
  ! -------------------------------------------------------------------------------------------
  function radiance_Name(var, number)
    character(len=*), intent(in) :: var
    integer,          intent(in) :: number
    character(len=64)            :: radiance_Name 

    radiance_Name = "radiance"
    if(len_trim(var) > 0) radiance_Name = trim(radiance_Name) // "_" // trim(var)
    if (number > 1) radiance_Name = trim(radiance_Name) // "." // trim(IntToChar(number))
  end function radiance_Name 
  ! -------------------------------------------------------------------------------------------
  function IntToChar(integerValue)
    integer, intent( in) :: integerValue
    
    character(len = 25)  :: IntToChar
    !
    !   Creates the character representation of an integer.  
    !  

    write(IntToChar, *) integerValue
    IntToChar = AdjustL(IntToChar)
  end function IntToChar
  ! -------------------------------------------------------------------------------------------
end module shdom_netcdf
