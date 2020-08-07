; IDL procedure  -  read_prp.pro
;
;****************************************************************************
;
; NAME:
;      read_prp
;
; PURPOSE:
;      Read an SHDOM input property file
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_prp, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
;                          temps, extinct, ssalbedo
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      temps:      temperature array [fltarr(nx,ny,nz)]
;      extinct:    extinction array  [fltarr(nx,ny,nz)]
;      ssalbedo:   single scattering albedo [fltarr(nx,ny,nz)]
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      Does not assume any particular order to the points.
;      
;      The tabulated phase functions must have 200 Legendre terms
;        per line as cloudprp outputs (because IDL can't handle more 
;        than 2048 characters per line).
;      The extinction only format phase function must be on one line.
;
;****************************************************************************

PRO read_prp, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
                        temps, extinct, ssalbedo

;****************************************************************************

;  Open file and determine the property file format
openr,lun, filename, /get_lun
exstr=" "
readf,lun, exstr
prptype = strmid(exstr,0,1)
nx=1 & ny=1 & nz=1  &  delx=0.0  &  dely=0.0

    ;  Read the grid size, Z levels, and make the X and Y grids
if (prptype ne 'E' and prptype ne 'T') then point_lun,lun,0
readf,lun, nx, ny, nz
zgrid = fltarr(nz)
readf,lun, delx, dely, zgrid
xgrid = delx*findgen(nx)
ygrid = dely*findgen(ny)

;  This section for extinction only format
if (prptype eq 'E') then begin
    ;  Read the temperature profile and single scattering albedo
    temp = fltarr(nz)
    readf,lun, temp
    ssalb = 0.0
    junk = strarr(1)
    readf,lun, ssalb, junk
    
    ;  Read in the extinction values one at a time
    temps=fltarr(nx,ny,nz)
    extinct=fltarr(nx,ny,nz)
    ssalbedo=fltarr(nx,ny,nz)
    ix=1  &  iy=1  &  iz=1  & ext=0.0
    while not EOF(lun) do begin
      ext=0.0
      if (ny eq 1) then readf,lun, ix,iz,ext  $
        else readf,lun, ix,iy,iz,ext
      extinct(ix-1,iy-1,iz-1) = ext
      ssalbedo(ix-1,iy-1,iz-1) = ssalb
      temps(ix-1,iy-1,iz-1) = temp(iz-1)
    endwhile   
endif $


;  This section for tabulated phase function format
else if (prptype eq 'T') then begin
    ;  Skip over the phase functions
    numphase=1  &  nleg=1
    readf,lun, numphase
    for i = 0, numphase-1 do begin
      readf,lun, nleg
      nl = fix((nleg-1)/200.)
      if (nl gt 0) then begin
        junk = strarr(nl)
        readf,lun, junk
      endif
    endfor
    
    ;  Read in the grid points one at a time
    temps=fltarr(nx,ny,nz)
    extinct=fltarr(nx,ny,nz)
    ssalbedo=fltarr(nx,ny,nz)
    ix=1  &  iy=1  &  iz=1  &  temp=0.0  &  ext=0.0  & ssalb=0.0
    while not EOF(lun) do begin
      temp=0.0  &  ext=0.0  &  ssalb=0.0
      readf,lun, ix,iy,iz, temp,ext,ssalb
      temps(ix-1,iy-1,iz-1) = temp
      extinct(ix-1,iy-1,iz-1) = ext
      ssalbedo(ix-1,iy-1,iz-1) = ssalb
    endwhile   
endif $


;  This section for standard format
else begin
    ;  Read in the grid points one at a time
    temps=fltarr(nx,ny,nz)
    extinct=fltarr(nx,ny,nz)
    ssalbedo=fltarr(nx,ny,nz)
    ix=1  &  iy=1  &  iz=1  &  temp=0.0  &  ext=0.0  & ssalb=0.0
    while not EOF(lun) do begin
      temp=0.0  &  ext=0.0  &  ssalb=0.0
      readf,lun, ix,iy,iz, temp,ext,ssalb
      temps(ix-1,iy-1,iz-1) = temp
      extinct(ix-1,iy-1,iz-1) = ext
      ssalbedo(ix-1,iy-1,iz-1) = ssalb
    endwhile   
endelse

free_lun,lun

END

