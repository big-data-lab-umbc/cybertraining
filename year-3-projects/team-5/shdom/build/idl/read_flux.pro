; IDL procedure  -  read_flux.pro
;
;****************************************************************************
;
; NAME:
;      read_flux
;
; PURPOSE:
;      Read an SHDOM-produced flux file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_flux, filename, filetype, nx, ny, nz, xgrid, ygrid, zgrid, $
;                                     fluxup, fluxdn, fluxdir
;
; INPUTS:
;      filename:  a string with the directory path and file name
;      filetype:  the type of SHDOM flux file (i.e. 1, 2, or 4)
;                 same as the output format given to SHDOM
;
; OUTPUTS:
;      nx, ny, nz: output grid size  (nz=1 for format 1 and 2)
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      fluxup:     upward flux array 
;                     format 1,2  [fltarr(ny,nx)]
;                     format 4    [fltarr(nz,ny,nx)]
;      fluxdn:     downward flux array (diffuse + direct flux)
;      fluxdir:    direct beam flux for solar case (optional for thermal source)
;
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      For use with files created by the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;      Currently, flux file types 1, 2, and 4 are implemented.
;
;****************************************************************************

PRO read_flux, filename, filetype,  nx, ny, nz, xgrid, ygrid, zgrid, $
                                    fluxup, fluxdn, fluxdir

;****************************************************************************

;  Open file and get the base grid size from the header
ss = ' '
openr,lun, filename, /get_lun
readf,lun, ss
readf,lun, ss
readf,lun, ss
nx = fix(strmid(ss,strpos(ss,'NX=')+3,4))
ny = fix(strmid(ss,strpos(ss,'NY=')+3,4))
nz = fix(strmid(ss,strpos(ss,'NZ=')+3,4))

;  Find the source type and skip rest of header
repeat begin
    readf,lun, ss
endrep until (strpos(ss,'SOURCE_TYPE') ne -1)
source = 's'
if (strpos(ss,'THERMAL') ne -1) then source = 't'


;  flux output type 1 - upward/downward/direct flux on base grid top and bottom
if (filetype eq 1) then begin
    repeat begin
      readf,lun, ss
    endrep until (strpos(ss,'UPWELLING FLUX') ne -1)
    ;  Read the top and bottom Z levels
    nz = 2
    zgrid = fltarr(nz)
    zgrid(1) = float(strmid(ss,strpos(ss,'Z=')+2,7))
    readf,lun, ss
    zgrid(0) = float(strmid(ss,strpos(ss,'Z=')+2,7))
    readf,lun, ss
    ;  Read in the flux grid and reformat    
    n = long(nx)*ny
    if (source eq 't') then buffer=fltarr(4,n) else buffer=fltarr(5,n)
    readf,lun, buffer
    temp=reform(buffer(0,*),ny,nx)
    xgrid=reform(temp(0,*),nx)
    temp=reform(buffer(1,*),ny,nx)
    ygrid=reform(temp(*,0),ny)
    fluxup = reform(buffer(2,*),ny,nx)
    fluxdn = reform(buffer(3,*),ny,nx)
    if (source eq 's') then begin
        fluxdir = reform(buffer(4,*),ny,nx)
        fluxdn = fluxdn + fluxdir
    endif
endif


;  flux output type 2 - upward/downward/direct flux at one level
if (filetype eq 2) then begin
    repeat begin
      readf,lun, ss
    endrep until (strpos(ss,'UPWELLING FLUX') ne -1)
    ;  Read the Z level and number of output locations
    nz = 1
    zgrid = fltarr(nz)
    zgrid(0) = float(strmid(ss,strpos(ss,'Z=')+2,7))
    nx = fix(strmid(ss,strpos(ss,'NXO=')+4,4))
    ny = fix(strmid(ss,strpos(ss,'NYO=')+4,4))
    readf,lun, ss
    readf,lun, ss

    ;  Read in the flux grid and reformat    
    n = long(nx)*long(ny)*nz
    if (source eq 't') then buffer=fltarr(4,n) else buffer=fltarr(5,n)
    readf,lun, buffer
    temp=reform(buffer(0,*),nx,ny)
    xgrid=reform(temp(*,0),nx)
    temp=reform(buffer(1,*),nx,ny)
    ygrid=reform(temp(0,*),ny)
    fluxup = transpose(reform(buffer(2,*),nx,ny))
    fluxdn = transpose(reform(buffer(3,*),nx,ny))
    if (source eq 's') then begin
        fluxdir = transpose(reform(buffer(4,*),nx,ny))
        fluxdn = fluxdn + fluxdir
    endif
endif


;  flux output type 4 - upward/downward/direct flux at every base grid point
if (filetype eq 4) then begin
    repeat begin
      readf,lun, ss
    endrep until (strpos(ss,'UP') ne -1) and (strpos(ss,'DOWN') ne -1)
    n = long(nx)*long(ny)*nz
    if (source eq 't') then buffer=fltarr(5,n) else buffer=fltarr(6,n)
    readf,lun, buffer
    temp=reform(buffer(0,*),nz,ny,nx)
    xgrid=reform(temp(0,0,*),nx)
    temp=reform(buffer(1,*),nz,ny,nx)
    ygrid=reform(temp(0,*,0),ny)
    temp=reform(buffer(2,*),nz,ny,nx)
    zgrid=reform(temp(*,0,0),nz)
    fluxup = reform(buffer(3,*),nz,ny,nx)
    fluxdn = reform(buffer(4,*),nz,ny,nx)
    if (source eq 's') then begin
        fluxdir = reform(buffer(5,*),nz,ny,nx)
        fluxdn = fluxdn + fluxdir
    endif
endif


free_lun, lun

END

