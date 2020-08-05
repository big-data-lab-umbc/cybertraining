; IDL procedure  -  read_broad.pro
;
;****************************************************************************
;
; NAME:
;      read_broad
;
; PURPOSE:
;      Read an SHDOM broadband flux and heating rate file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_broad, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
;                                     fluxup, fluxdn, heating
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      fluxup:     upward flux array (W/m^2) [fltarr(nz,ny,nx)]
;      fluxdn:     downward flux array (W/m^2)  [fltarr(nz,ny,nx)]
;      heating:    heating rate array (K/hr)  [fltarr(nz,ny,nx)]
;
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      Broadband files are created by summing SHDOM k-distribution output
;      files (e.g. with the runrtkdistles script).
;
;****************************************************************************

PRO read_broad, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
                          fluxup, fluxdn, heating

;****************************************************************************

;  Open file and get the base grid size from the header
ss = ' '
openr,lun, filename, /get_lun
readf,lun, ss
nx=1  &  ny=1  &  nz=1
readf,lun, nx,ny,nz
readf,lun, ss
readf,lun, ss

;  Read in all the data and reformat into the grid arrays
n = long(nx)*long(ny)*nz
buffer=fltarr(6,n)
readf,lun, buffer
temp=reform(buffer(0,*),nz,ny,nx)
xgrid=reform(temp(0,0,*),nx)
temp=reform(buffer(1,*),nz,ny,nx)
ygrid=reform(temp(0,*,0),ny)
temp=reform(buffer(2,*),nz,ny,nx)
zgrid=reform(temp(*,0,0),nz)

fluxup = reform(buffer(3,*),nz,ny,nx)
fluxdn = reform(buffer(4,*),nz,ny,nx)
heating = reform(buffer(5,*),nz,ny,nx)

free_lun, lun

END

