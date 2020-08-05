; IDL procedure  -  read_sh.pro
;
;****************************************************************************
;
; NAME:
;      read_sh
;
; PURPOSE:
;      Read an SHDOM-produced spherical harmonic (net fluxes) file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_sh, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
;                         meanrad, xnetflux, ynetflux, znetflux
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      meanrad:    mean radiance array  [fltarr(nz,ny,nx)]
;      xnetflux:   net flux in X direction  [fltarr(nz,ny,nx)]
;      ynetflux:   net flux in Y direction  [fltarr(nz,ny,nx)]
;      znetflux:   net flux in Z direction  [fltarr(nz,ny,nx)]
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      For use with files created by the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;****************************************************************************

PRO read_sh, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
                       meanrad, xnetflux, ynetflux, znetflux

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

;  Skip the rest of the header
repeat begin
  readf,lun, ss
endrep until (strpos(ss,'Imean') ne -1)

;  Read in the data and reformat to the output arrays
n = long(nx)*long(ny)*nz
buffer=fltarr(7,n)
readf,lun, buffer
temp=reform(buffer(0,*),nz,ny,nx)
xgrid=reform(temp(0,0,*),nx)
temp=reform(buffer(1,*),nz,ny,nx)
ygrid=reform(temp(0,*,0),ny)
temp=reform(buffer(2,*),nz,ny,nx)
zgrid=reform(temp(*,0,0),nz)
meanrad = reform(buffer(3,*),nz,ny,nx)
xnetflux = reform(buffer(4,*),nz,ny,nx)
ynetflux = reform(buffer(5,*),nz,ny,nx)
znetflux = reform(buffer(6,*),nz,ny,nx)

free_lun, lun

END

