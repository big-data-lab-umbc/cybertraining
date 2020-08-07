; IDL procedure  -  read_heat.pro
;
;****************************************************************************
;
; NAME:
;      read_heat
;
; PURPOSE:
;      Read an SHDOM-produced net flux convergence (heating rate) file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_heat, filename, nx, ny, nz, xgrid, ygrid, zgrid, fluxconv
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      fluxconv:   net flux convergence array [fltarr(nz,ny,nx)]
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      For use with files created by the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;      Reads heating rate file type 2 and 3.
;
;****************************************************************************

PRO read_heat, filename, nx, ny, nz, xgrid, ygrid, zgrid, fluxconv

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

repeat begin
  readf,lun, ss
endrep until (strpos(ss,'-DIV') ne -1)
n = long(nx)*long(ny)*nz
buffer=fltarr(4,n)
readf,lun, buffer
temp=reform(buffer(0,*),nz,ny,nx)
xgrid=reform(temp(0,0,*),nx)
temp=reform(buffer(1,*),nz,ny,nx)
ygrid=reform(temp(0,*,0),ny)
temp=reform(buffer(2,*),nz,ny,nx)
zgrid=reform(temp(*,0,0),nz)
fluxconv = reform(buffer(3,*),nz,ny,nx)

free_lun, lun

END

