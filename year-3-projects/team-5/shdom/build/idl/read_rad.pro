; IDL procedure  -  read_rad.pro
;
;****************************************************************************
;
; NAME:
;      read_rad
;
; PURPOSE:
;      Read an SHDOM-produced radiance file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_rad, filename, nx, ny, ndir, xgrid, ygrid, zlev, mu, phi, radiance
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny:     output grid size
;      ndir:       number of output radiance directions
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zlev:       height of Z level
;      mu:         array of mu directions (cosine zenith angle)
;      phi:        array of phi directions (azimuth in degrees)
;      radiance:   radiance array [fltarr(nx,ny,ndir)]
;
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      For use with files created by the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;
;****************************************************************************

PRO read_rad, filename, nx, ny, ndir, xgrid, ygrid, zlev, mu, phi, radiance

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

;  Find the Z level and number of radiance outputs
repeat begin
  readf,lun, ss
endrep until (strpos(ss,'RADIANCE AT') ne -1)
zlev = float(strmid(ss,strpos(ss,'Z=')+2,7))
nx = fix(strmid(ss,strpos(ss,'NXO=')+4,4))
ny = fix(strmid(ss,strpos(ss,'NYO=')+4,4))
ndir = fix(strmid(ss,strpos(ss,'NDIR=')+5,4))
readf,lun, ss

;  Read each block (direction) of radiances 
mu = fltarr(ndir)
phi = fltarr(ndir)
radiance = fltarr(nx,ny,ndir)
n = long(nx)*long(ny)
buffer = fltarr(3,n)
for i = 0, ndir-1 do begin
    readf,lun, ss
    mu(i) = float(strmid(ss,3,8))
    phi(i) = float(strmid(ss,12,6))
    readf,lun, buffer
    temp=reform(buffer(0,*),nx,ny)
    xgrid=reform(temp(*,0),nx)
    temp=reform(buffer(1,*),nx,ny)
    ygrid=reform(temp(0,*),ny)
    radiance(*,*,i) = reform(buffer(2,*),nx,ny)
endfor

free_lun, lun

END

