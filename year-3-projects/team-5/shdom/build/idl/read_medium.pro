; IDL procedure  -  read_medium.pro
;
;****************************************************************************
;
; NAME:
;      read_medium
;
; PURPOSE:
;      Read an SHDOM-produced medium property file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_medium, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
;                         extinct, ssalbedo, asymmetry, temperature
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      extinct:    extinction array  [fltarr(nz,ny,nx)]
;      ssalbedo:   single scattering albedo array  [fltarr(nz,ny,nx)]
;      asymmetry:  asymmetry parameter array  [fltarr(nz,ny,nx)]
;      temperature: temperature (K) array  [fltarr(nz,ny,nx)]
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      For use with files created by the Spherical Harmonics Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;      The SHDOM output medium properties may have been delta-M scaled
;      and hence are not the same as the input property values.
;
;****************************************************************************

PRO read_medium, filename, nx, ny, nz, xgrid, ygrid, zgrid, $
                           extinct, ssalbedo, asymmetry, temperature


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
endrep until (strpos(ss,'Albedo') ne -1)

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
extinct = reform(buffer(3,*),nz,ny,nx)
ssalbedo = reform(buffer(4,*),nz,ny,nx)
asymmetry = reform(buffer(5,*),nz,ny,nx)
temperature = reform(buffer(6,*),nz,ny,nx)

free_lun, lun

END

