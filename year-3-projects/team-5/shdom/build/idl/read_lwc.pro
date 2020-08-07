; IDL procedure  -  read_lwc.pro
;
;****************************************************************************
;
; NAME:
;      read_lwc
;
; PURPOSE:
;      Read a liquid water content file (input to cloudprp)
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_lwc, filename, nx, ny, nz, xgrid, ygrid, zgrid, temp, lwc, reff
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nx, ny, nz: output grid size
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zgrid:      Z grid locations  [fltarr(nz)]
;      temp:       1D array of temperature profile  [fltarr(nz)]
;      lwc:        liquid water content array  [fltarr(nx,ny,nz)]
;      reff:       effective radius array  [fltarr(nx,ny,nz)]
;                         (optional, used if LWC type is 2)
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, May 1999
;
; NOTES:
;      Does not assume any particular order to the LWC points.
;      See companion.doc file for LWC format input to cloudprp.
;
;****************************************************************************

PRO read_lwc, filename, nx, ny, nz, xgrid, ygrid, zgrid, temp, lwc, reff

;****************************************************************************

;  Open file and read the grid size and spacing and temperature profile
;    from the header
openr,lun, filename, /get_lun
lwctype=1
nx=1  &  ny=1  &  nz=1
delx=1.0  &  dely=1.0
readf,lun, lwctype
readf,lun, nx,ny,nz
readf,lun, delx, dely
xgrid=delx*findgen(nx)
ygrid=dely*findgen(ny)
zgrid=fltarr(nz)
temp=fltarr(nz)
readf,lun, zgrid
readf,lun, temp

;  Read in the LWC values one at a time
lwc=fltarr(nx,ny,nz)
if (lwctype eq 2) then reff=fltarr(nx,ny,nz)
ix=1  &  iy=1  &  iz=1
while not EOF(lun) do begin
  tlwc=0.0  &  treff=0.0
  if (lwctype eq 2) then readf,lun, ix,iy,iz, tlwc, treff $
    else readf,lun, ix,iy,iz, tlwc
  lwc(ix-1,iy-1,iz-1) = tlwc
  if (lwctype eq 2) then reff(ix-1,iy-1,iz-1) = treff
endwhile   

free_lun, lun

END

