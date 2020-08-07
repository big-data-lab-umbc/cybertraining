; IDL procedure  -  read_particle.pro
;
;****************************************************************************
;
; NAME:
;      read_particle
;
; PURPOSE:
;      Read a particle property file (input to propgen)
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_particle, filename, partype, nx, ny, nz, xgrid, ygrid, zgrid, 
;                     temp, lwc, reff
;
; INPUTS:
;      filename:  a string with the directory path and file name
;      partype:   particle type number of component to extract
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
;      Written by:   Frank Evans, May 2003
;
; NOTES:
;      Does not assume any particular order to the points.
;      See companion.doc file for particle file format input to propgen.
;
;****************************************************************************

PRO read_particle, filename,partype, nx,ny,nz, xgrid,ygrid,zgrid, temp,lwc,reff

;****************************************************************************

;  Open file and read the grid size and spacing and temperature profile
;    from the header
openr,lun, filename, /get_lun
lwctype=1
nx=1  &  ny=1  &  nz=1
delx=1.0  &  dely=1.0
readf,lun, lwctype
if (lwctype ne 3) then $
  print, 'read_particle: Only reads type 3 files. Use read_lwc for types 1 and 2'
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
reff=fltarr(nx,ny,nz)
ncomp=1
tempstr=' '
ix=1  &  iy=1  &  iz=1
while not EOF(lun) do begin
  readf,lun, tempstr
  reads,tempstr, ix,iy,iz,ncomp
  temp=fltarr(3*ncomp)
  reads,tempstr, ix,iy,iz,ncomp,temp
  for i = 0, ncomp-1 do begin
    if (FIX(temp[3*i]+0.5) eq partype) then begin
      lwc[ix-1,iy-1,iz-1] = temp[3*i+1]
      reff[ix-1,iy-1,iz-1] = temp[3*i+2]
    endif
  endfor
endwhile   

free_lun, lun

END

