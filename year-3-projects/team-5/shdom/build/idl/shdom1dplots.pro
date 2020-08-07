; Example IDL script for making 1D line graphs of 
;   SHDOM radiative transfer results.
;   Examples use solar results from script runrtmonoles.

;   The following plot types are available:
; plottype='lwptau'    liquid water path and optical depth vs X
; plottype='flux'      upwelling and downwelling flux vs X
; plottype='netflux'   X and Z net flux vs X
; plottype='radiance'  radiance for several directions vs X
plottype='lwptau'


; This script supports 3 output formats: 
;   X windows (out='x'), Postscript ('ps'), Encapsulated Postscript ('eps')
out='ps'
psfile='shdom1dplot.'+out

if (out eq 'ps') then begin
  ;  For Postscript use device fonts (pick Times Roman)
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /times, $
   	 xsize=18, ysize=24, xoffset=1.5, yoffset=2.0
  !p.font=0
  !p.charsize=1.2
endif
if (out eq 'eps') then begin
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /encapsulated, /times, $
   	xsize=9.0, ysize=12.0, xoffset=1.0, yoffset=1.0
  !p.font=0
  !p.charsize=0.7
endif
if (out eq 'x') then begin
  ;  For X windows use the Hershey stroked fonts
  set_plot, 'x'
  window, xsize=600, ysize=800, title='SHDOM Plots'
  !p.font=-1
  !p.charsize=1.2
endif

; Set other global plotting variables
!p.multi = [0,1,2]
!p.thick=1.0  &  !x.thick=1.0  &  !y.thick=1.0
!x.ticklen = 0.013  &  !y.ticklen = 0.010




;  Plot liquid water path and optical depth vs X
;      Reads LWC and property files
if (plottype eq 'lwptau') then begin
  lwcfile='les2y21.lwc'     ; LWC file 
  prpfile='les2y21w16.prp'  ; SHDOM property file to plot
  iy=0                      ; Y index slice (0 for 2D field)
  title1='LES FIRE Liquid Water Path'    ; title for first plot
  title2='LES FIRE Optical Depth'        ; title for second plot

  read_lwc, lwcfile, nx,ny,nz,xgrid1,ygrid,zgrid, temp,lwc,reff
  lwp = 1000*vertical_integral(zgrid, lwc, 3)

  read_prp, prpfile, nx,ny,nz, xgrid2,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)


  ;  Plot LWP vs X
                     ;  if X windows pick a nice Hershey font for title
  if (out eq 'x') then title1='!17'+title1
  plot, xgrid1, lwp(*,iy), /xstyle, /ystyle, $
	xtitle='X (km)', ytitle='LWP (g/m!u2!n)', $
	title=title1

  ;  Plot optical depth vs X
  plot, xgrid2, tau(*,iy), /xstyle, /ystyle, $
	xtitle='X (km)', ytitle='Optical Depth', $
	title=title2

endif






;  Plot upwelling and downwelling flux vs X
;       Reads SHDOM flux format 1, 2, or 4 files
if (plottype eq 'flux') then begin
  prpfile='les2y21w16.prp'     ; SHDOM property file
  fluxfile='les2y21w16af2.out' ; SHDOM flux output file
  fluxform=2                   ; flux file format (1,2,4)
  iz=22                        ; Z level index for format 4
  iy=0                         ; Y index slice (0 for 2D field)
  taumin=0.0                   ; minimum optical depth
  taumax=12.                   ; maximum optical depth
  fluxmin=0.00                 ; minimum flux
  fluxmax=1.20                 ; maximum flux
  title1='LES FIRE Optical Depth'    ; title for first plot
  title2='LES FIRE Fluxes'           ; title for second plot

  read_prp, prpfile, nx,ny,nz, xgrid1,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)

  read_flux, fluxfile,fluxform, nx,ny,nz,xgrid2,ygrid,zgrid, fluxup,fluxdn,fdir
  if (fluxform eq 4) then begin
    fup = reform(fluxup(iz,iy,*))
    fdn = reform(fluxdn(iz,iy,*))
    zlevdn = zgrid(iz)
    zlevup = zgrid(iz)
  endif else begin
    fup = reform(fluxup(iy,*))
    fdn = reform(fluxdn(iy,*))
    zlevdn = zgrid(0)
    zlevup = zgrid(0)
    if (fluxform eq 1) then zlevup = zgrid(1)
  endelse

  ;  Plot optical depth vs X
  if (out eq 'x') then title1='!17'+title1
  plot, xgrid1, tau(*,iy), xstyle=1, yrange=[taumin,taumax], $
	xtitle='X (km)', ytitle='Optical Depth', $
	title=title1

  ;  Plot upwelling and downwelling flux vs X
  plot, xgrid2, fup, /xstyle, yrange=[fluxmin,fluxmax], linestyle=0, $
	xtitle='X (km)', ytitle='Flux (relative)', $
	title=title2
  oplot, xgrid2, fdn, linestyle=1
  ;  Plot legend
  x0=xgrid2(0) &  dx=xgrid2(nx-1)-xgrid2(0)
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  & x3=x0+0.10*dx
  y0=fluxmin   &  dy=fluxmax-fluxmin  & ey=0.02*dy
  y1=y0+0.11*dy  &  y2=y0+0.05*dy
  oplot, [x1,x2], [y1,y1], linestyle=0
  xyouts, x3, y1-ey, 'Upwelling flux at Z=' $
	+string(zlevup,format='(F4.2)')+' km', charsize=1.0*!p.charsize
  oplot, [x1,x2], [y2,y2], linestyle=1
  xyouts, x3, y2-ey, 'Downwelling flux at Z=' $
	+string(zlevdn,format='(F4.2)')+' km', charsize=1.0*!p.charsize

endif






;  Plot X and Z net flux at one level vs X
;       Reads SHDOM spherical harmonic file
if (plottype eq 'netflux') then begin
  prpfile='les2y21w16.prp'     ; SHDOM property file
  shfile='les2y21w16as.out'    ; SHDOM SH (net flux) output file
  iz=18                        ; Z level index
  iy=0                         ; Y index slice (0 for 2D field)
  taumin=0.0                   ; minimum optical depth
  taumax=12.                   ; maximum optical depth
  netfluxmin=-2.00             ; minimum net flux
  netfluxmax=2.00              ; maximum net flux
  title1='LES FIRE Optical Depth'    ; title for first plot
  title2='LES FIRE Net Fluxes'       ; title for second plot

  read_prp, prpfile, nx,ny,nz, xgrid1,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)

  read_sh, shfile, nx,ny,nz,xgrid2,ygrid,zgrid, meanrad, fx,fy,fz
  fx = reform(fx(iz,iy,*))
  fz = reform(fz(iz,iy,*))
  zlevdn = zgrid(iz)
  zlevup = zgrid(iz)

  ;  Plot optical depth vs X
  if (out eq 'x') then title1='!17'+title1
  plot, xgrid1, tau(*,iy), xstyle=1, yrange=[taumin,taumax], $
	xtitle='X (km)', ytitle='Optical Depth', $
	title=title1

  ;  Plot X and Z net flux vs X
  plot, xgrid2, fz, /xstyle, yrange=[netfluxmin,netfluxmax], linestyle=0, $
	xtitle='X (km)', ytitle='Flux (relative)', $
	title=title2
  oplot, xgrid2, fx, linestyle=1
  ;  Plot legend
  x0=xgrid2(0) &  dx=xgrid2(nx-1)-xgrid2(0)
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  & x3=x0+0.10*dx
  y0=netfluxmin  &  dy=netfluxmax-netfluxmin  & ey=0.02*dy
  y1=y0+0.11*dy  &  y2=y0+0.05*dy
  oplot, [x1,x2], [y1,y1], linestyle=1
  xyouts, x3, y1-ey, 'X net flux at Z=' $
	+string(zlevup,format='(F4.2)')+' km', charsize=1.0*!p.charsize
  oplot, [x1,x2], [y2,y2], linestyle=0
  xyouts, x3, y2-ey, 'Z net flux at Z=' $
	+string(zlevdn,format='(F4.2)')+' km', charsize=1.0*!p.charsize

endif






;  Plot radiances for several directions vs X
;       Reads SHDOM radiance file
if (plottype eq 'radiance') then begin
  prpfile='les2y21w16.prp'     ; SHDOM property file
  radfile='les2y21w16ar.out'   ; SHDOM radiance output file
  iy=0                         ; Y index slice (0 for 2D field)
  idir=[6,0,3]                 ; indices for directions
  taumin=0.0                   ; minimum optical depth
  taumax=12.                   ; maximum optical depth
  radmin=-0.10                 ; minimum radiance
  radmax=0.30                  ; maximum radiance
  title1='LES FIRE Optical Depth'         ; title for first plot
  title2='LES FIRE Upwelling Radiances'   ; title for second plot

  read_prp, prpfile, nx,ny,nz, xgrid1,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)

  read_rad, radfile, nx,ny,ndir, xgrid2,ygrid,zlev, mu, phi, radiance
  nrad=n_elements(idir)
  rad=fltarr(nx,nrad)
  for i = 0, nrad-1 do $
    rad(*,i) = radiance(*,iy,idir(i))

  ;  Plot optical depth vs X
  if (out eq 'x') then title1='!17'+title1
  plot, xgrid1, tau(*,iy), xstyle=1, yrange=[taumin,taumax], $
	xtitle='X (km)', ytitle='Optical Depth', $
	title=title1

  ;  Plot for the radiances vs X
  plot, xgrid2, rad(*,0), /xstyle, yrange=[radmin,radmax], linestyle=0, $
	xtitle='X (km)', ytitle='Radiance', $
	title=title2
  for i = 1, nrad-1 do $
    oplot, xgrid2, rad(*,i), linestyle=i

  ;  Plot legend
  if (out eq 'x') then begin
    smu='!7l!17'  &  sphi='!7u!17'
  endif else begin
    smu='!9m!7'   &  sphi='!9f!7'
  endelse
  x0=xgrid2(0)   &  dx=xgrid2(nx-1)-xgrid2(0)
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  &  x3=x0+0.10*dx
  y0=radmin      &  dy=radmax-radmin  & ey=0.02*dy
  for i = 0, nrad-1 do begin
    y1=y0+(0.06*(nrad-i)-0.01)*dy
    oplot, [x1,x2], [y1,y1], linestyle=i
    xyouts, x3, y1-ey, 'Direction: ' $
	+smu+'='+string(mu(idir(i)),format='(F5.3)')+'  ' $
	+sphi+'='+string(fix(phi(idir(i))),format='(I3)'), $
	charsize=1.0*!p.charsize
  endfor

endif



if (out ne 'x') then device, /close

end
