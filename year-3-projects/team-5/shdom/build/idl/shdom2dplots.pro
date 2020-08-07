; Example IDL script for making 2D color contour plots of 
;   SHDOM radiative transfer results.
;   Examples use solar results from runrtmonoles, runrtkdistles, and run3dgaus.

;   The following plot types are available:
; plottype='xz_extalb'  X-Z  extinction and single scattering albedo
; plottype='xz_flux'    X-Z  upwelling and downwelling fluxes
; plottype='xz_lwcheat' X-Z  liquid water content and heating rate
; plottype='xz_extIm'   X-Z  extinction and mean radiance/flux vectors
; plottype='xy_lwp'     X-Y  liquid water path
; plottype='xy_tau'     X-Y  optical depth
; plottype='xy_rad'     X-Y  radiance at one angle
; plottype='xy_flux'    X-Y  up and down flux at one level
; plottype='xy_taufcnv' X-Y  tau and vertically integrated net flux convergence
plottype='xz_extalb'


; This script supports 3 output formats: 
;   X windows (out='x'), Postscript ('ps'), Encapsulated Postscript ('eps')
out='ps'
psfile='shdom2dplot.'+out
nlev=20  ; number of color levels

if (out eq 'ps') then begin
  ;  For Postscript use device fonts (pick Times Roman)
  set_plot, 'ps'            ; use 4:3 ratio in plot size (cm)
  device, filename=psfile, /color, /portrait, /times, $
   	 xsize=18, ysize=24, xoffset=1.5, yoffset=2.0
  cola=(240.0/nlev)*indgen(nlev)+1  &  !p.color = 0
  !p.font=0
  !p.charsize=1.2
endif
if (out eq 'eps') then begin
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /encapsulated, /times, $
   	xsize=9.0, ysize=12.0, xoffset=1.0, yoffset=1.0
  cola=(240.0/nlev)*indgen(nlev)+1  &  !p.color = 0
  !p.font=0
  !p.charsize=0.7
endif
if (out eq 'x') then begin
  ;  For X windows use the Hershey stroked fonts
  set_plot, 'x'
  window, xsize=600, ysize=800, title='SHDOM Plots'
  cola=(180.0/nlev)*indgen(nlev)+1  &  !p.color = 255
  !p.font=-1
  !p.charsize=1.2
endif

;  Select the color table
;loadct, 0    ; gray scale
loadct, 39    ; blue to red color

; Set other global plotting variables
!p.thick=1.0  &  !x.thick=1.0  &  !y.thick=1.0
!x.ticklen = -0.013  &  !y.ticklen = -0.010




;  Plot X-Z slice of extinction and single scattering albedo 
;      Reads SHDOM input property file
if (plottype eq 'xz_extalb') then begin
  prpfile='les2y21w16.prp'  ; SHDOM property file to plot
  zmin=0.0  &  zmax=1.0     ; height range to plot
  extinctmax=100.           ; maximum extinction to plot
  ssalbmin=0.98             ; minimum single scattering albedo
  ssalbmax=1.00             ; maximum single scattering albedo
  iy=0                      ; Y index slice (0 for 2D field)
  title='LES FIRE Optical Properties'    ; title for plot

  ;  Read in property file
  read_prp, prpfile, nx,ny,nz, xgrid,ygrid,zgrid, temps,extinct,ssalb
  plotdata1 = reform(extinct(*,iy,*))
  plotdata2 = reform(ssalb(*,iy,*))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot: X and Z are to scale if 4:3 aspect ratio output
  x1plot=0.15
  delxplot=0.70
  y1plot=0.50+0.18
  delyplot = (3./4.)*delxplot*(zmax-zmin)/(xgrid(nx-1)-xgrid(0))
  if (delyplot gt 0.25) then begin
    delxplot = delxplot*0.25/delyplot
    x1plot=0.50-delxplot/2
    delyplot = 0.25
  endif

  ;  Make color contour plot with color bar of extinction field
  contlev=0.0+((extinctmax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = 0.1 > (0.999*extinctmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.0)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Extinction (km!U-1!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of single scattering albedo
  y1plot=0.18
  contlev=ssalbmin+((ssalbmax-ssalbmin)/nlev)*findgen(nlev+1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.3)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Single scattering albedo', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif





;  Plot X-Z slice of upwelling and downwelling fluxes
;       Reads SHDOM flux format 4 or broadband file
if (plottype eq 'xz_flux') then begin
  fluxfile='les2y21sw.out'  ; flux output file 
  ;fluxfile='les2y21u1f.out' ; SHDOM flux output file 
  broadband=1               ; 1 for broadband file, 0 for SHDOM flux file
  zmin=0.0  &  zmax=1.0     ; height range to plot
  fluxupmin=0.              ; upward flux min
  fluxupmax=400.            ; upward flux max
  fluxdnmin=200.            ; downward flux min
  fluxdnmax=600.            ; downward flux max
  iy=0                      ; Y index slice (0 for 2D field)
  title='LES FIRE Broadband Flux'    ; title for plot

  if (broadband eq 1) then $
    read_broad, fluxfile, nx,ny,nz,xgrid,ygrid,zgrid, fluxup,fluxdn,heat $
  else $
    read_flux, fluxfile,4, nx,ny,nz,xgrid,ygrid,zgrid, fluxup,fluxdn,fdir

  plotdata1 = transpose(reform(fluxup(*,iy,*)))
  plotdata2 = transpose(reform(fluxdn(*,iy,*)))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.70
  y1plot=0.50+0.18
  delyplot = (3./4.)*delxplot*(zmax-zmin)/(xgrid(nx-1)-xgrid(0))
  if (delyplot gt 0.25) then begin
    delxplot = delxplot*0.25/delyplot
    x1plot=0.50-delxplot/2
    delyplot = 0.25
  endif

  ;  Make color contour plot with color bar of upward flux field
  contlev=fluxupmin+((fluxupmax-fluxupmin)/nlev)*findgen(nlev+1)
  plotdata1 = fluxupmin > (0.999*fluxupmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.0)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Upwelling Flux (W/m!U2!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of downward flux field
  y1plot=0.18
  contlev=fluxdnmin+((fluxdnmax-fluxdnmin)/nlev)*findgen(nlev+1)
  plotdata2 = fluxdnmin > (0.999*fluxdnmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.0)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Downwelling Flux (W/m!U2!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif






;  Plot X-Z slice of liquid water content and heating rate
;       Reads LWC file and broadband files
if (plottype eq 'xz_lwcheat') then begin
  lwcfile='les2y21.lwc'     ; LWC file 
  heatfile='les2y21sw.out'  ; broadband output file 
  zmin=0.0  &  zmax=1.0     ; height range to plot
  lwcmax=1.0                ; max liquid water content
  heatmax=2.0               ; max heating rate
  iy=0                      ; Y index slice (0 for 2D field)
  title='LES FIRE LWC and Broadband Solar Heating Rate'    ; title for plot

  read_lwc, lwcfile, nx,ny,nz,xgrid,ygrid,zgrid1, temp,lwc,reff
  read_broad, heatfile, nx,ny,nz,xgrid,ygrid,zgrid2, fluxup,fluxdn,heat

  plotdata1 = reform(lwc(*,iy,*))
  plotdata2 = transpose(reform(heat(*,iy,*)))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.70
  y1plot=0.50+0.18
  delyplot = (3./4.)*delxplot*(zmax-zmin)/(xgrid(nx-1)-xgrid(0))
  if (delyplot gt 0.25) then begin
    delxplot = delxplot*0.25/delyplot
    x1plot=0.50-delxplot/2
    delyplot = 0.25
  endif

  ;  Make color contour plot with color bar of LWC field
  contlev=0.0+((lwcmax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = (0.999*lwcmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid1, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Liquid Water Content (g/m!U3!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of heating rate field
  y1plot=0.18
  contlev=0.0+((heatmax-0.0)/nlev)*findgen(nlev+1)
  plotdata2 = (0.999*heatmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid2, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Heating Rate (K/hr)', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif







;  Plot X-Z slice of extinction and mean radiance with net flux vectors
;       Reads medium property and SH file
if (plottype eq 'xz_extIm') then begin
  medfile='les2y21w16am.out' ; medium property file 
  shfile='les2y21w16as.out'  ; spherical harmonic (net flux) output file 
  zmin=0.0  &  zmax=1.0      ; height range to plot
  extinctmax=60.             ; max extinction
  meanradmax=0.32            ; max mean radiance
  iy=0                       ; Y index slice (0 for 2D field)
  title='LES FIRE Extinction and Mean Radiance'    ; title for plot

  read_medium, medfile, nx,ny,nz,xgrid,ygrid,zgrid, extinct,ssalb,asym,temp
  read_sh, shfile, nx,ny,nz,xgrid,ygrid,zgrid, meanrad, fx,fy,fz

  plotdata1 = transpose(reform(extinct(*,iy,*)))
  plotdata2 = transpose(reform(meanrad(*,iy,*)))
  fx = transpose(reform(fx(*,iy,*)))
  fz = transpose(reform(fz(*,iy,*)))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.70
  y1plot=0.50+0.18
  delyplot = (3./4.)*delxplot*(zmax-zmin)/(xgrid(nx-1)-xgrid(0))
  if (delyplot gt 0.25) then begin
    delxplot = delxplot*0.25/delyplot
    x1plot=0.50-delxplot/2
    delyplot = 0.25
  endif

  ;  Make color contour plot with color bar of extinction field
  contlev=0.0+((extinctmax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = (0.999*extinctmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.1)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Scaled Extinction (km!U-1!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of mean radiance field
  y1plot=0.18
  contlev=0.0+((meanradmax-0.0)/nlev)*findgen(nlev+1)
  plotdata2 = (0.999*meanradmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Mean Radiance/Net Flux Vectors', size=1.2*!p.charsize, $
	/normal, alignment=0.5

  ; Do overplot of net flux vectors that have been subsampled
  isub=4  &  arrowlen=1.0 &  arrowhead=0.4
  nxo=isub*fix(nx/isub)  &  nzo=isub*fix(nz/isub)
  sxo=0                  &  szo=nz-nzo
  exo=nxo+sxo-1          &  ezo=nzo+szo-1
  ovector, rebin(fx(sxo:exo,szo:ezo),nxo/isub,nzo/isub), $
	 rebin(fz(sxo:exo,szo:ezo),nxo/isub,nzo/isub), $
	 rebin(xgrid(sxo:exo),nxo/isub), rebin(zgrid(szo:ezo),nzo/isub), $
	 length=arrowlen, hsize=arrowhead

endif







;  Plot X-Y map of liquid water path
;       Reads LWC file and does vertical integral
if (plottype eq 'xy_lwp') then begin
  lwcfile='gaussol.lwc'     ; LWC file 
  lwpmax=10.                ; max liquid water path
  title='3D Gaussian Liquid Water Path'    ; title for plot

  read_lwc, lwcfile, nx,ny,nz,xgrid,ygrid,zgrid, temp,lwc,reff
  lwp = 1000*vertical_integral(zgrid, lwc, 3)

  plotdata = lwp

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.90, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.75
  y1plot=0.25
  delyplot = (3./4.)*delxplot

  ;  Make color contour plot with color bar of LWP field
  contlev=0.0+((lwpmax-0.0)/nlev)*findgen(nlev+1)
  plotdata = (0.999*lwpmax < plotdata)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.1)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.16, 'Liquid Water Path (g/m!U2!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif






;  Plot X-Y map of optical depth
;       Reads property file and does vertical integral
if (plottype eq 'xy_tau') then begin
  prpfile='gaussol.prp'     ; property file 
  taumax=2.0                ; max optical depth
  title='3D Gaussian Optical Depth'    ; title for plot

  read_prp, prpfile, nx,ny,nz, xgrid,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)

  plotdata = tau

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.90, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.75
  y1plot=0.25
  delyplot = (3./4.)*delxplot

  ;  Make color contour plot with color bar of optical depth field
  contlev=0.0+((taumax-0.0)/nlev)*findgen(nlev+1)
  plotdata = (0.999*taumax < plotdata)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.1)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.16, 'Optical Depth', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif






;  Plot X-Y map of radiance for one direction
;       Reads SHDOM radiance output file
if (plottype eq 'xy_rad') then begin
  radfile='gaussolr.out'    ; property file 
  radmin=0.0                ; min radiance
  radmax=0.04               ; max radiance
  idir=2                    ; direction number (0,1,...)
  title='3D Gaussian Radiance'    ; title for plot

  read_rad, radfile, nx,ny,ndir, xgrid,ygrid,zlev, mu, phi, radiance
  plotdata = radiance(*,*,idir)
  if (out eq 'x') then begin
    direction = 'Viewing angle: !7l!17='+string(mu(idir),format='(F5.3)') $
               +'  !7u!17='+string(fix(phi(idir)),format='(I3)')
  endif else begin
    direction = 'Viewing angle: !9m!7='+string(mu(idir),format='(F5.3)') $
               +'  !9f!7='+string(fix(phi(idir)),format='(I3)')
  endelse

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.90, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.75
  y1plot=0.25
  delyplot = (3./4.)*delxplot

  ;  Make color contour plot with color bar of optical depth field
  contlev=radmin+((radmax-radmin)/nlev)*findgen(nlev+1)
  plotdata = radmin > (0.999*radmax < plotdata)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title=direction, xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.16, 'Radiance', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif





;  Plot X-Y slice of upwelling and downwelling fluxes
;       Reads SHDOM flux format 1, 2, or 4
if (plottype eq 'xy_flux') then begin
  fluxfile='gaussolf1.out'  ; SHDOM flux output file
  fluxform=1                ; flux file format (1,2,4)
  iz=6                      ; Z level index for format 4
  fluxupmin=0.              ; upward flux min
  fluxupmax=0.08            ; upward flux max
  fluxdnmin=0.6             ; downward flux min
  fluxdnmax=1.0             ; downward flux max
  title='3D Gaussian'    ; title for plot

  read_flux, fluxfile,fluxform, nx,ny,nz,xgrid,ygrid,zgrid, fluxup,fluxdn,fdir
  if (fluxform eq 4) then begin
    plotdata1 = transpose(reform(fluxup(iz,*,*)))
    plotdata2 = transpose(reform(fluxdn(iz,*,*)))
    zlevdn = zgrid(iz)
    zlevup = zgrid(iz)
  endif else begin
    plotdata1 = transpose(fluxup)
    plotdata2 = transpose(fluxdn)
    zlevdn = zgrid(0)
    zlevup = zgrid(0)
    if (fluxform eq 1) then zlevup = zgrid(1)
  endelse

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.97, title, /normal, size=1.2*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.25
  delxplot=0.50
  y1plot=0.50+0.05
  delyplot = (3./4.)*delxplot

  ;  Make color contour plot with color bar of upward flux field
  contlev=fluxupmin+((fluxupmax-fluxupmin)/nlev)*findgen(nlev+1)
  plotdata1 = fluxupmin > (0.999*fluxupmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='Upwelling Flux at Z='+string(zlevup,format='(F4.2)')+' km', $
	xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot+delxplot+0.05,y1plot, 0.03,delyplot, 0.070,-0.01, /orientation, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize


  ;  Make color contour plot with color bar of downward flux field
  y1plot=0.07
  contlev=fluxdnmin+((fluxdnmax-fluxdnmin)/nlev)*findgen(nlev+1)
  plotdata2 = fluxdnmin > (0.999*fluxdnmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='Downwelling Flux at Z='+string(zlevdn,format='(F4.2)')+' km', $
	xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot+delxplot+0.05,y1plot, 0.03,delyplot, 0.070,-0.01, /orientation, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize

endif







;  Plot X-Y maps of optical depth and vertically integrated net flux convergence
;       Reads property file and does vertical integral
if (plottype eq 'xy_taufcnv') then begin
  prpfile='gaussol.prp'     ; property file 
  heatfile='gaussolh.out'    ; SHDOM net flux convergence output file
  taumax=2.0                ; max optical depth
  fconvmin=0.0              ; min net flux convergence
  fconvmax=0.02             ; max net flux convergence
  title='3D Gaussian'    ; title for plot

  read_prp, prpfile, nx,ny,nz, xgrid,ygrid,zgrid, temps,extinct,ssalb
  tau = vertical_integral(zgrid, extinct, 3)
  read_heat, heatfile, nx,ny,nz, xgrid,ygrid,zgrid, fluxconv
  fconv = transpose(vertical_integral(zgrid, fluxconv, 1))
  plotdata1 = tau
  plotdata2 = fconv

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.97, title, /normal, size=1.2*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.25
  delxplot=0.50
  y1plot=0.50+0.05
  delyplot = (3./4.)*delxplot

  ;  Make color contour plot with color bar of optical depth field
  contlev=0.0+((taumax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = (0.999*taumax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='Optical Depth', $
	xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot+delxplot+0.05,y1plot, 0.03,delyplot, 0.070,-0.01, /orientation, $
        labels=string(lablev,format='(F5.1)'), colors=cola, charsize=!p.charsize


  ;  Make color contour plot with color bar of net flux convergence field
  y1plot=0.07
  contlev=fconvmin+((fconvmax-fconvmin)/nlev)*findgen(nlev+1)
  plotdata2 = fconvmin+1.0E-4 > (fconvmax-1.0E-4 < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, ygrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, /ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='Integrated Net Flux Convergence', $
	xtitle='X (km)', ytitle='Y (km)'
  makekey, x1plot+delxplot+0.05,y1plot, 0.03,delyplot, 0.070,-0.01, /orientation, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize

endif



if (out ne 'x') then device, /close

end
