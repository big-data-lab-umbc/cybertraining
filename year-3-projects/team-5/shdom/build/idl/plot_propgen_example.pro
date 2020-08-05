; IDL script for plotting "run_propgen_example" files.

;   The following plot types are available:
; plottype='xz_lwcre'   X-Z  liquid water content and effective radius
; plottype='xz_extalb'  X-Z  extinction and single scattering albedo
; plottype='radiance'  radiance for several directions vs X
plottype='xz_lwcre'


; This script supports 3 output formats: 
;   X windows (out='x'), Postscript ('ps'), Encapsulated Postscript ('eps')
out='ps'
psfile='propgen_example.'+out
nlev=20  ; number of color levels

if (out eq 'ps') then begin
  ;  For Postscript use device fonts (pick Times Roman)
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /times, $
   	 xsize=18, ysize=24, xoffset=0.0, yoffset=0.0
;   	 xsize=18, ysize=24, xoffset=1.5, yoffset=2.0
  cola=(!d.n_colors/nlev)*indgen(nlev)+1  &  !p.color = 0
  !p.font=0
  !p.charsize=1.2
endif
if (out eq 'eps') then begin
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /encapsulated, /times, $
   	xsize=9.0, ysize=12.0, xoffset=1.0, yoffset=1.0
  cola=(!d.n_colors/nlev)*indgen(nlev)+1  &  !p.color = 0
  !p.font=0
  !p.charsize=0.7
endif
if (out eq 'x') then begin
  ;  For X windows use the Hershey stroked fonts
  set_plot, 'x'
  window, xsize=600, ysize=800, title='SHDOM Plots'
  cola=(!d.n_colors/nlev)*indgen(nlev)+1  &  !p.color = 255
  !p.font=-1
  !p.charsize=1.2
endif

;  Select the color table
;loadct, 0    ; gray scale
loadct, 39    ; blue to red color

; Set other global plotting variables
!p.thick=1.0  &  !x.thick=1.0  &  !y.thick=1.0
!x.ticklen = -0.013  &  !y.ticklen = -0.010




;  Plot X-Z slice of liquid water content and effective radius for one
;   particle type (component).    Reads a particle properties file.
if (plottype eq 'xz_lwcre') then begin
  parfile='nauru19990707.part'   ; particle file 
  partype=2                      ; particle type to extract from file
  zmin=0.0  &  zmax=15.0         ; height range to plot
  if (partype eq 1) then begin
    parname='Aerosol'            ; particle type name
    lwcmax=2.0E-4                ; max liquid water content
    lwcform='(E8.1)'
    reffmin=0.5 & reffmax=2.5    ; min and max effective radius
  endif
  if (partype eq 2) then begin
    parname='Water droplet'            ; particle type name
    lwcmax=0.4                   ; max liquid water content
    lwcform='(F4.2)'
    reffmin=4.0 & reffmax=16.0   ; min and max effective radius
  endif
  if (partype eq 3) then begin
    parname='Plate ice crystal'  ; particle type name
    lwcmax=0.02                  ; max liquid water content
    lwcform='(F5.3)'
    reffmin=0.001 & reffmax=80.0    ; min and max effective radius
  endif
  if (partype eq 4) then begin
    parname='Aggregate crystal'  ; particle type name
    lwcmax=0.04                  ; max liquid water content
    lwcform='(F5.3)'
    reffmin=0.001 & reffmax=80.0    ; min and max effective radius
  endif
  iy=0                         ; Y index slice (0 for 2D field)
  title='Nauru LWC/r!Deff!N example: '+parname   ; title for plot

  read_particle, parfile, partype, nx,ny,nz,xgrid,ygrid,zgrid, temp,lwc,reff

  plotdata1 = reform(lwc(*,iy,*))
  plotdata2 = reform(reff(*,iy,*))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.95, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.70
  y1plot=0.47+0.18
  delyplot = 0.25

  ;  Make color contour plot with color bar of LWC field
  contlev=0.0+((lwcmax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = (0.999*lwcmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format=lwcform), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Liquid Water Content (g/m!U3!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of effective radius field
  y1plot=0.18
  contlev=reffmin+((reffmax-reffmin)/nlev)*findgen(nlev+1)
  plotdata2 = 0.999*reffmin > (0.999*reffmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Effective Radius (!9m!7m)', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif





;  Plot X-Z slice of extinction and single scattering albedo 
;      Reads SHDOM input property file
if (plottype eq 'xz_extalb') then begin
;  prpfile='nauru19990707_w065.prp'  ; SHDOM property file to plot
  prpfile='nauru19990707_w164.prp'  ; SHDOM property file to plot
  zmin=0.0  &  zmax=15.0            ; height range to plot
  extinctmin=-2.0                   ; minimum log10 extinction to plot
  extinctmax=2.0                    ; maximum log10 extinction to plot
  ssalbmin=0.80                     ; minimum single scattering albedo
  ssalbmax=1.00                     ; maximum single scattering albedo
  iy=0                              ; Y index slice (0 for 2D field)
;  title='Nauru MMCR example: 0.65 !9m!7m'   ; title for plot
  title='Nauru MMCR example: 1.64 !9m!7m'   ; title for plot

  ;  Read in property file
  read_prp, prpfile, nx,ny,nz, xgrid,ygrid,zgrid, temps,extinct,ssalb
  plotdata1 = alog10(1.0E-4 > reform(extinct(*,iy,*)))
  plotdata2 = reform(ssalb(*,iy,*))

  ;  if X windows pick a nice Hershey font and plot title
  if (out eq 'x') then title='!17'+title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot: X and Z are to scale if 4:3 aspect ratio output
  x1plot=0.15
  delxplot=0.70
  y1plot=0.47+0.18
  delyplot = 0.25

  ;  Make color contour plot with color bar of extinction field
  contlev=extinctmin+((extinctmax-extinctmin)/nlev)*findgen(nlev+1)
;  plotdata1 = 1.001*extinctmax > (0.999*extinctmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F4.1)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Log!D10!N Extinction (km!U-1!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of single scattering albedo
  y1plot=0.18
  contlev=ssalbmin+((ssalbmax-ssalbmin)/nlev)*findgen(nlev+1)
  plotdata2 = 1.001*ssalbmin > (0.999*ssalbmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F4.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Single scattering albedo', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif






;  Plot radiances for several directions vs X
;       Reads SHDOM radiance file
if (plottype eq 'radiance') then begin
  prpfile='nauru19990707_w065.prp'    ; SHDOM property file
  radfile='nauru19990707_w065.radout' ; SHDOM radiance output file
  iy=0                         ; Y index slice (0 for 2D field)
  idir=[0,2,4]                 ; indices for directions
  taumin=0.0                   ; minimum optical depth
  taumax=50.                   ; maximum optical depth
  radmin=0.0                   ; minimum radiance
  radmax=0.25                  ; maximum radiance
  title1='Nauru MMCR Example 0.65 !9m!7m Optical Depth '      ; title for plot
  title2='Nauru MMCR Example 0.65 !9m!7m Radiances at 17 km'  ; title for plot

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
	position=[0.12,0.57,0.97,0.95], title=title1

  ;  Plot for the radiances vs X
  cols=FIX([0.25,0.60,0.95]*!d.n_colors)
  plot, xgrid2, rad(*,0), /xstyle, yrange=[radmin,radmax], $
	xtitle='X (km)', ytitle='Radiance', /noerase, /nodata, $
	position=[0.12,0.07,0.97,0.45], title=title2
  for i = 0, nrad-1 do $
    oplot, xgrid2, rad(*,i), linestyle=0, color=cols[i]

  ;  Plot legend
  if (out eq 'x') then begin
    smu='!7l!17'  &  sphi='!7u!17'
  endif else begin
    smu='!9m!7'   &  sphi='!9f!7'
  endelse
  dx=xgrid2[nx-1]-xgrid2[0]
  x0=xgrid2[0]+0.5*dx
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  &  x3=x0+0.10*dx
  y0=radmax      &  dy=radmax-radmin  & ey=0.02*dy
  for i = 0, nrad-1 do begin
    y1=y0-(0.06*i+0.05)*dy
    oplot, [x1,x2], [y1,y1], linestyle=0, color=cols[i]
    xyouts, x3, y1-ey, 'Direction: ' $
	+smu+'='+string(mu(idir(i)),format='(F5.3)')+'  ' $
	+sphi+'='+string(fix(phi(idir(i))),format='(I3)'), $
	charsize=1.0*!p.charsize
  endfor
endif


if (out ne 'x') then device, /close

end
