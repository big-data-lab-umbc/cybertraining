
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

;;;;;;;;;;;;;;;;;;;


ff1=addfile("/home/hsong/zzbatmos_user/ncl/data/uwnd.mon.mean.nc","r")
u850=ff1->uwnd(:,2,:,:)
alon=ff1->lon
alat=ff1->lat

ff2=addfile("/home/hsong/zzbatmos_user/ncl/data/vwnd.mon.mean.nc","r")
v850=ff2->vwnd(:,2,:,:)

ff3=addfile("/home/hsong/zzbatmos_user/ncl/data/slp.mon.mean.nc","r")
slp=ff3->slp(:,:,:)


;; Data start from Jan 1948
;;ENSO positive years: 1982-1983, 1997-1998; 2015-2016
u8501=(dim_avg_n(u850(420-2:420,:,:),0)+dim_avg_n(u850(600-2:600,:,:),0)+dim_avg_n(u850(816-2:816,:,:),0))/3.
v8501=(dim_avg_n(v850(420-2:420,:,:),0)+dim_avg_n(v850(600-2:600,:,:),0)+dim_avg_n(v850(816-2:816,:,:),0))/3.
ps1=(dim_avg_n(slp(420-2:420,:,:),0)+dim_avg_n(slp(600-2:600,:,:),0)+dim_avg_n(slp(816-2:816,:,:),0))/3.


;;ENSO negative years: 1999-2000, 2007-2008; 2010-2011
u8502=(dim_avg_n(u850(624-2:624,:,:),0)+dim_avg_n(u850(720-2:720,:,:),0)+dim_avg_n(u850(756-2:756,:,:),0))/3.
v8502=(dim_avg_n(v850(624-2:624,:,:),0)+dim_avg_n(v850(720-2:720,:,:),0)+dim_avg_n(v850(756-2:756,:,:),0))/3.
ps2=(dim_avg_n(slp(624-2:624,:,:),0)+dim_avg_n(slp(720-2:720,:,:),0)+dim_avg_n(slp(756-2:756,:,:),0))/3.

;;;;;;
f0=addfile("/home/hsong/zzbatmos_common/Data/CAM5.3_CLUBB_Po-Lun_2016_July/for_zhibo/cam53_f19_clubb.cam.h0.2009-01.nc","r")
land=f0->LANDFRAC(0,:,:)
;swcf_obs=mask(swcf_obs,(land.lt.1),1)

u8501!0 = "lat"
u8501&lat = alat
u8501!1 = "lon"
u8501&lon = alon
u8501@units=" "

v8501!0 = "lat"
v8501&lat = alat
v8501!1 = "lon"
v8501&lon = alon
v8501@units=" "

ps1!0 = "lat"
ps1&lat = alat
ps1!1 = "lon"
ps1&lon = alon
ps1@units=" "

u8502!0 = "lat"
u8502&lat = alat
u8502!1 = "lon"
u8502&lon = alon
u8502@units=" "

v8502!0 = "lat"
v8502&lat = alat
v8502!1 = "lon"
v8502&lon = alon
v8502@units=" "

ps2!0 = "lat"
ps2&lat = alat
ps2!1 = "lon"
ps2&lon = alon
ps2@units=" "

;printVarSummary(ps2)
printVarSummary(u8502)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;plot the figures;;;;;

wks=gsn_open_wks("png","plot_uv850_NCEP_NCAR_ENSO")

plot=new(7,graphic)

res=True
;res@gsnDraw=False
res@gsnFrame=False

cmap = read_colormap_file("BlAqGrYeOrReVi200")

;res@mpLandFillColor="white"

res@tiMainOn=True
res@tiMainOffsetYF=-0.005

res@tmYLMinorOn=False
res@tmXBMinorOn=False
res@tmYROn=False
res@tmXTOn=False

res@cnFillOn=True
res@cnLinesOn=False
res@cnLineLabelsOn=False
res@gsnScalarContour     = True               ; contours desired
res@cnFillPalette=cmap(17:,:)

res@lbLabelBarOn=True
res@lbOrientation="vertical"

res@cnLevelSelectionMode="ManualLevels"
res@cnMinLevelValF=1000
res@cnMaxLevelValF=1020
res@cnLevelSpacingF=1
;gsn_define_colormap(wks,"matlab_jet")
;gsn_define_colormap(wks,"amwg")
;res@cnLabelBarEndStyle="ExcludeOuterBoxes"

;res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=4
;res@gsnSpreadColorEnd=136
res@lbLabelStride=2

res@mpMinLatF=-60
res@mpMaxLatF=60
;res@tmYLMode="Explicit"
;res@tmYLValues=(/30,35,40,45,50/)
;res@tmYLLabels=(/"30N","35N","40N","45N","50N"/)
;res@mpShapeMode = "FreeAspect"

res@mpMinLonF=0
res@mpMaxLonF=360
res@mpCenterLonF=180
;res@tmXBMode="Explicit"
;res@tmXBValues=(/320,325,330,335,340,345,350/)
;res@tmXBLabels=(/"40W","35W","30W","25W","20W","15W","10W"/)

res@vcRefMagnitudeF           = 10.0             ; define vector ref mag
;res@vcRefLengthF              = 0.04           ; define length of vec ref
res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
res@vcRefLengthF            = 0.045            ; define length of vec ref
res@vcLineArrowColor        = "black"           ; change vector color
;res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
;res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
res@vcMinDistanceF          = 0.017

;res@pmLabelBarWidthF      = 0.085               ; default is shorter
;res@pmLabelBarHeightF     = 0.15               ; default is taller

;;;plot1

res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25

;;;plot1
res@tiMainString="NCEP/NCAR 850hPa Wind in El Nino Years"
;plot1=gsn_csm_vector_scalar_map_ce(wks,u8501,v8501,ps1,res)
plot1 = gsn_csm_vector_map(wks,u8501,v8501,res)

frame(wks)

;;;plot2
res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25

res@tiMainString="NCEP/NCAR 850hPa Wind in La Nina Years"
;plot2=gsn_csm_vector_scalar_map_ce(wks,u8502,v8502,ps2,res)
plot2 = gsn_csm_vector_map(wks,u8502,v8502,res)

;-- polymarker resources
pmres                        =  True
pmres@gsMarkerColor          = "blue"       ;-- marker color
pmres@gsMarkerSizeF          =  0.025        ;-- set size of marker
pmres@gsLineThicknessF       =  5.          ;-- marker line thickness
pmres@gsMarkerThicknessF=2.5

pmres@gsMarkerIndex = 15
;gsn_polymarker(wks, plot1, 332, 39, pmres) 
;gsn_polymarker(wks, plot2, 332, 39, pmres)     
                                

frame(wks)

end
