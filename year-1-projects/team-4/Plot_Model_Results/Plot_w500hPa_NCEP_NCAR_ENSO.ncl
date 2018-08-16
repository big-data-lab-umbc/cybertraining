
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

;;;;;;;;;;;;;;;;;;;
;;;;; NCEP/NCAR ;;;;;;;;

ff1=addfile("/home/hsong/zzbatmos_user/ncl/data/omega.mon.mean.nc","r")
w500=ff1->omega(:,5,:,:)*864.0

alon=ff1->lon
alat=ff1->lat
;gw=ff1>gw

;; Data start from Jan 1948 
;;ENSO positive years: 1982-1983, 1997-1998; 2015-2016
w500_el=(dim_avg_n(w500(420-2:420,:,:),0)+dim_avg_n(w500(600-2:600,:,:),0)+dim_avg_n(w500(816-2:816,:,:),0))/3.

;;ENSO negative years: 1999-2000, 2007-2008; 2010-2011
w500_la=(dim_avg_n(w500(624-2:624,:,:),0)+dim_avg_n(w500(720-2:720,:,:),0)+dim_avg_n(w500(756-2:756,:,:),0))/3.

w500_mean=(dim_avg_n(w500(10:718:12,:,:),0)+dim_avg_n(w500(11:719:12,:,:),0)+dim_avg_n(w500(12:720:12,:,:),0))/3.
w500_el=w500_el-w500_mean
w500_la=w500_la-w500_mean

;;;;;;
f0=addfile("/home/hsong/zzbatmos_common/Data/CAM5.3_CLUBB_Po-Lun_2016_July/for_zhibo/cam53_f19_clubb.cam.h0.2009-01.nc","r")
land=f0->LANDFRAC(0,:,:)
;swcf_obs=mask(swcf_obs,(land.lt.1),1)

w500_el!0 = "lat"
w500_el&lat = alat
w500_el!1 = "lon"
w500_el&lon = alon
w500_el@units=" "

w500_la!0 = "lat"
w500_la&lat = alat
w500_la!1 = "lon"
w500_la&lon = alon
w500_la@units=" "

w500_mean!0 = "lat"
w500_mean&lat = alat
w500_mean!1 = "lon"
w500_mean&lon = alon
w500_mean@units=" "

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;plot the figures;;;;;

wks=gsn_open_wks("png","plot_w500_NCEP_NCAR_ENSO")

plot=new(7,graphic)

res=True
;res@gsnDraw=False
res@gsnFrame=False

cmap = read_colormap_file("BlAqGrYeOrReVi200")

res@mpLandFillColor="white"

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
;res@cnFillPalette=cmap(17:,:)
res@cnFillPalette="BlRe"
;res@cnFillPalette="ViBlGrWhYeOrRe"

res@lbLabelBarOn=True
res@lbOrientation="vertical"

res@cnLevelSelectionMode="ManualLevels"
res@cnMinLevelValF=-80
res@cnMaxLevelValF=80
res@cnLevelSpacingF=5
;gsn_define_colormap(wks,"matlab_jet")
;gsn_define_colormap(wks,"amwg")
res@cnLabelBarEndStyle="ExcludeOuterBoxes"

res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=4
;res@gsnSpreadColorEnd=136
res@lbLabelStride=4

res@mpMinLatF=-60
res@mpMaxLatF=60
;res@tmYLMode="Explicit"
;res@tmYLValues=(/30,35,40,45,50/)
;res@tmYLLabels=(/"30N","35N","40N","45N","50N"/)
;;res@mpShapeMode = "FreeAspect"

res@mpMinLonF=0
res@mpMaxLonF=360
res@mpCenterLonF=180

res@pmLabelBarWidthF      = 0.085               ; default is shorter
;res@pmLabelBarHeightF     = 0.15               ; default is taller

res@mpLandFillColor="white"

;;;plot1
res@vpXF=0.05
res@vpWidthF=0.325
res@vpYF=0.9
res@vpHeightF=0.3

res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25
res@tiMainString="NCEP/NCAR 500hPa Omega in El Nino Years"
plot1=gsn_csm_contour_map(wks,w500_el,res)

frame(wks)

;;;plot2
res@vpXF=0.505
res@vpWidthF=0.325
res@vpYF=0.9
res@vpHeightF=0.3

res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25
res@tiMainString="NCEP/NCAR 500hPa Omega in La Nina Years"
plot2=gsn_csm_contour_map(wks,w500_la,res)


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

;;;plot1
res@cnLevelSelectionMode="ManualLevels"
res@cnMinLevelValF=-100
res@cnMaxLevelValF=100
res@cnLevelSpacingF=10
;gsn_define_colormap(wks,"matlab_jet")
;gsn_define_colormap(wks,"amwg")
res@cnLabelBarEndStyle="ExcludeOuterBoxes"

;res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=4
;res@gsnSpreadColorEnd=136
res@lbLabelStride=2

res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25
res@tiMainString="NCEP/NCAR 500hPa Omega Climatology"
plot1=gsn_csm_contour_map(wks,w500_mean,res)

frame(wks)

end
