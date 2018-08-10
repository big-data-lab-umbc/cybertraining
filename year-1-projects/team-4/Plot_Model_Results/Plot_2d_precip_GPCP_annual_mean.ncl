
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

all_files=systemfunc("ls /home/hsong/zzbatmos_user/ncl/data/precip.mon.mean.nc")
ff=addfiles(all_files,"r")
ListSetType(ff,"cat")

PRECT=ff[:]->precip

printVarSummary(PRECT)

;;;;;
precip=dim_avg_n(PRECT(0:467,:,:),0)
precip@_FillValue=-9.96921e+36

;;;;;;;;;;;;
alon=ff[0]->lon
alat=ff[0]->lat

precip!0 = "lat"
precip&lat = alat
precip!1 = "lon"
precip&lon = alon


;;;;;;
;;;;;;plot the spatial distributions of precipitation

wks=gsn_open_wks("png","/home/hsong/zzbatmos_user/ncl/plot/cybertraining/plot_2d_precip_GPCP_annual_mean")

plot=new(3,graphic)

colors1=(/"white","black","White"/)
gsn_define_colormap(wks,colors1)
gsn_merge_colormaps(wks,colors1,"Matlab_jet")

res=True
;res@gsnDraw=False
res@gsnFrame=False

res@tiMainOn=True
res@tiMainOffsetYF=-0.005

res@tmYLMinorOn=False
res@tmXBMinorOn=False
res@tmYROn=False
res@tmXTOn=False

res@cnLinesOn=False
res@cnLineLabelsOn=False
res@cnFillOn=True
;res@cnFillPalette="MPL_coolwarm"

res@lbLabelBarOn=True
res@lbOrientation="vertical"

res@cnLabelBarEndStyle="ExcludeOuterBoxes"
res@cnLevelSelectionMode:="ManualLevels"
res@cnMinLevelValF:=0
res@cnMaxLevelValF:=12
res@cnLevelSpacingF:=0.5
;gsn_define_colormap(wks,"MPL_seismic")
res@cnLabelBarEndStyle="ExcludeOuterBoxes"
;gsn_define_colormap(wks,"amwg")

res@pmLabelBarWidthF      = 0.085               ; default is shorter
;res@pmLabelBarHeightF     = 0.15               ; default is taller

res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=4
;res@gsnSpreadColorEnd=136
res@lbLabelStride=4

res@mpMinLatF=-60
res@mpMaxLatF=60

;;;Plot1
res@vpXF=0.1
res@vpWidthF=0.775
res@vpYF=0.9
res@vpHeightF=0.25

res@mpMinLonF=0
res@mpMaxLonF=360
res@mpCenterLonF=180

res@tiMainString="GPCP Precipitation [mm/day]"
plot(0)=gsn_csm_contour_map(wks,precip,res)

;;;Plot2
res@vpXF=0.1
res@vpWidthF=0.85
res@vpYF=0.68
res@vpHeightF=0.2

res@tiMainString="Stratiform Precipitation [mm/day]"
;plot(1)=gsn_csm_contour_map(wks,precip_ls,res)

res@lbLabelBarOn=True
res@lbLabelBarOn=True
res@lbLabelBarOn=True
res@pmLabelBarOrthogonalPosF=0.125
;res@pmLabelBarParallelPosF=1.120
res@lbLabelFontHeightF=0.008
res@pmLabelBarWidthF=0.5
res@pmLabelBarHeightF=0.06

;;;Plot3
res@vpXF=0.1
res@vpWidthF=0.85
res@vpYF=0.41
res@vpHeightF=0.2

res@tiMainString="Total Precipitation [mm/day]"
;plot(2)=gsn_csm_contour_map(wks,precip_tt,res)

a1=avg(precip)
;a2=avg(precip_ls)
a3=avg(precip)

;;;;add some texts
txres=True

txres@txFontHeightF=0.008
txres@txJust="CenterLeft"
;gsn_text_ndc(wks,"Mean="+sprintf("%4.2f",a1),0.7,0.91,txres)
;gsn_text_ndc(wks,"Mean="+sprintf("%4.2f",a2),0.7,0.69,txres)
;gsn_text_ndc(wks,"Mean="+sprintf("%4.2f",a3),0.7,0.42,txres)

txres@txFontHeightF=0.015
txres@txJust="CenterCenter"
;gsn_text_ndc(wks,"Annual Mean",0.525,0.11,txres)

frame(wks)


end
