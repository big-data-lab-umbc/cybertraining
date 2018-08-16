
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

all_files=systemfunc("ls /home/hsong/zzbatmos_user/ncl/data/precip.mon.mean.nc")
ff=addfiles(all_files,"r")
ListSetType(ff,"cat")

PRECT=ff[:]->precip(0:467,:,:)
PRECT@_FillValue=-9.96921e+36
printVarSummary(PRECT)

ENSO1=asciiread("ENSO_index_1950_2017.txt",-1,"float")
ENSO=ENSO1(348:815)
printVarSummary(ENSO)

ccr=escorc(ENSO,PRECT(lat|:,lon|:,time|:))
printVarSummary(ccr)
;;;;;

;;;;;;;;;;;;
alon=ff[0]->lon
alat=ff[0]->lat

ccr!0="lat"
ccr&lat = alat
ccr!1 = "lon"
ccr&lon = alon


;;;;;;
;;;;;;plot the spatial distributions of precipitation

wks=gsn_open_wks("png","plot_2d_precip_ENSO_correlation_mean")

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
res@cnMinLevelValF:=-1.0
res@cnMaxLevelValF:=1.0
res@cnLevelSpacingF:=0.05
;gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
res@cnLabelBarEndStyle="ExcludeOuterBoxes"
gsn_define_colormap(wks,"BlWhRe")
gsn_define_colormap(wks,"BlRe")

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

res@tiMainString="Correlation between ENSO-Index and Precipitation"
plot(0)=gsn_csm_contour_map(wks,ccr,res)


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
