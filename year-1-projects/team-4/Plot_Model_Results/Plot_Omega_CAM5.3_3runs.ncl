
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

;;;;;;;;;;;
;;;;CAM5.3_ctrl
all_files01:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/ctrl_data/camrun.cam.h0.0001-*.nc")
;all_files02:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/ctrl_data/camrun.cam.h0.0000-06*.nc")
;all_files03:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/ctrl_data/camrun.cam.h0.0000-07*.nc")

;all_files1b:=array_append_record(all_files01,all_files02,0)
;all_files0=array_append_record(all_files1b,all_files03,0)
ff1=addfiles(all_files01,"r")
ListSetType(ff1,"cat")

land=ff1[0]->LANDFRAC(0,:,:)
lon=ff1[0]->lon
lat=ff1[0]->lat
ijk1=ind(lat.ge.-45.and.lat.le.45)
ijk2=ind(lon.ge.0.0.and.lon.le.360)

cltb=ff1[:]->OMEGA(:,18,:,:)*864.0
clt1=dim_avg_n(cltb,0)
print(avg(clt1))

;;;CAM5.3_p2K;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
all_files01:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/camrun.cam.h0.0001-*.nc")
;all_files02:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/camrun.cam.h0.0000-06*.nc")
;all_files03:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test/bld/camrun.cam.h0.0000-07*.nc")

;all_files1b:=array_append_record(all_files01,all_files02,0)
;all_files=array_append_record(all_files1b,all_files03,0)
ff=addfiles(all_files01,"r")
ListSetType(ff,"cat")

cltb:=ff[:]->OMEGA(:,18,:,:)*864.0
clt2=dim_avg_n(cltb,0)
print(avg(clt2))

;;;CAM5.3_n2K;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
all_files01:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test_n2K/bld/camrun.cam.h0.0001-*.nc")
;all_files02:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test_n2K/bld/camrun.cam.h0.0000-06*.nc")
;all_files03:=systemfunc("ls /home/hsong/zzbatmos_user/CAM5_test/cam5.3_60months_cosp_parallel_res0.9x1.25procs32_test_n2K/bld/camrun.cam.h0.0000-07*.nc")

;all_files1b:=array_append_record(all_files01,all_files02,0)
;all_files2:=array_append_record(all_files1b,all_files03,0)
ff2:=addfiles(all_files01,"r")
ListSetType(ff2,"cat")

cltb:=ff2[:]->OMEGA(:,18,:,:)*864.0
clt3=dim_avg_n(cltb,0)
print(avg(clt3))

;clt3=mask(clt3,(land.lt.1),1)
;clt2=mask(clt2,(land.lt.1),1)
;clt1=mask(clt1,(land.lt.1),1)

clt2=clt2-clt1
clt3=clt3-clt1

mclt1=avg(clt1(ijk1,ijk2))
mclt3=avg(clt3(ijk1,ijk2))
mclt2=avg(clt2(ijk1,ijk2))

;;;;;;;
alon=ff[0]->lon
alat=ff[0]->lat

clt3!0 = "lat"
clt3&lat = alat
clt3!1 = "lon"
clt3&lon = alon
clt3@units=" "

copy_VarAtts(clt3,clt1)
copy_VarAtts(clt3,clt2)

copy_VarCoords(clt3,clt1)
copy_VarCoords(clt3,clt2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;plot the figures;;;;;

wks=gsn_open_wks("png","plot_OMEGA_CAM5.3_3runs.png")

plot=new(7,graphic)

res=True
;res@gsnDraw=False
res@gsnFrame=False
res@mpLandFillColor="white"

res@tiMainOn=True
res@tiMainOffsetYF=-0.005

res@tmYLMinorOn=False
res@tmXBMinorOn=False
res@tmYROn=False
res@tmXTOn=False

res@cnLinesOn=False
res@cnLineLabelsOn=False
res@cnFillOn=True

res@lbLabelBarOn=True
res@lbOrientation="vertical"

res@cnLevelSelectionMode="ManualLevels"
res@cnMinLevelValF=-150
res@cnMaxLevelValF=150
res@cnLevelSpacingF=15
gsn_define_colormap(wks,"BlRe")
res@cnLabelBarEndStyle="ExcludeOuterBoxes"

res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=4
;res@gsnSpreadColorEnd=136
res@lbLabelStride=2

res@mpMinLatF=-45
res@mpMaxLatF=45
res@tmYLMode="Explicit"
res@tmYLValues=(/-45,-30,-15,0,15,30,45/)
res@tmYLLabels=(/"-45","-30","-15","0","15","30","45"/)

res@lbLabelBarOn=True
;;;plot2b
res@vpXF=0.325
res@vpWidthF=0.4
res@vpYF=0.95
res@vpHeightF=0.2
res@mpMinLonF=0
res@mpMaxLonF=360
res@mpCenterLonF=180
res@mpShapeMode = "FreeAspect"
res@tiMainString="CAM5.3-Control"
plot(1)=gsn_csm_contour_map(wks,clt1,res)


gsn_define_colormap(wks,"BlRe")
res@lbLabelBarOn=True
;;;;plot3b
res@vpXF=0.325
res@vpWidthF=0.4
res@vpYF=0.68
res@vpHeightF=0.2
res@lbLabelBarOn=True
res@cnLevels:=(/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)
res@tiMainString="CAM5.3-p2K minus CAM5.3-Ctrl"
plot(4)=gsn_csm_contour_map(wks,clt2,res)


res@lbLabelBarOn=True
;;;;plot3b
res@vpXF=0.325
res@vpWidthF=0.4
res@vpYF=0.41
res@vpHeightF=0.2
res@lbLabelBarOn=True
;res@cnLevels:=(/-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50/)
res@tiMainString="CAM5.3-n2K minus CAM5.3-Ctrl"
res@gsnDraw=False
plot(6)=gsn_csm_contour_map(wks,clt3,res)

draw(plot(6))


;;;;;;;;;;;;;;;;;;;;;
txres=True

txres@txFontHeightF=0.008
txres@txJust="CenterLeft"
;gsn_text_ndc(wks,"Mean = "+sprintf("%4.2f",mclt1)+"%",0.7,0.956,txres)
;gsn_text_ndc(wks,"Mean = "+sprintf("%4.2f",mclt2)+"%",0.7,0.689,txres)
;gsn_text_ndc(wks,"Mean = "+sprintf("%4.2f",mclt3)+"%",0.7,0.419,txres)

txres@txFontHeightF=0.015
txres@txJust="CenterCenter"
gsn_text_ndc(wks,"500hPa Omega [hPa/day]",0.525,0.15,txres)

txres@txFontHeightF=0.0145
;gsn_text_ndc(wks,"Annual Mean of 2009-2010",0.725,0.12,txres)

frame(wks)


end
