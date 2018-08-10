load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

ff1=addfile("land.nc","r")
lon=ff1->lon
lat=ff1->lat
land=ff1->land(0,:,:)
mn=dimsizes(land)

max_cor=asciiread("max_cor_w500_lagged_by_ENSO.txt",-1,"float")
max_cor2=reshape(max_cor,(/mn(0),mn(1)/))
;max_cor2=mask(max_cor2,(land.lt.1),1)

max_lag=asciiread("max_lag_w500_lagged_by_ENSO.txt",-1,"float")
max_lag2=reshape(max_lag,(/mn(0),mn(1)/))
;max_lag2=mask(max_lag2,(land.lt.1),1)

max_ind=asciiread("max_ind_w500_lagged_by_ENSO.txt",-1,"float")
max_ind2=floattointeger(reshape(max_ind,(/mn(0),mn(1)/)))

mon=(/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
printVarSummary(mon)

max_cor2!0 = "lat"
max_cor2&lat = lat
max_cor2!1 = "lon"
max_cor2&lon = lon
;max_cor2@_FillValue=-9999.99

max_lag2!0 = "lat"
max_lag2&lat = lat
max_lag2!1 = "lon"
max_lag2&lon = lon
;max_lag2@_FillValue=-9999.99

max_ind2!0 = "lat"
max_ind2&lat = lat
max_ind2!1 = "lon"
max_ind2&lon = lon

;;;;;;;Plot the results on the map
;;;;;;;Open png for graphics.
  
wks  = gsn_open_wks("pdf","ENSO_w500_Corr_Lag_Ind_global")

res                      = True
res@gsnMaximize          = True          ; Maximize plot in frame.
res@gsnFrame             = False         ; Turn off so we can add markers and text
res@cnFillOn             = True          ; Turn on contour fill
res@cnFillMode           = "RasterFill"  ; Turn on raster fill
res@cnLinesOn            = False         ; Turn off contour lines
;res@cnFillPalette        = "WhViBlGrYeOrRe"
res@cnFillPalette        = "wh-bl-gr-ye-re"
res@cnFillPalette = "BlWhRe"
  
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF       =  -1.0
res@cnMaxLevelValF       =  1.0
res@cnLevelSpacingF      =  0.0125
res@lbOrientation        = "vertical"
res@cnLabelBarEndStyle="ExcludeOuterBoxes"
res@pmLabelBarWidthF      = 0.045               ; default is shorter
;res@pmLabelBarHeightF     = 0.025               ; default is taller
res@lbBoxLinesOn=False
;res@gsnSpreadColorStart=0
;res@gsnSpreadColorEnd=198
res@lbLabelStride=16

res@mpLandFillColor="white"  
res@mpMinLatF=-90
res@mpMaxLatF=90
res@mpMinLonF=0
res@mpMaxLonF=360
res@mpCenterLonF=180
;res@tmYLMode="Explicit"
;res@tmYLValues=(/-60,-45,-30,-15,0,15/)
;res@tmYLLabels=(/"60S","45S","30S","15S","0","15N"/)
;res@tmXBMode="Explicit"
;res@tmXBValues=(/270,280,290,300,310,320,330/)
;res@tmXBLabels=(/"90W","","70W","","50W","","30W"/)
res@mpGridAndLimbOn        = True
res@mpGridAndLimbDrawOrder = "PreDraw"
res@mpGridLineDashPattern  = 2                  ; lat/lon lines dashed
res@mpGridLatSpacingF      = 15
res@mpGridLonSpacingF      = 15  

res@tmXBMajorLengthF     = 0.0
res@tmYLMajorLengthF     = 0.0
res@tmXBMinorOn          = False
res@tmYLMinorOn          = False

res@tiMainString         = "Lag-Correlation between ENSO-Index and 500hPa Vertical Velocity"
res@tiMainOffsetYF=-0.005
res@tiMainFontHeightF   = 0.0175 

plot = gsn_csm_contour_map(wks,max_cor2,res)

;---Draw text and markers at data locations
txres               = True
mkres               = True
txres@txFontHeightF = 0.0035
txres@txJust        = "TopCenter"  ; text will be drawn under the marker

do j=0,mn(0)-1
do i=0,mn(1)-1
  ;if(land(j,i).lt.1.and.max_ind2(j,i).gt.0)then
    ;gsn_text(wks,plot,mon(max_ind2(j,i)-1),lon(i),lat(j)+0.65,txres)
    ;gsn_text(wks,plot," ("+max_lag2(j,i)+")",lon(i),lat(j)-0.45,txres)
    gsn_text(wks,plot," "+max_lag2(j,i)+" ",lon(i),lat(j)-0.45,txres)
  ;end if
end do
end do

frame(wks)

end
