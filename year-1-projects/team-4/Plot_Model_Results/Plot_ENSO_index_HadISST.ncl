load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

f2=addfile("/umbc/xfs1/zzbatmos/users/hsong/ncl/data/HadISST_sst.nc","r")

lat=f2->latitude
lon=f2->longitude
time=f2->time

latS  = -5.0
latN  =  5.0
lonL  = -170.0
lonR  = -120.0

nrun  = 3             ; length of running average

yrStrt= 1950  
yrLast= 2017      

clStrt= 1950          ; climatology start    
clLast= 2000          ;             last

sst=f2->sst(:,{latS:latN},{lonL:lonR})

ymStrt = yrStrt*100 +  1            
ymLast = yrLast*100 + 12

clStrt = clStrt*100 +  1              ; redefine
clLast = clLast*100 + 12

YYYYMM = cd_calendar(sst&time, -1)      ; ALL dates assciated with X
print(YYYYMM(0))
tStrt  = ind(YYYYMM.eq.ymStrt)        ; indices of selected times
tLast  = ind(YYYYMM.eq.ymLast)
delete(YYYYMM)

sstn=sst(tStrt:tLast,:,:)
yyyymm = cd_calendar(sstn&time, -1) 
dimx   = dimsizes(sstn)
ntim   = dimx(0)
print(ntim)

print(yyyymm(0:5))


iClmStrt = ind(yyyymm.eq.clStrt)     
iClmLast = ind(yyyymm.eq.clLast) 

xClm     = clmMonTLL(sstn(iClmStrt:iClmLast,:,:))
xAnom    = calcMonAnomTLL (sstn,  xClm )

xAnom_avg = wgt_areaave_Wrap(xAnom, 1.0, 1.0, 1)
xAnom_avg@long_name = "areal avg anomalies"
printVarSummary(xAnom_avg)


asciiwrite("ENSO_index_1950_2017.txt",sprintf("%16.8f",xAnom_avg))

;*********************************
; plot
;*********************************

  yrfrac = yyyymm_to_yyyyfrac(yyyymm, 0.0)
  ;print(yrfrac)
  wks = gsn_open_wks("png","ENSO_index_HadISST")     
  plot = new(1, graphic )

  res                       = True       
 ; res@gsnDraw               = False       
  res@gsnFrame              = False     
  res@vpHeightF             = 0.3            ; change aspect ratio of plot
  res@vpWidthF              = 0.7            
  res@vpXF                  = 0.2            ; start plot at x ndc coord 
  res@gsnYRefLine           = 0.0            ; create a reference line 
  res@gsnAboveYRefLineColor = "red"              ; above ref line fill red
  res@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue
  
  res@trYMinF               = -3.0           ; min value on y-axis
  res@trYMaxF               =  3.0           ; max value on y-axis

  
  res@trXMinF               = 1950           ; min value on x-axis
  res@trXMaxF               =  2020           ; max value on x-axis 
  res@tiMainOn=True
  res@tiMainString="ENSO Index (1950-2017)"
  res@tiMainOffsetYF=-0.005
  res@tiYAxisString    = "Anomalies (C)"    ; y-axis label      
  plot(0) = gsn_csm_xy (wks,yrfrac,xAnom_avg,res)

  res@tiYAxisString    = "Standardized Anomalies"    ; y-axis label      
  ;plot(1) = gsn_csm_xy (wks,yrfrac,xAnom_std,res)

  frame(wks)

end
