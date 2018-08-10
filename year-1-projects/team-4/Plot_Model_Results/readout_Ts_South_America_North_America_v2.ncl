
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

;;;;;;;;;;;;;;;;;;;


ff1=addfile("/home/hsong/zzbatmos_user/ncl/data/air.mon.mean.nc","r")
T=ff1->air(:,:,:)
lon=ff1->lon
lat=ff1->lat

ii1=ind(lat.gt.56.and.lat.lt.60)
jj1=ind(lon.gt.288.and.lon.lt.292)

ii2=ind(lat.gt.-30.and.lat.lt.-26)
jj2=ind(lon.gt.298.and.lon.lt.302)
print(ii1)
print(jj1)
print(ii2)
print(jj2)

Ts1=T(24:839,ii1,jj1)
Ts2=T(24:839,ii2,jj2)

asciiwrite("Ts_North_America_70W_58N_195001_201712.txt",sprintf("%16.8f",Ts1))
asciiwrite("Ts_South_America_60W_28S_195001_201712.txt",sprintf("%16.8f",Ts2))

end
