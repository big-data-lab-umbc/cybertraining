from __future__ import print_function,division
from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import interpolate
from matplotlib import cm
from heapq import nlargest
from heapq import nsmallest
import itertools
import os,calendar,fnmatch,datetime
import math
import sys
import glob
import pandas as pd
from scipy import interpolate

#Created for CyberTrainning by Qianqian April 23, 2018
#Read MODIS_CALIPSO collated data & MODIS data & CALIPSO data

def Read_Collocated_File(file_CALIPSO_MODIS):
    
    CAL_MOD = SD(file_CALIPSO_MODIS, SDC.READ)
    #---READ MODIS granule info---
    MODIS_Year = CAL_MOD.select('MODIS_Year')
    year = MODIS_Year.get()
    MODIS_Day = CAL_MOD.select('MODIS_DOY')
    day = MODIS_Day.get()
    MODIS_Hour = CAL_MOD.select('MODIS_Hour')
    hour = MODIS_Hour.get()
    MODIS_Minute = CAL_MOD.select('MODIS_Minute')
    minute = MODIS_Minute.get()
    Lat = CAL_MOD.select('CALIPSO_Lat1km')
    lat = Lat.get()
    Lon = CAL_MOD.select('CALIPSO_Lon1km')
    lon = Lon.get()
    MODIS_granule_ind = CAL_MOD.select('MODIS_Granule_Ind')
    modis_granule_ind = MODIS_granule_ind.get()
    
    Ind_x = CAL_MOD.select('MODIS_Cross_Ind')
    ind_x = Ind_x.get()
    Ind_y = CAL_MOD.select('MODIS_Along_Ind')
    ind_y = Ind_y.get()
    
    return year,day,hour,minute,lat,lon,modis_granule_ind,ind_x,ind_y
    
def Read_CALIPSO(file_CALIPSO):
    
    
    hdf = SD(file_CALIPSO, SDC.READ)
    lat_calip = hdf.select('Latitude').get()
    lon_calip = hdf.select('Longitude').get()
    vfm = hdf.select('Feature_Classification_Flags').get()
    nly = hdf.select('Number_Layers_Found').get()
    clm_cod_532 = hdf.select('Column_Optical_Depth_Cloud_532').get()
    feature_od_532 = hdf.select('Feature_Optical_Depth_532').get()
    return lat_calip,lon_calip,vfm,nly,clm_cod_532,feature_od_532
  
def BT(I,wvl): 
    # unit of I is W/m2/m/str; unit of wvl is mi; unit of returned Tb is K 
    h=6.62e-34
    k=1.38e-23
    c=3e8
    k1 = 2*h*c**2*wvl**(-5)
    k2 = (h*c)/(k*wvl)
    ln_term = np.log((k1/I)+1)
    Tb = k2/ln_term
    return(Tb)


def Read_MODIS(file_MODIS):
    MYD021KM = SD(file_MODIS, SDC.READ)
    lat_mod = MYD021KM.select('Latitude').get()
    lon_mod = MYD021KM.select('Longitude').get()
    solar_zenith = MYD021KM.select('SolarZenith').get()*0.01
    solar_azimuth = MYD021KM.select('SolarAzimuth').get()*0.01
    sensor_zenith = MYD021KM.select('SensorZenith').get()*0.01
    sensor_azimuth = MYD021KM.select('SensorAzimuth').get()*0.01

    
    EV_250_Aggr1km_RefSB = MYD021KM.select('EV_250_Aggr1km_RefSB')
    EV_250_Aggr1km_RefSB_data = EV_250_Aggr1km_RefSB.get()
    EV_500_Aggr1km_RefSB = MYD021KM.select('EV_500_Aggr1km_RefSB')
    EV_500_Aggr1km_RefSB_data = EV_500_Aggr1km_RefSB.get()
    EV_1km_RefSB = MYD021KM.select('EV_1KM_RefSB')
    EV_1km_RefSB_data = EV_1km_RefSB.get()
    EV_1KM_Emissive = MYD021KM.select('EV_1KM_Emissive')
    EV_1KM_Emissive_data = EV_1KM_Emissive.get()
    
    #===read two important attributes 'offset' and 'scale_factor'
    for key, value in EV_250_Aggr1km_RefSB.attributes().iteritems():
        if key == 'radiance_scales':
            EV_250_Aggr1km_RefSB_scale_factor = value
        if key == 'radiance_offsets':
            EV_250_Aggr1km_RefSB_offset = value
    for key, value in EV_500_Aggr1km_RefSB.attributes().iteritems():
        if key == 'radiance_scales':
            EV_500_Aggr1km_RefSB_scale_factor = value
        if key == 'radiance_offsets':
            EV_500_Aggr1km_RefSB_offset = value
    for key, value in EV_1KM_Emissive.attributes().iteritems():
        if key== 'radiance_scales':
            EV_1KM_Emissive_scale = value
        if key == 'radiance_offsets':
            EV_1KM_Emissive_offset = value
            
    #---transform the data and get radiance
    EV_250_Aggr1km_RefSB_radiance = np.full_like(EV_250_Aggr1km_RefSB_data,0,dtype=float)
    EV_500_Aggr1km_RefSB_radiance = np.full_like(EV_500_Aggr1km_RefSB_data,0,dtype=float)
    EV_1KM_Emissive_radiance =  np.full_like(EV_1KM_Emissive_data,0,dtype=float)
    for i in range(len(EV_250_Aggr1km_RefSB_offset)):
        EV_250_Aggr1km_RefSB_radiance[i,:,:] = (EV_250_Aggr1km_RefSB_data[i,:,:]-EV_250_Aggr1km_RefSB_offset[i])*EV_250_Aggr1km_RefSB_scale_factor[i]
    for i in range(len(EV_500_Aggr1km_RefSB_offset)):
        EV_500_Aggr1km_RefSB_radiance[i,:,:] = (EV_500_Aggr1km_RefSB_data[i,:,:]-EV_500_Aggr1km_RefSB_offset[i])*EV_500_Aggr1km_RefSB_scale_factor[i]
    for i in range(len(EV_1KM_Emissive_offset)):
        EV_1KM_Emissive_radiance[i,:,:] = (EV_1KM_Emissive_data[i,:,:]-EV_1KM_Emissive_offset[i])*EV_1KM_Emissive_scale[i]
    # ===Do the same thing for reflectance
    for key, value in EV_250_Aggr1km_RefSB.attributes().iteritems():
        if key == 'reflectance_scales':
            EV_250_Aggr1km_RefSB_scale_factor1 = value
        if key == 'reflectance_offsets':
            EV_250_Aggr1km_RefSB_offset1 = value
    for key, value in EV_500_Aggr1km_RefSB.attributes().iteritems():
        if key == 'reflectance_scales':
            EV_500_Aggr1km_RefSB_scale_factor1 = value
        if key == 'reflectance_offsets':
            EV_500_Aggr1km_RefSB_offset1 = value
    for key, value in EV_1km_RefSB.attributes().iteritems():
        if key == 'reflectance_scales':
            EV_1km_RefSB_scale_factor = value
        if key == 'reflectance_offsets':
            EV_1km_RefSB_offset = value

    #---transform the data and get reflectance
    EV_250_Aggr1km_RefSB_reflectance = np.full_like(EV_250_Aggr1km_RefSB_data,0,dtype=float)
    EV_500_Aggr1km_RefSB_reflectance = np.full_like(EV_500_Aggr1km_RefSB_data,0,dtype=float)
    EV_1km_RefSB_reflectance =  np.full_like(EV_1km_RefSB_data,0,dtype=float)
    for i in range(len(EV_250_Aggr1km_RefSB_offset)):
        EV_250_Aggr1km_RefSB_reflectance[i,:,:] = (EV_250_Aggr1km_RefSB_data[i,:,:]-EV_250_Aggr1km_RefSB_offset1[i])*EV_250_Aggr1km_RefSB_scale_factor1[i]
    for i in range(len(EV_500_Aggr1km_RefSB_offset)):
        EV_500_Aggr1km_RefSB_reflectance[i,:,:] = (EV_500_Aggr1km_RefSB_data[i,:,:]-EV_500_Aggr1km_RefSB_offset1[i])*EV_500_Aggr1km_RefSB_scale_factor1[i]
    for i in range(len(EV_1km_RefSB_offset)):
        EV_1km_RefSB_reflectance[i,:,:] = (EV_1km_RefSB_data[i,:,:]-EV_1km_RefSB_offset[i])*EV_1km_RefSB_scale_factor[i]    

    #convert EV_1KM_Emissive_radiance to brightness temperature
    EV_1KM_BT = np.full_like(EV_1KM_Emissive_radiance,0,dtype=float) 
    EV_1KM_wvl = np.array([3750e-9,3959e-9,3959e-9,4050e-9,4465.5e-9,4515.5e-9,6715e-9,7325e-9,8550e-9,9730e-9,11030e-9,12020e-9,13335e-9,13635e-9,13935e-9,14235e-9])
    for i in range(np.size(EV_1KM_Emissive_radiance,0)):
	wvl = EV_1KM_wvl[i] 
	I = EV_1KM_Emissive_radiance[i,:,:]*1e6
        EV_1KM_BT[i,:,:] = BT(I, wvl)
    return (lat_mod,lon_mod,solar_zenith,solar_azimuth,sensor_zenith,sensor_azimuth,EV_250_Aggr1km_RefSB_radiance,EV_500_Aggr1km_RefSB_radiance,EV_1KM_BT,EV_250_Aggr1km_RefSB_reflectance,EV_500_Aggr1km_RefSB_reflectance,EV_1km_RefSB_reflectance)

    
       
   
if __name__ == '__main__':
    
    #=====================
    filename = 'MYD021KM.A2007196.1525.061.2018041024335.hdf'
    lat_modis,lon_modis,SolarZenith,SolarAzimuth,SensorZenith,SensorAzimuth,rad_250,rad_500,BT_1km,ref_250, ref_500, ref_1km= Read_MODIS(filename) 
    print(ref_250.shape)
    print(ref_500.shape)
    print(ref_1km.shape)
    print(BT_1km.shape)
    print(SolarZenith.shape)

    y = np.arange(0,2026,5)
    x = np.arange(0,1351,5)
    y_new = np.arange(0,2025,1)
    x_new = np.arange(0,1350,1)

    f_SensorZenith = interpolate.interp2d(x,y,SensorZenith,kind = 'linear')
    f_SensorAzimuth = interpolate.interp2d(x,y,SensorAzimuth, kind = 'linear')
    f_SolarZenith = interpolate.interp2d(x,y,SolarZenith, kind='linear')
    f_SolarAzimuth = interpolate.interp2d(x,y,SolarAzimuth, kind='linear')
    SensorZenith_new = f_SensorZenith(x_new,y_new)
    SensorAzimuth_new = f_SensorAzimuth(x_new,y_new)
    SolarZenith_new = f_SolarZenith(x_new,y_new)
    SolarAzimuth_new = f_SolarAzimuth(x_new,y_new)
    print(ref_250[:,:2025,:1350].shape)
    print(SensorZenith_new[12,12])
    SolarZenith_new = np.reshape(SensorZenith_new,(1,2025,1350))
    SolarAzimuth_new = np.reshape(SolarAzimuth_new,(1,2025,1350))
    SensorZenith_new = np.reshape(SensorZenith_new,(1,2025,1350))
    SensorAzimuth_new = np.reshape(SensorAzimuth_new,(1,2025,1350))
    print(SensorZenith_new[0,12,12])
    print(SolarZenith_new.shape)

    is898=np.vstack((ref_250[:,:2025,:1350],ref_500[:,:2025,:1350],ref_1km[:,:2025,:1350],BT_1km[:,:2025,:1350],SolarZenith_new,SolarAzimuth_new,SensorZenith_new,SensorAzimuth_new))
    print(is898.shape)
    #convert 3-d is898[38,2030, 1354] to 2-d is698[38,2030*1354]
    is698 = is898.reshape(42,2025*1350)
    print(is698.shape)
    print(is898[0,0,0],is698[0,0])
    print(is898[0,0,1349],is698[0,1349])
    print(is898[0,1,0],is698[0,1350])
    is698=np.transpose(is698)
    modis=pd.DataFrame(is698,columns = ["Band%d" % (i + 1) for i in range(len(is698[0]))])
    modis.to_csv("/home/cd11735/cybertraining_team3/Research/Progress3_DustStorm/15July2007/ModGranule_07152007_include_geometry.csv")
    
    
