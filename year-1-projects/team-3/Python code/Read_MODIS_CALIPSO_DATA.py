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
    #---Collocated Data---
    #here we only use one track of collocated data
    data_dir = '/home/cd11735/cybertraining_team3/Research/Progress3_DustStorm/15July2007/'

    years,Jdays,hs,ms,lat,lon,modis_granule_ind,x,y=Read_Collocated_File(data_dir+'CAL_LID_L2_01kmCLay-ValStage1-V3-01.2007-07-15T15-04-26ZD_IND.hdf')
    #years,Jdays,hs,ms,lat,lon,modis_granule_ind,x,y=Read_Collocated_File(data_dir+'CAL_LID_L2_01kmCLay-ValStage1-V3-01.2009-07-18T15-20-56ZD_IND.hdf') 
    # focus on the interested region by specifying the lat 
    # lat_interest gives us marker for those collocated pixels of our interest
    lat_interest = (lat >= 0)&(lat <= 30)
    ind_collocate = np.where((lat >= 0)&(lat <= 30))
    print('The index of the interested pixels in the collocated data file:',ind_collocate[0][0],ind_collocate[0][-1] )
    # find out the index of the granule in the interested region
    granule_ind = modis_granule_ind[lat_interest] 
    granule_ind_interest = np.unique(granule_ind)
    #print(granule_ind_interest)
    #find out the datetime of interest, thereby we could find MODIS granule file of interest
    years_interest = years[granule_ind_interest]
    Jdays_interest = Jdays[granule_ind_interest]
    hs_interest = hs[granule_ind_interest]
    ms_interest = ms[granule_ind_interest]
    #print(years_interest,Jdays_interest,hs_interest,ms_interest)
    #find out the indices of pixels we need in each granule 
    
    #===========
    #---MODIS---
    ref1_250 = np.zeros([2,len(x)])
    ref1_500 = np.zeros([5,len(y)])
    ref1_1km = np.zeros([15,len(x)])
    BT1_1km = np.zeros([16,len(x)])
    solar_zenith1 = np.zeros(len(x))
    solar_azimuth1 = np.zeros(len(x))
    sensor_zenith1 = np.zeros(len(x))
    sensor_azimuth1 = np.zeros(len(x))
 
    lat_modis_allgranules = np.zeros(len(x))
    lon_modis_allgranules = np.zeros(len(x))
    
    i_pix = 0 # indicate the index of each collocated pixel
    num_granule_interest = len(years_interest)
    for i_granule in range(num_granule_interest):
        #-----Find out the x_index, y_index of the interested pixel in each MODIS granule
        x_ind = x[(lat_interest)][(granule_ind == granule_ind_interest[i_granule])]
        y_ind = y[(lat_interest)][(granule_ind == granule_ind_interest[i_granule])]
  
    
        #-----Find out the date and granule time of each MODIS granule--------------------
        yc = years_interest[i_granule]
        Jday = Jdays_interest[i_granule]
        t = datetime.time(hs_interest[i_granule],ms_interest[i_granule])
        modis_file = data_dir+'MYD021KM.A{:04d}{:03d}.{:2d}{:02d}.???.?????????????.hdf'.format(yc,Jday,t.hour,t.minute)
        
        #---find all pathnames matching a specified path pattern
        for filename in glob.glob(modis_file):
            print(filename)
            lat_modis,lon_modis,solar_zenith,solar_azimuth,sensor_zenith,sensor_azimuth,rad_250,rad_500,BT_1km,ref_250, ref_500, ref_1km= Read_MODIS(filename)
	    
            for i_modis in range(len(x_ind)):
		solar_zenith1[ind_collocate[0][i_pix]] = solar_zenith[math.ceil(y_ind[i_modis]/5.0)-1,math.ceil(x_ind[i_modis]/5.0)-1]
		solar_azimuth1[ind_collocate[0][i_pix]] = solar_azimuth[math.ceil(y_ind[i_modis]/5.0)-1,math.ceil(x_ind[i_modis]/5.0)-1]
		sensor_zenith1[ind_collocate[0][i_pix]] = sensor_zenith[math.ceil(y_ind[i_modis]/5.0)-1,math.ceil(x_ind[i_modis]/5.0)-1]
		sensor_azimuth1[ind_collocate[0][i_pix]] = sensor_azimuth[math.ceil(y_ind[i_modis]/5.0)-1,math.ceil(x_ind[i_modis]/5.0)-1]
		for m in range(0,2):
                   ref1_250[m][ind_collocate[0][i_pix]] = ref_250[m,y_ind[i_modis],x_ind[i_modis]]
                for m in range(0,5):
                   ref1_500[m][ind_collocate[0][i_pix]] = ref_500[m,y_ind[i_modis],x_ind[i_modis]]
                for m in range(0,15):
                   ref1_1km[m][ind_collocate[0][i_pix]] = ref_1km[m,y_ind[i_modis],x_ind[i_modis]]
                for m in range(0,16):
                   BT1_1km[m][ind_collocate[0][i_pix]] = BT_1km[m,y_ind[i_modis],x_ind[i_modis]]
                
                #only for sanity check
                lat_modis_allgranules[ind_collocate[0][i_pix]] = lat_modis[int(y_ind[i_modis]/5),int(x_ind[i_modis]/5)]
                lon_modis_allgranules[ind_collocate[0][i_pix]] = lon_modis[int(y_ind[i_modis]/5),int(x_ind[i_modis]/5)]
                #print(lat_modis_allgranules[ind_collocate[0][i_pix]])
                i_pix+=1

    #print('i_pix:',i_pix)  # This variable should be equeal to print(ind_collocate[0].shape)
    #--------
    #Now we have MODIS data (ref_460 and ref_860) corresponding to each pixel in the collocated data (resolution is 1km)
    #Note only the MODIS pixels in the interested region have value, other pixels have value 0.
    #Because CALIPSO aerosol layer product is 5km, we need to groud each 5 pixels of MODIS data and then compare to CALIPSO data
    #Here, we plan to do MODIS dust detection for every 5 pixels (5km) in order to compare with CALIPSO data
    #-------
    # The index of pixels in the interested region start at:ind_collocate[0][0] ends at ind_collocate[0][-1]
    # The index of pixels for dust detection in MODIS will start at 'ind1_for_mod' ends at 'ind2_for_mod'
    # The corresponding pixels in CALIPSO(5KM) file starts at 'ind1_calipso' ends at 'ind2_calipso'
    
    ind1_calipso = int(ind_collocate[0][0]/5)+1
    ind2_calipso = int(ind_collocate[0][-1]/5)
    ind1_for_mod = ind1_calipso *5
    ind2_for_mod = ind2_calipso *5
    #Therefore, MODIS data we are interested is:
    ref_250_interest = ref1_250[:,ind1_for_mod:ind2_for_mod]
    ref_500_interest = ref1_500[:,ind1_for_mod:ind2_for_mod]
    ref_1km_interest=ref1_1km[:,ind1_for_mod:ind2_for_mod]
    BT_1km_interest=BT1_1km[:,ind1_for_mod:ind2_for_mod]
    solar_zenith_interest = solar_zenith1[ind1_for_mod:ind2_for_mod]
    solar_azimuth_interest = solar_azimuth1[ind1_for_mod:ind2_for_mod]
    sensor_zenith_interest = sensor_zenith1[ind1_for_mod:ind2_for_mod]
    sensor_azimuth_interest = sensor_azimuth1[ind1_for_mod:ind2_for_mod]
    
    print(ref_250_interest.shape,solar_zenith_interest.shape)
    is898=np.vstack((ref_250_interest,ref_500_interest,ref_1km_interest,BT_1km_interest,solar_zenith_interest,\
    solar_azimuth_interest,sensor_zenith_interest,sensor_azimuth_interest))

    lat_modis_interest = lat_modis_allgranules[ind1_for_mod:ind2_for_mod]
    lon_modis_interest = lon_modis_allgranules[ind1_for_mod:ind2_for_mod]
    nn=int(len(ref_250_interest[0,:])/5) 
    '''
    #Average reflectance for every 5 pixels
    ref_250_ave = np.zeros((2,nn))
    ref_500_ave = np.zeros((5,nn))
    ref_1km_ave = np.zeros((15,nn))
    BT_1km_ave = np.zeros((16,nn))

    for m in range(0,2):
        ref_250_ave[m,:] = np.mean(ref_250_interest[m,:].reshape(nn,5),axis=1)
    for m in range(0,5):
        ref_500_ave[m,:] = np.mean(ref_500_interest[m,:].reshape(nn,5),axis=1)
    for m in range(0,15):
        ref_1km_ave[m,:] = np.mean(ref_1km_interest[m,:].reshape(nn,5),axis=1)
    for m in range(0,16):
        BT_1km_ave[m,:] = np.mean(BT_1km_interest[m,:].reshape(nn,5),axis=1)
    '''
    lat_mod = np.mean(lat_modis_interest.reshape(nn,5),axis=1)
    lon_mod = np.mean(lon_modis_interest.reshape(nn,5),axis=1)
    #=============
    #---CALIPSO---
    from CALIPSO import Extract_Feature_Info
    lat_calip,lon_calip,vfm,nly,cod,feature_od_532 = Read_CALIPSO(data_dir+'CAL_LID_L2_05kmALay-Prov-V3-01.2007-07-15T15-04-26ZD.hdf') 
    #lat_calip,lon_calip,vfm,nly,cod,feature_od_532 = Read_CALIPSO(data_dir+'CAL_LID_L2_05kmALay-Prov-V3-01.2009-07-18T15-20-56ZD.hdf')
    feature_type, feature_type_qa,\
    ice_water_phase, ice_water_phase_qa,\
    feature_subtype, cloud_aerosol_psc_type_qa,\
    horizontal_averaging=Extract_Feature_Info(vfm,nly)
    
    nprofile = len(nly)
    dust_indicator = np.zeros(nprofile)
    
    for ip in range(nprofile):
        
            dust_aod = 0.0
            for ily in range(nly[ip]):
                if feature_type[ip,ily]==3 and feature_subtype[ip,ily]==2\
                                    and feature_type_qa[ip,ily] > 1:
                    dust_aod += feature_od_532[ip,ily]
                    

            #------------------------------
            #Identify dust w/o clouds cases
            if cod[ip] == 0.0 and dust_aod > 0.0:
                dust_indicator[ip] =1
            #Identify cloud w/o dust cases
            if dust_aod == 0.0 and cod[ip] > 0.0:
                dust_indicator[ip] =2
            #Identify both cloud and dust cases
            if dust_aod > 0.0 and cod[ip] > 0.0:
                dust_indicator[ip] =3
                
            
    #Calipso dust detection for the interested region (5km resolution)
    calipso_dust_detection = dust_indicator[ind1_calipso:ind2_calipso]
    ###################################################################
    is698=np.transpose(is898)
    cap=np.arange(len(calipso_dust_detection)).reshape((-1,1))
    caplipso=pd.DataFrame(np.column_stack((cap,calipso_dust_detection)),columns=["id","dust"])
    mod=np.arange(len(is698))
    mod=np.divide(mod,5)
    mod_id=pd.DataFrame(mod.reshape((-1,1)))
    modis=pd.DataFrame(is698,columns = ["Band%d" % (i + 1) for i in range(len(is698[0]))])
    print(len(is698))
    modis["id"]=mod_id
    modis.to_csv("/home/cd11735/cybertraining_team3/Research/Progress3_DustStorm/15July2007/modis.csv")
    caplipso.to_csv("/home/cd11735/cybertraining_team3/Research/Progress3_DustStorm/15July2007/caplipso.csv")
    print(modis.head(n=3))
    is698_final=modis.merge(caplipso,on='id',how='left')

    print(is698_final['dust'].value_counts())
    is698_final.to_csv("/home/cd11735/cybertraining_team3/Research/Progress3_DustStorm/15July2007/is698.csv") 
    

    
