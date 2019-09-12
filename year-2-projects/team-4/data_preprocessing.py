# -*- coding: utf-8 -*-
import numpy as np
import h5py   
import os
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


train_data_predictor=[]
train_data_target=[]
train_dust_idx=[]
test_data_predictor=[]
test_data_target=[]

predictor_keys=['VIIRS_M01','VIIRS_M02','VIIRS_M03','VIIRS_M04','VIIRS_M05','VIIRS_M06',
 'VIIRS_M07','VIIRS_M08','VIIRS_M09','VIIRS_M10','VIIRS_M11','VIIRS_M12','VIIRS_M13',
 'VIIRS_M14','VIIRS_M15','VIIRS_M16','VIIRS_SAA','VIIRS_SZA','VIIRS_VAA','VIIRS_VZA']
target_keys=['Pixel_Label']

atype_keys = ['CALIOP_Alay_Optical_Depth_532','CALIOP_Clay_Optical_Depth_532']

unmatched_length=[]
unmatched_location=[]

hdfpath = '/home/pengwu/cybertraining2019_team4/research/2019Team4/calipso-virrs-collocated/data_for_test/2014/'

listed_days=os.listdir(hdfpath)
for dd in range(len(listed_days)):
    d=listed_days[dd]
    if 59<int(d)<151:
        listed_time=os.listdir(hdfpath+d)
        for t in range(len(listed_time)):
            predictor_fromfile=np.zeros([len(predictor_keys),3248])
            target_fromfile=np.zeros(3248)
            aod_fromfile=np.zeros(3248)
            cod_fromfile=np.zeros(3248)
            dust_idx=np.zeros(3248)
            predictor_fromfile[:,:]=np.nan
            target_fromfile[:]=np.nan
            aod_fromfile[:]=np.nan
            cod_fromfile[:]=np.nan
#            dust_idx[:]=np.nan

            f=h5py.File(hdfpath+d+'/'+listed_time[t])
            
            data_length=len(f[predictor_keys[0]][:])
            if data_length<2600:
                pass
            else:
                for v in range(len(predictor_keys)):
                    predictor_fromfile[v,:data_length]=f[predictor_keys[v][:]]
                target_fromfile[:data_length]=f[target_keys[0][:]]
                aod_prof = f['CALIOP_Alay_Optical_Depth_532'[:]]
                cod_prof = f['CALIOP_Clay_Optical_Depth_532'[:]]
                aod_fromfile[:data_length] = np.nansum(aod_prof,axis=0)
                cod_fromfile[:data_length] = np.nansum(cod_prof,axis=0)

                target_fromfile[(target_fromfile<5)&(target_fromfile>=7)]=0
                target_fromfile[target_fromfile>0]=1

                dust_idx[(target_fromfile>0)&(cod_fromfile<0.05)&( aod_fromfile>0.05) ] = 1
                dust_idx[(target_fromfile>0)&(cod_fromfile>0.05)&( aod_fromfile<0.05) ] = 2
                dust_idx[(target_fromfile>0)&(cod_fromfile>0.05)&( aod_fromfile>0.05) ] = 3

            # Split data into train/test dataset
            if data_length<2600:
                pass
            else:
                if 59<int(d)<90: # replace with training range 59,121
                    train_data_predictor.append(predictor_fromfile)
                    train_data_target.append(target_fromfile)
                    train_dust_idx.append(dust_idx)
                else:
                    test_data_predictor.append(predictor_fromfile)
                    test_data_target.append(target_fromfile)
    else:
        pass
    print("Day of Year  "+str(d)+" was read")


train_data_predictor=np.array(train_data_predictor)
train_data_target=np.array(train_data_target)
train_dust_idx=np.array(train_dust_idx)
test_data_predictor=np.array(test_data_predictor)
test_data_target=np.array(test_data_target)
print("All days were read")

for v in range(11):
    train_data_predictor[:,v,:][train_data_predictor[:,v,:]>1000]=np.nan
    train_data_predictor[:,v,:][train_data_predictor[:,v,:]<-1000]=np.nan
    test_data_predictor[:,v,:][test_data_predictor[:,v,:]>1000]=np.nan
    test_data_predictor[:,v,:][test_data_predictor[:,v,:]<-1000]=np.nan


# Normalize predictor data
for v in range(20):
    train_data_predictor[:,v,:]=\
    (train_data_predictor[:,v,:]-np.nanmean(train_data_predictor[:,v,:]))/np.nanstd(train_data_predictor[:,v,:])
    test_data_predictor[:,v,:]=\
    (test_data_predictor[:,v,:]-np.nanmean(test_data_predictor[:,v,:]))/np.nanstd(test_data_predictor[:,v,:])

train_data_predictor_save=train_data_predictor.swapaxes(1,2)
test_data_predictor_save=test_data_predictor.swapaxes(1,2)
#-----------------------------------------------------------------------------------
# SAVE DATA
#-----------------------------------------------------------------------------------

train_prdt = np.zeros([train_dust_idx.ravel().shape[0],len(predictor_keys)])
#train_prdt01 = np.zeros([train_dust_idx.ravel().shape[0],len(predictor_keys)])
for i in range(len(predictor_keys)):
    train_prdt[:,i]=train_data_predictor_save[:,:,i].ravel()

param = np.hstack((train_prdt, train_dust_idx.ravel().reshape([len(train_dust_idx.ravel()),1])))
param01 = np.hstack((train_prdt, train_data_target.ravel().reshape([len(train_data_target.ravel()),1])))

yy=param
ff=param01
for i in range(param.shape[1]):
    zz=np.delete(yy, np.where(np.isnan(yy[:,i])==1), axis=0)
    gg=np.delete(ff, np.where(np.isnan(ff[:,i])==1), axis=0)
    yy=zz
    ff=gg

sample_num = 2500

loc0 = np.where(zz[:,-1] == 0)[0]
loc1 = np.where(zz[:,-1] == 1)[0]
loc2 = np.where(zz[:,-1] == 2)[0]
loc3 = np.where(zz[:,-1] == 3)[0]

zz[:,-1][loc0[sample_num:]] = 'NaN'
zz[:,-1][loc1[sample_num:]] = 'NaN'
zz[:,-1][loc2[sample_num:]] = 'NaN'
zz[:,-1][loc3[sample_num:]] = 'NaN'
zz_save = np.delete(yy, np.where(np.isnan(zz[:,-1])==1), axis=0)

column_name = predictor_keys + ['dust']

import pandas as pd
data_to_save=pd.DataFrame(zz,columns = [i for i in column_name])
data_to_save.to_csv('/home/pengwu/cybertraining2019_team4/research/2019Team4/march_april/march2014_new.csv')

np.save('/home/pengwu/cybertraining2019_team4/research/2019Team4/march_april/march2014_new.npy', zz_save)
np.save('/home/pengwu/cybertraining2019_team4/research/2019Team4/march_april/march2014_01_new.npy', gg)

print('Data saved in cybertraining2019_team4/research/2019Team4/march_april')
