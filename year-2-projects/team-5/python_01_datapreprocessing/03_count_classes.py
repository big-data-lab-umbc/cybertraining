# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 20:44:53 2019

@author: Jangho
"""


import numpy as np
import pickle
import h5py   
import os

#--------------------------------------------------------------------------------
# This program reads collocated VIIRS and surrounding data for each season
# After reading the data, data is then reshaped bu 5x5 moving window for 
# CNN-applicable format
#--------------------------------------------------------------------------------
#
# First we define each Season, as well as testing and training period
#
# Spring : 60 to 151
# Summer : 152 to 243
# Fall : 244 to 334
# Winter : 335 to 365 and 1 to 59

# Validation set will also be described below
#--------------------------------------------------------------------------------

# DEFINE TRAINING AND TESTING PERIOD FOR EACH SEASON
spring_train_fname=[]
summer_train_fname=[]
fall_train_fname=[]
winter_train_fname=[]

spring_test_fname = ['060', '063', '066', '069', '072', '075', '078', '081', '084',\
             '087', '090', '093', '096', '099', '102', '105', '108',\
             '111', '114', '117', '120', '123', '126', '129', '132',\
             '135', '138', '141', '144', '147', '150' ]
summer_test_fname = [ '152', '155', '158', '161', '164', '167', '170',\
              '173', '176', '179', '182', '185', '188', '191',\
              '194', '197', '200', '203', '206', '209', '212',\
              '215', '218', '221', '224', '227', '230', '233',\
              '236', '239', '242']
fall_test_fname = [ '244', '247', '250', '253', '256', '259', '262',\
               '265', '268', '271', '274', '277', '280', '283',\
               '286', '289', '292', '295', '298', '301', '304',\
               '307', '310', '313', '316', '319', '322', '325',\
               '328', '331']
winter_test_fname = [ '002', '004', '007', '010', '013', '016', '019', '022', '025',\
              '028', '031', '034', '037', '040', '043', '046', '049', '052', '055',\
              '058', '335', '338', '341', '344', '347', '350', '353',\
              '356', '359', '362' ]
for d in range(365):
    d=d+1
    fname=str(d)
    if len(fname)==1:
        fname='00'+fname
    elif len(fname)==2:
        fname='0'+fname
    if 60<=d<152:
        if fname in spring_test_fname:
            pass
        else:
            spring_train_fname.append(fname)
    elif 152<=d<244:
        if fname in summer_test_fname:
            pass
        else:
            summer_train_fname.append(fname)
    elif 244<=d<334:
        if fname in fall_test_fname:
            pass
        else:
            fall_train_fname.append(fname)
    else:
        if fname in winter_test_fname:
            pass
        else:
            winter_train_fname.append(fname)

# These are the Keys for the 1d data
predictor_keys=['VIIRS_M01','VIIRS_M02','VIIRS_M03','VIIRS_M04','VIIRS_M05','VIIRS_M06',
 'VIIRS_M07','VIIRS_M08','VIIRS_M09','VIIRS_M10','VIIRS_M11','VIIRS_M12','VIIRS_M13',
 'VIIRS_M14','VIIRS_M15','VIIRS_M16','VIIRS_SAA','VIIRS_SZA','VIIRS_VAA','VIIRS_VZA']
target_keys=['Pixel_Label', 'Latitude', 'Longitude']

print("LOADED ALL PRE-DESCRIPTION", flush=True)

#--------------------------------------------------------------------------------
# READ DATA and SPLIT
#--------------------------------------------------------------------------------
# There are currently 167 days total - wich will be added later
filepath1='../data_org/VIIRS/'
filepath1x='../data_org/Collocated/'
files_total=os.listdir(filepath1)

spring_predictor_train=[]
spring_predictor_test=[]
spring_target_train=[]
spring_target_test=[]

summer_predictor_train=[]
summer_predictor_test=[]
summer_target_train=[]
summer_target_test=[]

fall_predictor_train=[]
fall_predictor_test=[]
fall_target_train=[]
fall_target_test=[]

winter_predictor_train=[]
winter_predictor_test=[]
winter_target_train=[]
winter_target_test=[]



target_classes=[]
count1=1
for f1 in range(len(files_total)): #len(files_total)
    filepath2=filepath1+files_total[f1]+'/'
    filepath2x=filepath1x+files_total[f1]+'/'
    files_in_a_day=os.listdir(filepath2)
    files_in_a_dayx=os.listdir(filepath2x)
    
    # Initialize predictor and target array (per day)
    predictor_all=[]
    target_all=[]
    
    count2=1
    for f2 in range(len(files_in_a_day)):
        filepath3=filepath2+files_in_a_day[f2]
        
        time_key1=files_in_a_day[f2][-14:-10]
        time_key2=files_in_a_day[f2][-8:-4]
        # Read surroundings data
        file_2d=open(filepath3, 'rb')
        data_2d=pickle.load(file_2d, encoding='latin1')
        file_2d.close()
        
        # Read matching 1d Data
        
        for t1 in range(len(files_in_a_dayx)):
            if ('t'+time_key1+'.cal' in files_in_a_dayx[t1]) and ('t'+time_key2+'.h5' in files_in_a_dayx[t1]):
                filepath3x=filepath2x+files_in_a_dayx[t1]
            else:
                pass
        if filepath3x==filepath2x:
            print("NO MATCHING DATA FOUND FOR : "+files_in_a_day[f2], flush=True)
        
        else:
            data_1d=h5py.File(filepath3x)
            data_length_2d=np.shape(data_2d)[2]
            data_length_1d=len(data_1d[predictor_keys[0]][:])
            
            if data_length_1d==data_length_2d:
                
                data_predictor=np.zeros([5, 20+3+9, data_length_2d])
                data_target=np.zeros([data_length_2d, 4])
                
                data_predictor[:,:16,:]=data_2d
                data_predictor[:,16,:]=data_1d[predictor_keys[-4]]
                data_predictor[:,17,:]=data_1d[predictor_keys[-3]]
                data_predictor[:,18,:]=data_1d[predictor_keys[-2]]
                data_predictor[:,19,:]=data_1d[predictor_keys[-1]]
                
                data_predictor[:,20,:]=np.array(list(data_1d[target_keys[1]]))[:,0]
                data_predictor[:,21,:]=np.array(list(data_1d[target_keys[2]]))[:,0]
                data_predictor[:,22,:]=int(files_total[f1])
                
                data_target[:,0]=data_1d[target_keys[0]]
                data_target[:,1]=np.array(list(data_1d[target_keys[1]]))[:,0]
                data_target[:,2]=np.array(list(data_1d[target_keys[2]]))[:,0]
                data_target[:,3]=int(files_total[f1])
                
                # Replace aerosol data with 0 and 1
                #data_target[:,0][data_target[:,0]==54]=0
                #data_target[:,0][data_target[:,0]==64]=0
                #data_target[:,0][data_target[:,0]>=70]=0
                #data_target[:,0][data_target[:,0]<50]=0
                #data_target[:,0][data_target[:,0]>=60]=2
                #data_target[:,0][data_target[:,0]>=50]=1
                
                # Remove clouds
                if np.sum(data_target)==np.nansum(data_target):
                    pass
                else:
                    data_predictor[:,0,:]=9999
                
                # Replace overflow values with nan
                data_predictor[data_predictor>1000]=np.nan
                data_predictor[data_predictor<-1000]=np.nan
                
                # Extract data with moving window of 5X5, and interval of 5
                for cc in range((data_length_1d-4)//15):
                    c=cc*15+2 # This is the center location
                    predictor=data_predictor[:,:,c-2:c+3].swapaxes(1,2) 
                    # surrounding, window, channel
                    target=data_target[c,:]
                    targets1=data_target[c-2,:]
                    targets2=data_target[c-1,:]
                    targets3=data_target[c+1,:]
                    targets4=data_target[c+2,:]
                    
                    if (np.sum(np.isnan(predictor))==0) and (np.sum(np.isnan(target))==0) and (target[0] != 2)\
                    and (targets1[0]!=2) and (targets2[0]!=2) and (targets3[0]!=2) and (targets4[0]!=2):
                        predictor_all.append(predictor)
                        target_all.append(target)
                        target_classes.append(target[0])
    
            else:
                pass
                #print("DATA LENGTH UNMATCHED FOR : "+files_in_a_day[f2])
        
        print("DAY "+str(count1)+" of "+str(len(files_total)), "   TIME "+\
              str(count2)+" of "+str(len(files_in_a_day)), flush=True)
        count2=count2+1

np.save('../data/target_classes.npy', np.array(target_classes))