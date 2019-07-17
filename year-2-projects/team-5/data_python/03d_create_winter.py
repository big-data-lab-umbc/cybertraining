#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:51:52 2019

@author: jangho.lee.92
"""

import numpy as np


# READ DATAFRAME and INITIALIZE ARRAY
fn_train=np.load('../dataframe/train_fn.npy')
dust_train=np.load('../dataframe/train_dust.npy')
fn_test=np.load('../dataframe/test_fn.npy')
dust_test=np.load('../dataframe/test_dust.npy')


season_train=[]
season_test=[]
season_train_dust=[]
season_test_dust=[]

c1=60
c2=334
root1='/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/CNN/data/train_predictor_norm/'
root2='/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/CNN/data/test_predictor_norm/'

cou1=0
cou2=0
for i in range(len(fn_train)):
    if int(fn_train[i][:3])<=c1 or c2<int(fn_train[i][:3]): 
        if cou1/5==cou1//5:
            season_train.append(np.load(root1+fn_train[i]))
            season_train_dust.append(dust_train[i])
        else:
            pass
    if i//1000==i/1000:
        print('TRAIN: '+str(i*100/len(fn_train))[:5]+' %', i, flush=True)
    cou1=cou1+1
        
        
for i in range(len(fn_test)):
    if int(fn_train[i][:3])<=c1 or c2<int(fn_train[i][:3]): 
        if cou2/5==cou2//5:
            season_test.append(np.load(root2+fn_test[i]))
            season_test_dust.append(fn_test[i])
        else:
            pass
    if i//1000==i/1000:
        print('TRAIN: '+str(i*100/len(fn_train))[:5]+' %', i, flush=True)
    cou2=cou2+1

np.save('../data/seasonal_data/winter_train_predictor.npy', np.array(season_train))
np.save('../data/seasonal_data/winter_test_predictor.npy', np.array(season_test))
np.save('../data/seasonal_data/winter_train_target.npy', np.array(season_train_dust))
np.save('../data/seasonal_data/winter_test_target.npy', np.array(season_test_dust))