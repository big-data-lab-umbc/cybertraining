#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:56:42 2019

@author: jangho.lee.92
"""

import numpy as np


# READ DATAFRAME and INITIALIZE ARRAY
fn_train=np.load('../dataframe/train_fn.npy')
dust_train=np.load('../dataframe/train_dust.npy')
fn_test=np.load('../dataframe/test_fn.npy')
dust_test=np.load('../dataframe/test_dust.npy')

train_sample_org=len(fn_train)
test_sample_org=len(fn_test)
sample_ratio=5
train_sample=int(train_sample_org//sample_ratio)
test_sample=int(test_sample_org//sample_ratio)

fn_train_5p=np.zeros([train_sample, 5, 5, 32])
dust_train_5p=np.zeros([train_sample])
fn_test_5p=np.zeros([test_sample, 5, 5, 32])
dust_test_5p=np.zeros([test_sample])


# CREATE ARRAY WITH ONLY 5% OF TOTAL DATA
# THIS IS ONLY TO REDUCE TIME AND SEE THE RESULTS
root1='/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/CNN/data/train_predictor_norm/'
root2='/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/CNN/data/test_predictor_norm/'
"""
# CREATE TRAIN ARRAY
for i in range(train_sample):
    fn=fn_train[i*sample_ratio]
    dust=dust_train[i*sample_ratio]
    fn_train_5p[i,:,:,:]=np.load(root1+fn)
    dust_train_5p[i]=dust
    
    if i//1000==i/1000:
        print('TRAIN: '+str(i*100/train_sample)[:5]+' %', i, flush=True)
"""
# CREATE TEST ARRAY
for i in range(test_sample):
    fn=fn_test[i*sample_ratio]
    dust=dust_test[i*sample_ratio]
    fn_test_5p[i,:,:,:]=np.load(root2+fn)
    dust_test_5p[i]=dust
    
    if i//1000==i/1000:
        print('TEST: '+str(i*100/test_sample)[:5]+' %', i, flush=True)

#np.save('../data/20p_data/train_predictor_20p.npy', fn_train_5p)
np.save('../data/20p_data/test_predictor_20p.npy', fn_test_5p)
#np.save('../data/20p_data/train_target_20p.npy', dust_train_5p)
np.save('../data/20p_data/test_target_20p.npy', dust_test_5p)
