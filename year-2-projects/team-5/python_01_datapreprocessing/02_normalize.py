#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:20:00 2019

@author: jangho.lee.92
"""

import numpy as np


#---------------------------------------------------------------------------------------------------------
# -----  LOAD DATA  -----
# This data is for only 5% of total data due to overflowing memory problem
# 
# train_predictor, train_target, test_predictor, test_target
#---------------------------------------------------------------------------------------------------------

# Load training data
train_predictor=np.concatenate([np.load('../data/2gap/spring_predictor_train.npy'),\
                                np.load('../data/2gap/summer_predictor_train.npy'),\
                                np.load('../data/2gap/fall_predictor_train.npy'),\
                                np.load('../data/2gap/winter_predictor_train.npy')], axis=0)
print("Loaded Data")
train_target=np.concatenate([np.load('../data/2gap/spring_target_train.npy'),\
                                np.load('../data/2gap/summer_target_train.npy'),\
                                np.load('../data/2gap/fall_target_train.npy'),\
                                np.load('../data/2gap/winter_target_train.npy')], axis=0)
print("Loaded Data")
# Load test data
test_predictor=np.concatenate([np.load('../data/2gap/spring_predictor_test.npy'),\
                                np.load('../data/2gap/summer_predictor_test.npy'),\
                                np.load('../data/2gap/fall_predictor_test.npy'),\
                                np.load('../data/2gap/winter_predictor_test.npy')], axis=0)
print("Loaded Data")
test_target=np.concatenate([np.load('../data/2gap/spring_target_test.npy'),\
                                np.load('../data/2gap/summer_target_test.npy'),\
                                np.load('../data/2gap/fall_target_test.npy'),\
                                np.load('../data/2gap/winter_target_test.npy')], axis=0)
print("Loaded Data")

np.save('../data/train_predictor_nonorm.npy', train_predictor)
np.save('../data/test_predictor_nonorm.npy', test_predictor)
np.save('../data/train_target_nonorm.npy',train_target)
np.save('../data/test_target_nonorm.npy', test_target)

means1=np.load('../data/means/spring.npy')
means2=np.load('../data/means/summer.npy')
means3=np.load('../data/means/fall.npy')
means4=np.load('../data/means/winter.npy')
means=np.zeros([3,32])
means[0,:]=(means1[0,:]+means2[0,:]+means3[0,:]+means4[0,:])/4.
means[1,:]=np.max(np.concatenate([means1[1,:].reshape(1,32),\
                                 means2[1,:].reshape(1,32),\
                                 means3[1,:].reshape(1,32),\
                                 means4[1,:].reshape(1,32)], axis=0), axis=0)
means[2,:]=np.min(np.concatenate([means1[2,:].reshape(1,32),\
                                 means2[2,:].reshape(1,32),\
                                 means3[2,:].reshape(1,32),\
                                 means4[2,:].reshape(1,32)], axis=0), axis=0)

m=means[0,:]
r=means[1,:]-means[2,:]

np.save('../data/train_predictor.npy', (train_predictor-m)/r)
np.save('../data/test_predictor.npy', (test_predictor-m)/r)
np.save('../data/train_target.npy',train_target[:,0])
np.save('../data/test_target.npy', test_target[:,0])