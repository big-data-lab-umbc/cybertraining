#usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:12:15 2020

@author mkubacki and Ting

To access the correct version of Python with scikit-learn: 

module load Python/3.7.6-intel-2019a

This script creates and tests a Random Forest Classifier for the RUC sounding
data (without pressure) that has been interpolated over 37 fixed heights. It 
uses the oversamping method to address the class imbalance issue. 

This script outputs one of the following predictions

    1. nontornadic (class [0])
    2. weakly tornadic (class [1])
    3. significantly tornadic (class [2])

"""
import os
import sys
from datetime import datetime

import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

from read_data_RF import readRUCsoundings
from interpolate_avg_height import interpolate_avg_height
from interpolate_avg_height_no_pressure import interpolate_avg_height_no_pressure
import pandas as pd
import numpy as np

from collections import Counter
from sklearn.datasets import make_classification
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import SMOTE, ADASYN

################## Create Directory and output files #########################

directory_path = '/umbc/xfs1/cybertrn/cybertraining2020/team8/research/RandomForest/'

################### Load and Preprocess Data #################################

# Load Data into Pandas DataFrame
path='/umbc/xfs1/cybertrn/cybertraining2020/team8/research/RUCsoundings/files/'
RUCdata = readRUCsoundings(path)
# Interpolate Data over fixed height grid and update DataFrame
RUCdata = interpolate_avg_height_no_pressure(RUCdata)

###################  Split Training and Testing Data #########################

# Drop classifications from input data
X = RUCdata.drop('event',axis='columns')

# Create classification vector
events = pd.Categorical(RUCdata['event'],categories=['nontornadic','weakly\ntornadic','significantly\ntornadic'])
y, labels = pd.factorize(events, sort = True)

# Create Training and Testing Data
X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2)

print('    -Number of storms in the training data:', len(X_train))
print('    -Number of storms in the test data:',len(X_test))

################### Create and Train Random Forest Classifier ################

#print(Counter(y_train))

# oversampling
oversample = RandomOverSampler(sampling_strategy='minority')
#oversample = RandomOverSampler(random_state=0)
X_train, y_train = oversample.fit_resample(X_train, y_train)
#X_train, y_train = SMOTE().fit_resample(X_train, y_train)
#X_train, y_train = ADASYN().fit_resample(X_train, y_train)
print(Counter(y_train))

model = RandomForestClassifier(n_jobs=2,random_state=0,
                               class_weight="balanced",
                               n_estimators = 200,
                               criterion = "gini",
                               max_depth = 200,                                
                               )


# Train model
model.fit(X_train,y_train)

################### Evaluate model and output results ########################

# Apply model to testing data
y_pred = model.predict(X_test)

# Compute, print, and save  Accuracy Score
score= model.score(X_test,y_test)
print('Score:',score)

# Create and save Confusion Matrix
CM = pd.crosstab(y_test,y_pred,rownames=['Actual Event'],colnames=['Predicted Event'])
print('Confusion Matrix: \n',CM)

# Calculate the accuracy scores by class
y_test_l = list(y_test)
y_pred_l = list(y_pred)
res = list(zip(y_test_l, y_pred_l))
acc0 = res.count((0,0))/y_test_l.count(0)*100
acc1 = res.count((1,1))/y_test_l.count(1)*100
acc2 = res.count((2,2))/y_test_l.count(2)*100

print('Class 0 accuracy: ',acc0)
print('Class 1 accuracy: ',acc1)
print('Class 2 accuracy: ',acc2)


