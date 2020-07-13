#usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:12:15 2020

@author mkubacki

This script creates and tests a Random Forest Classifier for the RUC sounding
data WITHOUT THE PRESSURE VARIABLE that has been interpolated over 37 fixed heights.  
The Random Forest Classifier uses the 6x37 features

'temperature 1-37', 'dewpoint 1-37', 'humidity 1-37', 'u-wind 1-37',
'v-wind 1-37'

to create decision trees and output one of the following predictions

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
from interpolate_avg_height_no_pressure import interpolate_avg_height_no_pressure
import pandas as pd
import numpy as np

################## Create Directory and output files #########################

now = datetime.now()
now_string = now.strftime("%Y-%b%d-%H:%M:%S/")
print(now_string)
directory_path = '/umbc/xfs1/cybertrn/cybertraining2020/team8/research/RandomForest/'
output_directory = directory_path + now_string

os.mkdir(output_directory)
FI_file = open(output_directory + 'feature_importance.txt','w')
sum_file = open(output_directory + 'summary.txt','w')
 
sum_file.write('Summary of Random Forest Model output from ' + now_string + '\nModel information:\n')

################### Load and Preprocess Data #################################

# Load Data into Pandas DataFrame
path='/umbc/xfs1/cybertrn/cybertraining2020/team8/research/RUCsoundings/files/'
RUCdata = readRUCsoundings(path)
# Interpolate Data over fixed height grid and update DataFrame
RUCdata = interpolate_avg_height_no_pressure(RUCdata)
print('    -Number of height levels is 37 using height averages. \n',file = sum_file)
features = RUCdata.columns[0:5*37]
print('    -Number of features is ', len(features), '\n',file = sum_file)

###################  Split Training and Testing Data #########################

# Drop classifications from input data
X = RUCdata.drop('event',axis='columns')

# Create classification vector
events = pd.Categorical(RUCdata['event'],categories=['nontornadic','weakly\ntornadic','significantly\ntornadic'])
y, labels = pd.factorize(events, sort = True)

# Create Training and Testing Data
X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2)
print('    -Number of storms in the training data:', len(X_train),'\n',file = sum_file)
print('    -Number of storms in the test data:',len(X_test),'\n',file = sum_file)


################### Create and Train Random Forest Classifier ################

# Create model
model = RandomForestClassifier(n_jobs=2,random_state=0,
			 	class_weight="balanced",
				n_estimators = 200,
				criterion = "gini",
				max_depth = 200,				
				)

# Output model parameters into summary file
print(model,'\n', file = sum_file)

# Train model
model.fit(X_train,y_train)

################### Evaluate model and output results ########################

# Apply model to testing data
y_pred = model.predict(X_test)

# Compute, print, and save  Accuracy Score
score= model.score(X_test,y_test)
print('Percent Accuracy is ',score*100,'\n',file = sum_file)
print('Score:',score)

# Create and save Confusion Matrix
CM = pd.crosstab(y_pred,y_test,rownames=['Predicted Event'],colnames=['Actual Event'])
print('Confusion Matrix: \n',CM,file = sum_file)

################### Feature importance analysis ##############################

# Compute feature importance scores and save to file
feature_list = list(features)
importances = list(model.feature_importances_)
feature_importances = [(feature, round(importance, 4)) for feature, importance in zip(feature_list, importances)]
feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
[print('Variable: {:20} Importance: {}'.format(*pair),file = FI_file) for pair in feature_importances]
FI_file.close

# Compute number of features needed to preserve 95% importance and save to file
sorted_importances = [importance[1] for importance in feature_importances]
sorted_features = [importance[0] for importance in feature_importances]
cumulative_importances = np.cumsum(sorted_importances)
importance95 = np.where(cumulative_importances>0.95)[0][0]+1
print('-Number of features for 95% importance:',importance95,'\n',file = sum_file)
sum_file.close


