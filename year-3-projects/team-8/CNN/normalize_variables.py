# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:54:37 2020

@author: Brice Coffer

Function to normalize each variables from 0 to 1

input:  3D array of original RUC soundings (Array of float64)
output: 3D array of normalized RUC soundings (Array of float64)

"""

import numpy as np
from sklearn.preprocessing import StandardScaler

def normalize_data(X_train, X_test): 
       
    # normalize data
    for i in range(X_train.shape[2]):
         sc = StandardScaler()
         X_train[:, :, i] = sc.fit_transform(X_train[:,:,i])
         X_test[:, :, i] = sc.transform(X_test[:,:,i])

    return X_train, X_test 
