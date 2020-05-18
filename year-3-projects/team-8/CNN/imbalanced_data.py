# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:41:21 2020

@author: Brice Coffer

Functions to handle imbalanaced datasets using the package imblearn.

input:  training_data, training_labels (3D arrays of float 64)
output: modified training_data, training_labels (3D arrays of float 64)

"""

import numpy as np
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.combine import SMOTEENN

def smoteenn(X_train, y_train):
## DOES NOT WORK CORRECTLY 
    smoteenn = SMOTEENN(random_state=42)

    n_samples, n_levels, n_variables = X_train.shape[0], \
                                       X_train.shape[1], \
                                       X_train.shape[2]

    X_train = X_train.reshape((n_samples, -1), order='F')
    X_train, y_train = smoteenn.fit_sample(X_train, y_train)
    X_train = np.reshape(X_train, (-1, n_levels, n_variables))

    return X_train, y_train

def undersample(X_train, y_train):
## DOES NOT WORK CORRECTLY 
    rus = RandomUnderSampler(random_state=42)

    n_samples, n_levels, n_variables = X_train.shape[0], \
                                       X_train.shape[1], \
                                       X_train.shape[2]

    X_train = X_train.reshape((n_samples, -1), order='F')
    X_train, y_train = rus.fit_sample(X_train, y_train)
    X_train = np.reshape(X_train, (-1, n_levels, n_variables))

    return X_train, y_train

def oversample(X_train, y_train):
## DOES NOT WORK CORRECTLY 
    ros = RandomOverSampler(random_state=42)

    n_samples, n_levels, n_variables = X_train.shape[0], \
                                       X_train.shape[1], \
                                       X_train.shape[2]

    X_train = X_train.reshape((n_samples, -1), order='F')
    X_train, y_train = ros.fit_sample(X_train, y_train)
    X_train = np.reshape(X_train, (-1, n_levels, n_variables))

    return X_train, y_train

def oversample_old(training_data, training_labels):
    ros = RandomOverSampler(random_state=42)

    training_data_th, _ = \
        ros.fit_sample(training_data[:,:,0], training_labels)
    training_data_q,  _ = \
        ros.fit_sample(training_data[:,:,1], training_labels)
    training_data_rh, _ = \
        ros.fit_sample(training_data[:,:,2], training_labels)
    training_data_u,  _ = \
        ros.fit_sample(training_data[:,:,3], training_labels)
    training_data_v, training_labels  = \
        ros.fit_sample(training_data[:,:,4], training_labels)

    training_data = np.dstack((training_data_th, training_data_q,
                               training_data_rh, training_data_u,
                                                 training_data_v))
    return training_data, training_labels

def undersample_old(training_data, training_labels):
    rus = RandomUnderSampler(random_state=42)

    training_data_th, _ = \
        rus.fit_sample(training_data[:,:,0], training_labels)
    training_data_q,  _ = \
        rus.fit_sample(training_data[:,:,1], training_labels)
    training_data_rh, _ = \
        rus.fit_sample(training_data[:,:,2], training_labels)
    training_data_u,  _ = \
        rus.fit_sample(training_data[:,:,3], training_labels)
    training_data_v, training_labels  = \
        rus.fit_sample(training_data[:,:,4], training_labels)

    training_data = np.dstack((training_data_th, training_data_q,
                               training_data_rh, training_data_u,
                                                 training_data_v))
    return training_data, training_labels
