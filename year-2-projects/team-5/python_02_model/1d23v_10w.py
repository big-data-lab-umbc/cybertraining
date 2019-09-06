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
train_predictor=np.load('../data/train_predictor.npy')[:,2,:,:23]
train_target=np.load('../data/train_target.npy')
# Load test data
test_predictor=np.load('../data/test_predictor.npy')[:,2,:,:23]
test_target=np.load('../data/test_target.npy')



#---------------------------------------------------------------------------------------------------------
# ----- Create and Compile CNN Model -----
#---------------------------------------------------------------------------------------------------------

from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, MaxPooling1D, Conv1D,\
GlobalAveragePooling2D, GlobalAveragePooling1D, Dropout, Reshape
import keras
from keras.utils import to_categorical


surr=5
window=5
channel=23
ker1= 3
ker2= 2
ker3=2
fil1= 50
fil2= 100
mpl= 2
dr= 0.5
num_classes=2
class_weight = {0: 1.,
                1: 10.}

model=Sequential()
model.add(Reshape((window, channel), input_shape=(window, channel)))
model.add(Conv1D(fil1, ker1, activation='relu', input_shape=(window, channel)))
model.add(Conv1D(fil2, ker2, activation='relu'))
model.add(MaxPooling1D(mpl))
model.add(GlobalAveragePooling1D())
model.add(Dropout(dr))
model.add(Dense(1, activation='sigmoid'))
print(model.summary())

model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])


callbacks_list = [keras.callbacks.ModelCheckpoint(\
        filepath='../CNN_models/1d23v_10w/best_model.{epoch:02d}-{val_acc:.2f}.h5')]

#---------------------------------------------------------------------------------------------------------
# ----- FIT CNN Model -----
#---------------------------------------------------------------------------------------------------------

history=model.fit(train_predictor, train_target, batch_size=2**12, epochs=200,\
                  callbacks=callbacks_list, class_weight=class_weight,\
                  validation_data=[test_predictor, test_target], verbose=1)

from keras import backend as K 
K.clear_session()
