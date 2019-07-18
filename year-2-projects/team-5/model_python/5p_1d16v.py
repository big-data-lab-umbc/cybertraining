#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:32:29 2019

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
train_predictor=np.load('../data/5p_data/train_predictor_5p.npy')[:,2,:,:16]
train_target=np.load('../data/5p_data/train_target_5p.npy')
# Load test data
test_predictor=np.load('../data/5p_data/test_predictor_5p.npy')[:,2,:,:16]
test_target=np.load('../data/5p_data/test_target_5p.npy')



#---------------------------------------------------------------------------------------------------------
# ----- Create and Compile CNN Model -----
#---------------------------------------------------------------------------------------------------------

from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, GlobalAveragePooling2D, Dropout, Reshape
from keras.layers import Conv1D, MaxPooling1D, GlobalAveragePooling1D
import keras
from keras.utils import to_categorical

surr=5
window=5
channel=16
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

# Model
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
        filepath='../CNN_models/5p_1d16v/best_model.{epoch:02d}-{val_loss:.2f}.h5')]


#---------------------------------------------------------------------------------------------------------
# ----- FIT CNN Model -----
#---------------------------------------------------------------------------------------------------------
history=model.fit(train_predictor, train_target, batch_size=2000, epochs=100,\
                  callbacks=callbacks_list, class_weight=class_weight,\
                  validation_data=[test_predictor, test_target], verbose=1)

htry=np.array([history['acc'], history['loss'], history['val_acc'], history['vall_loss']])
np.save('../CNN_history/5p_1d16v.npy', htry)