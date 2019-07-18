#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 17:14:31 2019

@author: Jangho
"""

import numpy as np
import os 


#--------------------------------------------------------------------------------
# LOAD DATAFRAME OF FILENAMES
#--------------------------------------------------------------------------------
# Define base environment for training dataset
root='/umbc/xfs1/cybertrn/cybertraining2019/team5/Research/CNN/data/train_predictor_norm/'
df_train_fn=np.load('../dataframe/train_fn.npy', allow_pickle=True)
df_train_dust=np.load('../dataframe/train_dust.npy')
labels=dict(zip(df_train_fn, df_train_dust))

#--------------------------------------------------------------------------------
# PREPROCESS DATA AND DEFINE BATCH TRAINING
#--------------------------------------------------------------------------------
# This function is for image preprocessing
# This is used for slicing the data in this case, using only
# Some of the variables, of 1D CNN
def preprocess_image(img):
    # Some code for preprocessing
    img=img[2,:,:16]
    return(img)

# This function creates image_generator
def image_generator(input_ids, batch_size=100000):
    while True:
        batch_paths=np.random.choice(a=input_ids, size=batch_size)
        
        batch_input=[]
        batch_output=[]
        
        for input_id in batch_paths:
            inputx=np.load(root+input_id)
            outputx=labels[input_id]
            
            inputx=preprocess_image(inputx)
            
            batch_input += [inputx]
            batch_output += [outputx]
        
        batch_x=np.array(batch_input)
        batch_y=np.array(batch_output)
        
        yield(batch_x, batch_y)

img_ids=list(labels.keys())

# This is the generator object
train_generator=image_generator(img_ids, batch_size=100000)


# Creat 2D CNN model with 32 variables
from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, MaxPooling1D, Conv1D,\
GlobalAveragePooling2D, GlobalAveragePooling1D, Dropout, Reshape
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

model=Sequential()
model.add(Reshape((window, channel), input_shape=(window, channel)))
model.add(Conv1D(fil1, ker1, activation='relu', input_shape=(window, channel)))
model.add(Conv1D(fil2, ker2, activation='relu'))
model.add(MaxPooling1D(mpl))
model.add(GlobalAveragePooling1D())
model.add(Dropout(dr))
model.add(Dense(1, activation='sigmoid'))
print(model.summary())

model.compile(loss='binary_crossentropy',
                optimizer='adam', metrics=['accuracy'])



# FIT MODEL
batch_size=100000
class_weight = {0: 1.,
                1: 10.}
train_steps=len(img_ids)//batch_size


from keras.callbacks import Callback

class WeightsSaver(Callback):
    def __init__(self, N):
        self.N = N
        self.batch = 0

    def on_batch_end(self, batch, logs={}):
        if self.batch % self.N == 0:
            name = 'weights%08d.h5' % self.batch
            self.model.save_weights(name)
        self.batch += 1

callbacks_list = [
    keras.callbacks.ModelCheckpoint(
        filepath='../1D16V/best_model.{epoch:02d}-{val_loss:.2f}.h5')
]

history=model.fit_generator(train_generator, epochs=50, \
                    steps_per_epoch=train_steps, callbacks=callbacks_list,
                    class_weight=class_weight)
