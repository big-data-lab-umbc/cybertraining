# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 14:42:31 2019

@author: Jangho
"""

import numpy as np
from keras.models import load_model
from os import listdir

test_predictor=np.load('../data/test_predictor.npy')
test_target=np.load('../data/test_target.npy')

fn=['1w','5w','10w','15w','20w','25w']

for i in range(len(fn)):
    models=listdir('../CNN_models/2d32v_'+fn[i])
    for j in models:
        if 'best_model.200' in j:
            bm=j
    model=load_model('../CNN_models/2d32v_'+fn[i]+'/'+bm)
    
    res=model.predict(test_predictor)
    np.save('../data/prediction/pred_2d32v_'+fn[i]+'.npy', res[:,0])
    print(fn[i])