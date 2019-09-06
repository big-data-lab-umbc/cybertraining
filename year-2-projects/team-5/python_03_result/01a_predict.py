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

fn1=['1d16v','1d20v','1d23v','1d32v','2d16v','2d20v','2d23v','2d32v']
fn2=['10w','20w']

for i in range(len(fn1)):
    for j in range(len(fn2)):
        
        # Slice dataset
        if '1d' in fn1[i]:
            test_predictor_slice=test_predictor[:,2,:,:]
            
            if '16v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:16]
            elif '20v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:20]
            elif '23v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:23]
            else:
                pass
        else:
            test_predictor_slice=test_predictor[:,:,:,:]
            if '16v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:,:16]
            elif '20v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:,:20]
            elif '23v' in fn1[i]:
                test_predictor_slice=test_predictor_slice[:,:,:,:23]
            else:
                pass

        
        # Load model and predict
        models=listdir('../CNN_models/'+fn1[i]+'_'+fn2[j])
        for k in models:
            if 'best_model.200' in k:
                bm=k
        model=load_model('../CNN_models/'+fn1[i]+'_'+fn2[j]+'/'+bm)
        
        res=model.predict(test_predictor_slice)
        np.save('../data/prediction/pred_'+fn1[i]+'_'+fn2[j]+'.npy', res[:,0])
        print(fn1[i]+'_'+fn2[j])