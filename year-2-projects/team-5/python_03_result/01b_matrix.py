# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 15:08:26 2019

@author: Jangho
"""

import numpy as np

test_target=np.load('../data/test_target.npy')

fn1=['1d16v','1d20v','1d23v','1d32v','2d16v','2d20v','2d23v','2d32v']
fn2=['10w','20w']

for i in range(len(fn1)):
    for j in range(len(fn2)):
        pred=np.load('../data/prediction/pred_'+fn1[i]+'_'+fn2[j]+'.npy')
        pred[pred>=0.5]=1
        pred[pred<0.5]=0
        
        dd=(test_target+pred==2).sum()*100/np.sum(test_target)
        dnd=100-dd
        
        ndnd=(test_target+pred==0).sum()*100/(len(test_target)-np.sum(test_target))
        ndd=100-ndnd
        
        acc=((test_target+pred==2).sum()+(test_target+pred==0).sum())*100/len(test_target)
            
        print(fn1[i], fn2[j], dd, dnd, ndd, ndnd, acc)
