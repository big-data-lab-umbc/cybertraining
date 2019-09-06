# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 15:08:26 2019

@author: Jangho
"""

import numpy as np

test_target=np.load('../data/test_target.npy')

fn=['1w','5w','10w','15w','20w','25w']

for i in range(len(fn)):
    pred=np.load('../data/prediction/pred_2d32v_'+fn[i]+'.npy')
    pred[pred>=0.5]=1
    pred[pred<0.5]=0
    
    dd=(test_target+pred==2).sum()*100/np.sum(test_target)
    dnd=100-dd
    
    ndnd=(test_target+pred==0).sum()*100/(len(test_target)-np.sum(test_target))
    ndd=100-ndnd
    
    acc=((test_target+pred==2).sum()+(test_target+pred==0).sum())*100/len(test_target)
    
    print(dd, dnd, ndd, ndnd, acc)