#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:29:25 2019

@author: jangho.lee.92
"""

import numpy as np

dust=np.load('../data/target_classes.npy')

classes=[]
for i in range(np.shape(dust)[0]):
    if dust[i,0] in classes:
        pass
    else:
        classes.append(dust[i,0])
classes.sort()
count=np.zeros(len(classes))

for i in range(np.shape(dust)[0]):
    ind=classes.index(dust[i,0])
    count[ind]=count[ind]+1