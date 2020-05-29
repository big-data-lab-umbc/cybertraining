#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Oct 23 10:32:01 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
A script to downlod *.hdf5 files from https://opendap.larc.nasa.gov/opendap/
"""

import numpy as np
import os
from pyhdf.SD import SD, SDC

fname='CAL_LID_L2_01kmCLay-Standard-V4010.2007_data.txt'
fl=open(fname,'r')

str1=fl.readlines()[:]
fn=[]
for i in np.arange(0,np.size(str1)):
    fn.append(str1[i].split('\n',1)[0])
    
    
#wget files

linkf="https://opendap.larc.nasa.gov/opendap/CALIPSO/LID_L2_01kmCLay-Standard-V4-10/2007/"
outdr="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_01kmCLay/2007/"

missing_count=0
for i in np.arange(0,np.size(fn)-1,1):
    if i<2:
        link=linkf+fn[i].split('V4-10.2006',1)[1][1:3]+'/'+fn[i]
    else:
        link=linkf+fn[i].split('V4-10.2007',1)[1][1:3]+'/'+fn[i]
    print(link)
    if not(os.path.isfile(outdr+fn[i])):
        try:
            CAL1km=SD(outdr+fn[i])
        except Exception as e:
            print(e)
            print('missing: %d'%(missing_count))
            missing_count+=1
            os.system("wget "+link+' -P '+outdr)
    else:
        print('Already downloaded.')