"""
Extracting sample data from original data

Original: /data4/djin1/MODIS_C6/Gridded/aqua_d3_c6_tvp_pcl.12012002_4799dy.cf0
          [4799 days * 180 lats * 360 lons * 42 vars], float32

Output: Excluding missings and clear-sky (CF=0%)
        [nn * 42 vars]
"""

import numpy as np
import sys
import os.path
from datetime import timedelta, date


### Original File Parameter
ndy=4799; nlat=180; nlon=360; nvar=42
idate=date(2002,12,1) 
indir='/data4/djin1/MODIS_C6/Gridded/'
infn = indir+'aqua_d3_c6_tvp_pcl.12012002_4799dy.cf0'

### Target File Parameter
tgdate1=date(2005,1,1); nt1=(tgdate1-idate).days 
tgdate2=date(2005,12,31); nt2=(tgdate2-idate).days+1
tglat=(75,105); nlat2=tglat[1]-tglat[0]  ## -15S to 15N

### Read original file
chist = np.memmap(infn,dtype=np.float32,mode='r',shape=(ndy,nlat,nlon,nvar))
chist = np.array(chist[nt1:nt2,tglat[0]:tglat[1],:,:])
print('Dimension1=',chist.shape)
print('Min value= {}, Max value={}'.format(chist.min(),chist.max()))

### Filtering
tcf = chist.sum(axis=3)
chistidx = np.logical_and(tcf>0.,tcf<=1.)
chist=chist[chistidx,:]
print('Dimension2=',chist.shape)
print('Min value= {}, Max value={}'.format(chist.min(),chist.max()))

outdir = '/data1/djin1/Scratch/'
outfn = outdir+'aqua_d3_c6_tvp_pcl.noMissing.{}-{}_{}x{}.float32.dat'.format(tgdate1.strftime('%Y%m%d'),tgdate2.strftime('%Y%m%d'),*chist.shape)

with open(outfn,'wb') as fout:
    chist.tofile(fout)


