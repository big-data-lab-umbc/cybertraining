#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pyhdf.SD import SD, SDC
import os,datetime,sys,fnmatch
from jdcal import gcal2jd
from plot_global_map_no_colorbar import *
import math
from netCDF4 import Dataset

##########################
#### Main Program

with open('causality_ENSO_by_Ts_old.txt') as f:
    data = [map(float, line.split()) for line in f]
print(len(data))
 
nc_f= 'air.mon.mean.nc'
nc_fid=Dataset(nc_f,'r')
lat_bnd = nc_fid.variables['lat'][:]  
lon_bnd = nc_fid.variables['lon'][:]

nlon= len(lon_bnd)
nlat= len(lat_bnd)
causality_Ts = np.zeros((nlat,nlon))

nc_f2= 'land.nc'
nc_fid2=Dataset(nc_f2,'r')
land = nc_fid2.variables['land'][:]

for i in np.arange(nlat):
  for j in np.arange(nlon):
     print(i,j)
     causality_Ts[i,j]=data[i][j]

causality_Ts[causality_Ts!=1]=np.NaN

lon_bnd[nlon-1]=360.01
print('plot global map')
plot_global_map(lat_bnd,lon_bnd,causality_Ts, cmap= plt.get_cmap('rainbow'), \
            vmin=0.01,vmax=1.0001,title='Causality of ENDO by Ts', figure_name='ENSO_caused_by_Ts_wrong')
