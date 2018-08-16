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

with open('causality_SLP_by_ENSO.txt') as f:
    data = [map(float, line.split()) for line in f]
print(len(data))
 
nc_f= 'air.mon.mean.nc'
nc_fid=Dataset(nc_f,'r')
lat_bnd = nc_fid.variables['lat'][:]  # extract/copy the data
lon_bnd = nc_fid.variables['lon'][:]

nlon= len(lon_bnd)
nlat= len(lat_bnd)
causality_SLP = np.zeros((nlat,nlon))

nc_f2= 'land.nc'
nc_fid2=Dataset(nc_f2,'r')
land = nc_fid2.variables['land'][:]

for i in np.arange(nlat):
  for j in np.arange(nlon):
     print(i,j)
     causality_SLP[i,j]=data[i][j]*(1-land[0,i,j])

causality_SLP[causality_SLP!=1]=np.NaN

lon_bnd[nlon-1]=360.01
print('plot global map')
plot_global_map(lat_bnd,lon_bnd,causality_SLP, cmap= plt.get_cmap('rainbow'), \
            vmin=0.01,vmax=1.0001,title='Granger Causality, ENSO ->SLP', figure_name='SLP_caused_by_ENSO')
