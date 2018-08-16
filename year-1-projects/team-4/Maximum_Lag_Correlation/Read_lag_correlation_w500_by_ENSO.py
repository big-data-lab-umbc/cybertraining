#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.api import VAR, DynamicVAR
from netCDF4 import Dataset


def lag_correlation(var1,var2,nlag):

    ntime=len(var1)
    max_cor=0.0        # maximum correlation coefficient
    max_lag=0.0        # Lag number with maximum correlation coefficient

    for i in np.arange(nlag):
       cor = np.corrcoef(var1[0:ntime-i-1], var2[i:ntime-1])[0, 1]
       if (abs(cor) > abs(max_cor)):
          max_cor = cor
          max_lag = i
   
    max_ind = var2.argmax()    # Index of time with maximum value of Var2
    max_ind = (max_ind+1)%12   # Convert the number to 1-12 (Jan-Dec)
    if (max_ind ==0): max_ind=12
 
    return (max_cor,max_lag,max_ind)  


################################
# Main program starts here #

ENSO=np.loadtxt("ENSO_index_195001_201712.txt")

nc_f= 'omega.mon.mean.nc'
nc_fid=Dataset(nc_f,'r')
lats = nc_fid.variables['lat'][:] 
lons = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][:]
w500 = nc_fid.variables['omega'][:,5,:,:]*864.0

nlon= len(lons)
nlat= len(lats)
max_cor2 = np.zeros((nlat,nlon))
max_lag2 = np.zeros((nlat,nlon))
max_ind2 = np.zeros((nlat,nlon))

nc_f2= 'land.nc'
nc_fid2=Dataset(nc_f2,'r')
land = nc_fid2.variables['land'][:]

for i in np.arange(nlat):
  for j in np.arange(nlon):
   #if (land[0,i,j] < 1):
     var2=w500[24:624,i,j]
     var1=ENSO[0:600]
     results=lag_correlation(var1,var2,7)
     max_cor2[i,j]=results[0]
     max_lag2[i,j]=float(results[1])
     max_ind2[i,j]=float(results[2]) 

file = open('max_cor_w500_lagged_by_ENSO.txt','w')
np.savetxt(file,max_cor2,'%10.6f') 
file.close() 

file = open('max_lag_w500_lagged_by_ENSO.txt','w')
np.savetxt(file,max_lag2,'%6.2f')
file.close()

file = open('max_ind_w500_lagged_by_ENSO.txt','w')
np.savetxt(file,max_ind2,'%6.2f')
file.close()

