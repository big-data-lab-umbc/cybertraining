#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.api import VAR, DynamicVAR
from netCDF4 import Dataset

def causality_test(var1,var2):
    data1= pd.Series(var1, name='Var1')
    data2= pd.Series(var2, name='Var2')
    mdata=pd.concat([data1, data2], axis=1)
    mdata.index = pd.date_range('1950-01-01', periods=600, freq='M')
    
    model = VAR(mdata)
    results=model.fit(7)
   
    foo=crit=results.test_causality('Var2', ['Var1'], kind='f')
    crit=foo['crit_value']
    stat=foo['statistic']
    if(stat > crit):
       cause=1
    else:
       cause=0

    return cause


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
result = np.zeros((nlat,nlon))

for i in np.arange(nlat):
  for j in np.arange(nlon):
    var2=w500[24:624,i,j]
    var1=ENSO[0:600]
    result[i,j]=causality_test(var1,var2)

file = open('causality_w500_by_ENSO.txt','w')
np.savetxt(file,result,'%4.2f') 
file.close() 

