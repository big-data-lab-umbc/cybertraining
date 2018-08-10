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

ENSO=np.loadtxt("ENSO_index_195001_201712.txt")    # ENSO index

nc_f= 'slp.mon.mean.nc'
nc_fid=Dataset(nc_f,'r')
lats = nc_fid.variables['lat'][:]  
lons = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][:]
slp = nc_fid.variables['slp'][:,:,:] 

nlon= len(lons)
nlat= len(lats)
result1 = np.zeros((nlat,nlon))

for i in np.arange(nlat):
  for j in np.arange(nlon):
    var1=slp[24:624,i,j]
    var2=ENSO[0:600]
    print(i,j)
    result1[i,j]=causality_test(var1,var2)

file = open('causality_ENSO_by_SLP.txt','w')
np.savetxt(file,result1,'%4.2f') 
file.close() 

print(np.sum(result1))


