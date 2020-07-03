
#process ERA5 monthly data from 1979 to 2019
#Check missing values
#Average data in the Arctic (>60N) and save time series into csv file


import sys
import csv
import numpy as np
from netCDF4 import Dataset
import pandas as pd

###### read data from nc file #####

filename='ERA5_monthly_clouds_1979_2019_1deg.nc'
yy=filename
data=Dataset(yy)
time=np.array(data.variables['time'][:])
lat=np.array(data.variables['latitude'][:])	
lon=np.array(data.variables['longitude'][:])
data0=data.variables['tcc'][:,:,:]
data1=data.variables['tciw'][:,:,:]
data2=data.variables['tclw'][:,:,:]


#check missing values
missing=np.where(data0 == data0.fill_value)
if missing[0].size == 0:
	print('total cloud cover: No missing values')
else:
	sys.exit()

missing1=np.where(data1 == data1.fill_value)
if missing1[0].size == 0:
        print('cloud ice water: No missing values')
else:
        sys.exit()

missing2=np.where(data2 == data2.fill_value)
if missing2[0].size == 0:
        print('cloud liquid water: No missing values')
else:
        sys.exit()


nmons=time.size

#unit conversion
cf=np.array(data0)*100. #from 0-1 to percent
cwp=(np.array(data1)+np.array(data2))*1000. #from kg to g


###### calcuate Arctic domain average using weighted mean #####
arc=lat[lat>60]
arcr=np.deg2rad(arc)
weight=np.cos(arcr)
cf_zonal=np.mean(cf,axis=2)
cwp_zonal=np.mean(cwp,axis=2)

m0=[]
y0=[]
val0=[]
val1=[]

year=np.arange(1979,2020,1,dtype=int)

for i in range(0,nmons):
	cfmean=np.average(cf_zonal[i,lat>60],weights=weight)
	cwpmean=np.average(cwp_zonal[i,lat>60],weights=weight)
	yy0=i//12
	yy=year[yy0]
	mm=np.mod(i,12)+1	
	y0.append(yy)
	m0.append(mm)
	val0.append(cfmean)
	val1.append(cwpmean)

###### write data into csv file #####
dict={'Year':y0, 'Month':m0,'Total cloud cover (%)':val0, 'Total column cloud water (g m**-2)':val1}
data=pd.DataFrame(dict)
data.to_csv('ERA5_Arctic_cloud_variables_1979_2019.csv',header=True,index=True)




