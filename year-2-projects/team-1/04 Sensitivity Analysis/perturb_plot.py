#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: PLOT CV & IQR PERTURBATION FROM VIC OUPUT RESULTS

Created on 2019

@author: Jianyu Zheng
"""
import os 
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#lat,long,mean,sd,mean0,sd0,zero
# files store output
#perturb11，perturb12，perturb21，perturb22
data = np.array(pd.read_csv('output_data_from_VIC/perturb22.csv',header=0,delimiter=','))
month = np.zeros(6,dtype=object)
lat = data[:,0]
lon = data[:,1]

print(data[28:32,2])

numRow, numCol = 1,3 
month = np.zeros(numCol,dtype=object)
month[0] = 'Original IQR' #= data[:,2]
month[1] = '1.5 x IQR' #= data[:,3]
month[2] = '2.0 x IQR' #= data[:,4]

lat_uniq = np.unique(lat)
lon_uniq = np.unique(lon)
dim_lat  = len(lat_uniq)
dim_lon  = len(lon_uniq)

llcrnrlon = -80.0 #np.min(lon_uniq)
llcrnrlat = 37.0 #np.min(lat_uniq) 
urcrnrlon = -76.0 #np.max(lon_uniq)
urcrnrlat = 41.0 #np.max(lat_uniq)

Lon,Lat = np.meshgrid(lon_uniq,lat_uniq)
map_mean0 = np.zeros((dim_lat,dim_lon),dtype=np.float)

# Initiate figures
fig1, axe1 = plt.subplots(numRow, numCol, figsize=(16,12),sharex=True,sharey=True) 
fig1.subplots_adjust(wspace =0.2, hspace =0.2)
axlist = axe1[:]

m = Basemap(projection='cyl',fix_aspect=False,  llcrnrlon = llcrnrlon, \
					llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, \
					resolution = 'h')

	
data[np.where(np.isnan(data))]=0.0
val_min = np.min(data[:,2:5])
val_max = np.max(data[:,2:5])
data[np.where(data == 0.0)] = np.nan
print('data:',val_max,val_min)

data_cnt = 1		
for j in range(numCol):
	
	data_cnt += 1
	print(month[j])
	
	for k in range(len(lat)):
		map_mean0[np.where((Lat == lat[k]) & (Lon == lon[k]))] = data[k,data_cnt]
	map_mean0[np.where(map_mean0 == 0.0)] = np.nan
	
	levels = 100

	#print(map_mean0)
	m.drawrivers(ax=axe1[j])
	m.etopo(ax=axe1[j])
	m.drawcoastlines(linewidth=0.4)
	if j == 2:
		cset1 = axe1[j].contourf(Lon, Lat, map_mean0, levels, cmap=cm.jet)
		cset1.set_clim(val_min,val_max)
	else: 
		cset = axe1[j].contourf(Lon, Lat, map_mean0, levels, cmap=cm.jet)
		cset.set_clim(val_min,val_max)

	axe1[j].set_xticks(np.linspace(np.min(lon_uniq),np.max(lon_uniq),5))
	axe1[j].set_yticks(np.linspace(np.min(lat_uniq),np.max(lat_uniq),5))
	axe1[j].set_xlabel("Longitude (degree)",fontsize=15)
	if (j == 0): axe1[j].set_ylabel("Latitude (degree)",fontsize=15)
	axe1[j].set_title(month[j]+' (mm)',fontsize=15)

cg = plt.colorbar(cset1,ax=axe1[:],orientation='horizontal')
cg.set_label('Water Budget Sensitivity (mm)',fontsize=15)
			
plt.savefig('Perturb22_color.png',dpi=600)
plt.show()
plt.close()


