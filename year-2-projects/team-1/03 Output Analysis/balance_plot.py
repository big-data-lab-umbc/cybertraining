#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: PLOT MONTHLY WATER BALANCE FROM VIC OUPUT RESULTS

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
#data = np.loadtxt('heatmap.csv',skiprows=1)
#data = np.array(pd.read_csv('heatmap.csv',header=0,delimiter=','))
data = np.array(pd.read_csv('output_data_from_VIC/balance.csv',header=0,delimiter=','))
lat = data[:,0]
lon = data[:,1]

numRow, numCol = 2,3 
month = np.zeros((numRow,numCol),dtype=object)
month[0,0] = 'April' #= data[:,2]
month[0,1] = 'May' #= data[:,3]
month[0,2] = 'June' #= data[:,4]
month[1,0] = 'July' #= data[:,5]
month[1,1] = 'August' #= data[:,6]
month[1,2] = 'September' #= data[:,7]

#mean0   = data[:,4]
#zero    = data[:,6]

lat_uniq = np.unique(lat)
lon_uniq = np.unique(lon)
dim_lat  = len(lat_uniq)
dim_lon  = len(lon_uniq)
#print(lat_uniq,lon_uniq)

llcrnrlon = -80.0 #np.min(lon_uniq)
llcrnrlat = 37.0 #np.min(lat_uniq) 
urcrnrlon = -76.0 #np.max(lon_uniq)
urcrnrlat = 41.0 #np.max(lat_uniq)

Lon,Lat = np.meshgrid(lon_uniq,lat_uniq)
map_mean0 = np.zeros((dim_lat,dim_lon),dtype=np.float)

# Initiate figures
fig1, axe1 = plt.subplots(numRow, numCol, figsize=(16,10)) 
fig1.subplots_adjust(wspace =0.2, hspace =0.2)
axlist = axe1[:,:]

#plt.figure(figsize=(25,22)) 
#plt_num = np.array([321,322,323,324,325,326])

m = Basemap(projection='cyl',fix_aspect=False, llcrnrlon = llcrnrlon, \
				llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat, \
				resolution = 'h')

data[np.where(np.isnan(data))]=0.0
val_min = np.min(data[:,2:8])
val_max = np.max(data[:,2:8])
print('data:',val_max,val_min)

data_cnt = 1		
for i in range(numRow):
	for j in range(numCol):
		
		data_cnt += 1
		print(month[i,j])
		
		#print('data:',data[:,data_cnt])
		for k in range(len(lat)):
			map_mean0[np.where((Lat == lat[k]) & (Lon == lon[k]))] = data[k,data_cnt]
		map_mean0[np.where(map_mean0 == 0.0)] = np.nan
	
		levels = 100
		
		#plt.subplot(plt_num[i-2])

		m.drawrivers(ax=axe1[i,j])
		m.etopo(ax=axe1[i,j])
		m.drawcoastlines(linewidth=0.4,ax=axe1[i,j])
		cset = axe1[i,j].contourf(Lon, Lat, map_mean0, levels, cmap=cm.jet)
		cset.set_clim(val_min,val_max)
		cg = plt.colorbar(cset,orientation='horizontal',ax=axe1[i,j])
		cg.set_label('Water Balance (mm)')
		#ax = plt.gca()
		axe1[i,j].set_xticks(np.linspace(np.min(lon_uniq),np.max(lon_uniq),5))
		axe1[i,j].set_yticks(np.linspace(np.min(lat_uniq),np.max(lat_uniq),5))
		axe1[i,j].set_xlabel("Longitude (degree)",fontsize=12)
		if ((i == 0) & (j == 0)) | ((i == 1) & (j == 0)) : axe1[i,j].set_ylabel("Latitude (degree)",fontsize=12)
		axe1[i,j].set_title(month[i,j]+' (mm)',fontsize=15)
		axe1[i,j].tick_params(axis='both', which='major', labelsize=7)

#cg = plt.colorbar(cset,orientation='horizontal')
#cg = fig1.colorbar(cset,ax=axlist)
#cg.set_clim(val_min,val_max)
#cg.set_label('Water Budget (mm)',fontsize=10)
#cg.ax.tick_params(labelsize=8)

plt.savefig('Balance_color.png',dpi=600)
plt.show()
plt.close()
