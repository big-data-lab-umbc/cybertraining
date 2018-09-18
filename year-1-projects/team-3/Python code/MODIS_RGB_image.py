from __future__ import print_function,division
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
from pyhdf.SD import SD, SDC
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

#ModPred = np.loadtxt('ModPred.txt')
ModPred = np.loadtxt('622predict2.txt')
print(ModPred.shape)
ModPred = np.ravel(ModPred)
ModPred = np.reshape(ModPred,(2030,1354))
print(ModPred.shape)
'''
fig = plt.figure()
img = plt.imshow(ModPred,cmap=cm.PRGn,vmin=-1,vmax=1)
plt.colorbar()
plt.show()
'''

#read MODIS file
def Read_MODIS(filename):
    MYD = SD(filename, SDC.READ)
    lat = MYD.select('Latitude').get()
    lon = MYD.select('Longitude').get()
    return lat,lon

lat,lon = Read_MODIS('/home/pshi1/team3/Research/Progress3_DustStorm/22June2009_PhyAgm/MYD021KM.A2009173.1505.061.2018044055419.hdf')
print(lat.shape,lon.shape)

y = np.arange(0,2026,5)
x = np.arange(0,1351,5)
print(x.shape,y.shape)
y_new = np.arange(0,2025,1)
x_new = np.arange(0,1350,1)
print(x_new.shape,y_new.shape)
f_lat = interpolate.interp2d(x,y,lat,kind = 'linear')
f_lon = interpolate.interp2d(x,y,lon,kind = 'linear')

lat_new = f_lat(x_new,y_new)
lon_new = f_lon(x_new,y_new)
print(lat_new.shape,lon_new.shape)
   
ModPred_interest = ModPred[:2025,:1350]
print(ModPred_interest.shape)

#plot projected image
ll_lat = np.min(lat)
ll_lon = np.min(lon)
ur_lat = np.max(lat)
ur_lon = np.max(lon)

fig,ax = plt.subplots(figsize=(11,8))
map=Basemap(projection = 'cyl', llcrnrlat = ll_lat, urcrnrlat = ur_lat, llcrnrlon = ll_lon, urcrnrlon = ur_lon)
map.drawparallels(np.arange(ll_lat,ur_lat,10),labels=[1,0,0,0],fontsize=10)
map.drawmeridians(np.arange(ll_lon,ur_lon,20),labels=[0,0,0,0],fontsize=10)
map.pcolormesh(lon_new,lat_new,ModPred_interest,cmap='binary',vmin=0,vmax=1.0)
#divider = make_axes_locatable(ax)
plt.colorbar(orientation='vertical',pad=0.05,fraction = 0.1,aspect=20)
fig.savefig('/home/pshi1/test.png', dpi=300)

