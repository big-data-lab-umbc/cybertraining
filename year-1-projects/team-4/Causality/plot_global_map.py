#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np
from pyhdf.SD import SD, SDC
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as CS
import os,datetime,sys,fnmatch
import mpl_toolkits.basemap as bm

def plot_global_map(lat_grids, lon_grids, data,\
             proj = 'cyl',\
             cmap = plt.cm.gray_r,\
             vmin = None, vmax = None,\
             Log_Norm = False,\
             title=None, figure_name=None):

    mapproj = bm.Basemap(projection='cyl',
              llcrnrlat=lat_grids.min(),\
              llcrnrlon=lon_grids.min(),\
              urcrnrlat=lat_grids.max(),\
              urcrnrlon=lon_grids.max())

    nx, ny = lon_grids.size,lat_grids.size
    lonall, latall = mapproj.makegrid(nx, ny)
    lonproj, latproj = mapproj(lonall, latall)
    latlines = np.arange(-90,100,30.0)
    lonlines = np.arange(-180.0,180.0,60.0)

    fig,ax = plt.subplots()
    if vmin is None: vmin = np.nanmin(data)
    if vmax is None: vmax = np.nanmax(data)

    if not Log_Norm:
        norm = CS.Normalize(vmin=vmin,vmax=vmax)
    else:
        norm = CS.LogNorm(vmin=vmin, vmax=vmax)
    #check if the latitude direction is ascending
    if lat_grids[-1] < lat_grids[0]: data = data[::-1,:]
    masked_data= np.ma.masked_where(np.isnan(data),data)

    ctf = mapproj.pcolormesh(lonproj,latproj,masked_data,\
                        norm=norm,cmap=cmap)
    plt.title(title)
    mapproj.drawcoastlines()
    mapproj.drawparallels(latlines, labels=[1,0,0,0])
    mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
    #mapproj.fillcontinents(color='0.5')
    cb = mapproj.colorbar(ctf,"bottom", size="5%",pad = '10%')
    plt.savefig(figure_name+'.png', bbox_inches='tight',dpi=200)
    plt.close(fig)
