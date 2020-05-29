import os
from datetime import datetime
from datetime import timedelta

import matplotlib
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

days = [] 
day = datetime(2007, 1, 1)

path = "/umbc/xfs1/cybertrn/cybertraining2020/team2/research/SeaIce/NHData"
#path = os.getcwd()

while day <= datetime(2007, 1, 31):
    
    days.append(day)
    day += timedelta(days = 1)
    
for day in tqdm(days, total = len(days)):

    filename = os.path.join(path, f"asicd25e2_{day.strftime('%Y%m%d')}_v01r02.nc")

    dataset = nc.Dataset(filename)

    latitudes = dataset["latitude"][:].data
    longitudes = dataset["longitude"][:].data
    sea_ice = dataset["sea_ice_cover"][:].data.astype(np.float32)[0, :, :]

    latitudes[latitudes == dataset["latitude"]._FillValue] = np.nan
    longitudes[longitudes == dataset["longitude"]._FillValue] = np.nan
    sea_ice[sea_ice == dataset["sea_ice_cover"]._FillValue] = np.nan

    latbins = np.arange(0, 90 + 0.5, 0.5)
    lonbins = np.arange(-180, 180, 0.5)

    latbins2 = np.arange(0, 90 + 1, 1)
    lonbins2 = np.arange(-180, 180, 1)
    
    z, x, y = np.histogram2d(longitudes.flatten()[~np.isnan(longitudes.flatten())],
                             latitudes.flatten()[~np.isnan(latitudes.flatten())], (lonbins2, latbins2)) 
    z = z.T
    
    latitudes = latitudes[sea_ice == 30]
    longitudes = longitudes[sea_ice == 30]
    sea_ice = sea_ice[sea_ice == 30]
    #
    #Xm = ma.masked_invalid(longitudes)
    #Ym = ma.masked_invalid(latitudes)
    #Zm = ma.masked_invalid(sea_ice)

    
    z1, x1, y1 = np.histogram2d(longitudes, latitudes, (lonbins, latbins)) 
    z1 = z1.T
        
    df = pd.read_csv("2007-01_water_anomalies.csv")
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    df = df[(day <= df["timestamp"]) & (df["timestamp"] < day + timedelta(days = 1))]
    df["latitude"] = df["latitude"]
    df["longitude"] = df["longitude"]

    z2, x2, y2 = np.histogram2d(df["longitude"], df["latitude"], (lonbins, latbins))
    z2 = z2.T

    mask = np.where((z1 != 0) & (z2 != 0) & (~np.isnan(z1)))

    x3 = x1[mask[1]]
    y3 = y1[mask[0]]
    
    if x3.size > 0:
    
        cmap = plt.get_cmap("gray")
        cmap.set_bad('k', alpha = 1.)
        
        with np.errstate(divide='ignore', invalid='ignore'):
        
            z4, x4, y4 = np.histogram2d(longitudes, latitudes, (lonbins2, latbins2)) 
            z4 = z4.T / z

        fig, ax = plt.subplots(1, 1, figsize=(8,4),
                               subplot_kw={'projection': ccrs.PlateCarree()})
        ax.set_title("Anomaly / Sea Ice Collocations - " + day.strftime("%Y-%m-%d"))
        ax.set_extent([-180, 180, 0, 90], ccrs.PlateCarree())
        ax.coastlines(color = "b")
        im = ax.pcolormesh(x4, y4, z4, cmap = cmap, vmin = 0, vmax = 1)
        fig.colorbar(im, orientation="horizontal", pad=0.2)
        ax.scatter(x3, y3, c = "r", edgecolors = "k", zorder = 3)
        plt.savefig(f"{day.strftime('%Y-%m-%d')}.png")
        plt.close()
