# Standard library imports.
from pprint import pprint

# Third party imports.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from pyhdf.SD import SD
from pyhdf.SD import SDC
import my_lib as ml
import pandas as pd
import netCDF4 as nc
import math
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

SPATIAL_RESOLUTION = 1.25
def editTimestamp(x):

    head, sep, tail = str(x).partition(' ')
    head = head.replace("-","")
    return head

def findGridPosition(lat, lon):

    latRange = np.arange(-90,90 + SPATIAL_RESOLUTION, SPATIAL_RESOLUTION)
    latRange = np.append(latRange, lat)
    latRange.sort()
            
    lonRange = np.arange(0,360 + SPATIAL_RESOLUTION, SPATIAL_RESOLUTION)
    lonRange = np.append(lonRange, lon)
    lonRange.sort()

    # determine if lat and lon are along edge of grid box
    if lat%1 != 0 and lon%1 != 0:
        latPos = np.where(latRange == lat)[0][0] - 1
        lonPos = np.where(lonRange == lon)[0][0] - 1

    elif lat%1 != 0:
        latPos = np.where(latRange == lat)[0][0] - 1
        lonPos = np.where(lonRange == lon)[0][0]

    elif lon%1 != 0:
        latPos = np.where(latRange == lat)[0][0]
        lonPos = np.where(lonRange == lon)[0][0] - 1

    else:
        latPos = np.where(latRange == lat)[0][0]
        lonPos = np.where(lonRange == lon)[0][0]

    return latPos, lonPos

if __name__ == "__main__":
        
    # num rows and cols in NH data
    rowsNH = 896
    colsNH = 608

    jNH = np.arange(1,colsNH)
    iNH = np.arange(1,rowsNH)

    # NH coordinates
    NHlat = np.zeros((rowsNH, colsNH))
    NHlon = np.zeros((rowsNH, colsNH))

    # convert each i, j position in grid into lat, long
    for i in iNH:
        for j in jNH:

            # Call the fortran routine.
            # gtype (1 = 12.5km, 2 = 25km)
            # ihem (1 = NH, 2 = SH)
            # itrans(1 = i,j to lat, lon, 2 = lat,lon to i,j)
            lat, lon = ml.locate(gtype = 1, ihem = 1, itrans = 1, i = i, j = j)
            
            NHlat[i][j] = lat
            NHlon[i][j] = lon

    # num rows and cols in SH data
    rowsSH = 664
    colsSH = 632

    jSH = np.arange(1,colsSH)
    iSH = np.arange(1,rowsSH)

    # SH coordinates
    SHlat = np.zeros((rowsSH, colsSH))
    SHlon = np.zeros((rowsSH, colsSH))

    # convert each i, j position in grid into lat, long
    for i in iSH:
        for j in jSH:

            # Call the fortran routine.
            lat, lon = ml.locate(gtype = 1, ihem = 2, itrans = 1, i = i, j = j)
            
            SHlat[i][j] = lat
            SHlon[i][j] = lon    
        
    filename   = "data/AMSR_E_L3_SeaIce12km_V15_20070105.hdf"
    data_file  = SD(filename, SDC.READ)
    
    # https://nsidc.org/data/AE_SI12/versions/3
    
    NHconc = np.array(data_file.select("SI_12km_NH_ICECON_ASC")[:]).astype(np.float32)
    NHconc[(NHconc == 110) | (NHconc == 105) | (NHconc ==120) | (NHconc == 0)] = np.nan
    
    SHconc = np.array(data_file.select("SI_12km_SH_ICECON_ASC")[:]).astype(np.float32)
    SHconc[(SHconc == 110) | (SHconc == 105) | (SHconc == 120) | (SHconc == 0)] = np.nan
    
    gridXSize = int(180/SPATIAL_RESOLUTION)
    gridYSize = int(360/SPATIAL_RESOLUTION)
    gridZSize = 2
    grid = np.zeros((gridXSize, gridYSize, gridZSize))
    
    for i in range(len(NHlat)):
        for j in range(len(NHlat[i])):
            
            # uses lat and long values to find where on the 1x1 grid they find, returns index of grid
            latPos, lonPos = findGridPosition(NHlat[i][j], NHlon[i][j])
                            
            # add to sum of probablities of ice being present in the data grid
            grid[latPos][lonPos][0] += NHconc[i][j]

            # add one to the number of points inside the grid
            grid[latPos][lonPos][1] += 1


    for i in range(len(SHlat)):
        for j in range(len(SHlat[i])):
            
            # uses lat and long values to find where on the 1x1 grid they find, returns index of grid
            latPos, lonPos = findGridPosition(SHlat[i][j], SHlon[i][j])
                                
            grid[latPos][lonPos][0] += SHconc[i][j]
            grid[latPos][lonPos][1] += 1


    # array of avg probability of sea ice being in each grid box
    concProb = np.zeros((gridYSize, gridXSize))

    for i in range(len(grid)):
        for j in range(len(grid[i])):
                
            # avg probability = sum of probabilities / total data points in the grid box
            gridIceProb = grid[i][j][0]/grid[i][j][1]

            #if gridIceProb == 0:
            #    gridIceProb = np.nan

            concProb[j][i] = gridIceProb

    # create meshgrid
    lat = np.arange(-90,90 + SPATIAL_RESOLUTION, SPATIAL_RESOLUTION)
        
    lon = np.arange(0,360 + SPATIAL_RESOLUTION, SPATIAL_RESOLUTION)
    latv, longv = np.meshgrid(lat, lon)

    # plot
    cmap = plt.get_cmap('viridis')
    #cmap = plt.get_cmap('Blues_r')
    #cmap.set_bad('k', alpha = 1.)
    #levels = MaxNLocator(nbins=50).tick_values(0,np.nanmax(concProb))
    #norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True) 

    fig, ax = plt.subplots(figsize = (12,6))
    ax  = plt.axes(projection = ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()

    im = ax.pcolormesh(longv, latv,  concProb, cmap = cmap, transform =ccrs.PlateCarree())
    im.set_clim(0, 100)
    cbar = fig.colorbar(im, ax = ax)
    cbar.set_label('Average Concentration of Sea Ice')
    gridlines = ax.gridlines(draw_labels=True)
    
    ax.set_extent([145, 149, 51, 57])
    #ax.set_extent([-79.75, -79, 55, 55.75])
    ax.set_title("2007-01-05")


    # read in anomalies csv as pandas dataframe
    df = pd.read_csv("../../2007-01_water_anomalies.csv")
    
    # remove time
    df['timestamp'] = df['timestamp'].apply(editTimestamp)
    
    # go through each of the anomalies individually
    for index, row in df.iterrows():
        
        date = row['timestamp']
        anomalyLatitude = row['latitude']
        anomalyLongitude = row['longitude']

        if date == "20070105":

            # convert so lon goes from 0 to 360 instead of -180 to 180
            # then find corresponding pos in grid from lat, lon
            latPos, lonPos = findGridPosition(anomalyLatitude, anomalyLongitude + 180)

            if lonPos >= 180 / SPATIAL_RESOLUTION:
                lonPos -= int((180 / SPATIAL_RESOLUTION))
            else:
                lonPos += int((180 / SPATIAL_RESOLUTION))
    
            prob = concProb[lonPos][latPos]
            
            if math.isnan(prob):
                ax.scatter(anomalyLongitude + 180, anomalyLatitude, c = "b",edgecolors = "k", zorder = 3, marker = "x")
        
            else:
                ax.scatter(anomalyLongitude + 180, anomalyLatitude, c = "r",edgecolors = "k", zorder = 3, marker = "x")

    plt.savefig("20070105_zoom.png")
    plt.close()
