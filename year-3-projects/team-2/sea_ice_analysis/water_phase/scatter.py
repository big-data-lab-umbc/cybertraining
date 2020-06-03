import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs

df = pd.read_csv("2007_over-water_water_cloud_with_sea_ice.csv")

longitude = df["longitude"]
latitude  = df["latitude"]

ice = df["sea_ice_concentration"]

# this line cause a weird error(SettingWithCopyWarning)
#ice.loc[ice.isnull()] = -20


fig = plt.figure()
cmap = plt.get_cmap('jet')
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
im = ax.scatter(longitude, latitude, c = ice, cmap = cmap, vmin = 0, vmax = 100)
cbar = fig.colorbar(im, ax = ax, orientation="horizontal", pad=0.2)
cbar.set_label('Average Concentration of Sea Ice')
gridlines = ax.gridlines(draw_labels=True)
ax.set_extent([-180, 180, -90, 90])


plt.savefig("anomalies_locations+conc.png")
