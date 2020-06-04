import os
from operator import add

import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER
from cartopy.mpl.gridliner import LATITUDE_FORMATTER
from matplotlib import gridspec
from pyspark.sql import SparkSession
from tqdm import tqdm

df = pd.read_csv("2007-01_water_worldview_anomalies.csv")
df2 = pd.read_csv("2007-01_water_anomalies.csv")
df["t"] = pd.to_datetime(df2["timestamp"]).dt.strftime("%Y-%m-%d %H:%M:%S UTC")


def plot(row):

    index, row = row
    
    image = plt.imread(os.path.join("Images", row["image_filename"]))

    gridlines_kwargs = {
        "crs"         : ccrs.PlateCarree(),
        "draw_labels" : True,
        "linewidth"   : 0.5, 
        "color"       : "black",
        "linestyle"   : "--"
    }

    spec = gridspec.GridSpec(ncols=2, nrows=1,
                             width_ratios=[2, 1])

    fig = plt.figure(figsize = (10, 4))

    plt.subplots_adjust(wspace = 0.5)

    ax1 = fig.add_subplot(spec[0], projection = ccrs.PlateCarree())
    ax1.set_title(f"Anomaly at ({row['latitude']:.6f}, {row['longitude']:.6f})")

    ax1.stock_img()
    ax1.add_feature(cfeature.LAND)
    ax1.add_feature(cfeature.COASTLINE)

    #ax1.coastlines()
    ax1.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    ax1.scatter(row["longitude"], row["latitude"], zorder = 3, edgecolors = "k", color = "r", transform = ccrs.PlateCarree())

    gl = ax1.gridlines(**gridlines_kwargs)
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xformatter    = LONGITUDE_FORMATTER
    gl.yformatter    = LATITUDE_FORMATTER

    ax2 = fig.add_subplot(spec[1])
    ax2.set_title("NASA Worldview Snapshot")
    ax2.imshow(image)
    ax2.set_xticks([0, image.shape[0]])
    ax2.set_yticks([0, image.shape[1]])
    ax2.set_xticklabels([f"${x:.6f}^{{\circ}}$" for x in [row["minimum_bounding_longitude"], row["maximum_bounding_longitude"]]], rotation = 45)
    ax2.set_yticklabels([f"${x:.6f}^{{\circ}}$" for x in [row["maximum_bounding_latitude"], row["minimum_bounding_latitude"]]])

    plt.suptitle(f"{row['t']}", y = 0.9, fontweight = "bold")
    
    filename = f"Anomaly {row.name + 1:04d} ({row['latitude']:.6f}, {row['longitude']:.6f}, {row['timestamp'].replace(':', '-')}).png"
    
    plt.savefig(os.path.join("Figures", filename))
    plt.close()
    
    return 1
    

#for index, row in tqdm(df.iterrows(), total = df.index.size):
#        
#    plot(row)
    
if __name__ == "__main__":
        
    spark   = SparkSession.builder.appName("sandbox2").getOrCreate()
    sc      = spark.sparkContext
    rdd     = sc.parallelize(df.iterrows(), 100)
    results = rdd.map(plot).reduce(add)
    
    print("Mapping complete.")

    #spark.stop()
