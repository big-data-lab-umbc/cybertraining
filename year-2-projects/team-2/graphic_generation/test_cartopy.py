#! /home/sthussun/intelpython3/bin/python3.6

#Try matplotlib.patches.Polygon next!

#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import pdb

# - - - Parameters - - - #

#Boolean
show_to_user = True
save_to_file = True
color_map_background = False

#Plotting
text_offset = [2,2]
coast_and_grid_color = 'grey'


#Projection Choice
projection = ccrs.Robinson()
data_trans = ccrs.PlateCarree()

#For arrows
line_trans = data_trans
#line_trans = ccrs.Geodetic()

#This is the fraction of the line in (0,1) that determines the arrow direction
arrow_fraction = 0.001

#Title
plt.title("Robinson Projection")


# - - - Set up graph - - - #
#Establish axes
ax = plt.axes(projection=projection)
ax.stock_img()

#Set to show globe
ax.set_global()

#Draw coasts and parallels
ax.coastlines(color=coast_and_grid_color)
ax.gridlines(color=coast_and_grid_color)

#Draw "satellite" image of earth.
if color_map_background:
    ax.stock_img()

# - - - Plot Nodes - - - #

#Define nodes
fl_lon, fl_lat = -80.191, 25.7616,
ny_lon, ny_lat = -75, 43
delhi_lon, delhi_lat = 77.23, -80 

for lon, lat, name in [[fl_lon, fl_lat, "Miami"], [ny_lon, ny_lat, "NY"], [delhi_lon, delhi_lat, "Delhi"]]:
    plt.text(lon+text_offset[0], lat+text_offset[1], name, transform = ccrs.Geodetic())
    plt.plot(lon, lat, color="blue", marker='o', transform = data_trans)

# - - - Plot lines - - - #
#Plot line from ny to delhi
plt.plot([ny_lon, delhi_lon], [ny_lat, delhi_lat],
        color = "blue", linewidth=2, marker='o',
        #transform=ccrs.Geodetic(),
        transform=line_trans,
        )

#Plot arrow on top of line (only draw along last fraction to get correct direction)
alpha=arrow_fraction
a = plt.arrow(alpha*ny_lon + (1-alpha)*delhi_lon, alpha*ny_lat + (1-alpha)*delhi_lat, alpha*(delhi_lon-ny_lon), alpha*(delhi_lat-ny_lat),
        linewidth=2, head_width=6, head_length=8, fc='blue', ec='blue',
        length_includes_head=True,
        transform=line_trans,
        )
a.set_closed(False)

# - - - Display and Save - - - #

#Save to file
if save_to_file:
    plt.savefig('FloridaToBerlinPNG.png', bbox_inches="tight")
    plt.savefig('FloridaToBerlinPDF.pdf', bbox_inches="tight")

#Show to user
if show_to_user:
    plt.show()

