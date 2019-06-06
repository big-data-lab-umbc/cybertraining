#! /home/sthussun/intelpython3/bin/python3.6

#Standard Imports
import csv
import sys
import pdb

#Science Imports
import math as math
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#My imports
from graphic_supplemental import *

"""
This program will create a .png and .pdf file for inclusion in publications. It will take information about the global location of nodes and connections between several variables located at each node. (Or possibly one variable per node--doesn't matter.) The goal is to automatically visualize the causality graph given by any of several methods to discover global teleconnections among different climate phenomena.

Arguments: 0 or 2
    If 2, should be the node_location filename and then the connection location filename, e.g.,

    ```
    python generate_graph.py my_nodes.csv my_connections.csv
    ```


Input File Details: node_locations.cs connections.csv v

The connections.csv and node_locations.csv are provided as sample files to show the csv format. You can write directly into the csv format, or write to file with "," as a delimiter. (Make sure not to use "," in any cell files, or it will separate your data!)

The node_locations.csv file should be generated along with the timeseries files. The display name can be set manually, but should be short. To be honest, 3-5 letters would be optimal, but we could make do with SST-13, Si10-20, ENSO if need be.

The connections.csv file should only list the causal connections to plot, so you should sort and filter by whatever value your method uses. (If you list all the connections the plot will be unreadable.)

Make sure there are no extra spaces, and that the node names match exactly between the node_connections.csv and the connections.csv file, or the plot script will not be able to match the two.

Author: Steve Hussung
Date: May-June 2019
"""

# - - - Parameters - - - #

#Specify file defaults
node_location_file = "node_locations.csv"
conn_location_file = "connections.csv"

args = sys.argv
if len(args) == 3:
    node_location_file = args[1]
    conn_location_file = args[2]


#Boolean
show_to_user = False
save_to_file = True
color_map_background = False
straight_line = False

#Data Reading
#Make True to avoid graphing nodes not involved in connections, False to graph them.
remove_inactive_nodes = True

Colors = {"sst"  : (0.0, 0.0, 1.0),
          "t2m"  : (0.8, 0.0, 0.0),
          "msl"  : (0.9, 0.9, 0.0),
          "si1"  : (0.0, 0.8, 0.0)} #Keep it si1 to keep these three letters long


def getshort(node_name):
    for var in Colors.keys():
        if var in node_name:
            return var
    return "???"

Colors["???"] = (0.5, 0.5, 0.5)

#Plotting
text_offset = [2,2]
time_lag_offset = 4

if color_map_background: 
    coast_and_grid_color = 'grey'
else:
    coast_and_grid_color = (0.7, 0.7, 0.7)

#Projection Choice
projection = ccrs.Robinson()
data_trans = ccrs.PlateCarree()

#For arrows
#line_trans = data_trans
line_trans = ccrs.Geodetic()

#This is the fraction of the line in (0,1) that determines the arrow direction
arrow_fraction = 0.001

# - - - Read in Node Location and Connection Data - - - #
output_file_root = "causal_graph"
Nodes, Connections = read_data(node_location_file, conn_location_file)

#Title
plt.title("Robinson Projection")

# - - - Set up graph - - - #
#Establish axes
ax = plt.axes(projection=projection)

if color_map_background:
    ax.stock_img()

#Set to show globe
ax.set_global()

#Draw coasts and parallels
ax.coastlines(color=coast_and_grid_color)
ax.gridlines(color=coast_and_grid_color)

#Draw "satellite" image of earth.
if color_map_background:
    ax.stock_img()

# - - - Plot Connections - - - #

#Currently only plotting nodes involved in connections, so we iterate over connections
for connection in Connections:
    #Load node dictionaries
    node_a = Nodes[connection["cause"]]
    node_b = Nodes[connection["effect"]]

    cause_color =  Colors[getshort(connection["cause"])]
    effect_color = Colors[getshort(connection["effect"])]

    # - - - Plot Nodes - - - #
    arrows = []
    lines = []
    plt.text(node_a["longitude"]+text_offset[0], node_a["latitude"]+text_offset[1], node_a["display_name"], transform = ccrs.Geodetic(), zorder=10)
    plt.plot(node_a["longitude"], node_a["latitude"], color=cause_color, marker='o', transform = data_trans, zorder = 9)

    plt.text(node_b["longitude"]+text_offset[0], node_b["latitude"]+text_offset[1], node_b["display_name"], transform = ccrs.Geodetic(), zorder=10)
    plt.plot(node_b["longitude"], node_b["latitude"], color=effect_color, marker='o', transform = data_trans, zorder = 9)

    # - - - Plot lines and Arrows - - - #
    lon_a = node_a["longitude"]
    lon_b = node_b["longitude"]
    lat_a = node_a["latitude"]
    lat_b = node_b["latitude"]

    #Plot arrow (only draw along last fraction to get correct direction)
    alpha=arrow_fraction

    arrows.append(plt.arrow(alpha*lon_a + (1-alpha)*lon_b, 
                  alpha*lat_a + (1-alpha)*lat_b, 
                  alpha*(lon_b-lon_a), 
                  alpha*(lat_b-lat_a),
            linewidth=2, head_width=8, head_length=10, fc=effect_color, ec=effect_color,
            length_includes_head=True,
            zorder=5,
            transform=line_trans,
            ))
    arrows[-1].set_closed(False)

    #print("before transform", lon_a, lat_a)

    gradient_line([lon_a, lon_b], [lat_a, lat_b],
            [cause_color, effect_color], 
            line_trans,
            data_trans,
            projection,
            straight_line=straight_line,
            linewidth=2, marker='',
            zorder=4,
            )


    # - - - Plot timelag - - - #
    mid_point = np.array([0.5*(lon_a + lon_b), 0.5*(lat_a + lat_b)])
    line_angle = math.atan((lat_b - lat_a)/(lon_b - lon_a))
    offset_angle = line_angle + math.pi/2.0
    offset_mid_point = mid_point + \
        time_lag_offset*np.array([math.cos(offset_angle), math.sin(offset_angle)])

    plt.text(offset_mid_point[0], offset_mid_point[1], connection["timelag"], fontsize = 8, transform = line_trans, zorder=10)



# - - - Display and Save - - - #

#Save to file
if save_to_file:
    plt.savefig(output_file_root + '_png.png', bbox_inches="tight", dpi=300)
    plt.savefig(output_file_root + '_pdf.pdf', bbox_inches="tight")

#Show to user
if show_to_user:
    plt.show()

