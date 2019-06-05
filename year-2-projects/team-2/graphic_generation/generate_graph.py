#! /home/sthussun/intelpython3/bin/python3.6

#Standard Imports
import csv
import sys
import pdb

#Science Imports
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#My imports
from read_data import read_data

"""
Todo:
    [ ] What if a node is really in two places. Can I duplicate the graph connections in some smart way?
    [ ] How to display timelag? (Should be easy!)
"""

"""
Program Description:

This program will create a .png and .pdf file for inclusion in publications. It will take information about the global location of nodes and connections between several variables located at each node. (Or possibly one variable per node--doesn't matter.) The goal is to automatically visualize the causality graph given by any of several methods to discover global teleconnections among different climate phenomena.

Input File Details: connections.csv node_locations.csv

The connections.csv and node_locations.csv are provided as sample files to show the csv format. You can write directly into the csv format, or write to file with "," as a delimiter. (Make sure not to use "," in any cell files, or it will separate your data!)

The node_locations.csv file should be generated along with the timeseries files. The display name can be set manually, but should be short. To be honest, 3-5 letters would be optimal, but we could make do with SST-13, Si10-20, ENSO if need be.

The connections.csv file should only list the causal connections to plot, so you should sort and filter by whatever value your method uses. (If you list all the connections the plot will be unreadable.)

Make sure there are no extra spaces, and that the node names match exactly between the node_connections.csv and the connections.csv file, or the plot script will not be able to match the two.

Author: Steve Hussung
Date: May-June 2019
"""

# - - - Parameters - - - #

#Boolean
show_to_user = True
save_to_file = True
color_map_background = False

#Data Reading
#Make True to avoid graphing nodes not involved in connections, False to graph them.
remove_inactive_nodes = True

Colors = {"sst"  : "blue",
          "t2m"  : "red",
          "msl"  : "yellow",
          "si1" : "green"} #Keep it si1 to keep these three letters long

def getshort(node_name):
    return node_name[7:10]

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

# - - - Read in Node Location and Connection Data - - - #
output_file_root = "causal_graph"
Nodes, Connections = read_data("node_locations.csv", "connections.csv")

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

#Currently only plotting nodes involved in connections, so we iterate over connections
for connection in Connections:
    #Load node dictionaries
    node_a = Nodes[connection["cause"]]
    node_b = Nodes[connection["effect"]]

    # - - - Plot Nodes - - - #
    arrows = []
    for node in [node_a, node_b]:
        plt.text(node["longitude"]+text_offset[0], node["latitude"]+text_offset[1], node["display_name"], transform = ccrs.Geodetic(), zorder=10)
        plt.plot(node["longitude"], node["latitude"], color="blue", marker='o', transform = data_trans, zorder = 9)

    # - - - Plot lines and Arrows - - - #
    lon_a = node_a["longitude"]
    lon_b = node_b["longitude"]
    lat_a = node_a["latitude"]
    lat_b = node_b["latitude"]

    #Plot arrow (only draw along last fraction to get correct direction)
    alpha=arrow_fraction

    cause_color =  Colors[getshort(connection["cause"])]
    effect_color = Colors[getshort(connection["effect"])]

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

    #Plot line from node_a to node_b on top of arrow 
    #(well... not quite on top. change color to see)
    plt.plot([lon_a, lon_b], [lat_a, lat_b],
            color = cause_color, linewidth=2, marker='o',
            zorder=4,
            transform=line_trans,
            )

    # - - - Plot timelag - - - #
    plt.text(0.5*(lon_a + lon_b), 0.5*(lat_a + lat_b) + 4, connection["timelag"], fontsize = 8, transform = ccrs.Geodetic(), zorder=10)



# - - - Display and Save - - - #

#Save to file
if save_to_file:
    plt.savefig(output_file_root + '_png.png', bbox_inches="tight")
    plt.savefig(output_file_root + '_pdf.pdf', bbox_inches="tight")

#Show to user
if show_to_user:
    plt.show()

