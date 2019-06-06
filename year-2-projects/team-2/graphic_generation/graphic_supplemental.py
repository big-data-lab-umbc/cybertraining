import csv
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def read_data(node_filename, connection_filename):
    #Nodes - dictionary
    #The nodes are stored in node_locations.csv, and will be stored in a dictionary

    Nodes = {}
    print("Loading nodes from", node_filename)
    with open(node_filename,"r") as node_file:
        csv_reader = csv.reader(node_file)
        csv_reader.__next__()
        for row in csv_reader:
            if row != []:
                #Add new dictionary to Nodes
                Nodes[row[1]] = {}

                Nodes[row[1]]["display_name"] = row[0]
                Nodes[row[1]]["latitude"]  = float(row[2])
                Nodes[row[1]]["longitude"] = float(row[3])

    #Connections - list
    Connections = []
    print("Loading connections from", connection_filename)
    with open(connection_filename,"r") as node_file:
        csv_reader = csv.reader(node_file)
        csv_reader.__next__()
        for row in csv_reader:
            if row != []:
                #Read in Connections
                Connections.append({})
                #Need nodeA, nodeB, variable(1, 2, 3, 4)
                Connections[-1]["effect"]    = row[0]
                Connections[-1]["cause"]     = row[1]
                Connections[-1]["timelag"]   = row[2]
                Connections[-1]["strength"]  = row[3]
                Connections[-1]["method"]    = row[4]
     
    return Nodes, Connections


def gradient_line(x_vals, y_vals, color_vals, transform, data_trans, projection, straight_line=False, linewidth=2, marker='', zorder=4, divisions=100):
    """
    This method will return a *bunch* of plt.plot objects. Later we might combine into a line collection, but for now separate will do.

    Make sure to use extend when incorporating the return
    """

    #If doing a straight line, we transform *now*, and not later.
    if straight_line:
        x_vals = projection.transform_point(*x_vals, data_trans)
        y_vals = projection.transform_point(*y_vals, data_trans)
    
    #We interpolate each element.
    fractions = [float(i)/divisions for i in range(divisions+1)]
    X = [x_vals[0] + (x_vals[1] - x_vals[0])*alpha for alpha in fractions]
    Y = [y_vals[0] + (y_vals[1] - y_vals[0])*alpha for alpha in fractions]

    #First convert colors to np arrays, for doing arithmetic
    color_vals = list(map(lambda A : np.array(A), color_vals))

    #Interpolate, then recast to tuples again.
    Colors = [color_vals[0] + (color_vals[1] - color_vals[0])*alpha for alpha in fractions]
    Colors = list(map(lambda A : tuple(A), Colors))

    Lines = []
    for i in range(divisions):
        color = Colors[i]


        lon_a = X[i]
        lat_a = Y[i]
        lon_b = X[i+1]
        lat_b = Y[i+1]

        if straight_line:
            Lines.append(
                    plt.plot([lon_a, lon_b], [lat_a, lat_b],
                    color=color, 
                    linewidth=linewidth, marker=marker,
                    zorder=zorder,
                    ))
        else:
            Lines.append(
                    plt.plot([lon_a, lon_b], [lat_a, lat_b],
                    color=color, 
                    transform=transform,
                    linewidth=linewidth, marker=marker,
                    zorder=zorder,
                    ))

    return Lines
