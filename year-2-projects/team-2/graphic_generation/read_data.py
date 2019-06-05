import csv

def read_data(node_filename, connection_filename):
    #Nodes - dictionary
    #The nodes are stored in node_locations.csv, and will be stored in a dictionary
    Nodes = {}
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
