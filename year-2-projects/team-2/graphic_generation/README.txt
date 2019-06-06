| - - - - - - - - - - - - - - - - - - - |
| - - - READ THIS IF NOTHING ELSE - - - |
| - - - - - - - - - - - - - - - - - - - |

The main code is graphic_generation.py

To run it on your_node.csv and your_connections.csv, run
python graphic_generation.py your_node.csv your_connections.csv

It will output a png and a pdf file named causal_graph_png.png and causal_graph_pdf.pdf.

Let me know if you have questions / things that should change!

-Steve

| - - - - - - - - - - |
| - - - DETAILS - - - |
| - - - - - - - - - - |

the makefile is provided to automatically generate the sample.pdf file, which can be displayed using evince by "make display". 

The connections.csv and node_locations.csv are provided as sample files to show the csv format. You can write directly into the csv format, or write to file with "," as a delimiter. (Make sure not to use "," in any cell files, or it will separate your data!)

The node_locations.csv file should be generated along with the timeseries files. The display name can be set manually, but should be short. To be honest, 3-5 letters would be optimal, but we could make do with SST-13, Si10-20, ENSO if need be.

The connections.csv file should only list the causal connections to plot, so you should sort and filter by whatever value your method uses. (If you list all the connections the plot will be unreadable.)

Make sure there are no extra spaces, and that the node names match exactly between the node_connections.csv and the connections.csv file, or the plot script will not be able to match the two.
