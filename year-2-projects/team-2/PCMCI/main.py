#! /umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/python3.6m
# Imports
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

## use `%matplotlib notebook` for interactive figures
# plt.style.use('ggplot')
import sklearn

import tigramite
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb

#Mine
import pdb
from load_data import load_data
import sys
import csv

"""
This file will run the PCMCI method on a datafile.

To specify the csv file to be used, just write it as the first command line.
"""

#Run settings
#there is another tau_max in lagged dependencies that might be much longer!
tau_max = 3

#Verbosity:
# 0 - nothing
# 1 - final graph only
# 2 - everything
verbose_max = 2
verbose = 2

display_images = False
save_images = True

#File Load settings
data_random = False
data_file = "new.csv"
write_data = False

#True if time_column must be removed, False if no time column
time_column = False

delimiter = ","
quotechar = "|"

#Check command line settings.
args = sys.argv
print(args)
if len(args) > 1:
    data_file = args[1]

#The following requires data_random, data_file, write_data
#It returns data, headers, 
data, headers = load_data(data_file, data_random, write_data, delimiter, quotechar, time_column)
    
#pdb.set_trace()
    
#Data is now full

#Write to file if needed
if write_data:
    with open("test.csv", "w", newline = "") as csvfile:
        data_writer = csv.writer(csvfile, delimiter = delimiter, quotechar=quotechar)
        for line in data:
            #data_writer.writerow("|".join([str(s) for s in line]))
            data_writer.writerow(line)
    
    exit()

T, N = data.shape

# Initialize dataframe object, specify time axis and variable names
#var_names = [r'$X^0$', r'$X^1$', r'$X^2$', r'$X^3$']
dataframe = pp.DataFrame(data, datatime = np.arange(len(data)), var_names=headers)

if verbose > 0:
    plot = tp.plot_timeseries(dataframe)[0]
    if display_images:
        plot.show()
    if save_images:
        plot.savefig("timeseries.png")

parcorr = ParCorr(significance='analytic')
pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr, verbosity=1)

correlations = pcmci.get_lagged_dependencies(tau_max=3)
lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations, setup_args={'var_names':headers, 'x_base':5, 'y_base':.5})

if verbose > 1:
    if display_images:
        lag_func_matrix.savefig()
    if save_images:
        lag_func_matrix.savefig("lag_func.png")

pcmci.verbosity = 1
results = pcmci.run_pcmci(tau_max=tau_max, pc_alpha=None)

#Print results 
print("p-values")
print (results['p_matrix'].round(3))
print("MCI partial correlations")
print (results['val_matrix'].round(2))

#Save results to file
p_matrix = results['p_matrix']
with open("p-values.csv", "w") as csv_file:
    writer = csv.writer(csv_file, delimiter = ",", quotechar ="|", quoting = csv.QUOTE_MINIMAL) 
    #[[[1 2 3]]] Three brackets to get through.
    for sector in p_matrix:
        print("sector: ", sector)
        for row in sector:
            print("row: ", row)
            writer.writerow(row)
        writer.writerow([])

val_matrix = results['val_matrix']
with open("val-values.csv", "w") as csv_file:
    writer = csv.writer(csv_file, delimiter = ",", quotechar ="|", quoting = csv.QUOTE_MINIMAL) 
    #[[[1 2 3]]] Three brackets to get through.
    for sector in val_matrix:
        print("sector: ", sector)
        for row in sector:
            print("row: ", row)
            writer.writerow(row)
        writer.writerow([])

q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'], fdr_method='fdr_bh')
pcmci.print_significant_links( p_matrix = results['p_matrix'], q_matrix = q_matrix, val_matrix = results['val_matrix'], alpha_level = 0.01)

link_matrix = pcmci.return_significant_parents(pq_matrix=q_matrix, val_matrix=results['val_matrix'], alpha_level=0.01)['link_matrix']

graph = tp.plot_graph( val_matrix=results['val_matrix'], link_matrix=link_matrix, var_names=headers, link_colorbar_label='cross-MCI', node_colorbar_label='auto-MCI',)

if verbose > 1:
    if display_images:
        graph[0].show()
    if save_images:
        graph[0].savefig("causal.png")

if verbose > 0:
    # Plot time series graph
    if display_images:
        plot = tp.plot_time_series_graph(
                val_matrix=results['val_matrix'],
                link_matrix=link_matrix,
                var_names=headers,
                link_colorbar_label='MCI')
    if save_images:
        plot = tp.plot_time_series_graph(
                val_matrix=results['val_matrix'],
                link_matrix=link_matrix,
                var_names=headers,
                link_colorbar_label='MCI',
                save_name = "causal_time_series.png")


