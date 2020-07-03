from plot_graph import *
import numpy as np
import pandas as pd
import os


static_file = 'combined_decomposed_drop_temp_all_norm_1980_2018'
temporal_file = 'lagged_combined_decomposed_drop_temp_all_norm_1980_2018'
# epochs = 200 # not using anymore

lag_range=np.arange(1,13)
# threshold=[0.1,0.3,0.5,0.7]


short_feature_names=['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m',
           'sea_ice','CC','CW', 'GH']

df = pd.read_csv("../data/" + static_file + ".csv")
full_feature_names = list(df.columns)

df = pd.read_csv("../data/" + temporal_file + ".csv")
all_variable_names = list(df.columns)

static_root = os.path.join('../src/param-tuning2/')
folder_list = os.listdir(static_root)

for folder in folder_list:
    if ('.py' in folder) or ('__pycache__' in folder):
        continue
    os.chdir(os.path.join(static_root,folder))

    adja_mat = np.loadtxt('predG')
    thresh = float(folder[-3:])
    adja_mat[np.abs(adja_mat)<thresh] = 0
    if 'lagged' in folder:
        visualize_reduced_graph(np.around(adja_mat,2),short_feature_names,full_feature_names,lag_range,
                                    all_variable_names,"temporal_reduced.png")
    else:
        visualize_static(np.around(adja_mat,2),short_feature_names,"static.png")




