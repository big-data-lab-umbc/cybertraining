import numpy as np
import pandas
import sys
sys.path.append('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/code_Matt/notears-2.1/src')
import notearsORIGINAL, utils
from plot_graph import *
import matplotlib.pyplot as plt




file_names=['combined_decomposed_drop_temp_all_norm_1980_2018']

feature_names_dict={'combined_decomposed_drop_temp_all_norm_1980_2018': ['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']}


for file in file_names:
    
    df = pandas.read_csv('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/'+file+'.csv', delimiter=',',header=0)
    df.dropna(axis=0,inplace=True)
    X=df.to_numpy()
    feature_names=feature_names_dict[file]
    plotname="StaticModel__"+file

    #lambda1 used in example provided by authors is 0.1
    #default w_threshold=0.3

    for lambda1 in [0.1]:
        for w_threshold in [0.3]:
        
            W_est = notearsORIGINAL.notears_linear_l1(X.copy(), lambda1=lambda1, loss_type='l2',w_threshold=w_threshold)
            if not utils.is_dag(W_est):
                plotname_add='!!!NODAG!!!'
            else:
                plotname_add=''

            np.savetxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1)+'_Wthresh='+str(w_threshold)+plotname_add+'.csv', W_est, delimiter=',')
            visualize_static(np.around(W_est,2),feature_names=feature_names,plotname=plotname+'_lambda1='+str(lambda1)+'_Wthresh='+str(w_threshold)+plotname_add+'.pdf')


