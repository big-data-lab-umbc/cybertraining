import numpy as np
import pandas
import sys
sys.path.append('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/code_Matt/notears-2.1/src')
import notearsORIGINAL, utils
from plot_graph import *
import matplotlib.pyplot as plt



file_names=['lagged_combined_decomposed_drop_temp_all_norm_1980_2018']

feature_names_dict={'lagged_combined_decomposed_drop_temp_all_norm_1980_2018': ['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']}

full_feature_names_dict={'lagged_combined_decomposed_drop_temp_all_norm_1980_2018': ['Residual_heat_flux','Residual_shortwave','Residual_longwave','Residual_SLP',
                                                                                  'Residual_tot_precip','Residual_RH','Residual_u10m','Residual_v10m','Residual_sea_ice',
                                                                                   'Residual_cloud_cover','Residual_cloud_water','Residual_GH_mean']}



for file in file_names:

    df = pandas.read_csv('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/lagged_data_May_12/'+file+'.csv', delimiter=',',header=0)              
    df.dropna(axis=0,inplace=True)                                                                                                                                            
    X=df.to_numpy()                                                                                                                                                    
    feature_names=feature_names_dict[file] 
    full_feature_names=full_feature_names_dict[file]                                                     
    lag_range=range(1,12+1)                                                                                                                                            
    all_variable_names=list(df.columns)                                                                                                                                
    plotname="TemporalModel__"+file                                                                                                                        
    plotname_REDUCED="TemporalModel*REDUCED*__"+file


    #lambda1 used in example provided by authors is 0.1
    #default w_threshold=0.3
  
    for lambda1 in [0.1]:
        for w_threshold in [0.3]:
            
            W_est = notearsORIGINAL.notears_linear_l1(X.copy(), lambda1=lambda1, loss_type='l2',w_threshold=w_threshold)
            if not  utils.is_dag(W_est):
                plotname_add='!!!NODAG!!!'
            else:
                plotname_add=''
            
            np.savetxt('AdjMatrix_'+plotname+'_lambda1='+str(lambda1)+'_Wthresh='+str(w_threshold)+plotname_add+'.csv', W_est, delimiter=',')
            
            
            visualize_temporal(np.around(W_est,2),feature_names=feature_names,full_feature_names=full_feature_names,lag_range=lag_range,
                       all_variable_names=all_variable_names,plotname=plotname+'_lambda1='+str(lambda1)+'_Wthresh='+str(w_threshold)+plotname_add+'.pdf')
            visualize_reduced_graph(np.around(W_est,2),feature_names=feature_names,full_feature_names=full_feature_names,lag_range=lag_range,
                            all_variable_names=all_variable_names,plotname=plotname_REDUCED+'_lambda1='+str(lambda1)+'_Wthresh='+str(w_threshold)+plotname_add+'.pdf')

