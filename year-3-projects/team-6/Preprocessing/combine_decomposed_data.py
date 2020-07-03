
#Combined residual term of each variable
#Save it into a single CSV file

import numpy as np
import pandas as pd
import os

path = '/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/'
# list of all content in a directory, filtered so only directories are returned
dir_list = [directory for directory in os.listdir(path) if os.path.isdir(path+directory)]

#print(dir_list)

#for each_dir in dir_list:
 
file_path_list = [] 
file_name_list = [] 

from pathlib import Path

for i in dir_list:
    for file in Path(path+i).glob('decomposed*'):
        #print(file)
        #print(os.path.basename(file))
        file_path_list.append(str(file))
#        file_name_list.append(os.path.basename(file))
#file_name_list.remove('lagged*')

#print(file_path_list)


for file_path_1 in file_path_list:
    if 'lag' in file_path_1:
        file_path_list.remove(file_path_1)

#for file_path_2 in file_path_list:
#    if 'decompose' in file_path_2:
#        file_path_list.remove(file_path_2)


#print(file_path_list)
#exit()

#file_path_list.remove('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/temp_vertical/decomposed_ERA5_Arctic_temperature_1000-300hPa_1979_2018.csv')
#file_path_list.remove('/umbc/xfs1/cybertrn/cybertraining2020/team6/research/data/sea_ice/N_09_extent_v3.0.csv')

#print(file_path_list)
#exit()

#for file_name in file_name_list:
#    if file_name.startswith('lag'):
#        file_name_list.remove(file_name)

#for file_name in file_name_list:
#    if file_name.startswith('decomposed'):
#        file_name_list.remove(file_name)

#print(file_name_list)

#file_name_list.remove("N_09_extent_v3.0.csv")
#file_name_list.remove('decomposed_ERA5_Arctic_temperature_1000-300hPa_1979_2018.csv')
#print(file_name_list)

variable_list = []
 
#v1 = pd.read_csv(file_name_list[0]) 
#print(v1.head()) 

df_ = pd.DataFrame()

for v in file_path_list:
    df_v = pd.read_csv(v) 
    #df_ = df_.join(df_v)
    df_=pd.concat([df_,df_v],axis=1)

#print(df_.head())

df_ = df_.drop(df_.columns[0], axis=1)
df_ = df_.filter(like='Residual_',axis=1)
print(df_.head())


df_.to_csv("combined_decomposed.csv",index=False)
