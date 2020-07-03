
#Get lagged data for temporal model

import pandas as pd
import numpy as np
import sys

data_file_name = sys.argv[1]
maxlag = 12

data = pd.read_csv(data_file_name, header='infer')
#print(data.head())
#exit()

#data_noindex = data.drop(data.columns[0], axis=1)
df = data
#df = data_noindex.drop(['Year','Month'],axis=1)
print(df.head())

#df_list=[]

for x_name in list(df):
	for lag in range(1, maxlag + 1):
		df['{}-{}'.format(x_name, str(lag))] = df['{}'.format(x_name)].shift(lag)
		#df_list.append(df['{}-{}'.format(x_name, str(lag))])

print(df.head())
#print(df_list)


df.to_csv("lagged_"+data_file_name,index=False)

