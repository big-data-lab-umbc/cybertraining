
#Read monthly means of each variable in combined daset
#Normalize data using max-min method
#Drop all missing values in CSV file
#The dataset produced is used for model input

import pandas as pd
import numpy as np
import sys

df0=pd.read_csv('combined_decomposed.csv',sep=',')
print(df0.columns)

###### Average GH values at different vertival levels
a1=df0['Residual_GH_200hPa']
a2=df0['Residual_GH_500hPa']
a3=df0['Residual_GH_850hPa']
df0['Residual_GH_mean']=(a1+a2+a3)/3

df1=df0.drop(columns=['Residual_temp','Residual_temp_2m','Residual_windspeed_10m','Residual_GH_200hPa','Residual_GH_850hPa','Residual_GH_500hPa'])
df11=df1.dropna() #drop missing values

###### normalize data 
cols=len(df11.columns)
for j in range (0,cols):
        df11.iloc[:,j]=(df11.iloc[:,j]-df11.iloc[:,j].min())/(df11.iloc[:,j].max()-df11.iloc[:,j].min())


###### Save data into csv file
df11.to_csv('combined_decomposed_drop_temp_all_norm_1980_2018.csv',header=True,index=False)

