
#Read monthly means of each variable
#Decompose time series data into T_trend, T_seasonal and Residual by applying additive model
#One-sided moving averages is applied, so the values for the first few months are NaN

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from statsmodels.tsa.seasonal import seasonal_decompose
from scipy import signal

nyear=2018-1979+1
nmons=nyear*12

filename='ERA5_Arctic_clouds_1979_2019.csv'

df0=pd.read_csv(filename,sep=',')
var0_raw=df0['Total cloud cover (%)']
var1_raw=df0['Total column cloud water (g m**-2)']
var0=var0_raw[0:nmons]
var1=var1_raw[0:nmons]
y0_raw=df0['Year']
m0_raw=df0['Month']
y0=y0_raw[0:nmons]
m0=m0_raw[0:nmons]

######Decompose times series data
# Time Series Decomposition
result_mul0 = seasonal_decompose(var0.values, model='additive', freq=12, two_sided=False)
deseason0 = var0.values/result_mul0.seasonal
detrend0 = var0.values-result_mul0.trend
re0=result_mul0.resid

result_mul1 = seasonal_decompose(var1.values, model='additive', freq=12, two_sided=False)
deseason1 = var1.values/result_mul1.seasonal
detrend1 = var1.values-result_mul1.trend
re1=result_mul1.resid

######Save data into csv file
dict={'Year':y0, 'Month':m0,'Original_cloud_cover(%)':var0, 'Detrend_cloud_cover':detrend0,'Deseason_cloud_cover':deseason0,'Residual_cloud_cover':re0,'Original_cloud_water(gm**-2)':var1, 'Detrend_cloud_water':detrend1,'Deseason_cloud_water':deseason1,'Residual_cloud_water':re1}
data=pd.DataFrame(dict)

data.to_csv('decomposed_ERA5_Arctic_clouds_1979_2018.csv',header=True,index=True)

