import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from pyhdf.SD import SD, SDC
from matplotlib import colors as CS
from mpl_toolkits.axes_grid1 import make_axes_locatable

#BT_data_for_plot_dust_nighttime_v7.dat
#../screen_v4_result/BT_data_for_plot_dust_fenn_v6.dat
#--------------------Read IIR 3 Bands & CALIOP dust layers----------------------
#data4 = np.array(pd.read_csv('nigtime_order.txt', header=None, delim_whitespace=True))
#order = data4[:,0]

data2 = np.array(pd.read_csv('nigtime_simulaton.txt', header=None, delim_whitespace=True))
LBL_BT08  = data2[:,0]
LBL_BT10  = data2[:,1]
LBL_BT12  = data2[:,2]

#sort_idx=np.argsort(order)
##print(order,sort_idx,order[sort_idx])
#LBL_BT08=LBL_BT08[sort_idx]
#LBL_BT10=LBL_BT10[sort_idx]
#LBL_BT12=LBL_BT12[sort_idx]
#sys.exit()

data3 = np.array(pd.read_csv('BT_data_for_plot_dust_fenn_v7.dat', header=0, delim_whitespace=True))
lat       = data3[:,0]
lon       = data3[:,1]
IIR_BT08  = data3[:,2] 
IIR_BT10  = data3[:,3] 
IIR_BT12  = data3[:,4]
FAS_BT08  = data3[:,5]
FAS_BT10  = data3[:,6]
FAS_BT12  = data3[:,7]
aod_532   = data3[:,8]
FAS_AOD_fen = data3[:,9]
dust_flag = data3[:,10]
clen_flag = data3[:,11]
clou_flag = data3[:,12]
srta_flag = data3[:,13]

data1 = np.array(pd.read_csv('BT_data_for_plot_dust_nighttime_v7.dat', header=0, delim_whitespace=True))
surf_temp = data1[:,10]

#fout = open('surface_temp_nighttime.txt','w')
#for i in range(len(surf_temp)):
#    fout.write ('{0:<10.3f}\n'\
#                .format(surf_temp[i]))
#fout.close()
#sys.exit()

dbt_fas_08 = IIR_BT08 - FAS_BT08
dbt_fas_10 = IIR_BT10 - FAS_BT10
dbt_fas_12 = IIR_BT12 - FAS_BT12

dbt_0810 = IIR_BT08 - IIR_BT10
dbt_1012 = IIR_BT10 - IIR_BT12

aod_532 [np.where(dust_flag == 0)]=np.nan
FAS_AOD_fen[np.where(dust_flag == 0)]=np.nan

#Average AOD of FASDOM from 1km to 5km
uni_idx = np.unique(aod_532)
for i in uni_idx:
    address_idx = [x for x in range(len(aod_532)) if aod_532[x] == i]
    FAS_AOD_fen[address_idx] = np.mean(FAS_AOD_fen[address_idx])

# Specify Region
NTA_lats = [3,16]    
NTA_lons = [-180,180]

#--------------------Read CALIOP Attenuate Profile------------------------
f2 = SD("../../dust_case_night/CAL_LID_L1-Standard-V4-10.2008-07-18T03-30-04ZN.hdf", SDC.READ)

# Geolocation & Time
CALIOP_lon = f2.select('Longitude').get()
CALIOP_lat = f2.select('Latitude').get()

# Parameters for comparison
scat_prof_1064 = f2.select('Total_Attenuated_Backscatter_532').get()

# Use _FillValue to remove fill data
for key, value in f2.select('Total_Attenuated_Backscatter_532').attributes().items():
    if key == 'fillvalue':
        fill_value = value

scat_prof_1064 [np.where(scat_prof_1064 == fill_value)] = np.nan

# Refine CALIOP lat & lon at the case region
condition2 = np.where((CALIOP_lon >= NTA_lons[0]) & (CALIOP_lon <= NTA_lons[1]) & \
                      (CALIOP_lat >= NTA_lats[0]) & (CALIOP_lat <= NTA_lats[1]))

scat_prof_1064 = scat_prof_1064[condition2[0],:]
CALIOP_lon = CALIOP_lon[condition2]
CALIOP_lat = CALIOP_lat[condition2]

#--------------------Read MODIS IR Bands------------------------
f1 = SD("../2008-07-18T03-30-04ZN.hdf", SDC.READ)

# Geolocation & Time
MODIS_lon = f1.select('CALIPSO_Longitude_5km').get()
MODIS_lat = f1.select('CALIPSO_Latitude_5km').get()
cross_idx_ori = f1.select('MODIS_CrossTrack_Index').get()
along_idx_ori = f1.select('MODIS_AlongTrack_Index').get()

# Parameters for comparison
MODIS_BT   = f1.select('MYD02_1km_BrightnessTemperature').get()

# Use _FillValue to remove fill data
for key, value in f1.select('MYD02_1km_BrightnessTemperature').attributes().items():
    if key == '_FillValue':
        fill_value = value

MODIS_BT [np.where(MODIS_BT == fill_value)] = np.nan

BT_0375 = MODIS_BT[:,0]
BT_0855 = MODIS_BT[:,8]
BT_1102 = MODIS_BT[:,10]
BT_1203 = MODIS_BT[:,11]

condition1 = np.where((MODIS_lon[:,1] >= NTA_lons[0]) & (MODIS_lon[:,1] <= NTA_lons[1]) & \
                      (MODIS_lat[:,1] >= NTA_lats[0]) & (MODIS_lat[:,1] <= NTA_lats[1]))
MODIS_lon = MODIS_lon[condition1[0],1]
MODIS_lat = MODIS_lat[condition1[0],1]
cross_idx_ori = cross_idx_ori[condition1[0]]
along_idx_ori = along_idx_ori[condition1[0]]

BT_0375 = BT_0375[condition1[0]]
BT_0855 = BT_0855[condition1[0]]
BT_1102 = BT_1102[condition1[0]]
BT_1203 = BT_1203[condition1[0]]

dbt_MODIS_1 = BT_0375-BT_1102
dbt_MODIS_2 = BT_1102-BT_1203

jud_lat = np.round(lat[np.where((dust_flag != 1))],2) #& (clen_flag == 0))],2)
jud_lon = np.round(lon[np.where((dust_flag != 1))],2) #& (clen_flag == 0))],2)

# Matching MODIS 5km index to IIR profiles
cross_idx = np.zeros(lat.shape[0]) + np.nan
along_idx = np.zeros(lat.shape[0]) + np.nan
for i in range(MODIS_lat.shape[0]):    
    L2_idx = np.where(np.abs(np.round(lat,3) - np.round(MODIS_lat[i],3)) < 0.0001)
    if len(L2_idx[0]) == 1: 
        L2_idx = L2_idx[0][0]
        break
for j in range(i,MODIS_lat.shape[0]):
    cross_idx[L2_idx:L2_idx+5] = cross_idx_ori[j]
    along_idx[L2_idx:L2_idx+5] = along_idx_ori[j]
    L2_idx += 5

for i in range(len(jud_lat)):
    condition3 = np.where(np.round(MODIS_lat,2) == jud_lat[i])
    #print(condition3,MODIS_BT.shape)
    dbt_MODIS_1[condition3] = np.nan
    dbt_MODIS_2[condition3] = np.nan

plt.figure(figsize=(14,16)) 

#--------------------------First Plot: Dust attenuate backscatter profile---------------
# Initialize bin array    
binalt = np.zeros(583)
# 30.1km to 40.0km
for i in list (range (33)):
    ind=32-i
    binalt[ind] = 30.1 + (i+1)*0.3
# 20.2km to 30.1km
for i in list (range (55)):
    ind=87-i
    binalt[ind] = 20.2 + (i+1)*0.18
# 8.3km to 20.2km
for i in list (range (200)):
    ind=287-i
    binalt[ind] = 8.3 + (i+1)*0.06
# -0.5 to 8.3km
for i in list (range (290)):
    ind=577-i
    binalt[ind] = -0.5 + (i+1)*0.03
# -2.0 to -0.5km
for i in list (range (5)):
    ind=582-i
    binalt[ind] = -2.0 + (i+1)*0.3

levels  = 21
X = CALIOP_lat
Y = binalt[287:577]
Z = scat_prof_1064[:,287:577].T
Z[np.where((Z > 0.01) | (Z < 0.0001))] = np.nan
#print(X.shape,Y.shape,Z.shape)

# Define a discrete colormap
colors1 = plt.cm.jet(np.linspace(0., 1, 200))
colors2 = plt.cm.Greys_r(np.linspace(0.5, 1, 56))
colors = np.vstack((colors1, colors2))

# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', colors)#cmaplist, cmap.N)

# define the bins and normalize
bounds1 = np.logspace(-4,-1, 31)
bounds = np.linspace(1e-3,0.1,levels)
norm = mpl.colors.LogNorm(1e-3,0.1, cmap.N)

plt.subplot(411)
cset = plt.pcolor(X, Y, Z, cmap=cmap, norm=norm, vmin=0.001, vmax=0.1)

# Adjust the panel's size to agree with those without colarbar.
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.5)

m = plt.cm.ScalarMappable(cmap=cmap)
m.set_array(Z)
m.set_clim(0.001, 0.1)
cg = plt.colorbar(m, cmap=cmap,norm=norm,spacing='proportional',ticks=bounds,boundaries=bounds,cax=cax,format='%4.1e')
cg.set_label('532nm Total Attenuated Backscatter (km^-1 sr^-1)')
ax.set_xticks(np.linspace(min(lat),max(lat),9))
ax.set_ylabel("Altitude (km)",fontsize=12)
ax.set_title("Dust Case 2: 2008-07-18 03:30:04Z Nighttime",fontsize=18)

#--------------------------Second Plot: IIR BT & dBT of IIR-FASDOM---------------
ax1 = plt.subplot(412)
#IIR_BT08[np.where(dust_flag != 1)]=np.nan
#ax1.plot(lat[np.where(dust_flag == 1)],IIR_BT08[np.where(dust_flag == 1)],'-',color='C1',linewidth=1.2)
ax1.plot(lat,IIR_BT08,'-',color='C1',linewidth=1.5)
ax1.plot(lat,FAS_BT08,'-',color='C0',linewidth=1.5)
ax1.plot(lat,LBL_BT08,'-',color='C2',linewidth=1.5)
ax1.plot(MODIS_lat,BT_0855,'--',color='C5',linewidth=2.0)
ax1.plot(lat,surf_temp,'-',color='red',linewidth=1.5)
ax1.legend(("IIR BT (8.65 um)","FAS BT (8.65 um)","LBL BT (8.65 um)","MODIS BT (8.55 um)","Surf_Tmp"),bbox_to_anchor=(1.01,0.25), loc="lower left")

ax1.margins(0)
plt.xticks(np.linspace(min(lat),max(lat),9))
plt.yticks(color='C1')
ax1.set_ylim(275,303)
ax1.set_ylabel('IIR BT at 08.65um (K)',fontsize=12,color='C1')

# Adjust the panel's size to agree with those without colarbar.
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="3%", pad=0.5)
cax1.axis('off')

#--------------------------Third Plot: AOD at 532nm from CALIOP,MODIS,FASDOM---------------
ax1 = plt.subplot(413)
#IIR_BT10[np.where(dust_flag != 1)]=np.nan
ax1.plot(lat,IIR_BT10,'-',color='C1',linewidth=1.5)
ax1.plot(lat,FAS_BT10,'-',color='C0',linewidth=1.5)
ax1.plot(lat,LBL_BT10,'-',color='C2',linewidth=1.5)
ax1.plot(MODIS_lat,BT_1102,'--',color='C5',linewidth=2.0)
ax1.plot(lat,surf_temp,'-',color='red',linewidth=1.5)
ax1.legend(("IIR BT (10.60 um)","FAS BT (10.60 um)","LBL_BT (10.60 um)","MODIS BT (11.02 um)","Surf_Tmp"),bbox_to_anchor=(1.01,0.25), loc="lower left")

ax1.margins(0)
plt.xticks(np.linspace(min(lat),max(lat),9))
plt.yticks(color='C1')
ax1.set_ylim(275,303)
ax1.set_ylabel('IIR BT at 10.60um (K)',fontsize=12,color='C1')

# Adjust the panel's size to agree with those without colarbar.
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="3%", pad=0.5)
cax1.axis('off')

#--------------------------Fourth Plot: Simulated BT of FASDOM & dBT (IIR-FASDOM)---------------
ax1 = plt.subplot(414)
#IIR_BT12[np.where(dust_flag != 1)]=np.nan
ax1.plot(lat,IIR_BT12,'-',color='C1',linewidth=1.5)
ax1.plot(lat,FAS_BT12,'-',color='C0',linewidth=1.2)
ax1.plot(lat,LBL_BT12,'-',color='C2',linewidth=1.5)
ax1.plot(MODIS_lat,BT_1203,'--',color='C5',linewidth=2.0)
ax1.plot(lat,surf_temp,'-',color='red',linewidth=1.5)
ax1.legend(("IIR BT (12.05 um)","FAS BT (12.05 um)","LBL_BT (12.05 um)","MODIS BT (12.03 um)","Surf_Tmp"),bbox_to_anchor=(1.01,0.25), loc="lower left")

ax1.margins(0)
plt.xticks(np.linspace(min(lat),max(lat),9))
plt.yticks(color='C1')
ax1.set_ylim(275,303)
ax1.set_ylabel('IIR BT at 12.05um (K)',fontsize=12,color='C1')
ax1.set_xlabel('Latitude (degree)',fontsize=12)

# Adjust the panel's size to agree with those without colarbar.
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="3%", pad=0.5)
cax1.axis('off')

plt.savefig("nighttime_LBL_MODIS_BT.png",dpi=600)
#plt.show()
#plt.close()