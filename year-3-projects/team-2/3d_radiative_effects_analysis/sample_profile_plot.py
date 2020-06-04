from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime

from vertical_feature_mask import classify_vfm

# Constant definitions.
REFERENCE_DATE = datetime.datetime(1993, 1, 1)


def get_timestamp(seconds : int) -> datetime.datetime:
    
    '''
    
    Accepts as input the number of seconds since the reference date and returns
    the corresponding datetime.
    
    '''
    
    seconds      = np.array(seconds).astype(object)
    time_elapsed = np.vectorize(datetime.timedelta)(seconds = seconds)
    timestamp    = REFERENCE_DATE + time_elapsed
    
    return timestamp


if __name__ == "__main__":

    filename1 = "CAL_LID_L2_01kmCLay-ValStage1-V3-01.2007-01-01T01-09-14ZD_IND.hdf"
    filename2 = "CAL_LID_L2_01kmCLay-Standard-V4-10.2007-01-01T01-09-14ZD.hdf"
 
    # read in CALIPSO-MODIS collocated file
    f1 = SD(filename1, SDC.READ)

    # read in original CALIPSO file
    f2 = SD(filename2, SDC.READ)

    # Geolocation & Time
    lon = f2.select('Longitude').get()
    lat = f2.select('Latitude').get()
    time = f2.select('Profile_Time').get()
    clon = f1.select('CALIPSO_Lon1km').get()
    clat = f1.select('CALIPSO_Lat1km').get()

    time = get_timestamp(time)

    #time = get_timestamp(time)

    # MODIS Cloud Top Height
    MODIS_height = f1.select('MYD06_Cloud_Top_Height_1km').get()

    for key, value in f1.select('MYD06_Cloud_Top_Height_1km').attributes().items():
        if key == '_FillValue':
            fill_value = value
    
    # replace all the fill values with nan
    MODIS_height [np.where(MODIS_height == fill_value)] = np.nan
    MODIS_height[:] = MODIS_height[:]/1000

    # CALIPSO Cloud Top Height
    CAL_height = f2.select('Layer_Top_Altitude').get()

    # find the fill value being used
    for key, value in f2.select('Layer_Top_Altitude').attributes().items():
        if key == 'fillvalue':
            fill_value = value

    # replace all the fill values with nan
    CAL_height [np.where(CAL_height == fill_value)] = np.nan

    # heights are given as first value in each of the rows of cloud_top_alt
    CAL_height = [item[0] for item in CAL_height]


    feature_type = f2.select('Feature_Classification_Flags').get()

    # read in anomalies csv as pandas dataframe
    df = pd.read_csv("../2007-01_water_anomalies.csv")
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    
    #print(df.iloc[0])
    #print(time)
    
    feature_type = f2.select('Feature_Classification_Flags').get()

        
    loc =  np.where(lat == df['latitude'].iloc[0])[0][0]
    
    fig,ax = plt.subplots()
    CAL_height_zoom = CAL_height[loc-20:loc+21]
    MODIS_height_zoom = MODIS_height[loc-20:loc+21]
    lat_around_anomaly = np.around(lat[loc-20:loc+21],3)
    lon_around_anomaly = np.around(lon[loc-20:loc+21],3)
    print(len(lat_around_anomaly))

    for i in range(loc-20, loc+21):
        print(i)
        vfm = classify_vfm(feature_type[i])
        if vfm[0][0] != 2:
            CAL_height_zoom[i-(loc-20)] = np.nan


    labels = list(map(str, lat_around_anomaly[0]))

    ax.plot(MODIS_height_zoom, label = "MODIS")
    ax.plot(CAL_height_zoom, label = "CALIPSO")
    #ax.set_xticklabels(labels)
    ax.set_ylabel("Altitude (km)",fontsize=12)
    ax.axvline(x=20,c ='r', linestyle= '--', label = "Anomaly")
    ax.set_title("Jan 1, 2007 1:52:08 \n Latitude: " + str(lat_around_anomaly[0][0]) + " to " + str(lat_around_anomaly[-1][0]) + "\n Longitude: " + str(lon_around_anomaly[0][0]) + " to " + str(lon_around_anomaly[-1][0]))
    ax.set_xlim(0, 41)
    ax.legend()
    plt.savefig("cloud_top_alt_zoom.png")
    plt.close()

    """
    # go through each of the anomalies individually
    for index, row in df.iterrows():
    
        date = row['timestamp']
        anomalyLatitude = row['latitude']
        anomalyLongitude = row['longitude']
        
        # find index of where anomaly occurs using anomaly's time
        loc =  np.where(time == date)[0]
        
        fig,ax = plt.subplots()
        heights_around_anomaly = heights[loc-20:loc+21]
        lat_around_anomaly = lat[loc-20:loc+21]
        lon_around_anomaly = lon[loc-20:loc+21]
        
        ax.plot(heights_around_anomaly)
        ax.set_xticks(np.linspace(min(lat_around_anomaly),max(lat_around_anomaly),9))
        ax.set_ylabel("Altitude (km)",fontsize=12)
        plt.savefig("cloud_top_alt_zoom.png")
        plt.close()"""
    


#fig,ax = plt.subplots()
#ax.plot(heights)
#ax.set_xticks(np.linspace(min(lat),max(lat),9))
#ax.set_ylabel("Altitude (km)",fontsize=12)
#plt.savefig("cloud_top_alt.png")
#plt.close()


"""with file as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    data = list(f[a_group_key])"""
