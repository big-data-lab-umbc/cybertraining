# Standard library imports.
import datetime
import os
import re

# Third party imports.
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from tqdm import tqdm

# Local application imports.
import configuration
import file_loader
from collocate import lambert
from file_finder import get_filenames

# TODO: Parallelize with my imap function.

# Constant definitions.
DIRECTORIES = configuration.get_taki_directories()
START_DATE  = datetime.datetime(2007, 1, 1)
END_DATE    = datetime.datetime(2007, 12, 31)

# Constant definitions.
ANOMALY_PATH  = "2007_over-water_worldview_anomalies_with_sea_ice.csv"
CLOUD_TOP_KEY = "MYD06_Cloud_Top_Height_1km"
LATITUDE_KEY  = "CALIPSO_Lat1km"
LONGITUDE_KEY = "CALIPSO_Lon1km"
TIME_KEY      = "Profile_Time"
LOOKUP_RANGE  = 5


def truncate_datasets(datasets):
    
    for dataset in datasets:
        
        time = np.squeeze(dataset[TIME_KEY])
        time = np.array([datetime.timedelta(seconds = t) for t in time])
        time = datetime.datetime(1993, 1, 1) + time
        
        latitude         = dataset[LATITUDE_KEY]
        longitude        = dataset[LONGITUDE_KEY]
        cloud_top_height = dataset[CLOUD_TOP_KEY]
        
        try:
        
            truncated_dataset = pd.DataFrame({
                "time"             : time,
                "latitude"         : latitude,
                "longitude"        : longitude,
                "cloud_top_height" : cloud_top_height
            })
            
            yield truncated_dataset
        
        except ValueError:
        
            print("CALOPSO/COLLOCATION FILE LENGTH MISMATCH!")
            print(time.shape)
            print(latitude.shape)
            print(longitude.shape)
            print(cloud_top_height.shape)


if __name__ == "__main__":
    
    datasets = {
        "collocation_data" : [
            "Collocation_Flag",
            "CALIPSO_Lat1km",
            "CALIPSO_Lon1km",
            "MYD06_Cloud_Optical_Thickness",
            "MODIS_SZA1km",
            "MYD06_Cloud_Top_Height_1km"
        ],
        
        "collocation_indices" : [
        ],
        
        "CALIPSO_01km_data" : [
            "Opacity_Flag",
            "Feature_Classification_Flags",
            "Number_Layers_Found",
            "Solar_Zenith_Angle",
            "Profile_UTC_Time",
            "Profile_Time",
            "IGBP_Surface_Type",
            "Day_Night_Flag",
        ],
    }
    
    # Naive filename filtering.
    filenames = get_filenames(DIRECTORIES, START_DATE, END_DATE)
    datasets  = file_loader._preprocess_data(filenames, datasets)
    datasets  = truncate_datasets(datasets)
    
    # Loads anomalies.
    anomalies = pd.read_csv(ANOMALY_PATH)
    anomalies["timestamp"] = pd.to_datetime(anomalies["timestamp"])
    
    latitude_envelope_series   = pd.Series(np.empty(anomalies.size), dtype = object)
    longitude_envelope_series  = pd.Series(np.empty(anomalies.size), dtype = object)
    distance_series            = pd.Series(np.empty(anomalies.size), dtype = object)
    cloud_top_height_series    = pd.Series(np.empty(anomalies.size), dtype = object)
    cloud_top_slope_series     = pd.Series(np.empty(anomalies.size), dtype = np.float32)
    cloud_top_intercept_series = pd.Series(np.empty(anomalies.size), dtype = np.float32)
    
    latitude_envelope_series[:]   = np.nan
    longitude_envelope_series[:]  = np.nan
    distance_series[:]            = np.nan
    cloud_top_height_series[:]    = np.nan
    cloud_top_slope_series[:]     = np.nan
    cloud_top_intercept_series[:] = np.nan
    
    previous_dataset = pd.DataFrame(None)
    current_dataset  = next(datasets)
    next_dataset     = next(datasets)
    
    for i, anomaly in tqdm(anomalies.iterrows(), total = anomalies.index.size):
        
        anomaly_in_file = False
        
        while not anomaly_in_file:
            
            dataset = pd.concat([previous_dataset,
                                 current_dataset,
                                 next_dataset]).reset_index(drop = True)
            
            anomaly_time = anomaly.timestamp
            dataset_time = current_dataset.time
            
            anomaly_index   = np.where(dataset_time.isin([anomaly_time]))[0]
            anomaly_in_file = bool(anomaly_index)
                        
            if anomaly_in_file:
            
                # Aliasing for brevity.
                x  = previous_dataset.index.size + anomaly_index[0]
                dx = LOOKUP_RANGE

                surrounding_data = dataset.iloc[x - dx: x + dx + 1]
                
                latitudes  = surrounding_data.latitude
                longitudes = surrounding_data.longitude
                points     = list(zip(latitudes, longitudes))
                reference  = (latitudes.iloc[dx], longitudes.iloc[dx])
                
                distances = np.array(list(map(lambda p: lambert(reference, p), points)))
                distances[:LOOKUP_RANGE] *= -1
                cloud_top_height = surrounding_data["cloud_top_height"].values
                
                nan_mask = np.isnan(distances) | np.isnan(cloud_top_height)
                
                m, b = np.polyfit(distances[~nan_mask],
                                  cloud_top_height[~nan_mask],
                                  1)
                
                latitude_envelope_series[i]   = latitudes
                longitude_envelope_series[i]  = longitudes
                distance_series[i]            = distances
                cloud_top_height_series[i]    = cloud_top_height
                cloud_top_slope_series[i]     = m
                cloud_top_intercept_series[i] = b
                
            else:

                previous_dataset = current_dataset
                current_dataset  = next_dataset
                next_dataset     = next(datasets)
    
    anomalies["latitude_envelope"]   = latitude_envelope_series
    anomalies["longitude_envelope"]  = longitude_envelope_series
    anomalies["distance"]            = distance_series
    anomalies["cloud_top_height"]    = cloud_top_height_series
    anomalies["cloud_top_slope"]     = cloud_top_slope_series
    anomalies["cloud_top_intercept"] = cloud_top_intercept_series
    
    anomalies.to_csv("2007_over-water_worldview_water_cloud_anomalies_with_sea_ice_and_slopes.csv", index = False)
                