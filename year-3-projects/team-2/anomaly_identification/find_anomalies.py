# Standard library imports.
import datetime
import platform
import re
from multiprocessing import Pool

# Third party imports.
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from tqdm import tqdm

# Local application imports.
import configuration
import file_loader
from file_finder import get_filenames

# TODO: Parallelize with my imap function.

# Constant definitions.
if re.match(r"taki-usr\d+", platform.node()):
    
    DIRECTORIES = configuration.get_taki_directories()
    START_DATE  = datetime.datetime(2007, 1, 1)
    END_DATE    = datetime.datetime(2007, 12, 31)
    PROCESSES   = 4
    
else:
    
    DIRECTORIES = configuration.get_local_directories()
    START_DATE  = datetime.datetime(2007, 1, 1)
    END_DATE    = datetime.datetime(2007, 12, 31)
    PROCESSES   = 2

    
if __name__ == "__main__":
    
    # Defines the datasets to extract from each respective file.
    #
    # Note:
    #
    #   Cloud layer product "Number_Layers_Found" flag only reports the number
    #   of cloud layers detected, and vice versa for aerosol layers.
    #
    # Note:
    #
    #   MODIS has a maximum cloud optical depth of 150; this causes a histogram
    #   of COD values to have a stacked tail at the bin corresponding to a COD
    #   of 150.
    #
    #       Source:
    #
    #           MODIS Cloud Optical Properties
    #
    #           User Guide for the Collection 6/6.1 Level-2 MOD06/MYD06 Product
    #           and Associated Level-3 Datasets
    #
    #           2.10.10. Maximum retrievable cloud optical thickness extended
    #           to 150
    #
    #           Page 75
    
    datasets = {
        "collocation_data" : [
            "Collocation_Flag",
            "CALIPSO_Lat1km",
            "CALIPSO_Lon1km",
            "MYD06_Cloud_Optical_Thickness",
            "MODIS_SZA1km",
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
        
#        "CALIPSO_05km_data" : [
#            "Opacity_Flag",
#            "Day_Night_Flag",
#            "Feature_Classification_Flags",
#            "Profile_UTC_Time",
#            "Solar_Zenith_Angle",
#            "Number_Layers_Found",
#            "Column_Optical_Depth_Cloud_532"
#        ]
    }
    
    # Initializes the dataset generator.
    filenames = get_filenames(DIRECTORIES, START_DATE, END_DATE)
    datasets  = file_loader._preprocess_data(filenames, datasets)
    
    nc_outfile = nc.Dataset("2007_water_histograms.nc", "w", format = "NETCDF4_CLASSIC")
    
#    x = next(datasets)
#    print(x)
#    print(x["timestamp"].strftime("%Y-%m"))
#    print("Exiting...")
#    quit()
    
    # Defines logscale bin edges and centers.
    cod_bin_edges   = np.logspace(start = -2, stop = 3, num = 151, base = 10)
    cod_bin_centers = cod_bin_edges[:-1] + np.diff(cod_bin_edges) / 2
    
    sza_bin_edges   = np.linspace(0, 90, 91)
    sza_bin_centers = sza_bin_edges[:-1] + np.diff(sza_bin_edges) / 2
    
    cod_grid, sza_grid = np.meshgrid(cod_bin_centers, sza_bin_centers)
    
    nc_outfile.createDimension("cod_bin_edges", cod_bin_edges.size)
    nc_outfile.createDimension("sza_bin_edges", sza_bin_edges.size)
    nc_outfile.createDimension("cod_bin_centers", cod_bin_centers.size)
    nc_outfile.createDimension("sza_bin_centers", sza_bin_centers.size)
    nc_outfile.createVariable("cod_bin_edges", np.float64, ("cod_bin_edges",))
    nc_outfile.createVariable("sza_bin_edges", np.float64, ("sza_bin_edges",))
    nc_outfile.createVariable("cod_bin_centers", np.float64, ("cod_bin_centers",))
    nc_outfile.createVariable("sza_bin_centers", np.float64, ("sza_bin_centers",))
    
    nc_outfile["cod_bin_edges"][:] = cod_bin_edges
    nc_outfile["sza_bin_edges"][:] = sza_bin_edges
    nc_outfile["cod_bin_centers"][:] = cod_bin_centers
    nc_outfile["sza_bin_centers"][:] = sza_bin_centers
        
    for product in ["transparent", "opaque"]:
        
        for month in range(1, 13):
            
            key = f"2007-{month:02d}_{product}"

            nc_outfile.createVariable(key, np.float64, ("cod_bin_centers", "sza_bin_centers"))
            nc_outfile[key][:] = np.zeros((cod_bin_centers.size, sza_bin_centers.size))
    
    # Sets up empty logscale histograms to add to.
    opaque_histogram      = np.histogram([], cod_bin_edges)[0]
    transparent_histogram = np.histogram([], cod_bin_edges)[0]
    
    opaque_2d_histogram      = np.histogram2d([], [], (cod_bin_edges, sza_bin_edges))[0]
    transparent_2d_histogram = np.histogram2d([], [], (cod_bin_edges, sza_bin_edges))[0]
    
    anomalies = []
    
    i = 0
    
    # TODO: "continue" lines in loader seem to make tqdm inaccurate.
    
    print(filenames)
    
    # Iterates through datasets.
    for dataset in tqdm(datasets, total = filenames.index.size):
                
        i += 1
        
        base_key        = dataset["timestamp"].strftime("%Y-%m")
        transparent_key = f"{base_key}_transparent"
        opaque_key      = f"{base_key}_opaque"
                
        # Aliases arrays from the dataset for brevity.
        vfm     = dataset["Feature_Classification_Flags"]
        opacity = dataset["Opacity_Flag"]
        cod     = dataset["MYD06_Cloud_Optical_Thickness"]
        slay    = dataset["Number_Layers_Found"]
        sza     = dataset["MODIS_SZA1km"]  # dataset["Solar_Zenith_Angle"]

        try:
            
            slay = np.squeeze(slay)
            sza  = np.squeeze(sza)

            # Masks water clouds.
            water_clouds = ((vfm.T[0] == 2) & (vfm.T[2] == 2)).T
            
            # TEMP: Masks ice clouds.
            #water_clouds = ((vfm.T[0] == 2) & ((vfm.T[2] == 1) | (vfm.T[2] == 3))).T

            # Masks single-layer water clouds.
            single_layer_water_clouds = water_clouds[:, 0] & (slay == 1)

            # Identifies single-layer opaque and transparent clouds.
            single_layer_opacity      = single_layer_water_clouds & (opacity[:, 0] == 1)
            single_layer_transparency = single_layer_water_clouds & (opacity[:, 0] == 0)
        
            # Separates transparent and opaque cloud optical depth measurements.
            opaque_cod      = cod[single_layer_opacity]
            transparent_cod = cod[single_layer_transparency]

            opaque_sza      = sza[single_layer_opacity]
            transparent_sza = sza[single_layer_transparency]
                        
            transparent_time = np.squeeze(dataset["Profile_Time"][single_layer_transparency])
            transparent_time = np.array([datetime.timedelta(seconds = t) for t in transparent_time])
            transparent_time = datetime.datetime(1993, 1, 1) + transparent_time
            
            transparent_lat = dataset["CALIPSO_Lat1km"][single_layer_transparency]
            transparent_lon = dataset["CALIPSO_Lon1km"][single_layer_transparency]
            
            anomaly_mask = transparent_cod == 150
            
            file_anomalies = pd.DataFrame()
            file_anomalies["timestamp"] = transparent_time[anomaly_mask]
            file_anomalies["daylight"] = np.squeeze(dataset["Day_Night_Flag"][single_layer_transparency][anomaly_mask]) == 0
            file_anomalies["latitude"] = transparent_lat[anomaly_mask]
            file_anomalies["longitude"] = transparent_lon[anomaly_mask]
            file_anomalies["over_water"] = np.squeeze(dataset["IGBP_Surface_Type"][single_layer_transparency][anomaly_mask]) == 17
            file_anomalies["cloud"] = single_layer_transparency[single_layer_transparency][anomaly_mask]
            file_anomalies["single_layer"] = single_layer_transparency[single_layer_transparency][anomaly_mask]
            file_anomalies["transparent"] = single_layer_transparency[single_layer_transparency][anomaly_mask]
            file_anomalies["cod"] = transparent_cod[anomaly_mask]
            file_anomalies["sza"] = transparent_sza[anomaly_mask]
            
            anomalies.append(file_anomalies)
            
            #print(pd.Series.dt.strftime(transparent_time))
            
        except:  # REMOVE MODULENOTFOUNDERROR FOR SERVER> LOCAL TESTING ONLY
            
            print(filenames.iloc[i])
            print(vfm.shape)
            print(opacity.shape)
            print(cod.shape)
            print(sza.shape)
            print(slay.shape)
        
        # Adds cloud optical depth results to the histogram.
        opaque_histogram      += np.histogram(opaque_cod, cod_bin_edges)[0]
        transparent_histogram += np.histogram(transparent_cod, cod_bin_edges)[0]
        
        opaque_2d_histogram      += np.histogram2d(opaque_cod, opaque_sza, (cod_bin_edges, sza_bin_edges))[0]
        transparent_2d_histogram += np.histogram2d(transparent_cod, transparent_sza, (cod_bin_edges, sza_bin_edges))[0]
        
        nc_outfile[transparent_key][:] += np.histogram2d(transparent_cod, transparent_sza, (cod_bin_edges, sza_bin_edges))[0]
        nc_outfile[opaque_key][:]      += np.histogram2d(opaque_cod, opaque_sza, (cod_bin_edges, sza_bin_edges))[0]
        
    # Sums the transparent and opaque histogram as the observation count.
    total_histogram = \
        transparent_histogram + opaque_histogram
        
    total_2d_histogram = \
        transparent_2d_histogram + opaque_2d_histogram
    
    anomalies = pd.concat(anomalies)
    anomalies.to_csv("Water_Anomalies.csv", index = False)
    
    nc_outfile.close()
    
    print("Script execution completed.")
    
#    print(transparent_2d_histogram.shape)
#    print(total_2d_histogram.shape)
#    print(cod_bin_centers.shape)
#    print(sza_bin_centers.shape)
#    
#    # Plots the histogram.
#    plt.figure()
#    
#    plt.suptitle("MODIS Cloud Optical Depth Histogram")
#    plt.title(f"{START_DATE} to {END_DATE}")
#    
#    plt.step(cod_bin_centers,
#             opaque_histogram,
#             c     = "b",
#             label = "MODIS Opaque")
#    
#    plt.step(cod_bin_centers,
#             transparent_histogram,
#             c     = "r",
#             label = "MODIS Transparent")
#    
#    plt.step(cod_bin_centers,
#             total_histogram,
#             c     = "k",
#             label = "MODIS Total")
#    
#    plt.xscale("log")
#    plt.xlim(cod_bin_edges.min(), cod_bin_edges.max())
#    plt.ylim(bottom = 0)
#    plt.xlabel(r"Cloud Optical Depth ($COD$)")
#    plt.ylabel(r"Observations ($N$)")
#    plt.legend()
#    
#    plt.savefig("COD_Histogram.png")
#    
#    """
#    __main__.py:240: RuntimeWarning: divide by zero encountered in true_divide
#      plt.pcolormesh(cod_grid, sza_grid, (transparent_2d_histogram / total_2d_histogram).T, cmap = "jet", vmin = 0, vmax = 1)
#    __main__.py:240: RuntimeWarning: invalid value encountered in true_divide
#      plt.pcolormesh(cod_grid, sza_grid, (transparent_2d_histogram / total_2d_histogram).T, cmap = "jet", vmin = 0, vmax = 1)
#    """
#    
#    plt.figure()
#    plt.pcolormesh(cod_grid, sza_grid, (transparent_2d_histogram / total_2d_histogram).T, cmap = "jet", vmin = 0, vmax = 1)
#    cbar = plt.colorbar()
#    cbar.set_label("Transparent Cloud Fraction")
#    cs = plt.contour(cod_grid, sza_grid, total_2d_histogram.T / total_2d_histogram.max(), levels = np.arange(0.15, 1, 0.15), cmap = "jet", vmin = 0, vmax = 1)
#    plt.clabel(cs)
#    plt.xscale("log")
#    plt.xlim(10 ** -1, 150)
#    plt.ylim(20, 80)
#    plt.xlabel(r"Cloud Optical Depth ($COD$)")
#    plt.ylabel(r"Solar Zenith Angle ($SZA$)")
#    plt.savefig("COD_2D_Histogram.png")
