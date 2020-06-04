import pandas as pd
import netCDF4 as nc
import numpy as np

if __name__ == "__main__":

    csv_filename   = "Water_Anomalies.csv"
    ncin_filename  = "2007_water_histograms.nc"
    ncout_filename = "2007_MYD_CALIPSO_water_cloud_collocation_discrepancies.nc"

    df   = pd.read_csv(csv_filename)
    ncin = nc.Dataset(ncin_filename)

    ncout = nc.Dataset(ncout_filename, "w")
    
    ncout.createGroup("histograms")
    ncout.createGroup("histograms/1d")
    ncout.createGroup("histograms/2d")
    ncout.createGroup("histograms/1d/transparent")
    ncout.createGroup("histograms/1d/opaque")
    ncout.createGroup("histograms/1d/transparent_fraction")
    ncout.createGroup("histograms/1d/opaque_fraction")
    ncout.createGroup("histograms/1d/total")
    ncout.createGroup("histograms/2d/transparent")
    ncout.createGroup("histograms/2d/opaque")
    ncout.createGroup("histograms/2d/transparent_fraction")
    ncout.createGroup("histograms/2d/opaque_fraction")
    ncout.createGroup("histograms/2d/total")
    
    ncout["histograms"].createDimension("cod_bin_edges",   ncin["cod_bin_edges"][:].size)
    ncout["histograms"].createDimension("sza_bin_edges",   ncin["sza_bin_edges"][:].size)
    ncout["histograms"].createDimension("cod_bin_centers", ncin["cod_bin_centers"][:].size)
    ncout["histograms"].createDimension("sza_bin_centers", ncin["sza_bin_centers"][:].size)
    
    ncout["histograms"].createVariable("cod_bin_edges",   np.float64, ("cod_bin_edges",))
    ncout["histograms"].createVariable("sza_bin_edges",   np.float64, ("sza_bin_edges",))
    ncout["histograms"].createVariable("cod_bin_centers", np.float64, ("cod_bin_centers",))
    ncout["histograms"].createVariable("sza_bin_centers", np.float64, ("sza_bin_centers",))
    
    ncout["histograms/cod_bin_edges"][:]   = ncin["cod_bin_edges"][:]
    ncout["histograms/sza_bin_edges"][:]   = ncin["sza_bin_edges"][:]
    ncout["histograms/cod_bin_centers"][:] = ncin["cod_bin_centers"][:]
    ncout["histograms/sza_bin_centers"][:] = ncin["sza_bin_centers"][:]
    
    ncout.createGroup("anomalies")
    
    ncout["anomalies"].createDimension("time", df.timestamp.size)
    
    # np.unicode_ type generates incorrect timestamp? Use np.datetime64?
    ncout["anomalies"].createVariable("calipso_timestamp",    np.unicode_, ("time",))
    ncout["anomalies"].createVariable("calipso_latitude",     np.float64,  ("time",))
    ncout["anomalies"].createVariable("calipso_longitude",    np.float64,  ("time",))
    ncout["anomalies"].createVariable("ceres_over_water",     np.int16,    ("time",))
    ncout["anomalies"].createVariable("modis_sza",            np.float64,  ("time",))
    ncout["anomalies"].createVariable("modis_cod",            np.float64,  ("time",))
    ncout["anomalies"].createVariable("calipso_single_layer", np.int16,    ("time",))
    ncout["anomalies"].createVariable("calipso_daylight",     np.int16,    ("time",))
    ncout["anomalies"].createVariable("calipso_cloud",        np.int16,    ("time",))
    ncout["anomalies"].createVariable("calipso_transparent",  np.int16,    ("time",))
    
    ncout["anomalies/calipso_timestamp"][:]    = df.timestamp.values
    ncout["anomalies/calipso_latitude"][:]     = df.latitude.values
    ncout["anomalies/calipso_longitude"][:]    = df.longitude.values
    ncout["anomalies/ceres_over_water"][:]     = df.over_water.values.astype(np.int16)
    ncout["anomalies/modis_sza"][:]            = df.sza.values
    ncout["anomalies/modis_cod"][:]            = df.cod.values
    ncout["anomalies/calipso_single_layer"][:] = df.single_layer.values.astype(np.int16)
    ncout["anomalies/calipso_daylight"][:]     = df.daylight.values.astype(np.int16)
    ncout["anomalies/calipso_cloud"][:]        = df.cloud.values.astype(np.int16)
    ncout["anomalies/calipso_transparent"][:]  = df.transparent.values.astype(np.int16)
    
    for product in ["transparent", "opaque"]:
    
        for month in range(1, 13):
            
            opposite = ["transparent", "opaque"]
            opposite.remove(product)
            opposite = opposite[0]
            
            old_key = f"2007-{month:02d}_{opposite}"
            new_key = f"2007-{month:02d}"

            ncout[f"histograms/2d/{product}"].createVariable(new_key, np.float64, ("cod_bin_centers", "sza_bin_centers"))
            ncout[f"histograms/2d/{product}/{new_key}"][:] = ncin[old_key][:]
            
    for month in range(1, 13):
        
        key = f"2007-{month:02d}"
    
        ncout["histograms/2d/total"].createVariable(key, np.float64, ("cod_bin_centers", "sza_bin_centers"))        
        ncout["histograms/2d/transparent_fraction"].createVariable(key, np.float64, ("cod_bin_centers", "sza_bin_centers"))
        ncout["histograms/2d/opaque_fraction"].createVariable(key, np.float64, ("cod_bin_centers", "sza_bin_centers"))
        
        ncout[f"histograms/2d/total/{key}"][:]                = np.zeros((ncin["cod_bin_centers"][:].size, ncin["sza_bin_centers"][:].size))
        ncout[f"histograms/2d/transparent_fraction/{key}"][:] = np.zeros((ncin["cod_bin_centers"][:].size, ncin["sza_bin_centers"][:].size))
        ncout[f"histograms/2d/opaque_fraction/{key}"][:]      = np.zeros((ncin["cod_bin_centers"][:].size, ncin["sza_bin_centers"][:].size))
        
        for product in ["transparent", "opaque"]:
            
            ncout[f"histograms/2d/total/{key}"][:] += ncout[f"histograms/2d/{product}/{key}"]
            
        for product in ["transparent", "opaque"]:
            
            ncout[f"histograms/2d/{product}_fraction/{key}"][:] = ncout[f"histograms/2d/{product}/{key}"][:] / ncout[f"histograms/2d/total/{key}"][:]

    for product in ["transparent", "opaque", "total"]:
    
        for month in range(1, 13):
            
            key = f"2007-{month:02d}"
            
            ncout[f"histograms/1d/{product}"].createVariable(key, np.float64, ("cod_bin_centers",))
            ncout[f"histograms/1d/{product}/{key}"][:] = np.sum(ncout[f"histograms/2d/{product}/{key}"][:], axis = -1)
            