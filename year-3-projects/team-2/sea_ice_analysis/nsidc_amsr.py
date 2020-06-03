# Standard library imports.
import datetime
import os
import re
from typing import Dict
from typing import Tuple

# Third party imports.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

# Constant definitions.
AMSR_12KM_REGEX = (r"AMSR_E_L3_SeaIce12km_(?P<version>V\d+15)_"
                   r"(?P<year>\d\d\d\d)(?P<month>\d\d)(?P<day>\d\d).hdf")


def get_geodetic_crs() -> Dict:
    
    '''
    
    https://nsidc.org/data/polar-stereo/tools_geo_pixel.html
    
    '''
    
    def read_file(file : str, shape : Tuple[int, int]) -> np.ndarray:
        
        data = np.fromfile(file, dtype = np.int32).reshape(shape) / 1e5
    
        return data
    
    nh_shape = (896, 608)  # "nh" represents "northern hemisphere".
    sh_shape = (664, 632)  # "sh" represents "southern hemisphere".
    
    coordinates = {
        "northern_hemisphere" : {
            "latitudes"  : read_file("psn12lats_v3.dat", nh_shape),
            "longitudes" : read_file("psn12lons_v3.dat", nh_shape),
            "shape"      : nh_shape,
        },
        
        "southern_hemisphere" : {
            "latitudes"  : read_file("pss12lats_v3.dat", sh_shape),
            "longitudes" : read_file("pss12lons_v3.dat", sh_shape),
            "shape"      : sh_shape,
        }
    }
    
    return coordinates


class AMSR12kmFile:
    
    coordinates = get_geodetic_crs()
    
    def __init__(self, path : str) -> None:
        
        self.path     = os.path.abspath(path)
        self.filename = os.path.basename(self.path)
        
        self._parse_filename()
        
        self.dataset = nc.Dataset(self.path)
    
    def __getitem__(self, name) -> nc.Variable:
        
        return self.dataset[name]
        
    def _parse_filename(self) -> None:
        
        regex_match = re.match(AMSR_12KM_REGEX, self.filename)
        
        if regex_match:
            
            filename_metadata = regex_match.groupdict()
            
        else:
            
            error_message = (
                "the given file is not a valid NSIDC AMSR level 3 sea ice "
                "12km file."
            )
            
            raise ValueError(error_message)
            
        year      = filename_metadata["year"]
        month     = filename_metadata["month"]
        day       = filename_metadata["day"]
        version   = filename_metadata["version"]
        timestamp = datetime.date(year, month, day)
        
        self.version   = version
        self.timestamp = timestamp
    
    
if __name__ == "__main__":
    
    # https://nsidc.org/data/AE_SI12/versions/3
    dataset = nc.Dataset("AMSR_E_L3_SeaIce12km_V15_20070101.hdf")    
        
    nh_ice_concentration = dataset["SI_12km_NH_ICECON_ASC"][:].data.astype(np.float)
    sh_ice_concentration = dataset["SI_12km_SH_ICECON_ASC"][:].data.astype(np.float)
    
    nh_water_mask = np.where(nh_ice_concentration == 0)
    nh_ice_mask   = np.where((1 <= nh_ice_concentration) & (nh_ice_concentration <= 100))
    nh_nan_mask   = np.where(nh_ice_concentration == 110)
    nh_land_mask  = np.where(nh_ice_concentration == 120)
    
    sh_water_mask = np.where(sh_ice_concentration == 0)
    sh_ice_mask   = np.where((1 <= sh_ice_concentration) & (sh_ice_concentration <= 100))
    sh_nan_mask   = np.where(sh_ice_concentration == 110)
    sh_land_mask  = np.where(sh_ice_concentration == 120)
    
    nh_water = nh_ice_concentration[nh_water_mask]
    nh_ice   = nh_ice_concentration[nh_ice_mask]
    nh_land  = nh_ice_concentration[nh_land_mask]
    
    sh_water = sh_ice_concentration[sh_water_mask]
    sh_ice   = sh_ice_concentration[sh_ice_mask]
    sh_land  = sh_ice_concentration[sh_land_mask]
    
    plt.figure()
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90])
    
    img = np.tile(np.array([[[169, 169, 169]]], dtype = np.uint8), [2, 2, 1])
    
    ax.imshow(img,
              transform=ccrs.PlateCarree(),
              extent=[-180, 180, -90, 90])
    
    ax.coastlines("50m")
    
    nhshape = (896, 608)
    shshape = (664, 632)

    nh_lons = np.fromfile("psn12lons_v3.dat", dtype = np.int32).reshape(nhshape) / 1e5
    nh_lats = np.fromfile("psn12lats_v3.dat", dtype = np.int32).reshape(nhshape) / 1e5
    sh_lons = np.fromfile("pss12lons_v3.dat", dtype = np.int32).reshape(shshape) / 1e5
    sh_lats = np.fromfile("pss12lats_v3.dat", dtype = np.int32).reshape(shshape) / 1e5
    
    ax.scatter(nh_lons[nh_water_mask],
               nh_lats[nh_water_mask],
               c = "b", s = 0.5,
               transform = ccrs.PlateCarree())
    
    ax.scatter(nh_lons[nh_ice_mask],
               nh_lats[nh_ice_mask],
               c = nh_ice_concentration[nh_ice_mask], cmap = "gray", vmin = -200, s = 0.5,
               transform = ccrs.PlateCarree())
    
    ax.scatter(nh_lons[nh_land_mask],
               nh_lats[nh_land_mask],
               c = "g", s = 0.5,
               transform = ccrs.PlateCarree())
        
    ax.scatter(sh_lons[sh_water_mask],
               sh_lats[sh_water_mask],
               c = "b", s = 0.5,
               transform = ccrs.PlateCarree())
    
    ax.scatter(sh_lons[sh_ice_mask],
               sh_lats[sh_ice_mask],
               c = sh_ice_concentration[sh_ice_mask], cmap = "gray", vmin = -200, s = 0.5,
               transform = ccrs.PlateCarree())
    
    ax.scatter(sh_lons[sh_land_mask],
               sh_lats[sh_land_mask],
               c = "g", s = 0.5,
               transform = ccrs.PlateCarree())
    
    resolution = 0.25
    latbins    = np.arange(-90,  90  + resolution, resolution)
    lonbins    = np.arange(-180, 180 + resolution, resolution)
    
    from scipy.stats import binned_statistic_2d
    
    nh_ice_mask = np.where((0 <= nh_ice_concentration) & (nh_ice_concentration <= 100))
    sh_ice_mask = np.where((0 <= sh_ice_concentration) & (sh_ice_concentration <= 100))
    
    nh_ice = nh_ice_concentration[nh_ice_mask]
    sh_ice = sh_ice_concentration[sh_ice_mask]
        
    nh_sums, _, _, _   = binned_statistic_2d(nh_lons[nh_ice_mask], nh_lats[nh_ice_mask], nh_ice, "sum",   bins = [lonbins, latbins])
    nh_counts, _, _, _ = binned_statistic_2d(nh_lons[nh_ice_mask], nh_lats[nh_ice_mask], nh_ice, "count", bins = [lonbins, latbins])
    
    sh_sums, _, _, _   = binned_statistic_2d(sh_lons[sh_ice_mask], sh_lats[sh_ice_mask], sh_ice, "sum",   bins = [lonbins, latbins])
    sh_counts, _, _, _ = binned_statistic_2d(sh_lons[sh_ice_mask], sh_lats[sh_ice_mask], sh_ice, "count", bins = [lonbins, latbins])
    
    with np.errstate(divide='ignore', invalid='ignore'):
    
        nh_avg = nh_sums / nh_counts
        sh_avg = sh_sums / sh_counts
            
    lats = latbins[:-1]
    lons = lonbins[:-1]
    
    fig  = plt.figure()
    ax   = plt.axes(projection = ccrs.PlateCarree())
    img1 = ax.pcolormesh(lons, lats, nh_avg.T, vmin = 0, vmax = 100, cmap = "jet")
    img2 = ax.pcolormesh(lons, lats, sh_avg.T, vmin = 0, vmax = 100, cmap = "jet")
    ax.coastlines("50m")
    plt.colorbar(img1, orientation = "horizontal")
    
    plt.show()
