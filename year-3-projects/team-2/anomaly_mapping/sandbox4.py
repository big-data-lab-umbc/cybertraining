from skyfield.api import Topos, load
from skyfield.api import EarthSatellite
from skyfield.api import utc
import re
import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cartopy.feature.nightshade import Nightshade


stations_url = 'http://celestrak.com/NORAD/elements/science.txt'
satellites = load.tle_file(stations_url)
print('Loaded', len(satellites), 'satellites')

ts = load.timescale()
satellite = [s for s in satellites if s.name == "CALIPSO"][0]
print(satellite)

t = []
tt = []
x = datetime.datetime(2007, 1, 1,0, 0, tzinfo = utc)

for i in range(48 * 60 + 1):
    
    tt.append(x)
    t.append(ts.utc(x))
    x += datetime.timedelta(minutes = 1)

lats = []
lons = []
    
for i, ti in enumerate(t):
    
    geocentric = satellite.at(ti)
    subpoint = geocentric.subpoint()
    #print('Time:', tt[i])
    #print('Latitude:', subpoint.latitude.degrees)
    #print('Longitude:', subpoint.longitude.degrees)
    
    lats.append(subpoint.latitude.degrees)
    lons.append(subpoint.longitude.degrees)

# https://www.ncdc.noaa.gov/gridsat/docs/Angle_Calculations.pdf

def sza(latitude, longitude, time):
    
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    
    # Compute the fractional Julian day.
    reference_time  = datetime.datetime(time.year, time.month, time.day)
    fractional_day  = (time - reference_time).total_seconds() / (24 * 60 * 60)
    julian_day      = int(time.strftime("%j")) + fractional_day - 1
    
    # Aliasing for brevity.
    N = julian_day
    
    epsilon           = np.radians(23.44)
    orbital_period    = 365.256
    angular_frequency = 2 * np.pi / orbital_period
    
    declination = -epsilon * np.cos(2 * np.pi / 365 * (N + 10))
    
    print(declination)
    
#    C1 = np.sin(-epsilon)
#    C2 = angular_frequency
#    C3 = 2 * 0.0167
#    
    hours = longitude / np.pi * 12
#    
#    declination = np.arcsin(C1 * np.cos(C2 * (N + 10) + C3 * np.sin(C2 * (N - 2))))
    hour_angle  = np.vectorize(lambda a: (reference_time - time + datetime.timedelta(hours = a)).total_seconds() / (60 * 60))(hours)
#    
#    print(declination, declination.min(), declination.max())
#    print(hour_angle, hour_angle.min(), hour_angle.max())
#    
    sza = np.arccos(np.sin(latitude) * np.sin(declination) + np.cos(latitude) * np.cos(declination) * np.cos(hour_angle))

    print(sza)
    
    ##    sza = np.degrees(sza)
    
    return sza



import numpy as np  # numerics & matrix algebra
import time  # for measuring time
import matplotlib as mpl  # plotting
import matplotlib.pyplot as plt  # plotting
import cartopy.crs as ccrs


#################################################################################
# USER GIVEN PARAMETERS

# time
t = datetime.datetime(2017, 4, 6, 9, 30, 0)  # 6th April 2017 09:30:00 UTC

# grid dimensions
Nlats = 90
Nlons = 180

import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u


def get_sza(latitude, longitude, time = Time.now()):

    loc = coord.EarthLocation(lat = latitude * u.deg, lon = longitude * u.deg)

    altaz = coord.AltAz(location = loc, obstime = time)
    sun = coord.get_sun(time)

    return sun.transform_to(altaz).zen.degree


lat, lon = np.linspace(-90.0, 90.0, Nlats + 1), np.linspace(-180.0, 180.0, Nlons + 1)  # lat and lon vectors for grid boundaries
latC, lonC = 0.5 * (lat[:-1] + lat[1:]), 0.5 * (lon[:-1] + lon[1:])  # center points

# make grid
latgrid, longrid = np.meshgrid(latC, lonC)

t0 = time.time()  # measure time to compute szas
# compute solar zenith angle and azimuth (be careful with azimuth: I haven't checked this at all)
sza = get_sza(latgrid.ravel(), longrid.ravel())
print('Computed {} solar zenith angles and azimuths and it took {:.04f} seconds'.format(len(longrid.ravel()), time.time() - t0))

print(sza.min(), sza.max())

cmap = mpl.cm.viridis

fig = plt.figure()
ax = plt.axes(projection = ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
im = ax.contourf(longrid, latgrid, sza.reshape(longrid.shape), levels = 50, alpha = 0.3)
fig.colorbar(im)
ax.stock_img()
ax.add_feature(Nightshade(datetime.datetime.utcnow(), alpha=0.2))
plt.show()

