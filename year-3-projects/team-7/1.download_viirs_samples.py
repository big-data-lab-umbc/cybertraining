'''
This step is to use collocated calipso-viirs dataset to download VIIRS granule based on filename.
Substeps:
1) get_bbox()
    input/requires: path of calipso-viirs dataset on taki
    output: a viirs_cal_bbox.log file that records the collocated file name and their bounding box: minlon, maxlon, minlat, and maxlat
2) downsize()
    input/requires: viirs_cal_bbox.log, a desired bounding box: minlon, maxlon, minlat, and maxlat, a desired time period
    output: files_downsized.csv a file with the downsized set of collocated files and their bboxes
3) generateSh()
    input/requires: files_downsized.csv file
    output: a download_VIIRS.sh that has the wget commands of all files related to this spatiotemporal range
    This sh file should be ran after this python execution to download the viirs granule data into a specific (viirs_data) folder
4) cal_pixel_label()
    Run on taki, read sampled collocated h5 file and save lat, lon, pixel_label, and CALIOP_Alay_Aerosol_Type_Mode into a npy file
    input: viirs_cal_bbox.log
    output: npy files into the save_folder = '/umbc/xfs1/cybertrn/cybertraining2020/team7/research/VIIRS-SIPS/cal_label/'
'''

from netCDF4 import Dataset
import pandas as pd
import numpy as np
import os, datetime
from jdcal import gcal2jd
from datetime import datetime, timedelta
import requests
from retrying import retry

def get_bbox():
    # Loop through all data in the collocated folder and get their bounding boxes
    # done once, no need to call again
    root_path = '/umbc/xfs1/cybertrn/common/Data/calipso-virrs-collocated/data_for_testv2/2014/'

    listOfFiles = []
    for (dirpath, dirnames, filenames) in os.walk(root_path):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames if file.endswith('h5')]

    list = []
    for h5file in listOfFiles:
        VNP_CAL = Dataset(h5file)
        lat = VNP_CAL['Latitude'][:]
        lon = VNP_CAL['Longitude'][:]
        minlon, maxlon, minlat, maxlat = lon.min(), lon.max(), lat.min(), lat.max()
        list.append([h5file, minlon, maxlon, minlat, maxlat])
        print(h5file, minlon, maxlon, minlat, maxlat)

    df = pd.DataFrame(list, index=False, columns=['filename', 'minlon', 'maxlon', 'minlat', 'maxlat'])
    df.to_csv('viirs_cal_bbox.log')

def cal_pixel_label():
    # Run on taki, read sampled collocated h5 file and save lat, lon, pixel_label, and CALIOP_Alay_Aerosol_Type_Mode into a npy file
    # done once, no need to call again
    save_folder = '/umbc/xfs1/cybertrn/cybertraining2020/team7/research/VIIRS-SIPS/cal_label/'
    df = pd.read_csv('viirs_cal_bbox.log', index_col=False, header=None, delimiter=' ')
    for i in range(df.shape[0]):
        fname = df.loc[i, 0]
        VNP_CAL = Dataset(fname)
        lat = VNP_CAL['Latitude'][:]
        lon = VNP_CAL['Longitude'][:]
        Pixel_Label = VNP_CAL['Pixel_Label'][:]
        CALIOP_Alay_Aerosol_Type_Mode = VNP_CAL['CALIOP_Alay_Aerosol_Type_Mode'][:]
        CALIOP_N_Clay_5km = VNP_CAL['CALIOP_N_Clay_5km'][:]
        CALIOP_N_Clay_1km = VNP_CAL['CALIOP_N_Clay_1km'][:]
        IGBP_SurfaceType = VNP_CAL['IGBP_SurfaceType'][:]

        np.savez_compressed(save_folder + fname.split('.')[-3][3:], lon=lon, lat=lat, Pixel_Label=Pixel_Label,
                            CALIOP_Alay_Aerosol_Type_Mode=CALIOP_Alay_Aerosol_Type_Mode, CALIOP_N_Clay_5km=CALIOP_N_Clay_5km,
                            CALIOP_N_Clay_1km = CALIOP_N_Clay_1km, IGBP_SurfaceType = IGBP_SurfaceType)

def downsize(lon_begin, lon_end, lat_begin, lat_end, start_month, end_month, vnp_fname):
    # downsize the viirs_cal_bbox.log to a specific spatiotemporal range
    file_object = open("viirs_cal_bbox.log", "r")
    lines = file_object.readlines()
    file_object.close()

    y = 2014
    JD01, JD02 = gcal2jd(y, 1, 1)
    JD1, JD2 = gcal2jd(y, start_month, 1)
    julian_begin = np.int((JD2 + JD1) - (JD01 + JD02) + 1)
    JD3, JD4 = gcal2jd(y, end_month, 31)
    julian_end = np.int((JD4 + JD3) - (JD01 + JD02) + 1)

    keep_line = []
    for lineno in range(len(lines)):
        line_inter = lines[lineno].rstrip()
        line_split = line_inter.split(' ')
        min_lon = float(line_split[1])
        max_lon = float(line_split[2])
        min_lat = float(line_split[3])
        max_lat = float(line_split[4])
        line_split2 = line_split[0].split('.')
        julian_day = int(line_split2[1][7:10])

        if ((julian_day >= julian_begin) & (julian_day <= julian_end)):
            if ((min_lon >= lon_begin) & (max_lon <= lon_end)):
                if ((min_lat >= lat_begin) & (max_lat <= lat_end)):
                    keep_line.append(lines[lineno].rstrip())

    df = pd.DataFrame(keep_line)
    df.to_csv(vnp_fname, index=False, header=False)

def retry_if_connection_error(exception):
    print(datetime.now())
    print(exception)
    return True
@retry(retry_on_exception=retry_if_connection_error, wait_fixed=10000)
def safe_request(url, **kwargs):
    return requests.get(url, **kwargs)

def generateSh(vnp_fname, sh_fname):
    # from files_downsized.csv, extract all vnp file name identifiers, e.g., vnp2014203t0230
    df = pd.read_csv(vnp_fname, index_col=False, header=None, delimiter=' ')
    urls = []

    for h5fname in df[0].values:
        fname = h5fname.split('.')[-3]
        julian_day = int(fname[7:10])
        time = fname[-4:]
        startTime = datetime(2014, 1, 1) + timedelta(julian_day - 1) + timedelta(hours=int(time[:2])) + timedelta(minutes=int(time[2:]))
        endTime = datetime(2014, 1, 1) + timedelta(julian_day - 1) + timedelta(hours=int(time[:2])) + timedelta(minutes=int(time[2:]) + 3)
        s = startTime.strftime("%Y-%m-%dT%H:%M:%SZ")
        e = endTime.strftime("%Y-%m-%dT%H:%M:%SZ")
        # Use the API to request filenames
        url = "https://sips.ssec.wisc.edu/api/v1/products/search.sh?products=VNP02MOD|VNP03MOD&satellite=snpp&start={0}&end={1}".format(s,e)
        r = safe_request(url)
        allLines = r.text.split("\n")
        for l in allLines:
            if l.startswith("https://"):
                urls.append(l)

    with open(sh_fname, mode='w', encoding='UTF-8') as file:
        file.write('# !/usr/bin/env sh\n')
        for url in urls:
            file.write('wget '+url+'\n')

if __name__ == '__main__':
    get_bbox()
    cal_pixel_label()
    # north atlantic = [-74, -20, 13, 43]
    downsize(-74, -20, 13, 43, 1, 12, 'files_north_atlantic.csv')
    generateSh('files_north_atlantic.csv', 'download_VIIRS_north_atlantic.sh')
    # asian dust = [105, 140, 20, 50] spring month: March, April, May
    downsize(105, 140, 20, 50, 3, 5, 'files_asian_spring.csv')
    generateSh('files_asian_spring.csv', 'download_VIIRS_asian_spring.sh')
    # north africa = [-20, 60, 0, 60] summer months: June, July, August
    downsize(-20, 60, 0, 60, 6, 8, 'files_north_africa_summer.csv')
    generateSh('files_north_africa_summer.csv', 'download_VIIRS_north_africa_summer.sh')