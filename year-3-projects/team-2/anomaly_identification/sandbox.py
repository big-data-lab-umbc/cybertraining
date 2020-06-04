# Standard library imports.
import os
from operator import add
from pprint import pprint

# Third party imports.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from pyhdf.SD import SD
from pyhdf.SD import SDC
from pyspark.sql import SparkSession
from tqdm import tqdm

# https://worldview.earthdata.nasa.gov/?v=-448.0833737998156,-262.1407233297318,297.2291262001844,90.54677667026814&t=2007-01-01-T14%3A34%3A44Z
# https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&TIME=2007-01-01T00:00:00Z&BBOX=-77.36858054226036,-32.2674641148307,-77.30478468899082,-32.20366826156117&CRS=EPSG:4326&LAYERS=MODIS_Terra_CorrectedReflectance_TrueColor,Coastlines&WRAP=day,x&FORMAT=image/jpeg&WIDTH=1&HEIGHT=1&ts=1587656674112
# https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&TIME=2020-04-23T00:00:00Z&BBOX=38.03065191387561,-20.881130382775126,38.665819377990445,-20.404754784689004&CRS=EPSG:4326&LAYERS=MODIS_Terra_CorrectedReflectance_TrueColor,Coastlines&WRAP=day,x&FORMAT=image/jpeg&WIDTH=1734&HEIGHT=2313&ts=1587657363454

# TODO: Timestamp must be in ISO format!!!

def get_image(row):
    
    index, row = row
    url, path = row["snapshot_url"], row["image_path"]
        
    retry = True
    
    while retry:
    
        try:

            r     = requests.get(url, timeout = 10)
            retry = False
            
        except requests.exceptions.ConnectTimeout:
            
            pass
    
    with open(path, "wb+") as file:
        
        file.write(r.content)
        
    return 1


def generate_url(row):
    
    filename = f"Anomaly {row.name + 1} - ({row['latitude']:.6f}, {row['longitude']:.6f}, {row['timestamp']}).png"
    filename = filename.replace(":", "-")
    path = os.path.join("Images", filename)
    
    url = (f"https://wvs.earthdata.nasa.gov/api/v1/snapshot?"
           f"REQUEST=GetSnapshot&"
           f"TIME={row['timestamp']}&"
           f"BBOX={row['min_lat']},{row['min_lon']},{row['max_lat']},{row['max_lon']}"
           f"&CRS=EPSG:4326&LAYERS=MODIS_Aqua_CorrectedReflectance_TrueColor,Coastlines&"
           f"WRAP=day,x&"
           f"FORMAT=image/png&"
           f"WIDTH=800&HEIGHT=800")
    
    url2 = (f"https://worldview.earthdata.nasa.gov/?v="
            f"{row['min_lon']},{row['min_lat']},{row['max_lon']},{row['max_lat']}&"
            f"t={row['timestamp'].replace(':', '%3A')}"
            f"&l=Coastlines,MODIS_Aqua_CorrectedReflectance_TrueColor")
        
    row["worldview_url"]   = url2
    row["snapshot_url"]   = url
    row["image_filename"] = filename
    row["image_path"]     = path
    
    return row


def bytes_to_string(data):
    
    data   = data.astype(np.uint8).view("S3").squeeze()
    string = "".join([datum.decode("utf-8") for datum in data])

    return string

    
if __name__ == "__main__":
    
    df = pd.read_csv("2007_over-water_anomalies.csv")
    
    df["iso_timestamp"] = pd.to_datetime(df["timestamp"]).dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    df["min_lon"]   = df["longitude"] - 0.1
    df["max_lon"]   = df["longitude"] + 0.1
    df["min_lat"]   = df["latitude"] - 0.1
    df["max_lat"]   = df["latitude"] + 0.1
    df["worldview_url"]       = pd.Series([], dtype = object)
    df["snapshot_url"]   = pd.Series([], dtype = object)
    df["image_filename"] = pd.Series([], dtype = object)
    df["image_path"]      = pd.Series([], dtype = object)
    
    df = df.apply(generate_url, axis = 1)
    
    columns = {
        "min_lat" : "minimum_bounding_latitude",
        "max_lat" : "maximum_bounding_latitude",
        "min_lon" : "minimum_bounding_longitude",
        "max_lon" : "maximum_bounding_longitude",
    }
    
    df.rename(columns = columns, inplace = True)
    df.to_csv("2007_water_worldview_anomalies.csv", index = False)
    
#    for index, row in tqdm(df.iterrows(), total = df.index.size):
#        
#        get_image(row["snapshot_url"], row["image_path"])
        
    spark   = SparkSession.builder.appName("sandbox").getOrCreate()
    sc      = spark.sparkContext
    rdd     = sc.parallelize(df.iterrows(), 100)
    results = rdd.map(get_image).reduce(add)
    print("Mapping complete.")
    #spark.stop()
    
#    filename   = "NISE_SSMIF13_20061231.HDFEOS"
#    data_file  = SD(filename, SDC.READ)
#    attributes = data_file.attributes()
#    datasets   = data_file.datasets()
#    
#    pprint(attributes)
#    pprint(datasets)
#    
#    # https://nsidc.org/data/NISE/versions/5
#    
#    age    = np.array(data_file.select("Age")[:]).astype(np.float32)
#    extent = np.array(data_file.select("Extent")[:]).astype(np.float32)
#    
#    age[age == 255]       = np.nan
#    extent[(extent == 254) | (extent == 252)] = np.nan
#    
#    '''
#    
#    The file format consists of a 300-byte descriptive header followed by a
#    two-dimensional array of one-byte values containing the data. The file
#    header is composed of:
#    
#    a 21-element array of 6-byte character strings that contain information
#    such as polar stereographic grid characteristics
#    
#    a 24-byte character string containing the file name
#    
#    a 80-character string containing an optional image title
#    
#    a 70-byte character string containing ancillary information such as data
#    origin, data set creation date, etc.
#    
#    The data can be read with image processing software by specifying a
#    300-byte header, with an image size of 304 columns x 448 rows for Arctic
#    data and 316 columns x 332 rows for Antarctic data. For example, in a
#    high-level programming language or image processing software, declare a
#    300-byte array for the header and a 304 x 448 Arctic image array. Read the
#    300-byte header array first, then read the image array.
#    
#    https://nsidc.org/support/21680984-How-do-I-import-the-0051-sea-ice-concentration-data-into-ArcGIS-
#    
#    https://nsidc.org/support/how/how-do-i-convert-nsidc-0051-sea-ice-concentration-data-binary-geotiff
#    
#    '''
#    
#    a = np.fromfile("nt_19781026_n07_v1.1_n.bin", dtype = np.uint8)
#    header = bytes_to_string(a[:300])
#    data   = a[300:].reshape((448, 304))
#    
#    plt.figure()
#    plt.pcolormesh(data)
#    plt.colorbar()
#    plt.show()
#    
#    plt.figure()
#    plt.title("Age")
#    img  = plt.contourf(age)
#    cbar = plt.colorbar()
#    
#    plt.figure()
#    plt.title("Extent")
#    img  = plt.contourf(extent)
#    cbar = plt.colorbar()
#    
#    plt.show()
