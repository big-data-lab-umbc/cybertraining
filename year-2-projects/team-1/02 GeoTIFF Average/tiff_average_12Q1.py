#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Annual Mean Program: Average MCD12Q1 (Land Cover) GeoTIFF data converted by HEG monthly

Created on 2019

@author: Jianyu Zheng
"""

import os
import math
import itertools
import numpy as np
from osgeo import gdal
from netCDF4 import Dataset

def writeTiff(im_data,im_width,im_height,im_bands,im_geotrans,im_proj,path):
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'uint8' in im_data.dtype.name:
        datatype = gdal.GDT_UInt8
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_data = np.array([im_data])
    else:
        im_bands, (im_height, im_width) = 1,im_data.shape
        
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, im_width, im_height, im_bands, datatype)
    if(dataset!= None):
        dataset.SetGeoTransform(im_geotrans) 
        dataset.SetProjection(im_proj) 
    for i in range(im_bands):
        dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset

def readTif(fileName):
    dataset = gdal.Open(fileName)
    if dataset == None:
        print(fileName+"can't open the file")
        return
    im_width  = dataset.RasterXSize  
    im_height = dataset.RasterYSize 
    im_bands  = dataset.RasterCount  
    im_data   = dataset.ReadAsArray(0,0,im_width,im_height)
    im_geos   = dataset.GetGeoTransform()
    im_proj   = dataset.GetProjection()  

    return im_width,im_height,im_bands,im_data,im_geos,im_proj 
 
#---------------------STEP 1: Set up data path and time interval

MCD_dir= '/umbc/xfs1/cybertrn/cybertraining2019/team1/final_project/data/MODIS/MCD12Q1/tif_files/heg'
MCD_prefix = '/MCD12Q1.A'
regions = np.array(['.h11v04','.h11v05','.h12v04','.h12v05'])

for reg in range(len(regions)):
    reg_prefix=regions[reg]

    path = 'MCD12Q1.A2017'+reg_prefix+'_yearly.tif'

    str = os.popen("ls "+ MCD_dir + reg_prefix + MCD_prefix + "*.tif").read()
    filename = np.array(str.split("\n"))

    im_width,im_height,im_bands,im_data,im_geos,im_proj=readTif(filename[0])

    #----------For MCD12Q1: Shortwave Albedo-------------
    LC = np.array(im_data,dtype=np.int)
    LC[np.where(LC == np.float(255))] = 0.0 #np.nan 
    print(LC)
    writeTiff(LC,im_width,im_height,im_bands,im_geos,im_proj,path)
    #print(reg_prefix,numfile)



