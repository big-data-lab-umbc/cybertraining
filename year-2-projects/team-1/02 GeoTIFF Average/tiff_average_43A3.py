#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Annual Mean Program: Average MCD43A3 (Surface Albedo) GeoTIFF data converted by HEG monthly

Created on 2019

@author: Jianyu Zheng
"""

import os
import math
import itertools
import numpy as np
from osgeo import gdal
from netCDF4 import Dataset

#gdal_merge.py

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

MCD_dir= '/umbc/xfs1/cybertrn/cybertraining2019/team1/final_project/data/MODIS/MCD43A3/tif_files/heg'
MCD_prefix = '/MCD43A3.A'
regions = np.array(['.h11v04','.h11v05','.h12v04','.h12v05'])

numfile = 0

years = np.arange(2016,2018)
month = np.arange(1,13)
mc_days1 = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
mc_days2 = np.array([31,29,31,30,31,30,31,31,30,31,30,31])

for reg in range(len(regions)):
	reg_prefix=regions[reg]

	for m in month:
		mc = '%02i' % m
		path = 'MCD43A3.A2016-17'+reg_prefix+'_'+mc+'.tif'

		for y in years:
			if y == 2016: 
				mc_days = mc_days2
			else:
				mc_days = mc_days1

			yc = '%04i' % y
			if (m-1) == 0:
				days = np.arange(1,mc_days[m-1]+1)
			else: 
				days = np.arange(1,mc_days[m-1]+1) + np.sum(mc_days[0:m-1]) 
			print(yc,m,mc_days[m-1],reg_prefix)

			for day in days:
				dc = '%03i' % day
				#print(yc,dc)
				str = os.popen("ls "+ MCD_dir + reg_prefix + MCD_prefix + yc + dc + "*.tif").read()
				filename = np.array(str.split("\n"))
				#print(filename[0])
				#filename = np.delete(filename,len(filename)-1)

				im_width,im_height,im_bands,im_data,im_geos,im_proj=readTif(filename[0])
				
				#----------For MCD43A3: Shortwave Albedo-------------
				albedo = np.array(im_data,dtype=np.float)
				albedo[np.where(albedo == np.float(32767))] = 0.0 #np.nan 
				albedo = albedo*0.001
				#print(albedo)
				if (day == days[0]) & (y == years[0]):
					num_var= np.zeros((albedo.shape[0],albedo.shape[1]))
					num_var[np.where(albedo != 0.0)] += 1
					var = albedo
				else:
					num_var[np.where(albedo != 0.0)] += 1
					var += albedo
				
		ave_var = var / num_var #(len(days)*len(years))
		ave_var[np.where(np.isnan(ave_var))]=0.0
		writeTiff(ave_var,im_width,im_height,im_bands,im_geos,im_proj,path)
		#print(reg_prefix,numfile)


