#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Annual Mean Program: Average MCD15A2 (Leaf Area Index) GeoTIFF data stitched by HEG monthly

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

    return

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

MCD_dir= '/umbc/xfs1/cybertrn/cybertraining2019/team1/final_project/data/MODIS/MCD15A2/tif_files/stitched'
MCD_prefix = '/MCD15A2H.A'
regions = np.array(['.h11v04','.h11v05','.h12v04','.h12v05'])

numfile = 0

years = np.arange(2016,2018)
month = np.arange(1,13)
mc_days1 = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
mc_days2 = np.array([31,29,31,30,31,30,31,31,30,31,30,31])

for m in month:
	mc = '%02i' % m
	#path = '../data/MODIS/MCD15A2-Potomac/MCD15A2.A2016-17_'+mc+'.tif'
	path = '../data/MODIS/MCD15A2-Potomac-no255/MCD15A2.A2016-17_'+mc+'.tif'
	lai_num = mc_days1[m-1]+mc_days2[m-1]

	for y in years:
		if y == 2016: 
			mc_days = mc_days2
		else:
			mc_days = mc_days1

		yc = '%04i' % y
		tot  = np.sum(mc_days[0:m])
		tot_pre = np.sum(mc_days[0:m-1])

		if (m-1) == 0:
			days = np.arange(1,8*3+2,8)
			day_sta = 32*m+1
		else: 
			n = 3 
			day_end = 32*(m-1)+n*8+2
			if day_end > (tot+1):
				day_end = 32*(m-1)+(n-1)*8+2				
				if day_end > (tot):
					day_end = 32*(m-1)+(n-2)*8+2
			#print(day_sta,day_end)
			days = np.arange(day_sta , day_end, 8)#tot-(tot%8)+2,8)
			if y == 2017: day_sta = day_end+7
			#print(m,days,32*(m-1)+3*8+2,tot)

		for day in days:
			dc = '%03i' % day
			str = os.popen("ls "+ MCD_dir + MCD_prefix + yc + dc + "*.tif").read()
			filename = np.array(str.split("\n"))
			if len(filename[0]) == 0: 
				lai_num -= 8
				continue
			
			im_width,im_height,im_bands,im_data,im_geos,im_proj=readTif(filename[0])
			
			#----------For MCD15A2: Leaf Area Index-------------
			lai_ori = np.array(im_data,dtype=np.float)	
			lai     = np.array(im_data,dtype=np.float)	
			#print(np.where(lai >= np.float(248)))			
			lai[np.where(lai >= np.float(248))] = 0.0 #np.nan 
			lai = lai*0.1
			#print('days:',day)
			if ((y == 2016) | (y == 2017)) & (m ==1):
				if day > (tot-9):
					#if y == 2016: print('>(tot-8):',day,m,tot)
					weg_lai = lai*(tot-day+1)
					if y == 2016: 
						pre_lai1 = lai*(8-(tot-day+1))
					else:
						pre_lai2 = lai*(8-(tot-day+1))
					var += weg_lai
					#lai_num += tot-day
				else:
					#if y == 2016: print('entire:',day,m,tot)
					if y == 2016: 
						var = lai*8
					else:
						var += lai*8
					#lai_num += 8
			else:
				
				if day > (tot-9):
					#if y == 2016: print('>(tot-8):',day,m,tot)
					weg_lai = lai*(tot-day+1)
					if y == 2016: 
						pre_lai1 = lai*(8-(tot-day+1))
					else:
						pre_lai2 = lai*(8-(tot-day+1))
					var += weg_lai
					#lai_num += tot-day
				elif ((day - tot_pre) <= 9) & (day < (tot-9)):
					#if y == 2016: print('<(tot-8)~entire:',day,m,tot)
					if y == 2016: 
						var = pre_lai1 + lai*8
						#lai_num = day-mc_days[m-2]-1 + 8
					else:
						var += pre_lai2 + lai*8
						#lai_num += day-mc_days[m-2]-1 + 8
				else:
					#if y == 2016: print('entire:',day,m,tot)
					var += lai*8
					#lai_num += 8
			
	#print(mc_days[m-1]*len(years))
	#print(lai_num)
	ave_var = var / lai_num #(mc_days1[m-1]+mc_days2[m-1])
	#ave_var[np.where(lai_ori >= np.float(248))] = lai_ori[np.where(lai_ori >= np.float(248))]
	ave_var[np.where(ave_var == 0.0)]= 0.0001
	print(ave_var)
	#ave_var[np.where(np.isnan(ave_var))]=0.0
	writeTiff(ave_var,im_width,im_height,im_bands,im_geos,im_proj,path)
	#print(reg_prefix,numfile)

