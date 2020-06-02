# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:51:21 2020

@author: Brice Coffer

Function to interpolate varying height levels to a fixed vertical grid above
ground level (AGL).

input:  3D array of original RUC soundings (Array of float64)
output: 3D array of interpolated RUC soundings (Array of float64)

"""

import numpy as np

def interpolate_height(RUCsounding): 
       
    ## create a fixed height grid up to 13010 m, spaced every 200 m
    zlevs = np.linspace(10, 13010, 66)

    ## fixed 37 level height grid for interpolation
    zlevs = np.genfromtxt('/home/becoffer/research/RandomForest/DataInfo/avg_heights.txt')

#   zlevs =  np.array([1.00000000e+01, 1.58563842e+02, 3.84481862e+02, 
#      6.15252704e+02, 8.51109038e+02, 1.09229353e+03, 1.33908910e+03, 
#      1.59183214e+03, 1.85088787e+03, 2.11659829e+03, 2.38938211e+03, 
#      2.66966766e+03, 2.95778574e+03, 3.25452807e+03, 3.56008802e+03, 
#      3.87540070e+03, 4.20102989e+03, 4.53768426e+03, 4.88650312e+03, 
#      5.24808010e+03, 5.62423990e+03, 6.01540546e+03, 6.42405460e+03, 
#      6.85136080e+03, 7.29916298e+03, 7.77010236e+03, 8.26674070e+03, 
#      8.79218028e+03, 9.35030267e+03, 9.94611583e+03, 1.05856779e+04, 
#      1.12778015e+04, 1.20348609e+04, 1.28763354e+04, 1.38332626e+04, 
#      1.49549642e+04, 1.63174206e+04])   
 
    ## initialize 3D array for interpolated RUC soundings
    RUCsounding_interp = np.empty([RUCsounding.shape[0],len(zlevs), 
                                  RUCsounding.shape[2]])
    
    ## loop through each sounding variable in each file for interpolation 
    for i in range(RUCsounding.shape[0]):
        for n in range(RUCsounding.shape[2]):
            # trim zeros from end of shorter sounding files
            height = np.trim_zeros(RUCsounding[i, :, RUCsounding.shape[2]-1], 
                                   trim='b')
            interp_array = np.trim_zeros(RUCsounding[i, :, n], 
                                         trim='b')
            # interpolate to fixed height grid
            RUCsounding_interp[i, :, n] = np.interp(zlevs, height, 
                                                      interp_array)
            
    return RUCsounding_interp 
