# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:29:40 2020

@author: Brice Coffer

Function to read in RUC soundings (including metadata) and return 
a dictionary of metadata and an 3D array of sounding data, 
structured as a 20194 x 37 x 7 numpy array 
(i.e., number of soundings, x max size of soundings x sounding variables)

input: local path to files (str)
output: metadata from each RUC sounding (dict)
        3D array of RUC soundings (Array of float64)

"""

import os
import numpy as np

def readRUCsounding(path):

    ## list of metadata in filenames
    filename_vars=['datetime', 'eventtype', 'mag', 'lat', 'lon', 'st', 'radar', 
                   'radardist', 'mode', 'meso', 'tropical', 'mucape', 'mucin', 
                   'mulcl', 'shr6km', 'shr3km', 'shr1km', 'srh3km', 'srh1km', 
                   'tsfc', 'tdsfc', 'lr3km', 'lr700to500mb', 'pw', 'stp', 'ship']
    ## list of atmospheric variables in sounding files
    sounding_vars=['prs', 'tc', 'tdc', 'rh', 'u', 'v', 'z']
    
    # number of files
    numfiles = len(os.listdir(path))
    
    ## initialize dictionary for metadata
    RUCdata = {key:[] for key in filename_vars}
    
    ## initialize 3D array for RUC soundings
    ## hardwired value of 37 for maximum size of sounding
    RUCsounding = np.empty([numfiles, 37, len(sounding_vars)])
    
    ## initialize a list of event types (weak, significant, nontornadic) 
    eventtype = [" " for x in range(numfiles)]
    
    ## loop through each sounding file
    if path[-1] != '/': path = path + '/' 
    filenames=sorted(os.listdir(path))
    for i in range(numfiles):
        
        ## read in metadata from filename
        foo=filenames[i].split('_')
      
        ## loop through metadata from filenames and assign each to dictionary  
        for var in filename_vars:
            ## integer metadata
            if filename_vars.index(var) in [2, 10]:
                RUCdata[var].append(int(foo[filename_vars.index(var)]))
            ## float metadata
            elif filename_vars.index(var) in ([3, 4, 7] + list(range(11,26))):
                ## check for bad data
                if foo[filename_vars.index(var)] in (-9999.0, 'NA'):
                    RUCdata[var].append(float('NaN'))
                else:
                    RUCdata[var].append(float(foo[filename_vars.index(var)]))
            ## string metadata
            else:
                RUCdata[var].append(foo[filename_vars.index(var)])
             
        ## separate  reports into significantly tornadic, weakly tornadic, and nontornadic categories
        if 0 <= RUCdata['mag'][i] <= 1:
            RUCdata['eventtype'][i] = 1 # weakly tornadic
        elif 2 <= RUCdata['mag'][i] <= 5:
            RUCdata['eventtype'][i] = 2 # significantly tornadic
        else:
            RUCdata['eventtype'][i] = 0 # nontornadic
               
        ## read in sounding data for each file
        sounding = np.genfromtxt(path+filenames[i],names=sounding_vars)
        
        # put sounding into final array format
        for n in range(len(sounding_vars)):
            RUCsounding[i, 0:len(sounding), n] = sounding[sounding_vars[n]]
    
#   RUCsounding[216,1,2]   = 
#   RUCsounding[10777,1,2] =
#   RUCsounding[13247,1,2] =
#   RUCsounding[13369,1,2] =
#   RUCsounding[16592,1,2] =
#   RUCsounding[19250,1,2] = 

    return RUCdata, RUCsounding
