"""
Created on Thu Apr 16 15:29:40 2020

@author: Brice Coffer

Function to read in RUC soundings (including metadata) and return 
Used in RF_model.py.

input: local path to RUCdata files (string)

output: pandas dataframe
"""

import os
import numpy as np
import pandas as pd

def readRUCsoundings(path):

    ## list of metadata in filenames
    filename_vars=['datetime', 'eventtype', 'mag', 'lat', 'lon', 'st', 'radar', 
                   'radardist', 'mode', 'meso', 'tropical', 'mucape', 'mucin', 
                   'mulcl', 'shr6km', 'shr3km', 'shr1km', 'srh3km', 'srh1km', 
                   'tsfc', 'tdsfc', 'lr3km', 'lr700to500mb', 'pw', 'stp', 'ship']
    ## list of atmospheric variables in sounding files
    sounding_vars=['prs', 'tc', 'tdc', 'rh', 'u', 'v', 'z']
    
    ## initialize dictionaries for metadata and soundings
    RUCdata = {key:[] for key in filename_vars}
    RUCsounding = {key:[] for key in sounding_vars} 
    
    ## initialize a list of event types (weak, significant, nontornadic) 
    eventtype = [" " for x in range(len(sorted(os.listdir(path))))]
    
    ## loop through each sounding file
    i=0
    if path[-1] != '/': path = path + '/' 
    for filename in sorted(os.listdir(path)):
        
        ## read in metadata from filename
        foo=filename.split('_')
      
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
             
        ## separate  reports into weakly tornadic, significantly tornadic, and 
        ## nontoradic categories
        if 0 <= RUCdata['mag'][i] <= 1:
            eventtype[i] = "weakly\ntornadic"
        elif 2 <= RUCdata['mag'][i] <= 5:
            eventtype[i] = "significantly\ntornadic"
        else:
            eventtype[i] = "nontornadic"
        i+=1
               
        ## read in sounding data for each file
        sounding = np.genfromtxt(path+filename,names=sounding_vars)
        
        RUCsounding['prs'].append(sounding['prs'])
        RUCsounding['tc'].append(sounding['tc'])
        RUCsounding['tdc'].append(sounding['tdc'])
        RUCsounding['rh'].append(sounding['rh'])
        RUCsounding['u'].append(sounding['u'])
        RUCsounding['v'].append(sounding['v'])
        RUCsounding['z'].append(sounding['z'])    
    
    # create dataFrame        
    df = pd.DataFrame({ 'event'  : eventtype, 
                        'temperature' : RUCsounding['tc'],
                        'dewpoint'    : RUCsounding['tdc'],
                        'humidity'    : RUCsounding['rh'],
                        'u-wind'      : RUCsounding['u'],
                        'v-wind'      : RUCsounding['v'],
                        'pressure'    : RUCsounding['prs'],
                        'height'      : RUCsounding['z']})
    return df
