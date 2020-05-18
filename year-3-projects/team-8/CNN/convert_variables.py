# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 15:23:04 2020

@author: Brice Coffer

Function to convert base state variables in RUC soundings to other 
possible representations.

input:  3D array of original RUC soundings (Array of float64)
output: 3D array RUC soundings with converted variables (Array of float64)

"""

import numpy as np
import metpy.calc as mpcalc
from metpy.units import units 

def convert_variables(input_data): 

    output_data = input_data
 
    for i in range(input_data.shape[0]):

        ## tdc to specific humidity
        prs = input_data[i, :, 0]
        tdc = input_data[i, :, 2]
        q = mpcalc.specific_humidity_from_dewpoint(tdc * units.degC,
                                                   prs * units.hPa)

        output_data[i, :, 2] = np.array(q)

        ## tc to potential temperature
        tc  = input_data[i, :, 1]
#       theta = mpcalc.potential_temperature(prs * units.hPa, 
#                                            tc * units.degC)
        thv = mpcalc.virtual_potential_temperature(prs * units.hPa,
                                                   tc * units.degC,
                                                   q)
        output_data[i, :, 1] = np.array(thv)

    return output_data
