'''
Created by Qianqian 03/22/2018
To get familiar with CALIPSO data
'''
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import os, calendar, sys, fnmatch,datetime
from pyhdf.SD import SD, SDC
import netCDF4
from matplotlib import colors as CS
import mpl_toolkits.basemap as bm
from scipy.ndimage.filters import gaussian_filter
import itertools
import glob
import sys

#=======================================================
#Two functions below are copied directly from Chamara
def Extract_Feature_Info(vfm_array,nlay):
    npro = nlay.size
    #print(type(vfm_array))
    feature_type    = np.full_like(vfm_array,-1)
    feature_type_qa = np.full_like(vfm_array,-1)
    ice_water_phase = np.full_like(vfm_array,-1)
    ice_water_phase_qa = np.full_like(vfm_array,-1)
    feature_subtype = np.full_like(vfm_array,-1)
    cloud_aerosol_psc_type_qa = np.full_like(vfm_array,-1)
    horizontal_averaging = np.full_like(vfm_array,-1)
    for i in range(npro):
        for l in range(nlay[i]):
            feature_type[i,l],feature_type_qa[i,l],\
            ice_water_phase[i,l],ice_water_phase_qa[i,l],\
            feature_subtype[i,l],cloud_aerosol_psc_type_qa[i,l],\
            horizontal_averaging[i,l]=vfm_feature_flags(vfm_array[i,l])

    return feature_type, feature_type_qa,\
           ice_water_phase, ice_water_phase_qa,\
           feature_subtype, cloud_aerosol_psc_type_qa,\
           horizontal_averaging
def vfm_feature_flags(val,verbose=True):
    #  this routine demonstrates how to read and extract values from a feature
    #  classification flag 16-bit integer value in CALIPSO Level 2 Vertical
    #  Feature Mask files
    #
    #  INPUT:
    #  val - the feature classification flag value to be decoded
    #
    #  OUTPUT:
    #  all information is printed into the IDL log window
    val_bit = np.binary_repr(val)
    feature_type              = int(val_bit[-3:],2)
    feature_type_qa           = int(val_bit[-5:-3],2)
    ice_water_phase           = int(val_bit[-7:-5],2)
    ice_water_phase_qa        = int(val_bit[-9:-7],2)
    feature_subtype           = int(val_bit[-12:-9],2)
    cloud_aerosol_psc_type_qa = int(val_bit[-13],2)
    horizontal_averaging      = int(val_bit[-16:-13],2)
    if verbose:
           if   feature_type ==0 :
               print("Feature Type : invalid (bad or missing data)")
           elif feature_type ==1 :
               print("Feature Type : clear air")
           elif feature_type ==2 :
               print("Feature Type : cloud")
               if feature_subtype ==  0 :
                   print( "Feature Subtype : low overcast, transparent")
               elif feature_subtype ==  1 :
                   print( "Feature Subtype : low overcast, opaque")
               elif feature_subtype ==  2 :
                   print( "Feature Subtype : transition stratocumulus")
               elif feature_subtype ==  3 :
                   print( "Feature Subtype : low, broken cumulus")
               elif feature_subtype ==  4 :
                   print( "Feature Subtype : altocumulus (transparent)")
               elif feature_subtype ==  5 :
                   print( "Feature Subtype : altostratus (opaque)")
               elif feature_subtype ==  6 :
                   print( "Feature Subtype : cirrus (transparent)")
               elif feature_subtype ==  7 :
                   print( "Feature Subtype : deep convective (opaque)")
               else :
                   print("*** error getting Feature Subtype")
           elif feature_type == 3 :
               print("Feature Type : aerosol")
               if feature_subtype ==  0 : print( "Feature Subtype : not determined")
               elif feature_subtype ==  1 : print( "Feature Subtype : clean marine")
               elif feature_subtype ==  2 : print( "Feature Subtype : dust")
               elif feature_subtype ==  3 : print( "Feature Subtype : polluted continental")
               elif feature_subtype ==  4 : print( "Feature Subtype : clean continental")
               elif feature_subtype ==  5 : print( "Feature Subtype : polluted dust")
               elif feature_subtype ==  6 : print( "Feature Subtype : smoke")
               elif feature_subtype ==  7 : print( "Feature Subtype : other")
               else : print("*** error getting Feature Subtype")
           elif feature_type ==4:
             print("Feature Type : stratospheric feature--PSC or stratospheric aerosol")
             if feature_subtype ==0 : print( "Feature Subtype : not determined")
             elif feature_subtype ==1 : print( "Feature Subtype : non-depolarizing PSC")
             elif feature_subtype ==2 : print( "Feature Subtype : depolarizing PSC")
             elif feature_subtype ==3 : print( "Feature Subtype : non-depolarizing aerosol")
             elif feature_subtype ==4 : print( "Feature Subtype : depolarizing aerosol")
             elif feature_subtype ==5 : print( "Feature Subtype : spare")
             elif feature_subtype ==6 : print( "Feature Subtype : spare")
             elif feature_subtype ==7 : print( "Feature Subtype : other")
             else : print("*** error getting Feature Subtype")
           elif feature_type == 5 : print("Feature Type : surface")
           elif feature_type == 6 : print("Feature Type : subsurface")
           elif feature_type == 7 : print("Feature Type : no signal (totally attenuated)")
           else : print("*** error getting Feature Type")

           if feature_type_qa == 0 : print("Feature Type QA : none")
           elif feature_type_qa ==1 : print("Feature Type QA : low")
           elif feature_type_qa ==2 : print("Feature Type QA : medium")
           elif feature_type_qa == 3 : print("Feature Type QA : high")
           else : print("*** error getting Feature Type QA")

           if ice_water_phase == 0 : print("Ice/Water Phase : unknown/not determined")
           elif ice_water_phase == 1 : print("Ice/Water Phase : ice")
           elif ice_water_phase == 2 : print("Ice/Water Phase : water")
           elif ice_water_phase == 3 : print("Ice/Water Phase : mixed phase")
           else : print("*** error getting Ice/Water Phase")
           
           if ice_water_phase_qa ==0 : print("Ice/Water Phase QA: none")
           elif ice_water_phase_qa == 1 : print("Ice/Water Phase QA: low")
           elif ice_water_phase_qa ==2 : print("Ice/Water Phase QA: medium")
           elif ice_water_phase_qa ==3 : print("Ice/Water Phase QA: high")
           else : print("*** error getting Ice/Water Phase QA")

           if (cloud_aerosol_psc_type_qa == 0):
               print("Cloud/Aerosol/PSC Type QA : not confident")
           else:
               print("Cloud/Aerosol/PSC Type QA : confident")

           if horizontal_averaging == 0:
               print("Horizontal averaging required for detection: not applicable")
           elif horizontal_averaging == 1 :
               print("Horizontal averaging required for detection: 1/3 km")
           elif horizontal_averaging ==2 :
               print("Horizontal averaging required for detection: 1 km")
           elif horizontal_averaging ==3 :
               print("Horizontal averaging required for detection: 5 km")
           elif horizontal_averaging ==4 :
               print("Horizontal averaging required for detection: 20 km")
           elif horizontal_averaging ==5 :
               print("Horizontal averaging required for detection: 80 km")
           else :
               print("*** error getting Horizontal averaging")
    return feature_type, feature_type_qa,\
           ice_water_phase, ice_water_phase_qa,\
           feature_subtype, cloud_aerosol_psc_type_qa,\
           horizontal_averaging

