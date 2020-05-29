#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Dec 12 10:06:46 2016
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

-A collection of classes and functions to read and analyze CALIPSO and CATS data.
-Used for my GRL_paper.
-22 Mar 2017: So far, only used level 2 "Layer" data. Going to add "profile" data.
 - some "profile" data added to see aca extinction profile
- 23 Apr 2017: dtcollect() changed to handle different seasons.
- 11 May 2017: Changed to handle calipso level 2 version 4 data.
03/22/2018:
    CALIOP data reading routings added to the end.
    Attenuated_Backscatter_profile added
 
"""

import scipy.io as sio
import numpy as np
import pickle, os
import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt
from pyhdf.SD import SD,SDC
from textwrap import wrap

class dsetCALCATS(object):
    def __init__(self,dflt=True): 
        
        if dflt:
            self.initdfltlayer()
        
            self.fstamp=[];self.figstamp=[]
            self.lon_mn=[];self.lon_mx=[];self.lat_mn=[];self.lat_mx=[]
            self.yrs=[];self.dnt=[];self.sat=[];self.loct=[]
            self.pro_data=False
    def read_CALIOPL2(self,filename,path=None):
        '''
        Read CALIOP lvele 2 layer data (CAL_LID_L2_05kmMLay-Standard-V4-10)
        '''
        self.filename=filename
        if path==None:
            hfL2 = SD(filename,SDC.READ)
        else:
            hfL2 = SD(path+filename,SDC.READ)
            self.path=path
        self.layer_IAB_532=hfL2.select('Integrated_Attenuated_Backscatter_532')[:][:]
        self.layer_DPR_532=hfL2.select('Integrated_Volume_Depolarization_Ratio')[:][:]
        self.layer_IAB_1064=hfL2.select('Integrated_Attenuated_Backscatter_1064')[:][:]
        self.Feature_Optical_Depth_532=hfL2.select('Feature_Optical_Depth_532')[:][:]
        self.Feature_Optical_Depth_1064=hfL2.select('Feature_Optical_Depth_1064')[:][:]
        self.layer_base=hfL2.select('Layer_Base_Altitude')[:][:]
        self.layer_top=hfL2.select('Layer_Top_Altitude')[:][:]
        self.horizontal_averaging=hfL2.select('Horizontal_Averaging')[:][:]  

        self.lat = hfL2.select('Latitude')[:][:]   
        self.lon = hfL2.select('Longitude')[:][:]
        self.vfm = hfL2.select('Feature_Classification_Flags')[:][:]
        self.nly = np.squeeze(hfL2.select('Number_Layers_Found')[:][:])
        self.layer_base=hfL2.select('Layer_Base_Altitude')[:][:]
        self.layer_top=hfL2.select('Layer_Top_Altitude')[:][:]

    def read_CALIPL1(self,filename,path=None):
        '''
        Read CALIOP level 1 data (CAL_LID_L1-Standard-V4-10)
        '''
        self.filename=filename
        if path==None:
            hfL1 = SD(filename,SDC.READ)
        else:
            hfL1 = SD(path+filename,SDC.READ)
            self.path=path          
        self.Attenuated_Backscatter_1064=hfL1.select('Attenuated_Backscatter_1064')[:][:]
        self.Attenuated_Backscatter_532=hfL1.select('Total_Attenuated_Backscatter_532')[:][:]
        self.time=hfL1.select('Profile_Time')[:][:]
        self.lat=hfL1.select('Latitude')[:][:]
        self.lon=hfL1.select('Longitude')[:][:]
        self.profile_id=hfL1.select('Profile_ID')[:][:]
        self.Bin_altitude289_578=np.flipud(np.linspace(-.5,8.3,289))#From the ATBD

    def initdfltlayer(self,):
        self.pnt_IAB_cld_532=[]
        self.pnt_DPR_cld_532=[]
        self.pnt_IAB_cld_1064=[]
        self.pnt_DPR_cld_1064=[]
        self.pnt_IAB_wcld_1064=[]
        self.pnt_DPR_wcld_1064=[]
        self.pnt_IAB_icld_1064=[]
        self.pnt_DPR_icld_1064=[]
        self.aer_base=[]
        self.wcld_top=[]
        self.icld_top=[]
        self.cba_top=[]
        self.cba_base=[]
        self.wcba_top=[]
        self.wcba_base=[]
        self.uphcba_top=[]
        self.uphcba_base=[]
        self.aac_base=[]
        self.aac_top=[]
        self.aawc_base=[]
        self.aawc_top=[]
        self.aaic_base=[]
        self.aaic_top=[]
        self.aauphc_base=[]
        self.aauphc_top=[]
        self.set_basemintop=[]
        self.set_wbasemintop=[]
        self.set_ibasemintop=[]
        self.set_ubasemintop=[]
        #location
        self.tot_lat=[]
        self.tot_lon=[]
                
        self.cba_lat=[]
        self.cba_lon=[]
        self.wcba_lat=[]
        self.wcba_lon=[]
        self.aac_lat=[]
        self.aac_lon=[]
        self.aawc_lat=[]
        self.aawc_lon=[]
        self.aaic_lat=[]
        self.aaic_lon=[]
        self.aauphc_lat=[]
        self.aauphc_lon=[]
        self.ld_lat=[]
        self.cld_lat=[]
        self.cld_lon=[]
        self.aer_lat=[]
        self.aer_lon=[]

        self.cldy_lon=[]
        self.cldy_lat=[]
        self.wcldy_lon=[]
        self.wcldy_lat=[]
        self.icldy_lon=[]
        self.icldy_lat=[]
        self.uphcldy_lat=[]
        self.uphcldy_lon=[]
        self.aerly_lon=[]
        self.aerly_lat=[]

    def initprodata(self,):
        self.pro_data=True
        self.norm_aca_ex_pro=[]
        self.norm_aca_thick=[]
    def initWyoming(self,):
        self.aer_cnt=[]
        self.cld_cnt=[]
        self.ttl_cnt=[]
        self.bin_alt=[]
        self.aawc_cnt=[]
        self.wcba_cnt=[]
        self.aer_ext=[]
    def appendprodata(self,totds,ds):
        totds.norm_aca_ex_pro=np.append(totds.norm_aca_ex_pro,ds.norm_aca_ex_pro)
        totds.norm_aca_thick=np.append(totds.norm_aca_thick,ds.norm_aca_thick)

    def appendset(self,totds,ds):
        totds.pnt_IAB_cld_532=np.append(totds.pnt_IAB_cld_532,ds.pnt_IAB_cld_532)
        totds.pnt_DPR_cld_532=np.append(totds.pnt_DPR_cld_532,ds.pnt_DPR_cld_532)
        totds.aer_base=np.append(totds.aer_base,ds.aer_base)
        totds.wcld_top=np.append(totds.wcld_top,ds.wcld_top)
        totds.icld_top=np.append(totds.icld_top,ds.icld_top)
        totds.cba_top=np.append(totds.cba_top,ds.cba_top)
        totds.cba_base=np.append(totds.cba_base,ds.cba_base)    
        totds.wcba_top=np.append(totds.wcba_top,ds.wcba_top)
        totds.wcba_base=np.append(totds.wcba_base,ds.wcba_base)
        totds.uphcba_top=np.append(totds.uphcba_top,ds.uphcba_top)
        totds.uphcba_base=np.append(totds.uphcba_base,ds.uphcba_base)
        totds.aac_base=np.append(totds.aac_base,ds.aac_base)
        totds.aac_top=np.append(totds.aac_top,ds.aac_top)
        totds.aawc_base=np.append(totds.aawc_base,ds.aawc_base)
        totds.aawc_top=np.append(totds.aawc_top,ds.aawc_top)
        totds.aaic_base=np.append(totds.aaic_base,ds.aaic_base)
        totds.aaic_top=np.append(totds.aaic_top,ds.aaic_top)
        totds.aauphc_base=np.append(totds.aauphc_base,ds.aauphc_base)
        totds.aauphc_top=np.append(totds.aauphc_top,ds.aauphc_top)
        totds.set_basemintop=np.append(totds.set_basemintop,ds.set_basemintop)
        totds.set_wbasemintop=np.append(totds.set_wbasemintop,ds.set_wbasemintop)
        totds.set_ibasemintop=np.append(totds.set_ibasemintop,ds.set_ibasemintop)
        totds.set_ubasemintop=np.append(totds.set_ubasemintop,ds.set_ubasemintop)
        #locations
        totds.tot_lat=np.append(totds.tot_lat,ds.tot_lat)
        totds.tot_lon=np.append(totds.tot_lon,ds.tot_lon)
        totds.cba_lat=np.append(totds.cba_lat,ds.cba_lat)
        totds.cba_lon=np.append(totds.cba_lon,ds.cba_lon)
        totds.wcba_lat=np.append(totds.wcba_lat,ds.wcba_lat)
        totds.wcba_lon=np.append(totds.wcba_lon,ds.wcba_lon)
        totds.aac_lat=np.append(totds.aac_lat,ds.aac_lat)
        totds.aac_lon=np.append(totds.aac_lon,ds.aac_lon)
        totds.aawc_lat=np.append(totds.aawc_lat,ds.aawc_lat)
        totds.aawc_lon=np.append(totds.aawc_lon,ds.aawc_lon)
        totds.aaic_lat=np.append(totds.aaic_lat,ds.aaic_lat)
        totds.aaic_lon=np.append(totds.aaic_lon,ds.aaic_lon)
        totds.aauphc_lat=np.append(totds.aauphc_lat,ds.aauphc_lat)
        totds.aauphc_lon=np.append(totds.aauphc_lon,ds.aauphc_lon)
        totds.cld_lat=np.append(totds.cld_lat,ds.cld_lat)
        totds.cld_lon=np.append(totds.cld_lon,ds.cld_lon)
        totds.aer_lat=np.append(totds.aer_lat,ds.aer_lat)
        totds.aer_lon=np.append(totds.aer_lon,ds.aer_lon)
        totds.cldy_lon=np.append(totds.cldy_lon,ds.cldy_lon)
        totds.cldy_lat=np.append(totds.cldy_lat,ds.cldy_lat)
        totds.wcldy_lon=np.append(totds.wcldy_lon,ds.wcldy_lon)
        totds.wcldy_lat=np.append(totds.wcldy_lat,ds.wcldy_lat)
        totds.icldy_lon=np.append(totds.icldy_lon,ds.icldy_lon)
        totds.icldy_lat=np.append(totds.icldy_lat,ds.icldy_lat)
        totds.uphcldy_lat=np.append(totds.uphcldy_lat,ds.uphcldy_lat)
        totds.uphcldy_lon=np.append(totds.uphcldy_lon,ds.uphcldy_lon)
        totds.aerly_lon=np.append(totds.aerly_lon,ds.aerly_lon)
        totds.aerly_lat=np.append(totds.aerly_lat,ds.aerly_lat)
    
        totds.fstamp=np.append(totds.fstamp,ds.fstamp)
        totds.lon_mn=ds.lon_mn
        totds.lon_mx=ds.lon_mx
        totds.lat_mn=ds.lat_mn
        totds.lat_mx=ds.lat_mx
        
    def dtcollect(self,totds,yrs,loct,dnt,sat,fileinfo=None,version=None):
        '''
        #Merging data ("Layer" data)
        loct: Location string
            sea: SE Atlantic region (for GRL paper)
            cNEA: central NE Atlantic dust study
        '''
        totds.yrs=yrs;totds.dnt=dnt;totds.sat=sat;totds.loct=loct
        if version==None:
            version=''
        else:
            version=version+'_'
        if fileinfo==None:
            st_date='_07_01to'
            ed_date='_10_30_'
            lats='lat-30.0_10.0_'
            lons='lon-20.0_20.0_'
        else:
            st_date=fileinfo['st_date']
            ed_date=fileinfo['ed_date']
            lats   =fileinfo['lats']
            lons   =fileinfo['lons']
        self.figstamp=str(yrs.min())+'to'+str(yrs.max())+st_date+ed_date+lats+lons
        for yr in range(np.size(yrs)):
            print('Reading year:'+str(yrs[yr]))
            if loct=='sea' or 'cNEA':
                with open('../'+sat+'/'+str(yrs[yr])+st_date+str(yrs[yr])+\
                ed_date+lats+lons+version+dnt+'.pkl', 'rb') as input:
                    ds = pickle.load(input)
#                    ds.lat_mn=-30.0;ds.lat_mx=10.0
#                    ds.lon_mn=-20.0;ds.lon_mx=20.0
            if loct=='glb':
                with open('../'+sat+'/'+str(yrs[yr])+'_07_01to'+str(yrs[yr])+\
                '_10_30_lat-60.0_60.0_lon-180.0_180.0_'+version+dnt+'.pkl', 'rb') as input:
                    ds = pickle.load(input)
#                    ds.lat_mn=-60.0;ds.lat_mx=60.0
#                    ds.lon_mn=-180.0;ds.lon_mx=180.0
            self.appendset(totds,ds)
            if self.pro_data:
                self.appendprodata(totds,ds)
        return totds
    
class setCal(object):
    def __init__(self,yrs,dnt,loct):
        self.dnt=dnt
        self.yrs=yrs
        self.loct=loct
        self.tot_lon=[];self.tot_lat=[];self.aac_lat=[];self.aac_lon=[];
        self.cld_lat=[];self.cld_lon=[];self.aer_lat=[];self.aer_lon=[];
        self.aac_base=[];self.cba_top=[];self.set_basemintop=[]
        self.wcld_lat=[];self.wcld_lon=[];self.icld_lat=[];self.icld_lon=[]
        self.pwcld_lat=[];self.pwcld_lon=[];
        self.cba_lat=[];self.cba_lon=[];self.aac_top=[];self.cba_base=[]
        self.cldy_lon=[];self.cldy_lat=[];self.wcldy_lon=[];self.wcldy_lat=[];
        self.icldy_lon=[];self.icldy_lat=[];self.aerly_lat=[];self.aerly_lon=[]
        
        self.aawc_lat=[];self.aawc_lon=[];self.aawc_base=[];self.wcba_top=[];
        self.set_wbasemintop=[];self.wcba_lat=[];self.wcba_lon=[];self.aawc_top=[];
        self.wcba_base=[]
        for yr in range(np.size(self.yrs)):
            print('Reading year:'+str(self.yrs[yr]))
            if self.loct=='sea':
                self.data=sio.loadmat(str(self.yrs[yr])+'_06_01to'+str(self.yrs[yr])+\
                     '_08_30_lat-30.0_10.0_lon-20.0_20.0_'+self.dnt+'.mat')
            elif self.loct=='glb':
                self.data=sio.loadmat(str(self.yrs[yr])+'_06_01to'+str(self.yrs[yr])+\
                     '_08_30_lat-60.0_60.0_lon-180.0_180.0_'+self.dnt+'.mat')
            self.tot_lon=np.append(self.tot_lon,self.data['tot_lon'].T)
            self.tot_lat=np.append(self.tot_lat,self.data['tot_lat'].T)
            self.aac_lat=np.append(self.aac_lat,self.data['aac_lat'].T)
            self.aac_lon=np.append(self.aac_lon,self.data['aac_lon'].T)
            self.cba_lat=np.append(self.cba_lat,self.data['cba_lat'].T)
            self.cba_lon=np.append(self.cba_lon,self.data['cba_lon'].T)
            self.cld_lat=np.append(self.cld_lat,self.data['cld_lat'].T)
            self.cld_lon=np.append(self.cld_lon,self.data['cld_lon'].T)
            self.aer_lat=np.append(self.aer_lat,self.data['aer_lat'].T)
            self.aer_lon=np.append(self.aer_lon,self.data['aer_lon'].T)
            self.wcld_lat=np.append(self.wcld_lat,self.data['wcld_lat'].T)
            self.wcld_lon=np.append(self.wcld_lon,self.data['wcld_lon'].T)
            self.icld_lat=np.append(self.icld_lat,self.data['icld_lat'].T)
            self.icld_lon=np.append(self.icld_lon,self.data['icld_lon'].T)
            self.pwcld_lat=np.append(self.pwcld_lat,self.data['pwcld_lat'].T)
            self.pwcld_lon=np.append(self.pwcld_lon,self.data['pwcld_lon'].T)
            self.aac_base=np.append(self.aac_base,self.data['aac_base'].T)
            self.aac_top=np.append(self.aac_top,self.data['aac_top'].T)
            self.cba_top=np.append(self.cba_top,self.data['cba_top'].T)
            self.cba_base=np.append(self.cba_base,self.data['cba_base'].T)
            self.set_basemintop=np.append(self.set_basemintop,self.data['set_basemintop'].T)
            
            self.aawc_lat=np.append(self.aawc_lat,self.data['aawc_lat'].T)
            self.aawc_lon=np.append(self.aawc_lon,self.data['aawc_lon'].T)
            self.wcba_lat=np.append(self.wcba_lat,self.data['wcba_lat'].T)
            self.wcba_lon=np.append(self.wcba_lon,self.data['wcba_lon'].T)
            self.aawc_base=np.append(self.aawc_base,self.data['aawc_base'].T)
            self.aawc_top=np.append(self.aawc_top,self.data['aawc_top'].T)
            self.wcba_top=np.append(self.wcba_top,self.data['wcba_top'].T)
            self.wcba_base=np.append(self.wcba_base,self.data['wcba_base'].T)
            self.set_wbasemintop=np.append(self.set_wbasemintop,self.data['set_wbasemintop'].T)
            
            self.cldy_lon=np.append(self.cldy_lon,self.data['cldy_lon'])
            self.cldy_lat=np.append(self.cldy_lat,self.data['cldy_lat'])
            self.wcldy_lon=np.append(self.wcldy_lon,self.data['wcldy_lon'])
            self.wcldy_lat=np.append(self.wcldy_lat,self.data['wcldy_lat'])
            self.icldy_lon=np.append(self.icldy_lon,self.data['icldy_lon'])
            self.icldy_lat=np.append(self.icldy_lat,self.data['icldy_lat'])
            self.aerly_lon=np.append(self.aerly_lon,self.data['aerly_lon'])
            self.aerly_lat=np.append(self.aerly_lat,self.data['aerly_lat'])
        self.lat_mn=self.data['lat_mn']
        self.lon_mn=self.data['lon_mn']
        self.lat_mx=self.data['lat_mx']
        self.lon_mx=self.data['lon_mx']
        self.fstamp=self.data['fstamp']
        self.figstamp=self.data['figstamp']
        self.figstamp=str(self.figstamp)

class setCat(object):
    def __init__(self,yrs,dnt,loct):
        self.dnt=dnt
        self.yrs=yrs
        self.tot_lon=[];self.tot_lat=[];self.aac_lat=[];self.aac_lon=[];
        self.cld_lat=[];self.cld_lon=[];self.aer_lat=[];self.aer_lon=[];
        self.aac_base=[];self.cba_top=[];self.set_basemintop=[]
        self.wcld_lat=[];self.wcld_lon=[];self.icld_lat=[];self.icld_lon=[]
        self.pwcld_lat=[];self.pwcld_lon=[];self.paer_lat=[];self.paer_lon=[]
        self.cba_lat=[];self.cba_lon=[];self.aac_top=[];self.cba_base=[]
        self.cldy_lon=[];self.cldy_lat=[];self.wcldy_lon=[];self.wcldy_lat=[];
        self.icldy_lon=[];self.icldy_lat=[];self.aerly_lat=[];self.aerly_lon=[]
        
        self.aawc_lat=[];self.aawc_lon=[];self.aawc_base=[];self.wcba_top=[];
        self.set_wbasemintop=[];self.wcba_lat=[];self.wcba_lon=[];self.aawc_top=[];
        self.wcba_base=[]
        for yr in range(np.size(self.yrs)):
            print('Reading year:'+str(self.yrs[yr]))
            if loct=='sea':
                self.data=sio.loadmat(str(self.yrs[yr])+'_06_01to'+str(self.yrs[yr])+\
                     '_08_30_lat-30.0_10.0_lon-20.0_20.0_'+self.dnt+'.mat')
            elif loct=='glb':
                self.data=sio.loadmat(str(self.yrs[yr])+'_06_01to'+str(self.yrs[yr])+\
                     '_08_30_lat-60.0_60.0_lon-180.0_180.0_'+self.dnt+'.mat')
            self.tot_lon=np.append(self.tot_lon,self.data['tot_lon'].T)
            self.tot_lat=np.append(self.tot_lat,self.data['tot_lat'].T)
            self.aac_lat=np.append(self.aac_lat,self.data['aac_lat'].T)
            self.aac_lon=np.append(self.aac_lon,self.data['aac_lon'].T)
            self.cba_lat=np.append(self.cba_lat,self.data['cba_lat'].T)
            self.cba_lon=np.append(self.cba_lon,self.data['cba_lon'].T)
#            self.cld_lat=np.append(self.cld_lat,self.data['cld_lat'].T)
#            self.cld_lon=np.append(self.cld_lon,self.data['cld_lon'].T)
#            self.aer_lat=np.append(self.aer_lat,self.data['aer_lat'].T)
#            self.aer_lon=np.append(self.aer_lon,self.data['aer_lon'].T)
#            self.wcld_lat=np.append(self.wcld_lat,self.data['wcld_lat'].T)
#            self.wcld_lon=np.append(self.wcld_lon,self.data['wcld_lon'].T)
#            self.icld_lat=np.append(self.icld_lat,self.data['icld_lat'].T)
#            self.icld_lon=np.append(self.icld_lon,self.data['icld_lon'].T)
#            self.pwcld_lat=np.append(self.pwcld_lat,self.data['pwcld_lat'].T)
#            self.pwcld_lon=np.append(self.pwcld_lon,self.data['pwcld_lon'].T)
#            self.paer_lat=np.append(self.paer_lat,self.data['paer_lat'].T)
#            self.paer_lon=np.append(self.paer_lon,self.data['paer_lon'].T)
            self.aac_base=np.append(self.aac_base,self.data['aac_base'].T)
            self.aac_top=np.append(self.aac_top,self.data['aac_top'].T)
            self.cba_top=np.append(self.cba_top,self.data['cba_top'].T)
            self.cba_base=np.append(self.cba_base,self.data['cba_base'].T)
            self.set_basemintop=np.append(self.set_basemintop,self.data['set_basemintop'].T)
            
            self.aawc_lat=np.append(self.aawc_lat,self.data['aawc_lat'].T)
            self.aawc_lon=np.append(self.aawc_lon,self.data['aawc_lon'].T)
            self.wcba_lat=np.append(self.wcba_lat,self.data['wcba_lat'].T)
            self.wcba_lon=np.append(self.wcba_lon,self.data['wcba_lon'].T)
            self.aawc_base=np.append(self.aawc_base,self.data['aawc_base'].T)
            self.aawc_top=np.append(self.aawc_top,self.data['aawc_top'].T)
            self.wcba_top=np.append(self.wcba_top,self.data['wcba_top'].T)
            self.wcba_base=np.append(self.wcba_base,self.data['wcba_base'].T)
            self.set_wbasemintop=np.append(self.set_wbasemintop,self.data['set_wbasemintop'].T)
            
            self.cldy_lon=np.append(self.cldy_lon,self.data['cldy_lon'])
            self.cldy_lat=np.append(self.cldy_lat,self.data['cldy_lat'])
            self.wcldy_lon=np.append(self.wcldy_lon,self.data['wcldy_lon'])
            self.wcldy_lat=np.append(self.wcldy_lat,self.data['wcldy_lat'])
            self.icldy_lon=np.append(self.icldy_lon,self.data['icldy_lon'])
            self.icldy_lat=np.append(self.icldy_lat,self.data['icldy_lat'])
            self.aerly_lon=np.append(self.aerly_lon,self.data['aerly_lon'])
            self.aerly_lat=np.append(self.aerly_lat,self.data['aerly_lat'])
        self.lat_mn=self.data['lat_mn']
        self.lon_mn=self.data['lon_mn']
        self.lat_mx=self.data['lat_mx']
        self.lon_mx=self.data['lon_mx']
        self.fstamp=self.data['fstamp']
        self.figstamp=self.data['figstamp']
        self.figstamp=str(self.figstamp)

def rmvxtrms(dM):
    if dM.size!=0:
        q3=np.percentile(dM,75);q1=np.percentile(dM,25)
        dMmin=q1-1.5*(q3-q1);dMmax=q3+1.5*(q3-q1)
        dM=dM[dM>dMmin];dM=dM[dM<dMmax]
#    else:
#        print('Zero size dM!!!')
    return dM

def drawmap(mapproj,latlines,lonlines):
    mapproj.drawcoastlines()
    mapproj.drawparallels(latlines, labels=[1,0,0,0])
    mapproj.drawmeridians(lonlines, labels=[0,0,0,1])

def draw1map(ax,seacatN,dlat,dlon,fcon):
    #plotting one row gridded figures in SEA region. 
    cmap = plt.cm.jet;
    lat_bnds=np.arange(seacatN.lat_mn,seacatN.lat_mx+dlat,dlat)
    lon_bnds=np.arange(seacatN.lon_mn,seacatN.lon_mx+dlon,dlon)
    lat_grids=(lat_bnds[0:-1]+lat_bnds[1:])/2.0
    lon_grids=(lon_bnds[0:-1]+lon_bnds[1:])/2.0
    for i in range(0,np.shape(ax)[0]):
        mapproj = bm.Basemap(ax=ax[i],projection='cyl',llcrnrlat= seacatN.lat_mn, \
                                 llcrnrlon= seacatN.lon_mn,urcrnrlat= seacatN.lat_mx, urcrnrlon= seacatN.lon_mx)
        [lonall, latall] = np.meshgrid(lon_grids, lat_grids)
        lonproj, latproj = mapproj(lonall, latall)
        latlines = np.arange(seacatN.lat_mn,seacatN.lat_mx,10.0)
        lonlines = np.arange(seacatN.lon_mn,seacatN.lon_mx,10.0)
        
        mapproj.drawcoastlines()
        mapproj.drawparallels(latlines, labels=[1,0,0,0])
        mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
        ax[i].patch.set_color('gray')
#            ax[i,j].patch.set_fill(True)
            
        if fcon=='y':
            mapproj.fillcontinents(color='0.8',alpha=0.5) 
    return cmap,lat_bnds,lon_bnds,lonall,latall,lonproj,latproj, lat_grids,\
            lon_grids
def draw6maps(ax,seacatN,dlat,dlon,fcon):
    #plotting 3*3 gridded figures in SEA region. 
    cmap = plt.cm.jet;
    lat_bnds=np.arange(seacatN.lat_mn,seacatN.lat_mx+dlat,dlat)
    lon_bnds=np.arange(seacatN.lon_mn,seacatN.lon_mx+dlon,dlon)
    lat_grids=(lat_bnds[0:-1]+lat_bnds[1:])/2.0
    lon_grids=(lon_bnds[0:-1]+lon_bnds[1:])/2.0
    for i in range(0,np.shape(ax)[0]):
        for j in range(0,np.shape(ax)[1]):
            mapproj = bm.Basemap(ax=ax[i,j],projection='cyl',llcrnrlat= seacatN.lat_mn, \
                                 llcrnrlon= seacatN.lon_mn,urcrnrlat= seacatN.lat_mx, urcrnrlon= seacatN.lon_mx)
            [lonall, latall] = np.meshgrid(lon_grids, lat_grids)
            lonproj, latproj = mapproj(lonall, latall)
            latlines = np.arange(seacatN.lat_mn,seacatN.lat_mx,10.0)
            lonlines = np.arange(seacatN.lon_mn,seacatN.lon_mx,10.0)
        
            mapproj.drawcoastlines()
            mapproj.drawparallels(latlines, labels=[1,0,0,0])
            mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
            ax[i,j].patch.set_color('gray')
#            ax[i,j].patch.set_fill(True)
            
            if fcon=='y':
                mapproj.fillcontinents(color='0.8',alpha=0.5) 
            
    return cmap,lat_bnds,lon_bnds,lonall,latall,lonproj,latproj, lat_grids,\
            lon_grids

def mask_non_awca(clp,lon_bnds,lat_bnds,dataMesh,stdM,op):
    #Ignoring the region based on awca_fq value
    val=0.80
    Hcldy, xedges, yedges = np.histogram2d(np.squeeze(clp.cldy_lon),\
                            np.squeeze(clp.cldy_lat), bins=(lon_bnds,lat_bnds))        
    Haawc,xedges,yedges=np.histogram2d(np.squeeze(clp.aawc_lon),\
                            np.squeeze(clp.aawc_lat), bins=(lon_bnds,lat_bnds))
    awca_fq=(Haawc/Hcldy*100).T
#    dataMesh[awca_fq<np.mean(awca_fq)*val]=np.nan
#    stdM[awca_fq<np.mean(awca_fq)*val]=np.nan
    rgn=awca_fq<np.mean(awca_fq)*val
    mask=np.empty_like(dataMesh)
    mask[rgn]=0
    mask[~rgn]=1

    sio.savemat('mask'+clp.sat+clp.dnt+'.mat',{'awca_fq':awca_fq,'mask':mask})
    
    if op==1:
        outf=np.sum(Haawc[rgn])/np.sum(Haawc)
        print(clp.sat+' - '+clp.dnt+': '+str(outf))
    
    return mask

def fig15_mer_avg(clp,lon_grids,lon_bnds,ax,legt,prob_thresh=None,lat_range=None):
    if prob_thresh==None:
        prob_thresh==0.250
    if lat_range==None:
        lat_range={'bot':-30.0,'top':10.0}
    data=np.zeros((4,np.size(lon_grids)))
    std=np.zeros((3,np.size(lon_grids)))
    for i in range(0,np.size(lon_grids)):
        dM=clp.aawc_base[(lon_bnds[i]<=clp.aawc_lon)*(lon_bnds[i+1]>clp.aawc_lon)*\
                    (lat_range['bot']<=clp.aawc_lat)*(lat_range['top']>clp.aawc_lat)]
        dM_base=dM;
        dM=rmvxtrms(dM)
        data[0,i]=np.mean(dM)
        std[0,i]=np.std(dM)
        dM=clp.aawc_top[(lon_bnds[i]<=clp.aawc_lon)*(lon_bnds[i+1]>clp.aawc_lon)*\
                    (lat_range['bot']<=clp.aawc_lat)*(lat_range['top']>clp.aawc_lat)]
        dM=rmvxtrms(dM)
        data[1,i]=np.mean(dM)
        std[1,i]=np.std(dM)
        dM=clp.wcba_top[(lon_bnds[i]<=clp.wcba_lon)*(lon_bnds[i+1]>clp.wcba_lon)*\
                    (lat_range['bot']<=clp.wcba_lat)*(lat_range['top']>clp.wcba_lat)]   
        dM_top=dM
        dM=rmvxtrms(dM)
        data[2,i]=np.mean(dM)
        std[2,i]=np.std(dM)
        
        AB2CT=dM_base-dM_top
        AB2CT=rmvxtrms(AB2CT)
        data[3,i]=np.array(AB2CT[AB2CT<prob_thresh].size,dtype=float)/AB2CT.size
      
    ax.plot(lon_grids,data[0,:],'r.-')
    ax.fill_between(lon_grids,data[0,:]-std[0,:],data[0,:]+std[0,:],\
                    facecolor='#ff9999',linewidth=0.0)
    ax.errorbar(lon_grids+0.1,data[1,:],yerr=std[1,:],fmt='r.--')
    ax.plot(lon_grids,data[2,:],'b.-')
    ax.fill_between(lon_grids,data[2,:]-std[2,:],data[2,:]+std[2,:],alpha=0.5,\
                    facecolor='#809fff',linewidth=0.0)
    if legt=='y':
        leg=ax.legend(['Aerosol base','Aerosol top','Cloud top'],prop={'size':10},loc='upper left')
        leg.get_frame().set_alpha(0.5)
    ax.set_ylim([0,7.0])
    ax.set_ylabel('Altitude(km)')
    
    Hcldy,edges=np.histogram(np.squeeze(clp.cldy_lon),bins=lon_bnds)
    Haawc,edges=np.histogram(np.squeeze(clp.aawc_lon),bins=lon_bnds)
    tH,edges=np.histogram(np.squeeze(clp.tot_lon),bins=lon_bnds)
    Hcldy=np.array(Hcldy,dtype=float)
    ax15_2 = ax.twinx()
#    x=lon_grids;y=Haawc/Hcldy
    x=lon_grids;y=data[3,:]
    x = np.ravel(zip(x,x+2));y = np.ravel(zip(y,y))
    ax15_2.plot(x,y,'k-', linewidth=1.5);
    ax15_2.set_ylim([0,1.0])
    if legt=='y':
        ax15_2.legend([r'Fraction of $AB2CT<360m$'],prop={'size':10},loc='upper right')
    ax15_2.set_ylabel(r'Fraction of $AB2CT<360m$',fontsize=9)
#    ax.set_xlim([-20.0,20.0])
    
    print(clp.sat+' - '+clp.dnt+' aerosol top: '+str(np.mean(data[1,:]))+' var. '+str(np.mean(std[1,:])))
    print(clp.sat+' - '+clp.dnt+' aerosol base: '+str(np.mean(data[0,:]))+' var. '+str(np.mean(std[0,:])))
    print(clp.sat+' - '+clp.dnt+' cloud top: '+str(np.mean(data[2,:]))+' var. '+str(np.mean(std[2,:])))
    
def fig3_aawc_fq(clp,lon_bnds,lat_bnds,lonproj,latproj,cmap,fig,ax0,ax1,ttl0,ttl1,txtclr=None,ignan=False):
    #ignan=True : ignore Nan values
    tH,xedges,yedges=np.histogram2d(np.squeeze(clp.tot_lon),\
                            np.squeeze(clp.tot_lat), bins=(lon_bnds,lat_bnds))
    Hcldy, xedges, yedges = np.histogram2d(np.squeeze(clp.cldy_lon),\
                            np.squeeze(clp.cldy_lat), bins=(lon_bnds,lat_bnds))        
    Haawc,xedges,yedges=np.histogram2d(np.squeeze(clp.aawc_lon),\
                            np.squeeze(clp.aawc_lat), bins=(lon_bnds,lat_bnds))
    ax1.set_title(ttl1,fontsize=10)
    v = np.linspace(0.0,100,num=50)
    ctf = ax1.contourf(lonproj, latproj,(Haawc/Hcldy*100).T,v,cmap=cmap)
    fig.colorbar(ctf,ax=ax1,orientation='horizontal',ticks=np.arange(0,101,20))
    ax0.set_title(ttl0,fontsize=10)
    v = np.linspace(0.0,100,num=50)
    ctf = ax0.contourf(lonproj, latproj,(Hcldy/tH*100).T,v,cmap=cmap)
    fig.colorbar(ctf,ax=ax0,orientation='horizontal',ticks=np.arange(0,101,20))
    
    if ignan:
        print('Domain avg CF: '+ str(np.nanmean(Hcldy/tH)))
        print('Domain avg ACA_F: '+str(np.nanmean(Haawc/Hcldy)))
        if txtclr==None:
            ax1.text(0.2,0.9,'Avg ACA_F: %0.2f '%np.nanmean(Haawc/Hcldy),\
             transform = ax1.transAxes,color='w')
        else:
            ax1.text(0.2,0.9,'Avg ACA_F: %0.2f '%np.nanmean(Haawc/Hcldy),\
             transform = ax1.transAxes,color=txtclr)
    else:
        print('Domain avg CF: '+ str(np.mean(Hcldy/tH)))
        print('Domain avg ACA_F: '+str(np.mean(Haawc/Hcldy)))
        if txtclr==None:
            ax1.text(0.2,0.9,'Avg ACA_F: %0.2f '%np.mean(Haawc/Hcldy),\
             transform = ax1.transAxes,color='w')
        else:
            ax1.text(0.2,0.9,'Avg ACA_F: %0.2f '%np.mean(Haawc/Hcldy),\
             transform = ax1.transAxes,color=txtclr)

'''
=================================================================================================
CALIOP data reading routings.
=================================================================================================
'''
def Extract_Feature_Info(vfm_array,nlay):
    npro = nlay.size
#    print(type(vfm_array))
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

def vfm_feature_flags(val,verbose=False):
    '''
  this routine demonstrates how to read and extract values from a feature
  classification flag 16-bit integer value in CALIPSO Level 2 Vertical
  Feature Mask files

  INPUT:
  val - the feature classification flag value to be decoded

  OUTPUT:
  all information is printed into the IDL log window
   '''
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
class Aeros(object):
    d_loc='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmALay/'
    fnd=0
    def getPat(self,dy,dm,date_char,version=None,dn=None):
#        fn_pattern='CAL_LID_L2_05kmALay-Standard-V4-10.'+date_char+'**'+dnt+'*.hdf'#-10-29T08-40-13ZN.hdf'
        if not(dn==None):
            dnt=dn
        if version==None:
            if dy>2012 and dy<2017:
                fn_pattern='CAL_LID_L2_05kmALay-Prov-V3-30.'+date_char+'**'+dnt+'*.hdf'#-10-29T08-40-13ZN.hdf'
            elif dy==2011 or dy==2012:
                if dy==2011 and dm<11:
                    fn_pattern='CAL_LID_L2_05kmALay-Prov-V3-01.'+date_char+'**'+dnt+'*.hdf'
                else:
                    fn_pattern='CAL_LID_L2_05kmALay-Prov-V3-02.'+date_char+'**'+dnt+'*.hdf'
    
            elif dy<2011 and dy>2005:
                fn_pattern='CAL_LID_L2_05kmALay-Prov-V3-01.'+date_char+'**'+dnt+'*.hdf'
        elif version=='V4':
            fn_pattern='CAL_LID_L2_05kmALay-Standard-V4-10.'+date_char+'**'+dnt+'*.hdf'
        else:
            print('Invalid version!!')
        return fn_pattern
        
class Cloud(object):
    d_loc='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/'
    fnd=0
    def getPat(self,dy,dm,date_char,version=None,dn=None):
#        fn_pattern='CAL_LID_L2_05kmCLay-Standard-V4-10.'+date_char+'**'+dnt+'*.hdf'#-10-29T08-40-13ZN.hdf'
        if not(dn==None):
            dnt=dn
        if version==None:
            if dy>2012 and dy<2017:
                fn_pattern='CAL_LID_L2_05kmCLay-Prov-V3-30.'+date_char+'**'+dnt+'*.hdf'#-10-29T08-40-13ZN.hdf'
            elif dy==2011 or dy==2012:
                if dy==2011 and dm<11:
                    fn_pattern='CAL_LID_L2_05kmCLay-Prov-V3-01.'+date_char+'**'+dnt+'*.hdf'
                else:
                    fn_pattern='CAL_LID_L2_05kmCLay-Prov-V3-02.'+date_char+'**'+dnt+'*.hdf'
            elif dy<2011 and dy>2005:
                fn_pattern='CAL_LID_L2_05kmCLay-Prov-V3-01.'+date_char+'**'+dnt+'*.hdf'
        elif version=='V4':
            fn_pattern='CAL_LID_L2_05kmCLay-Standard-V4-10.'+date_char+'**'+dnt+'*.hdf'
        else:
            print('Invalid version !!')
        return fn_pattern

def findCld(cdata_path,f):
    #dependencies: Extract_Feature_Info
    fnd=1;
    print(f)
    hf = SD(cdata_path+f,SDC.READ)
    lat = hf.select('Latitude')[:][:]   
    lon = hf.select('Longitude')[:][:]
    vfm = hf.select('Feature_Classification_Flags')[:][:]
    nly = hf.select('Number_Layers_Found')[:][:]
    layer_base=hf.select('Layer_Base_Altitude')[:][:]
    layer_top=hf.select('Layer_Top_Altitude')[:][:]
    nly=np.squeeze(nly)
                
    feature_type, feature_type_qa,\
    ice_water_phase, ice_water_phase_qa,\
    feature_subtype, cloud_aerosol_psc_type_qa,\
    horizontal_averaging=Extract_Feature_Info(vfm,nly)

    ix1=lat[:,0]>=lat_mn;ix2=lat[:,0]<=lat_mx
    lat_ix=ix1*ix2
    ix1=lon[:,0]>=lon_mn;ix2=lon[:,0]<=lon_mx
    lon_ix=ix1*ix2
    loc_ix=np.asarray([lat_ix*lon_ix]).T
    loc_ix=np.tile(loc_ix,vfm.shape[1])
    
    layer_IAB_532=hf.select('Integrated_Attenuated_Backscatter_532')[:][:]
    layer_DPR_532=hf.select('Integrated_Volume_Depolarization_Ratio')[:][:]
    layer_IAB_1064=hf.select('Integrated_Attenuated_Backscatter_1064')[:][:]
    opacity_flag=hf.select('Opacity_Flag')[:][:]      
    horizontal_averaging=hf.select('Horizontal_Averaging')[:][:]      
    
    cld_ix=feature_type==2#cloud
    aer_ix=feature_type==3#aerosols
    wcld_ix=(ice_water_phase==2)+(ice_water_phase==0)# water cloud
    icld_ix=(ice_water_phase==1) + (ice_water_phase==3)#ice cloud
    uphase_ix=ice_water_phase==0#unknown phase
    opq_ix=opacity_flag==1#opaque clouds (Try to use Opacity_Flag)
    op_wcld=wcld_ix*opq_ix#*opq_ix#opaque water cloud     
    hoz_av80=horizontal_averaging==5# 80km horizontal averaging        
                              
#                aac=np.zeros(feature_type.shape,dtype=bool)#aerosols above cloud
#                cba=np.zeros(feature_type.shape,dtype=bool)#cloud below aerosols
#             
    return fnd,layer_base, lat,lon,loc_ix,cld_ix,aer_ix,wcld_ix,\
        icld_ix,uphase_ix,op_wcld,layer_DPR_532,layer_IAB_532,\
        layer_IAB_1064,feature_type,layer_top,hoz_av80
        
def findAC(cdata_path,f,adata_path,g):
    #dependencies: Extract_Feature_Info, findCld
    cfnd,c_layer_base, c_lat,c_lon,c_loc_ix,cld_ix,_,wcld_ix,icld_ix,uphase_ix,\
    op_wcld,c_layer_DPR_532,c_layer_IAB_532,c_layer_IAB_1064,\
    c_feature_type,c_layer_top,c_hoz_av80\
    =findCld(cdata_path,f)
    afnd=1;
    print(g)
    hf = SD(adata_path+g,SDC.READ)
    a_lat = hf.select('Latitude')[:][:]   
    a_lon = hf.select('Longitude')[:][:]
    vfm = hf.select('Feature_Classification_Flags')[:][:]
    nly = hf.select('Number_Layers_Found')[:][:]
    a_layer_base=hf.select('Layer_Base_Altitude')[:][:]
    a_layer_top=hf.select('Layer_Top_Altitude')[:][:]
    nly=np.squeeze(nly)
    
    a_feature_type, a_feature_type_qa,\
    _, _,a_feature_sub_type, _,_=Extract_Feature_Info(vfm,nly)
    
    ix1=a_lat[:,0]>=lat_mn;ix2=a_lat[:,0]<=lat_mx
    lat_ix=ix1*ix2
    ix1=a_lon[:,0]>=lon_mn;ix2=a_lon[:,0]<=lon_mx
    lon_ix=ix1*ix2
    a_loc_ix=np.asarray([lat_ix*lon_ix]).T
    a_loc_ix=np.tile(a_loc_ix,vfm.shape[1])
    a_layer_IAB_532=hf.select('Integrated_Attenuated_Backscatter_532')[:][:]
    a_layer_DPR_532=hf.select('Integrated_Volume_Depolarization_Ratio')[:][:]
    a_layer_IAB_1064=hf.select('Integrated_Attenuated_Backscatter_1064')[:][:]
    a_horizontal_averaging=hf.select('Horizontal_Averaging')[:][:]  

    aer_ix=a_feature_type==3#aerosols
    smoke_ix=(a_feature_sub_type==6)+(a_feature_sub_type==5)# smoke(6) + polluted dust(5)
    aer_ix=aer_ix*smoke_ix
    a_hoz_av80=a_horizontal_averaging==5# 80km horizontal averaging 
#    cldy_ix=wcld_ix + icld_ix 
    if np.size(np.where(a_lat!=c_lat))!=0 or np.size(np.where(a_lon!=c_lon))!=0:
        print('Aerosol and cloud data files have different locations!!!')
    aac=np.zeros(a_feature_type.shape,dtype=bool)#aerosols above cloud
    aawc=np.zeros(a_feature_type.shape,dtype=bool)#aerosols above wcloud
    aaic=np.zeros(a_feature_type.shape,dtype=bool)#aerosols above icloud
    aauphc=np.zeros(a_feature_type.shape,dtype=bool)#aerosols above uphcloud
    cba=np.zeros(c_feature_type.shape,dtype=bool)#cloud below aerosols
    wcba=np.zeros(c_feature_type.shape,dtype=bool)#wcloud below aerosols
    icba=np.zeros(c_feature_type.shape,dtype=bool)#icloud below aerosols
    uphcba=np.zeros(c_feature_type.shape,dtype=bool)#uphcloud below aerosols
    
    cldy_ix=np.zeros(c_feature_type.shape,dtype=bool)#cloud
    aerbs_ix=np.zeros(a_feature_type.shape,dtype=bool)#aerosol
    wcldy_ix=np.zeros(c_feature_type.shape,dtype=bool)#water cloud
    icldy_ix=np.zeros(c_feature_type.shape,dtype=bool)#ice cloud
    uphcldy_ix=np.zeros(c_feature_type.shape,dtype=bool)#unknown phase clouds
    basemintop=[];wbasemintop=[];ibasemintop=[];ubasemintop=[]
    for i in range(0,a_feature_type.shape[0]):
        #selecting aerly
        if np.sum(aer_ix[i,:]*1)>0:
            botaer_ix=np.max(np.where(aer_ix[i,:]==True))
            aerbs_ix[i,botaer_ix]=True
        #selecting cldy
        if np.sum(cld_ix[i,:]*1)>0:
            topcld_ix=np.min(np.where(cld_ix[i,:]==True))
            cldy_ix[i,topcld_ix]=True
            #selecting wcldy
            if np.sum(wcld_ix[i,:]*1)>0:
                topwcld_ix=np.min(np.where(wcld_ix[i,:]==True))
                if topwcld_ix==np.min(np.where(cld_ix[i,:]==True)):
                        wcldy_ix[i,topwcld_ix]=True
            #selecting icldy   
            if np.sum(icld_ix[i,:]*1)>0:
                topicld_ix=np.min(np.where(icld_ix[i,:]==True))
                if topicld_ix==np.min(np.where(cld_ix[i,:]==True)):
                        icldy_ix[i,topicld_ix]=True
            #selecting unknown phase clouds   
            if np.sum(uphase_ix[i,:]*1)>0:
                #found "unknown phase cloud" layers without "cloud" 
                topupcld_ix=np.min(np.where(uphase_ix[i,:]==True))
                if topupcld_ix==np.min(np.where(cld_ix[i,:]==True)):
                    uphcldy_ix[i,topupcld_ix]=True
        #selecting aca       
        if (np.sum(cldy_ix[i,:]*1)>0) and \
        (np.sum(aer_ix[i,:]*1)>0):#having both aerosols and water cld
#                bmin=a_layer_base[i,aer_ix[i,:]].min()
                bs=a_layer_base[i,aer_ix[i,:]]
                cmax=c_layer_top[i,cldy_ix[i,:]].max()
                
                if np.sum((bs>cmax)*1)>0:
#                if bmin > cmax:#aerosol above op water cloud
                    bmin=bs[bs>cmax].min()
                    bmin_ix=np.where(a_layer_base[i,:]==bmin)
                    ctop_ix=np.where(c_layer_top[i,:]==cmax)
                    topaer=aer_ix[i,:]*0;topaer[np.max(bmin_ix)]=1;
                    blcld=cldy_ix[i,:]*0;blcld[np.min(ctop_ix)]=1;
                    aac[i,:]=topaer*a_loc_ix[i,:]
                    cba[i,:]=blcld*c_loc_ix[i,:]
#                    if bmin-cmax<0:
#                        print("MALAMAGULAY>>>>>>>>>>..")
#                    if np.sum(cba[i,:]*1)>0:
#                        basemintop=np.append(basemintop,bmin-cmax)
        #selecting awca
        if (np.sum(wcldy_ix[i,:]*1)>0) and \
        (np.sum(aer_ix[i,:]*1)>0):#having both aerosols and water wcld
#                bmin=a_layer_base[i,aer_ix[i,:]].min()
                bs=a_layer_base[i,aer_ix[i,:]]
                cmax=c_layer_top[i,wcldy_ix[i,:]]
                
                if np.sum((bs>cmax)*1)>0:#if any aerosol above cloud
#                if bmin > cmax:#aerosol above op water cloud
                    bmin=bs[bs>cmax].min()
                    bmin_ix=np.where(a_layer_base[i,:]==bmin)
                    ctop_ix=np.where(c_layer_top[i,:]==cmax)
                    topaer=aer_ix[i,:]*0;topaer[np.max(bmin_ix)]=1;
                    blwcld=wcldy_ix[i,:]*0;blwcld[np.min(ctop_ix)]=1;
                    aawc[i,:]=topaer*a_loc_ix[i,:]
                    wcba[i,:]=blwcld*c_loc_ix[i,:]
#                    if wcba[i,wcba[i,:]==True].shape[0]!=aawc[i,aawc[i,:]==True].shape[0]:
#                        print('wcba and aawc not same!!')
#                    if bmin-cmax<0:
#                        print("MALAMAGULAY>>>>>>>>>>>>..")
                    if np.sum(wcba[i,:]*1)>0:
                        wbasemintop=np.append(wbasemintop,bmin-cmax)
        #selecting aica
        if (np.sum(icldy_ix[i,:]*1)>0) and \
        (np.sum(aer_ix[i,:]*1)>0):#having both aerosols and water wcld
#                bmin=a_layer_base[i,aer_ix[i,:]].min()
                bs=a_layer_base[i,aer_ix[i,:]]
                cmax=c_layer_top[i,icldy_ix[i,:]]
                
                if np.sum((bs>cmax)*1)>0:
#                if bmin > cmax:#aerosol above op water cloud
                    bmin=bs[bs>cmax].min()
                    bmin_ix=np.where(a_layer_base[i,:]==bmin)
                    ctop_ix=np.where(c_layer_top[i,:]==cmax)
                    topaer=aer_ix[i,:]*0;topaer[np.max(bmin_ix)]=1;
                    blicld=icldy_ix[i,:]*0;blicld[np.min(ctop_ix)]=1;
                    aaic[i,:]=topaer*a_loc_ix[i,:]
                    icba[i,:]=blicld*c_loc_ix[i,:]
#                    if bmin-cmax<0:
#                        print("MALAMAGULAY>>>>>>>>>>>>.")
                    if np.sum(icba[i,:]*1)>0:
                        ibasemintop=np.append(ibasemintop,bmin-cmax)
        #selecting auphca
        if (np.sum(uphcldy_ix[i,:]*1)>0) and \
        (np.sum(aer_ix[i,:]*1)>0):#having both aerosols and uphase cld
#                bmin=a_layer_base[i,aer_ix[i,:]].min()
                bs=a_layer_base[i,aer_ix[i,:]]
                cmax=c_layer_top[i,uphcldy_ix[i,:]]
                
                if np.sum((bs>cmax)*1)>0:
#                if bmin > cmax:#aerosol above op water cloud
                    bmin=bs[bs>cmax].min()
                    bmin_ix=np.where(a_layer_base[i,:]==bmin)
                    ctop_ix=np.where(c_layer_top[i,:]==cmax)
                    topaer=aer_ix[i,:]*0;topaer[np.max(bmin_ix)]=1;
                    bluphcld=uphcldy_ix[i,:]*0;bluphcld[np.min(ctop_ix)]=1;
                    aauphc[i,:]=topaer*a_loc_ix[i,:]
                    uphcba[i,:]=bluphcld*c_loc_ix[i,:]
#                    if bmin-cmax<0:
#                        print("MALAMAGULAY>>>>>>>>>>>>.")
                    if np.sum(uphcba[i,:]*1)>0:
                        ubasemintop=np.append(ubasemintop,bmin-cmax)
#    print(str(cld_ix.size))
#    print(str(wcldy_ix.size+icldy_ix+uphcldy_ix.size))
    return cfnd,afnd,c_lat,c_lon,a_lat,a_lon,a_loc_ix,c_loc_ix,cld_ix,uphase_ix,\
        op_wcld,c_layer_DPR_532,c_layer_IAB_532,c_layer_IAB_1064,c_feature_type,\
        c_layer_base,c_layer_top,a_layer_base,a_layer_top,a_feature_type,\
        a_layer_IAB_532,a_layer_DPR_532,a_layer_IAB_1064,aer_ix,aac,aawc,aaic,\
        cba,wcba,icba,basemintop,wbasemintop,ibasemintop,wcldy_ix,icldy_ix,\
        uphcldy_ix,cldy_ix,aerbs_ix,aauphc,uphcba,ubasemintop
        
        
def findACA(obj,latlons=None):
    '''
    Return aac and cba indices
    The simplified function for Mumin
    obj: dsetCALCATS object
    latlons:[lon_mn,lon_mx,lat_mn,lat_mx]
    '''
    if latlons==None:
        lon_mn=obj.lon.min();lon_mx=obj.lon.max()
        lat_mn=obj.lat.min();lat_mx=obj.lat.max()
    else:
        lon_mn=latlons[0];lon_mx=latlons[1]
        lat_mn=latlons[2];lat_mx=latlons[3]
    
    feature_type, feature_type_qa,\
    ice_water_phase, ice_water_phase_qa,\
    feature_sub_type, cloud_aerosol_psc_type_qa,\
    horizontal_averaging=Extract_Feature_Info(obj.vfm,obj.nly)
    
    ix1=obj.lat[:,0]>=lat_mn;ix2=obj.lat[:,0]<=lat_mx
    lat_ix=ix1*ix2
    ix1=obj.lon[:,0]>=lon_mn;ix2=obj.lon[:,0]<=lon_mx
    lon_ix=ix1*ix2
    loc_ix=np.asarray([lat_ix*lon_ix]).T
    loc_ix=np.tile(loc_ix,obj.vfm.shape[1])

    aer_ix=feature_type==3#aerosols
    smoke_ix=(feature_sub_type==6)+(feature_sub_type==5)# smoke(6) + polluted dust(5)
    aer_ix=aer_ix*smoke_ix
    hoz_av80=horizontal_averaging==5# 80km horizontal averaging 
    cld_ix=feature_type==2#cloud
    wcld_ix=(ice_water_phase==2)+(ice_water_phase==0)# water cloud
    icld_ix=(ice_water_phase==1) + (ice_water_phase==3)#ice cloud
    uphase_ix=ice_water_phase==0#unknown phase
        
    aawc=np.zeros(feature_type.shape,dtype=bool)#aerosols above wcloud
    wcba=np.zeros(feature_type.shape,dtype=bool)#wcloud below aerosols
    wcldy_ix=np.zeros(feature_type.shape,dtype=bool)#water cloud
    cldy_ix=np.zeros(feature_type.shape,dtype=bool)#cloud
    uphcldy_ix=np.zeros(feature_type.shape,dtype=bool)#unknown phase clouds
    icldy_ix=np.zeros(feature_type.shape,dtype=bool)#ice cloud
    aac=np.zeros(feature_type.shape,dtype=bool)#aerosols above cloud
    cba=np.zeros(feature_type.shape,dtype=bool)#cloud below aerosols

    for i in range(0,feature_type.shape[0]):
        if np.sum(cld_ix[i,:]*1)>0:
            topcld_ix=np.min(np.where(cld_ix[i,:]==True))
            cldy_ix[i,topcld_ix]=True
            #selecting wcldy
            if np.sum(wcld_ix[i,:]*1)>0:
                topwcld_ix=np.min(np.where(wcld_ix[i,:]==True))
                if topwcld_ix==np.min(np.where(cld_ix[i,:]==True)):
                    wcldy_ix[i,topwcld_ix]=True
            #selecting icldy   
            if np.sum(icld_ix[i,:]*1)>0:
                topicld_ix=np.min(np.where(icld_ix[i,:]==True))
                if topicld_ix==np.min(np.where(cld_ix[i,:]==True)):
                    icldy_ix[i,topicld_ix]=True
            #selecting unknown phase clouds   
            if np.sum(uphase_ix[i,:]*1)>0:
                #found "unknown phase cloud" layers without "cloud" 
                topupcld_ix=np.min(np.where(uphase_ix[i,:]==True))
                if topupcld_ix==np.min(np.where(cld_ix[i,:]==True)):
                    uphcldy_ix[i,topupcld_ix]=True
        #selecting aca       
        if (np.sum(cldy_ix[i,:]*1)>0) and (np.sum(aer_ix[i,:]*1)>0):#having both aerosols and water cld
#       bmin=a_layer_base[i,aer_ix[i,:]].min()
            bs=obj.layer_base[i,aer_ix[i,:]]
            cmax=obj.layer_top[i,cldy_ix[i,:]].max()
                
            if np.sum((bs>cmax)*1)>0:
#               if bmin > cmax:#aerosol above op cloud
                bmin=bs[bs>cmax].min()
                bmin_ix=np.where(obj.layer_base[i,:]==bmin)
                ctop_ix=np.where(obj.layer_top[i,:]==cmax)
                topaer=aer_ix[i,:]*0;topaer[np.max(bmin_ix)]=1;
                blcld=cldy_ix[i,:]*0;blcld[np.min(ctop_ix)]=1;
                aac[i,:]=topaer*loc_ix[i,:]
                cba[i,:]=blcld*loc_ix[i,:]
    return aac,cba

def Attenuated_Backscatter_profile(chan,ds,ax,fig,ll_lim=None,v=None,cb_ticks=None,ll_ticks=500):
    '''
    To plot Attenuated_Backscatter profiles form CALIOP level1 data (2D contourf plot)
    ds: dsetCALCATS object (level1 data)
    ax: axis
    fig: figure 
    arg*
    ll_lim:[x1,y1] lat or lon limit
    v: v=np.linspace(a,b,N)
    cb_ticks: cb_ticks=np.arange(a,b,step)
    '''

    if chan=='1064':
        Attenuated_Backscatter=ds.Attenuated_Backscatter_1064[:,289:578].T
    elif chan=='532':
        Attenuated_Backscatter=ds.Attenuated_Backscatter_532[:,289:578].T
    else:
        print('Only 1064 and 532 channels!!! Error!!')
        
    if ll_lim==None:
        if v==None:
            ctf=ax.contourf(ds.lat[:,0],ds.Bin_altitude289_578,\
                 Attenuated_Backscatter,extend='both')
        else:
            ctf=ax.contourf(ds.lat[:,0],ds.Bin_altitude289_578,\
                 Attenuated_Backscatter,v,extend='both')
        lat=ds.lat[0::500]
        lon=ds.lon[0::500]
    else:
        if v==None:
            ctf=ax.contourf(ds.lat[ll_lim[0]:ll_lim[1],0],ds.Bin_altitude289_578,\
                 Attenuated_Backscatter[:,ll_lim[0]:ll_lim[1]],extend='both')
        else:
            ctf=ax.contourf(ds.lat[ll_lim[0]:ll_lim[1],0],ds.Bin_altitude289_578,\
                 Attenuated_Backscatter[:,ll_lim[0]:ll_lim[1]],v,extend='both')
        lat=ds.lat[ll_lim[0]:ll_lim[1]:ll_ticks]
        lon=ds.lon[ll_lim[0]:ll_lim[1]:ll_ticks]
    
    tick_labels=[]
    for i in np.arange(0,np.size(lat),1):
        tick_labels=np.append(tick_labels,'%0.1f,%0.1f'%(lon[i],lat[i]))
    ax.set_xlabel('lon,lat')
    ax.set_ylabel('Altitude (km)')
    ax.set_xticks(lat)
    ax.set_xticklabels(tick_labels,size=10)
    if cb_ticks==None:
        fig.colorbar(ctf,ax=ax,label="\n".join(wrap('Attenuated_Backscatter_'+chan+' (1/km/sr)',29)))
    else:
        fig.colorbar(ctf,ax=ax,ticks=cb_ticks,label="\n".join(wrap('Attenuated_Backscatter_'+chan+' (1/km/sr)',29)))

def check_CAl_avail(fltxtN,pdir):
    '''
    To check whether the all  CALIOP products exist(available) in the given directory. 
    Get the file list as a *.txt file form onine and use this function.
    fltxtN: file contains the list of cloud products name(ex.CALIOP_product_names/CAL_LID_L2_05kmCLay-Prov-V3-30.2016-06-to-08_cNEA_dust_study.txt)
    pdir:data product directory to check
    '''
    nlist=np.loadtxt(fltxtN,list)
    m=0#missing  count
    for i in np.arange(0,nlist.size):
#        print(nlist[i])
        if not(os.path.isfile(pdir+nlist[i])):
            m+=1
            if 'mc' in locals():
                mc.append(nlist[i])
            else:
                mc=[nlist[i]]
    print('%d missing files'%(m))
    if not(m==0):
        print('\n'.join(mc))
