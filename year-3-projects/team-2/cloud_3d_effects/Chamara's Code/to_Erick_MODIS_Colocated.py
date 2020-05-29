#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 10:32:20 2019

@author: cpnhere

The functions to share with Erik (MODIS-CALIPSO colocated data - spot 3D effects study)
"""

def getCOslwcld_from_CA1km(sd,colFlag):
    '''
    To get colocated single-layer water cloud
    For 1km data
    ------------
    cFlag: from MOD-CAL colocation data
    '''
    vfm = sd.select('Feature_Classification_Flags').get()
    nly = sd.select('Number_Layers_Found').get()
    feature_type, feature_type_qa,ice_water_phase, ice_water_phase_qa,\
    feature_subtype, cloud_aerosol_psc_type_qa,horizontal_averaging\
    =Extract_Feature_Info(vfm,nly)
    
    nlayers=feature_type.shape[1]
    colFlagEQ1km_feature=np.repeat(colFlag,nlayers).reshape(colFlag.size,nlayers)#successful collocations
    wcld_ix=(feature_type==2)*(ice_water_phase==2)
    single_layers=np.repeat(nly.squeeze()==1,nlayers).reshape(nly.size,nlayers)#single layer detections
    if wcld_ix.shape[0]==colFlagEQ1km_feature.shape[0]:
        COslwcld1km_feature=single_layers*wcld_ix*colFlagEQ1km_feature#Colocated single layer water cloud
        COslwcld1km=COslwcld1km_feature.any(True)
        col_mismatch=0
    else:
        col_mismatch=1
        COslwcld1km,COslwcld1km_feature=np.nan,np.nan
        print('stop')
    
    return COslwcld1km,COslwcld1km_feature,col_mismatch

def opacity_from1km(CAL1km,COslwcld1km_feature):
    opflag1km=readSDS(CAL1km,'Opacity_Flag',fill=99.0)
    op_1km=np.nansum(opflag1km*COslwcld1km_feature,axis=1)
    tra1km=np.nansum(-(opflag1km-1)*COslwcld1km_feature,axis=1)
    op1km=op_1km.astype(bool)
    tra1k=tra1km.astype(bool)
    return op1km,tra1k

    path='/umbc/xfs1/zzbatmos/common/Data/CALIPSO-MODIS-CloudSat/'
    indexP=path+'Index/'+year+'/'
    dataPa=path+'Data/'+year+'/'
    CALD1kmPa='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_01kmCLay/'+year+'/' # If 1km CALIOP data being used
    CALD5kmPa='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/'+year+'/'
    CA333p='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_333mMLay-Standard-V4-10/'
    
    for day in np.arange(total_days):
        fn_pattern='CAL_LID_L2_01kmCLay-ValStage1-V3-01.'+date_char+'*.hdf'
        for f in os.listdir(indexP):
            if fnmatch.fnmatch(f,fn_pattern):
                print(f)
                #reading data files
                #==================
                fI=SD(indexP+f)
                #from colocation index-----------------------------------------
                colFlag=fI.select('Collocation_Flag').get()
                sza=readSDS(fI,'MODIS_SZA1km')
                Solar Azimuth angle, viewing angle
                mod_lon=readSDS(fI,'MODIS_Lon1km')
                mod_lat=readSDS(fI,'MODIS_Lat1km')

                #from colocated MODIS data-------------------------------------
                fM=SD(dataPa+f)
                COD=readSDS(fM,'MYD06_Cloud_Optical_Thickness')
