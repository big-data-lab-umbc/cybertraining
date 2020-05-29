"""

********************************************
Created on Wed Sep  5 14:05:12 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

To spot 3D effects in MODIS COD retrievals by using MODIS-CALIPSO-CloudSat
co-located data. Check my 090518 log for more details.

11/04/2018:

    *_CA333m.pkl:
    
        Use 333mMLay V-4.10 data with MODIS colocated data. Other *.pkl s use
        5km CALIOP data.
    
    *_CA333m_v2.pkl:
    
        Use 5km Layer data to detect "single layer water clouds" and use 333m
        opacity flag.
    
11/09/2018:

    *_CA333m_v3.pkl:
    
        Use 1km layer data to identify "single layer water clouds" and use 333m
        opacity flag. Opacity qualtiy is defined. See the table in the
        comments.
    
11/13/2018:

    *_CA333m_v4.pkl:
    
        Uses 5km Layer data to identify "single layer water clouds" and use
        333m opacity flag
    
11/16/2018:

    *_CA1km_v1.pkl: Both opacity and water cloud detections from 1km

"""

# Standard library imports.
import datetime
import fnmatch
import glob
import os
import sys
import time
from textwrap import wrap

# Third party imports.
import numpy as np
import matplotlib.pyplot as plt

# Local application imports.
from cpnMODISlib import dsetMOD021K
from cpnMODISlib import drawrgb
from cpnCALCATSlib import Extract_Feature_Info
import cpnCommonlib as cpn


def readSDS(sd, label, fill = None):
    
    val = sd.select(label).get()
    
    if fill is None:
        
        fill = sd.select(label).attributes().get("_FillValue")
        
    val = val.astype(float)
    
    val[val == fill] = np.nan
    
    return np.squeeze(val)


def savefig(fig, fig_ttl):
    
    cpn.savefig(fig, fig_ttl, "figures/")
    
    
def getCOslwcld_from_CA1km(sd, colFlag):
    
    '''
    
    To get colocated single-layer water cloud
    For 1km data
    ------------
    cFlag: from MOD-CAL colocation data
    
    '''
    
    vfm = sd.select("Feature_Classification_Flags").get()
    nly = sd.select("Number_Layers_Found").get()
    
    (feature_type,
     feature_type_qa,
     ice_water_phase,
     ice_water_phase_qa,
     feature_subtype,
     cloud_aerosol_psc_type_qa,
     horizontal_averaging) = \
        Extract_Feature_Info(vfm, nly)
    
    nlayers = feature_type.shape[1]
    
    # successful collocations
    colFlagEQ1km_feature = np.repeat(colFlag, nlayers)
    colFlagEQ1km_feature = colFlagEQ1km_feature.reshape(colFlag.size, nlayers)
    
    wcld_ix = (feature_type == 2) * (ice_water_phase == 2)
    
    # single layer detections
    single_layers = np.repeat(nly.squeeze() == 1, nlayers)
    single_layers = single_layers.reshape(nly.size, nlayers)
    
    if wcld_ix.shape[0] == colFlagEQ1km_feature.shape[0]:
        
        # Colocated single layer water cloud
        COslwcld1km_feature = single_layers * wcld_ix * colFlagEQ1km_feature
        
        COslwcld1km  = COslwcld1km_feature.any(True)
        col_mismatch = 0
        
    else:
        
        col_mismatch = 1
        COslwcld1km, COslwcld1km_feature = np.nan, np.nan
        
        print("stop")
    
    return COslwcld1km, COslwcld1km_feature, col_mismatch


def getCOslwcld_from_CA5km(sd, colFlag):
    
    '''
    
    To get colocated single-layer water cloud
    For 5km data
    ------------
    colFlag: colFlag from MOD-CAL colocated data
    
    '''
    
    vfm = sd.select("Feature_Classification_Flags").get()
    nly = sd.select("Number_Layers_Found").get()
    
    (feature_type,
     feature_type_qa,
     ice_water_phase,
     ice_water_phase_qa,
     feature_subtype,
     cloud_aerosol_psc_type_qa,
     horizontal_averaging) = \
        Extract_Feature_Info(vfm, nly)
    
    nlayers = feature_type.shape[1]
    
    # all 5 1km MODIS-CALIOP successful collocation
    colFlagEQ5km = colFlag.reshape((5, colFlag.size / 5)).sum(axis = 0) > 4
    
    # successful collocations
    colFlagEQ5km_feature = np.repeat(colFlagEQ5km, nlayers)
    
    colFlagEQ5km_feature = colFlagEQ5km_feature.reshape(colFlagEQ5km.size,
                                                        nlayers)
    
    wcld_ix = (feature_type == 2) * (ice_water_phase == 2)
    
    # single layer detections
    single_layers = np.repeat(nly.squeeze() == 1, nlayers)
    single_layers = single_layers.reshape(nly.size, nlayers)
    
    # Colocated single layer water cloud
    COslwcld5km_feature = single_layers * wcld_ix * colFlagEQ5km_feature
    
    COslwcld5km = COslwcld5km_feature.any(True)
    
    return COslwcld5km, COslwcld5km_feature


def getCOslwcld_from_CA333(CA333m, cFlag):
    
    '''
    
    To get colocated single-layer water clouds from 333m CALIOP data and
    aggregate into 1km.
    -------------
    cFlag:boolean array
    
    '''
    
    vfm = CA333m.select("Feature_Classification_Flags").get()
    nly = CA333m.select("Number_Layers_Found").get()
    
    (feature_type,
     feature_type_qa,
     ice_water_phase,
     ice_water_phase_qa,
     feature_subtype,
     cloud_aerosol_psc_type_qa,
     horizontal_averaging) = \
        Extract_Feature_Info(vfm, nly)
    
    wcld_ix = (feature_type == 2) * (ice_water_phase == 2)
    
    # A successful colocated data point equivalent to 3 333m points
    colFlagEQ333m = cFlag.repeat(3)
    
    # successful collocations
    colFlagEQ333m_feature = np.repeat(colFlagEQ333m, feature_type.shape[1])
    colFlagEQ333m_feature = \
        colFlagEQ333m_feature.reshape(colFlagEQ333m.size,
                                      feature_type.shape[1])
    
    # single layer detections
    single_layers = np.repeat(nly.squeeze() == 1, feature_type.shape[1])
    single_layers = single_layers.reshape(nly.size, feature_type.shape[1])
    
    # Colocated single layer water cloud
    COslwcld = single_layers * wcld_ix * colFlagEQ333m_feature
    COslwcld_column = COslwcld.any(True)
    
    # Select if all neighboring 3 333m observations are single layer water
    # clouds
    COslwcld1km = COslwcld_column.reshape((3, COslwcld_column.size / 3))
    COslwcld1km = COslwcld1km.sum(axis = 0) > 2
    
    return COslwcld1km, COslwcld


def opacity_from333m(CA333m, COslwcld333m_feature):
    
    '''
    
    #Opacity qualtiy definition--------------------------------
     _________________________________________________________
    |  Label  |  x/3 opaque 333m shots  |Opacity qual.(1 best)|
    |_________|_________________________|_____________________|
    | 'Hig'   |           3             |         1.00        |
    | 'Med'   |           2             |         0.66        |
    | 'Low'   |           1             |         0.33        |
    |_________|_________________________|_____________________|
    
    '''
    
    opflag333m = readSDS(CA333m, "Opacity_Flag", fill = 99.0)
    
    op_333m = np.nansum(opflag333m * COslwcld333m_feature, axis = 1)
    tra333m = np.nansum(-(opflag333m - 1) * COslwcld333m_feature, axis = 1)
    
    # All 3 333m lidar shots see opaque cloud layer
    op1km1 = op_333m.reshape((3, op_333m.size / 3)).sum(axis = 0) > 2
    
    # 2/3 333m lidar shots see opaque cloud layer
    op1km2 = op_333m.reshape((3, op_333m.size / 3)).sum(axis = 0) > 1
    
    # 1/3 333m lidar shots see opaque cloud layer
    op1km3 = op_333m.reshape((3, op_333m.size / 3)).sum(axis = 0) > 0
    
    # All 3 333m lidar shots see non-opaque cloud layers
    tra1k1 = tra333m.reshape((3, tra333m.size / 3)).sum(axis = 0) > 2
    
    # 2/3 333m lidar shots see non-opaque cloud layers
    tra1k2 = tra333m.reshape((3, tra333m.size / 3)).sum(axis = 0) > 1
    
    # 1/3 333m lidar shots see non-opaque cloud layers
    tra1k3 = tra333m.reshape((3, tra333m.size / 3)).sum(axis = 0) > 0
    
    return op1km1, op1km2, op1km3, tra1k1, tra1k2, tra1k3


def opacity_from1km(CAL1km, COslwcld1km_feature):
    
    opflag1km = readSDS(CAL1km, "Opacity_Flag", fill = 99.0)
    
    op_1km = np.nansum(opflag1km * COslwcld1km_feature, axis = 1)
    tra1km = np.nansum(-(opflag1km - 1) * COslwcld1km_feature, axis = 1)
    
    op1km = op_1km.astype(bool)
    tra1k = tra1km.astype(bool)
    
    return op1km, tra1k


def find_mode(f, xe):
    
    x = (xe[0:-1] + xe[1:]) / 2
    
    return x[np.argmax(f)]


def plot_CODHist(modCOD, temp_stamp, ax1, tps, normed = True, bins = None):
    
    '''
    
    To plot COD histograms
    modCOD: Dictionary of different data sets
    tps: string array to select data sets from modCOD
    temp_stamp: Prefix to the file name
    
    '''
    
    if bins is None:
        
        bins = np.logspace(np.log10(0.01), np.log10(151.0), 100)
        
    fig1_ttl = temp_stamp + "_single_layer_wcld_CoLCALIPSO-MODIS_COD"
    
    kwargs = {"normed"    : normed,
              "histtype"  : "step",
              "linewidth" : 2.0,
              "bins"      : bins}
    
    f0, xe0, _ = ax1.hist(modCOD[tps[0]][~np.isnan(modCOD[tps[0]])],
                          color = "b",
                          label = tps[0],
                          **kwargs)
    
    f1, xe1, _ = ax1.hist(modCOD[tps[1]][~np.isnan(modCOD[tps[1]])],
                          color = "k",
                          label = tps[1],
                          **kwargs)
    
    f2, xe2, _ = ax1.hist(modCOD[tps[2]][~np.isnan(modCOD[tps[2]])],
                          color = "r",
                          label = tps[2],
                          **kwargs)
    
    # ax1.hist(modCOD["lowSun"][~np.isnan(modCOD["lowSun"])],
    #                           color = "g",
    #                           label = "Low Sun",
    #                           **kwargs)
    
    # ax1.hist(modCOD["highSn"][~np.isnan(modCOD["highSn"])],
    #                           color = "c",
    #                           label = "High Sun",
    #                           **kwargs)
    
    ax1.legend(loc="upper left")
    ax1.set_xlabel("COD")
    ax1.set_xscale("log")
    ax1.set_xlabel("COD")
    
    if normed:
        
        ax1.set_ylabel("PDF")
        
        ax1.text(0.6,
                 0.75,
                 (
                     "md, mn, std: {0.1f}, {0.1f}, {0.1f}"
                 ).format(find_mode(f2, xe2),
                          np.nanmean(modCOD[tps[2]]),
                          np.nanstd(modCOD[tps[2]])),
                 transform = ax1.transAxes,
                 color     = "r")
        
        ax1.text(0.6,0.80,"md,mn,std: %0.1f, %0.1f, %0.1f"%(find_mode(f1,xe1),np.nanmean(modCOD[tps[1]]),\
                                                   np.nanstd(modCOD[tps[1]])),transform = ax1.transAxes,color="k")
        
        ax1.text(0.6,0.85,"md,mn,std: %0.1f, %0.1f, %0.1f"%(find_mode(f0,xe0),np.nanmean(modCOD[tps[0]]),\
                                            np.nanstd(modCOD[tps[0]])),transform = ax1.transAxes,color="b")
    
    return fig1_ttl


def plot_CODCDF(modCOD,temp_stamp,tps=["calopq","allang","calTra","TraLow","TraHig"],bins=np.logspace(np.log10(0.01),np.log10(151.0), 100)):
    
    '''
    
    COD CDF
    
    '''
    
    b0,c0=cpn.find_CDF(modCOD[tps[0]][~np.isnan(modCOD[tps[0]])],bins=bins)
    b1,c1=cpn.find_CDF(modCOD[tps[1]][~np.isnan(modCOD[tps[1]])],bins=bins)
    b2,c2=cpn.find_CDF(modCOD[tps[2]][~np.isnan(modCOD[tps[2]])],bins=bins)
    b3,c3=cpn.find_CDF(modCOD[tps[3]][~np.isnan(modCOD[tps[3]])],bins=bins)
    b4,c4=cpn.find_CDF(modCOD[tps[4]][~np.isnan(modCOD[tps[4]])],bins=bins)
    

    fig3,ax3=plt.subplots()
    fig3_ttl=temp_stamp+'_COD_CDF'
    fig3.suptitle(fig3_ttl)
    ax3.plot(b0,c0,'b',label=tps[0])
    ax3.plot(b1,c1,'k',label=tps[1])
    ax3.plot(b2,c2,'r',label=tps[2])
    ax3.plot(b3,c3,'r--',label=tps[3])
    ax3.plot(b4,c4,'r:',label=tps[4])
    ax3.set_xscale('log')
    ax3.set_ylabel('CDF')
    ax3.set_xlabel('COD')
    ax3.legend(loc='best')
    fig3.tight_layout(rect=[0,0,1,0.98])
    fig3.show()
    
    return ax3, fig3, fig3_ttl


def Use_CALIOP5km(year,f,calCOD,modCOD,calOpa,pt_ill,pilmax,pe_names,COslwcld1km,lowSun,highSn,colFlag,COD,show_3Df,mod_lon,mod_lat):
    #Reading CALIOP 5km data
    CA5kmP='/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/'+year+'/'
    C5km=SD(CA5kmP+f.replace('01km','05km').replace('_IND','').replace('ValStage1','Prov'))
    colFlagEQ5km=colFlag.reshape((5,colFlag.size/5)).sum(axis=0)>4#all 5 1km MODIS-CALIOP successful collocation
    lowSunEQ_5km=lowSun.reshape((5,lowSun.size/5)).sum(axis=0)==5#all 5 1km MODIS observations have low sun geometry
    highSnEQ_5km=highSn.reshape((5,highSn.size/5)).sum(axis=0)==5#all 5 1km MODIS observations have high sun geometry
    colFlagEQ5km_feature=np.repeat(colFlagEQ5km,10).reshape(colFlagEQ5km.size,10)#successful collocations
    lowSunEQ_5km_feature=np.repeat(lowSunEQ_5km,10).reshape(lowSunEQ_5km.size,10)#MODIS low sun
    highSnEQ_5km_feature=np.repeat(highSnEQ_5km,10).reshape(highSnEQ_5km.size,10)#MODIS high sun
    vfm = C5km.select('Feature_Classification_Flags').get()
    nly = C5km.select('Number_Layers_Found').get()
    
    (feature_type,
     feature_type_qa,
     ice_water_phase,
     ice_water_phase_qa,
     feature_subtype,
     cloud_aerosol_psc_type_qa,
     horizontal_averaging) = \
        Extract_Feature_Info(vfm, nly)
        
    wcld_ix=(feature_type==2)*(ice_water_phase==2)
    single_layers=np.repeat(nly.squeeze()==1,10).reshape(nly.size,10)#single layer detections
    COslwcld=single_layers*wcld_ix*colFlagEQ5km_feature#Colocated single layer water cloud
    
    opflag5km=readSDS(C5km,'Opacity_Flag',fill=99.0)
    FOD_532=readSDS(C5km,'Feature_Optical_Depth_532',fill=-9999.0)
    
    calCOD['allang']=np.append(calCOD['allang'],np.nansum(FOD_532*(COslwcld*1),axis=1))#single layer water cloud MODIS-CALIOP colocated
    calCOD['lowSun']=np.append(calCOD['lowSun'],np.nansum(FOD_532*(COslwcld*lowSunEQ_5km_feature*1),axis=1))#single layer water cloud MODIS-CALIOP colocated low Sun
    calCOD['highSn']=np.append(calCOD['highSn'],np.nansum(FOD_532*(COslwcld*highSnEQ_5km_feature*1),axis=1))#single layer water cloud MODIS-CALIOP colocated high Sun
    
    op5km=np.nansum(opflag5km*COslwcld,axis=1)
    tra5k=np.nansum(-(opflag5km-1)*COslwcld,axis=1)
    op1km=np.repeat(op5km,5)
    tra1k=np.repeat(tra5k,5)
    calOpa['allang']=np.append(calOpa['allang'],op1km[COslwcld1km])
    modCOD['calopq']=np.append(modCOD['calopq'],COD[op1km.astype(bool)*COslwcld1km])
    modCOD['calTra']=np.append(modCOD['calTra'],COD[tra1k.astype(bool)*COslwcld1km])
    x=COD[tra1k.astype(bool)*lowSun]
    modCOD['TraLow']=np.append(modCOD['TraLow'],x)
    modCOD['TraHig']=np.append(modCOD['TraHig'],COD[tra1k.astype(bool)*highSn])
    pil=(x>20).sum()#potentila illumination effects (ie. low sun cases that have COD>20)
    if pilmax<pil and show_3Df:
        pe_names=np.append(pe_names,f)
        print('%d found: '%(pil))
        print('\t'+f)
        print('\t'+f.replace('01km','05km').replace('_IND','').replace('ValStage1','Prov'))
        pilmax=pil
    pot_il=tra1k.astype(bool)*lowSun*(COD>20)#potential illuminating effects
    pt_ill['lat']=np.append(pt_ill['lat'],mod_lat[pot_il])
    pt_ill['lon']=np.append(pt_ill['lon'],mod_lon[pot_il])
    C5km.end()
    return pilmax


def find_pot_3Deffects(show_3Df = True, year = "2010", end_m = 12, end_d = 31):
    
    '''
    
    Find potential 3D effects in CALIOP-MODIS co-located data.
    Save following *.pkl objects separately
    modCOD: MODIS COD for different criteria 
        allang: All data/geometries
        lowSun: Low sun
        highSn: High sun
        calopq: clouds that apear opaque to CALIOP
        calTra: clouds that apear transparent to CALIOP
        TraLow: non-opaque to CAL with lower solar geometry
        TraHig: non-opaque to CAL with higher solar geometry
    calCOD: CALIOP Layer optical thickness for different criteria
        allang: All data/geometries
        lowSun: Low sun
        highSn: High sun
    -----------------------------------------------------------------
    show_3Df: Show files with potention 3D effects to a file
    11/05/2018:
        Either 1km or 333m CALIOP data can be used to select single layer water clouds.
        Similarlly, either 5km or 333m CALIOP data can be used to get opacity flag.
        (Hard coded. Check and edit the code)
    11/09/2018:
        3 separate opacity quality flags are defined based on 333km opacity flag.
        #Opacity qualtiy definition--------------------------------
         _________________________________________________________
        |  Label  |  x/3 opaque 333m shots  |Opacity qual. factor |
        |_________|_________________________|____(1 best)_________|
        | 'Hig'   |           3             |         1.00        |
        | 'Med'   |           2             |         0.66        |
        | 'Low'   |           1             |         0.33        |
        |_________|_________________________|_____________________|
        
    '''
    
    version_append = "_CA1km_v1"
    
    start_date = datetime.date(int(year), 1, 1)
    end_date   = datetime.date(int(year), end_m, end_d)
    
    total_days = (end_date - start_date).days + 1
    
    temp_stamp="%d-%02d-%02d-to-%d-%02d-%02d"%(start_date.year,start_date.month,start_date.day,end_date.year,end_date.month,end_date.day)
    path="/umbc/xfs1/zzbatmos/common/Data/CALIPSO-MODIS-CloudSat/"
    
    indexP = path + "Index/" + year + "/"
    dataPa = path + "Data/"  + year + "/"
    
    op_from = "1km"  # or "333m"
    wc_from = "1km"  # or "5km"
    
    CALD1kmPa="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_01kmCLay/"+year+"/" # If 1km CALIOP data being used
    CALD5kmPa="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/"+year+"/"
    CA333p="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_333mMLay-Standard-V4-10/"
    
    modCOD = {"allang" : []}
    calCOD = {"allang" : []}
    
    calOpa = {"Hig"    : [],
              "Med"    : [],
              "Low"    : []}
    
    calTra = {"Hig"    : [],
              "Med"    : [],
              "Low"    : []}
    
    modSZA = {"allang" : []}
    
    locs={"lat":[],"lon":[]}#locations of the potential illuminating effects
    counts={"missing":0,"agg_incomp":0,"col_fail":0,"col_misMtch":0}#missing files(333m) and aggregation incompatible file counts
    
    t1=time.time()
    #Putting output into a file
    orig_stdout=sys.stdout
    outf = open(temp_stamp+version_append+'.dat','w')
    sys.stdout=outf
    for day in np.arange(total_days):
        d=start_date+datetime.timedelta(day)
        yc='%04i'%d.year
        mc='%02i'%d.month
        dc='%02i'%d.day
        date_char=yc+'-'+mc+'-'+dc
        fn_pattern='CAL_LID_L2_01kmCLay-ValStage1-V3-01.'+date_char+'*.hdf'
        for f in os.listdir(indexP):
            if fnmatch.fnmatch(f,fn_pattern):
                print(f)
                #reading data files
                #==================
                fI=SD(indexP+f)
                fM=SD(dataPa+f)
                if wc_from=='1km' or op_from=='1km':
                    try:
                        CAL1km=SD(CALD1kmPa+f.replace('_IND','').replace('ValStage1','Standard').replace('V3-01','V4-10'))# if 1km CAL data being used
                    except Exception as e:
                        print(e)
                        break
                if wc_from=='5km':
                    CAL5km=SD(CALD5kmPa+f.replace('01km','05km').replace('_IND','').replace('ValStage1','Prov'))
                if op_from=='333m':
                    if os.path.isfile(CA333p+f.replace('01km','333m').replace('CLay','MLay').replace('_IND','').replace('ValStage1','Standard').replace('V3-01','V4-10')):
                        CA333m=SD(CA333p+f.replace('01km','333m').replace('CLay','MLay').replace('_IND','').replace('ValStage1','Standard').replace('V3-01','V4-10'))
                    else:
                        print(f.replace('01km','333m').replace('CLay','MLay').replace('_IND','').replace('ValStage1','Standard').replace('V3-01','V4-10')+' MISSING!!')    
                        fI.end();fM.end()
                        counts['missing']+=1
                        break
                #from colocation index-----------------------------------------
                colFlag=fI.select('Collocation_Flag').get()
                sza=readSDS(fI,'MODIS_SZA1km')
                mod_lon=readSDS(fI,'MODIS_Lon1km')
                mod_lat=readSDS(fI,'MODIS_Lat1km')
                #from colocated MODIS data-------------------------------------
                COD=readSDS(fM,'MYD06_Cloud_Optical_Thickness')
                
                #Selecting single layer water clouds
                #===================================
                if wc_from=='1km':
                    COslwcld1km,COslwcld1km_feature,cmt=getCOslwcld_from_CA1km(CAL1km,colFlag)#from CALIPSO 1km
                    if cmt==1:
                        counts['col_misMtch']+=1
                        break
                elif wc_from=='5km':
                    colFlag5km=colFlag.reshape((5,colFlag.size/5)).sum(axis=0)>4
                    COslwcld5km,COslwcld5km_feature=getCOslwcld_from_CA5km(CAL5km,colFlag) #from CALIPSO 5km
                    COslwcld1km=COslwcld5km.repeat(5)
                    COslwcld1km_feature=COslwcld5km_feature.repeat(5,axis=0)
                    #Checking Colocation quality
                    if colFlag5km.sum()/colFlag5km.size<0.9 and colFlag.sum()/colFlag.size>0.9:
                        print('Colocation quality reduction threshold achieved when aggregate')
                        counts['col_fail']+=1
                        break  
                elif wc_from=='333m':
                    COslwcld1km,COslwcld=getCOslwcld_from_CA333(CA333m,colFlag==1)#from CALIPSO 333m(to select single layer water clouds)
                if wc_from=='5km' and op_from=='333m':
                    #from CALIOP333m data (to find single layers layer identifications)
                    vfm = CA333m.select('Feature_Classification_Flags').get()
                    nly = CA333m.select('Number_Layers_Found').get()
                    feature_type, _,_, _,_, _,_=Extract_Feature_Info(vfm,nly)
    
                    if COslwcld1km.size*3!=feature_type.shape[0]:
                        
                        print('Aggregated 333m CALIOP does not match with 1km.')
                        counts['agg_incomp']+=1
                        break
                        
                    if feature_type.shape[1]*2==COslwcld1km_feature.shape[1]:
                        
                        cfac = feature_type.shape[0] / COslwcld1km_feature.shape[0]
                        
                        COslwcld333m_feature=np.repeat(COslwcld1km_feature[:,0::2] + COslwcld1km_feature[:,1::2],cfac,axis=0)
                    
                    else:
                        
                        print('333m n feature is not twice as 5km')
                        counts['agg_incomp'] += 1
                        break
                        
                # ----------------------------------------------------------------
                # SZA and COD from MODIS
                # ======================
                modSZA["allang"] = np.append(modSZA["allang"],
                                             sza[COslwcld1km])
                
                modCOD["allang"] = np.append(modCOD["allang"],
                                             COD[COslwcld1km])
                
                locs["lat"] = np.append(locs["lat"], mod_lat[COslwcld1km])
                locs["lon"] = np.append(locs["lon"], mod_lon[COslwcld1km])
                
                # Opacity Flag
                # =============
                
                if op_from == "1km":
                    
                    op1km, tra1k = opacity_from1km(CAL1km,
                                                   COslwcld1km_feature)
                    
                    calOpa["Hig"] = np.append(calOpa["Hig"],
                                              op1km[COslwcld1km])
                    
                    calTra["Hig"] = np.append(calTra["Hig"],
                                              tra1k[COslwcld1km])
                    
                elif op_from == "333m":
                    
                    (op1km1,
                     op1km2,
                     op1km3,
                     tra1k1,
                     tra1k2,
                     tra1k3) = \
                        opacity_from333m(CA333m, COslwcld333m_feature)
                    
                    CA333m.end()
                    
                    calOpa["Hig"] = np.append(calOpa["Hig"],
                                              op1km1[COslwcld1km])
                    
                    calOpa["Med"] = np.append(calOpa["Med"],
                                              op1km2[COslwcld1km])
                    
                    calOpa["Low"] = np.append(calOpa["Low"],
                                              op1km3[COslwcld1km])
                    
                    calTra["Hig"] = np.append(calTra["Hig"],
                                              tra1k1[COslwcld1km])
                    
                    calTra["Med"] = np.append(calTra["Med"],
                                              tra1k2[COslwcld1km])
                    
                    calTra["Low"] = np.append(calTra["Low"],
                                              tra1k3[COslwcld1km])

#                pilmax=Use_CALIOP5km(year,f,calCOD,modCOD,calOpa,pt_ill,pilmax,pe_names,COslwcld1km,lowSun,highSn,colFlag,COD,show_3Df,mod_lon,mod_lat)
                fI.end();fM.end();
    print('Missing 333m CALIOP files: %d'%counts['missing'])
    print('333m and 1km CALIOP incompatibility issues: %d'%counts['agg_incomp'])
    print('Colocation quality reduction cases: %d'%counts['col_fail'])
    print('Colocation file 1km CALIOP incompatibility issues: %d'%counts['col_misMtch'])
    t2=time.time()
    print('Time: %0.2f mins'%((t2-t1)/60))
    obj={'modCOD':modCOD,'calCOD':calCOD,'locs':locs,'modCOD':modCOD,'calOpa':calOpa,'calTra':calTra,'modSZA':modSZA}
    cpn.save_obj(obj,temp_stamp+version_append)
    sys.stdout=orig_stdout
 
def MODIS_CALIOP_single():
    
    '''
    
    Single case of MODIS-CALIOP colocated
    Plot MODIS RGB and CALIOP track.
    Overlay coresponding colocated data of potential 3D effects (low sun)
    
    '''
    
    year = "2007"    
    # hS_thr = 30
    lS_thr=70
    path="/umbc/xfs1/zzbatmos/common/Data/CALIPSO-MODIS-CloudSat/"
    indexP=path+"Index/"+year+"/"
    dataPa=path+"Data/"+year+"/"
#    CALDPa="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_01kmCLay/"+year+"/"
    CA5kmP="/umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/"+year+"/"
    f="CAL_LID_L2_01kmCLay-ValStage1-V3-01.2007-01-11T01-46-36ZD_IND.hdf"
    #Co-located index file
    fI=SD(indexP+f)
    colFlag=fI.select("Collocation_Flag").get()#Collocation flag
#    N_MODIS=fI.select("N_modis").get()#Number of MODIS granules
    gran_ID=readSDS(fI,"MODIS_Granule_Ind")#MODIS granuele ID
    mod_DOY=readSDS(fI,"MODIS_DOY")#MODIS day of the year
    mod_min=readSDS(fI,"MODIS_Minute")#minute of the MODIS granule
    modhour=readSDS(fI,"MODIS_Hour")# hour MODIS
    mod_lon=readSDS(fI,"MODIS_Lon1km")
    mod_lat=readSDS(fI,"MODIS_Lat1km")
    cal_lon=readSDS(fI,"CALIPSO_Lon1km")
    cal_lat=readSDS(fI,"CALIPSO_Lat1km")
    sza=np.deg2rad(readSDS(fI,"MODIS_SZA1km"))
    #reading CALIPSO 1km
#    C1km=SD(CALDPa+f.replace("_IND",""))
#    COslwcld1km=getCOslwcld_from_CA1km(C1km,colFlag==1)
                    
    #reading colocated MODIS data
    fM=SD(dataPa+f)
    COD=readSDS(fM,"MYD06_Cloud_Optical_Thickness")
#    highSn=sza<np.deg2rad(hS_thr)
    lowSun=sza>np.deg2rad(lS_thr)
    
    #Reading CALIOP 5km data
    C5km=SD(CA5kmP+f.replace("01km","05km").replace("_IND","").replace("ValStage1","Prov"))
    colFlagEQ5km=colFlag.reshape((5,colFlag.size/5)).sum(axis=0)>4#all 5 1km MODIS-CALIOP successful collocation
    # lowSunEQ_5km=lowSun.reshape((5,lowSun.size/5)).sum(axis=0)==5#all 5 1km MODIS observations have low sun geometry
    # highSnEQ_5km=highSn.reshape((5,highSn.size/5)).sum(axis=0)==5#all 5 1km MODIS observations have high sun geometry
    colFlagEQ5km_feature=np.repeat(colFlagEQ5km,10).reshape(colFlagEQ5km.size,10)#successful collocations
    # lowSunEQ_5km_feature=np.repeat(lowSunEQ_5km,10).reshape(lowSunEQ_5km.size,10)#MODIS low sun
    # highSnEQ_5km_feature=np.repeat(highSnEQ_5km,10).reshape(highSnEQ_5km.size,10)#MODIS high sun
    
    vfm = C5km.select("Feature_Classification_Flags").get()
    nly = C5km.select("Number_Layers_Found").get()
    
    (feature_type,
     feature_type_qa,
     ice_water_phase,
     ice_water_phase_qa,
     feature_subtype,
     cloud_aerosol_psc_type_qa,
     horizontal_averaging) = \
        Extract_Feature_Info(vfm, nly)
    
    wcld_ix = (feature_type == 2) * (ice_water_phase == 2)
    
    # single layer detections
    single_layers = np.repeat(nly.squeeze() == 1, 10).reshape(nly.size, 10)
    
    # Colocated single layer water cloud
    COslwcld = single_layers * wcld_ix * colFlagEQ5km_feature
    
    opflag5km = readSDS(C5km, "Opacity_Flag", fill = 99.0)
    # FOD_532=readSDS(C5km,"Feature_Optical_Depth_532",fill=-9999.0)

#    op5km=np.nansum(opflag5km*COslwcld,axis=1)
    tra5k=np.nansum(-(opflag5km-1)*COslwcld,axis=1)
#    op1km=np.repeat(op5km,5)
    tra1k=np.repeat(tra5k,5)
    pot_il=tra1k.astype(bool)*lowSun*(COD>20)#potential illuminating effects
    print(gran_ID[pot_il])

    #MODIS RGB
#    mod_paths={"MYD021KM":"/umbc/xfs1/zzbatmos/common/Data/MODIS/6/MYD021KM/2012/","MYD06":}
    data_path="/umbc/xfs1/zzbatmos/common/Data/MODIS/MYD021KM/2007/"
    gid=1
#    fname="MYD021KM.A2012129.1750.006.2012186173545.hdf"
#    "MYD06_L2.A2015227.1520.006.2015228174651.hdf"

    fname = glob.glob(data_path+"MYD021KM.A%d%03d.%02d%02d.061.*.hdf"%(int(year),mod_DOY[gid],modhour[gid],mod_min[gid]))[0]
    hdf = SD(fname)
    
    res250 = dsetMOD021K("EV_250_Aggr1km_RefSB")  # [11]EV_250_Aggr1km_RefSB
    res250.readres(hdf)

    res500 = dsetMOD021K("EV_500_Aggr1km_RefSB")  # [17]EV_500_Aggr1km_RefSB
    res500.readres(hdf)

    fig1, ax1 = plt.subplots()
    fig1_ttl = fname.split("/", -1)[-1] + "_RGB"
    
    m = drawrgb(res250, res500, ax1, bfac = 0.5, rtn = True)
    
    x0, y0 = m(mod_lon[pot_il], mod_lat[pot_il])
    x2, y2 = m(cal_lon, cal_lat)
    
    ax1.plot(x2, y2, "y.", markersize = 0.5)
    ax1.plot(x0, y0, "r.", markeredgewidth = 0)
    
    fig1.suptitle(fig1_ttl)
    fig1.show()
    

if __name__ == "__main__":

    # end_m = int(sys.argv[1])
    # end_d = int(sys.argv[2])
    # find_pot_3Deffects(year = "2007", end_m = end_m, end_d = end_d)
    
    temp_stamp="2007-01-01-to-2007-06-30_CA1km_v1"
    obj=cpn.load_obj(temp_stamp)    
    modCOD=obj["modCOD"]
    modSZA=obj["modSZA"]["allang"] 
    lowSun=modSZA>70
    higSun=modSZA<30
    opQ="Hig" #Opacity quality
    modSZA=obj["modSZA"]["allang"] 
    calTra=obj["calTra"][opQ]
    calOpa=obj["calOpa"][opQ]
    modCOD["calTra"]=modCOD["allang"][calTra.astype(bool)]
    modCOD["TraLow"]=modCOD["allang"][calTra.astype(bool)*lowSun]
    modCOD["TraHig"]=modCOD["allang"][calTra.astype(bool)*higSun]
    modCOD["calopq"]=modCOD["allang"][calOpa.astype(bool)]
    cpn.setup_figures(plt)
    fig1,ax1=plt.subplots(2,1,figsize=(8,12))
    fig1_ttl=plot_CODHist(modCOD,temp_stamp+"_opQ"+opQ,ax1[0],tps=["calopq","allang","calTra"],normed=True) 
    _=plot_CODHist(modCOD,temp_stamp+"_opQ"+opQ,ax1[1],tps=["calopq","allang","calTra"],normed=False) 
    fig1.suptitle("\n".join(wrap(fig1_ttl,70)))
    fig1.tight_layout(rect=[0,0,1,0.97])
    fig1.show()
    fig2,ax2=plt.subplots(2,1,figsize=(8,12))
    fig2_ttl=plot_CODHist(modCOD,temp_stamp+"_Tran_opQ"+opQ,ax2[0],tps=["TraLow","calTra","TraHig"],normed=True)  
    _       =plot_CODHist(modCOD,temp_stamp+"_Tran_opQ"+opQ,ax2[1],tps=["TraLow","calTra","TraHig"],normed=False)
    fig2.suptitle("\n".join(wrap(fig2_ttl,70)))
    fig2.tight_layout(rect=[0,0,1,0.97])
    fig2.show()
    _,fig3,fig3_ttl=plot_CODCDF(modCOD,temp_stamp+"_opQ"+opQ)    


    #COT-SZA vs CALIOP transparency denisty plot
    dSZA=1;
    SZA_bnds=np.arange(22,80+dSZA,dSZA)
    COT_bnds=np.logspace(np.log10(0.2),np.log10(200.0), 100)
    SZA_grids=(SZA_bnds[0:-1]+SZA_bnds[1:])/2.0
    COT_grids=(COT_bnds[0:-1]+COT_bnds[1:])/2.0
    Ntot,xedges,yedges=np.histogram2d(modCOD["allang"],modSZA,bins=(COT_bnds,SZA_bnds))
#    Ntra,xedges,yedges=np.histogram2d(modCOD["allang"]*np.invert(calOpa.astype(bool)),modSZA*np.invert(calOpa.astype(bool)),bins=(COT_bnds,SZA_bnds))
    Ntra,xedges,yedges=np.histogram2d(modCOD["allang"]*calTra.astype(bool),modSZA*calTra.astype(bool),bins=(COT_bnds,SZA_bnds))
    [COT_mesh, SZA_mesh] = np.meshgrid(COT_grids, SZA_grids)
    fig5,ax5=plt.subplots()
    fig5_ttl=temp_stamp+"_CALIOP_single_layer_Tran_wcld_frac"
    v = np.linspace(0.0,0.25,num=50)
    #fraction
    ctf1 = ax5.contourf(COT_mesh,SZA_mesh,(Ntra/Ntot).T,v,cmap=plt.cm.jet,extend="both")
#    ctf2 = ax5.contour(COT_mesh,SZA_mesh,Ntot.T/Ntot.max(),8,linewidths=1.0,colors="k")
    
    ax5.set_xscale("log")
    ax5.set_xlabel("MODIS COT")
    ax5.set_ylabel("MODIS SZA")
    cpn.add_cb(fig5,ctf1,ax5,ticks=np.arange(0,0.251,0.05),label="CALIOP Transparent Cloud Fraction",orientation="vertical")
    # ax5.clabel(ctf2, inline=1, fontsize=10)
    fig5.suptitle(fig5_ttl)
    fig5.tight_layout(rect = [0, 0, 1, 0.97])
    fig5.show()
    
    # Calculations (*_CA333m_v4")
    # print("Single layer wcld opacity:"+temp_stamp)
    # Nt=modCOD["allang"].size
    # print("Total single layer water cloud cases: %0.4E"%Nt)
    # print("\t High quality opaque CF: %0.4f"%(obj["calOpa"]["Hig"].sum()/Nt))
    # print("\t Medm quality opaque CF: %0.4f"%(obj["calOpa"]["Med"].sum()/Nt))
    # print("\t Low  quality opaque CF: %0.4f"%(obj["calOpa"]["Low"].sum()/Nt))
    # print("\t High quality Transp. CF: %0.4f"%(obj["calTra"]["Hig"].sum()/Nt))
    # print("\t Medm quality Transp. CF: %0.4f"%(obj["calTra"]["Med"].sum()/Nt))
    # print("\t Low  quality Transp. CF: %0.4f"%(obj["calTra"]["Low"].sum()/Nt))
    
    # Calculations (*_CA1km_v1)
    print("Single layer wcld opacity:" + temp_stamp)
    
    Nt = modCOD["allang"].size
    
    print("Total single layer water cloud cases: %0.4E" % Nt)
    print("\t Opaque CF: %0.4f"  % (obj["calOpa"]["Hig"].sum() / Nt))
    print("\t Transp. CF: %0.4f" % (obj["calTra"]["Hig"].sum() / Nt))

    # Locations of the potential 3D effects projected on the world map
    # fig4,ax4=plt.subplots(figsize=(10,5))
    # fig4_ttl=temp_stamp+"_Potential_3d_effects"
    # m=cpn.mapProject(ax4)
    # lon_ill,lat_ill=m(pt_ill["lon"],pt_ill["lat"])
    # ax4.scatter(lon_ill,lat_ill,s=1.0)
    # fig4.suptitle(fig4_ttl)
    # fig4.tight_layout(rect=[0.03,0.03,1,0.98])
    # fig4.show()
