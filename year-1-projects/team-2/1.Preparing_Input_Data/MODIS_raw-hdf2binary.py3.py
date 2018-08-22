"""
Convert MODIS HDF to binary file using h5dump

By Daeho Jin, 2018.03.27

* Variable Name
: Cloud_Optical_Thickness_ISCCP_JHisto_vs_Pressure  <== histogram of count, not CF
:: Initial dimension=[lon(360),lat(180),COT(8),CTP(7)], fill_value=-9999
: Cloud_Retrieval_Fraction_Combined
:: Initial dimension=[lon(360),lat(180)], fill_value=-9999, scale_factor=0.0001

: Cloud_Optical_Thickness_PCL_ISCCP_JHisto_vs_Pressure  <== histogram of count, not CF
: Cloud_Retrieval_Fraction_PCL_Combined

* Processing for joint histogram
: (bin_count/tot_count) * grid_Cloud_Fraction

* How to treat missing
: 1. if cf=0, set as clear sky (histogram counts=0)
: 2. if all bins of histogram are missing, set missing
: 3. if there are partial missings, set it as zero.

* Processing for combining Normal and PCL
: if normal count>0 and pcl count>0 then sum two.
: if normal cound=0 and pcl count>0 then use pcl.
: if normal count>0 and pcl count=0 then use normal.
: if both counts=0 then mark missing.

*  Histo plot bin boundaries:
CTP – 0, 180, 310, 440, 560, 680, 800, 1100
COT– 0,  (0.3), 1.3, 3.6, 9.4, 23.0, 60.0, (100.0), 150.0

* Lons and Lats
Lons (XDim): -179.5 to 179.5
Lats (YDim): 89.5 to -89.5 (Should be flipped)

** mod: Terra (10:30am)
** myd: Aqua  (01:30pm)

"""

import sys
import numpy as np
import os  #.path
import glob
from subprocess import call
from datetime import timedelta, date
from pyhdf.SD import SD, SDC

def daterange(start_date, end_date):
    ### Including end date
    for n in range(int((end_date - start_date).days)+1):
        yield start_date + timedelta(n)

def get_hdf_data(hf,vnm):
    sds_obj = hf.select(vnm)
    vdata = sds_obj.get()
    att = sds_obj.attributes()
    offset = att['add_offset']
    scf = att['scale_factor']
    vdata = (vdata-offset)*scf
    return vdata


###------ Start ----
nx=360; ny=180
undef = -99.9

### upto end_date (both including start and end dates)
start_date = date(2002,12,1)  ### 2002-12-01
end_date = date(2016,11,30)   ### 2017-06-30

dset = ['MOD','MYD']; dset_id = 1   ### Set index for Aqua or Terra
dset2 = ['terra','aqua']
dsetnm = '{}08_D3'.format(dset[dset_id])

vnames = ['Cloud_Optical_Thickness_ISCCP_JHisto_vs_Pressure',
          'Cloud_Retrieval_Fraction_Combined',
          'Cloud_Optical_Thickness_PCL_ISCCP_JHisto_vs_Pressure',
          'Cloud_Retrieval_Fraction_PCL_Combined',]

indir = '/your_pass/MODIS_C6/Raw_hdf/{}/'.format(dsetnm)  ### Directory for MODIS raw hdf files
outdir = '/your_output_pass/'
outfn = '{}_d3-c6_tvp_pcl.{}_{}dy.cf0'.format(dset2[dset_id],start_date.strftime('%m%d%Y'),(end_date-start_date).days+1)


for oneday in daterange(start_date,end_date):
    dd=oneday.strftime('%Y%m%d')
    jday=oneday.timetuple().tm_yday

    infn = glob.glob(indir+'{}/{:03d}/'.format(dd[:4],jday)+'{}.A{}{:03d}.006.*.hdf'.format(dsetnm,dd[:4],jday))
    if len(infn) != 1:
        print(infn)
        sys.exit('# of file should be 1.')

    hfile = SD(infn[0],SDC.READ)
#    dset_dic=hfile.datasets()
#    print([x for x in dset_dic if 'ISCCP' in x])
#    sds_obj = hfile.select('YDim')
#    print(sds_obj.get())
    chist = get_hdf_data(hfile,vnames[0])  ### Should be [7,8,180,360]
    cf    = get_hdf_data(hfile,vnames[1])  ### Should be [180,360]
    chist_pcl = get_hdf_data(hfile,vnames[2])
    cf_pcl    = get_hdf_data(hfile,vnames[3])

    p,q,n,m = chist.shape
    print(chist[:,:,156,357],chist_pcl[:,:,156,357]);
    ### Screening clear sky first
    cf0idx = cf==0
    chist[:,:,cf0idx]=0.
    ### Screening missings
    idx_chist = chist<0.
    idxtsum = idx_chist.sum(axis=0).sum(axis=0)
    mask_chist = np.logical_or(idxtsum==p*q,cf<0.)  ### Missings, Now [180,360], Only for whole missings
    chist[idx_chist]=0.   ### Necessary for filling partial missing with 0
    cf[mask_chist]=0.
    ### Transform from count to cloud fraction
    tsum = chist.sum(axis=0).sum(axis=0)
    tsum_0idx = tsum>0.
    chist[:,:,tsum_0idx]=chist[:,:,tsum_0idx]/tsum[None,None,tsum_0idx]*cf[None,None,tsum_0idx]

    ### Repeat for PCL clouds
    cf0idx_pcl = cf_pcl==0
    chist_pcl[:,:,cf0idx_pcl]=0.
    idx_chist_pcl = chist_pcl<0.
    idxtsum_pcl = idx_chist_pcl.sum(axis=0).sum(axis=0)
    mask_chist_pcl = np.logical_or(idxtsum_pcl==p*q,cf_pcl<0.)  ### Missings, Now [180,360], Only for whole missings
    chist_pcl[idx_chist_pcl]=0.  ### Necessary for filling partial missing with 0
    cf_pcl[mask_chist_pcl]=0.
    tsum_pcl = chist_pcl.sum(axis=0).sum(axis=0)
    tsum_pcl_0idx = tsum_pcl>0.
    chist_pcl[:,:,tsum_pcl_0idx]=chist_pcl[:,:,tsum_pcl_0idx]/tsum_pcl[None,None,tsum_pcl_0idx]*cf_pcl[None,None,tsum_pcl_0idx]

    ### Sum with PCL
    chist+=chist_pcl
    ### Mask for both missings
    dualmask = np.logical_and(mask_chist,mask_chist_pcl)  ### Should be [180,360]

    ### Change dimension to ISCCP format
    chist[:,1,:,:]+=chist[:,0,:,:]
    chist[:,-2,:,:]+=chist[:,-1,:,:]
    chist=chist[:,1:-1,:,:]  ### Should be [7,6,180,360]
    ## Flipping along latitude
    chist=chist[:,:,::-1,:]; dualmask=dualmask[::-1,:]
    chist=np.swapaxes(np.reshape(chist,[42,ny,nx,1]),0,3).reshape([ny,nx,42]) ### Should be [180,360,42]

    ### Setting missing values
    chist[dualmask]=undef
    print(infn[0],dualmask.sum())

    ### Writing to binary file
    with open(outdir+outfn,'ab') as f:
        chist.astype(np.float32).tofile(f)
