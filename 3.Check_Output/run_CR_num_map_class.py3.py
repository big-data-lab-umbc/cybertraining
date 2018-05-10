from CR_num_map_class_py3 import get_crnums as crnum
from datetime import date
import numpy as np

idate=date(2005,1,1)  ## Start date
edate=date(2005,12,31)  ## End date
nday=(edate-idate).days+1

nx=360; ny=30
nelem=42
ncl=8

indir = '/input_data_location/'
infn = indir+'aqua_d3_c6_tvp_pcl.{}-{}_{}x{}x{}.float32.dat'.format(idate.strftime('%Y%m%d'),edate.strftime('%Y%m%d'),ny,nx,nelem)

indir = './CTD/'
inctd = indir+'MODIS_T+A_b42.cent_km{:02d}_sid01.dpdat'.format(ncl)

outdir = '/output_data_location/'
outfn = outdir+'aqua_CRnum_map.MODISc6_b42_CR{:02d}.{}-{}_{}x{}.int16dat'.format(ncl,idate.strftime('%Y%m%d'),edate.strftime('%Y%m%d'),ny,nx)

## Define object
crn=crnum(nday=nday,nlat=ny,nlon=nx,ncl=ncl,nelem=nelem)

## Check existence of out file
crn.file_existence(outfn,False)  ### Outfile shouldn't exist yet.

## Read centroid
ctd=crn.read_bin_data(inctd,dtp=np.float64)
ctd=ctd.reshape([ncl,nelem])

## Read input data
indata=crn.read_bin_data(infn)
indata=indata.reshape([nday,ny,nx,nelem])

## Call Main Calculation unit
outcrs=crn.crnums_main(ctd,indata)

## Save the result
crn.savefile(outfn,outcrs)
