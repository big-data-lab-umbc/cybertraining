from rfo_map_display_class_py3 import rfo_map_display as rmdisp
from datetime import date
import numpy as np

ncl=8   ### Number of clusters or K number
nelem=42

### CR_num file info
idate=date(2005,1,1)  ## Start date
edate=date(2005,12,31)  ## End date
nday=(edate-idate).days+1

lons_new = np.arange(-179.5,180.,1.); nlon=len(lons_new)
lats_new = np.arange(-14.5,15.,1.); nlat=len(lats_new)
lons2d, lats2d = np.meshgrid(lons_new,lats_new)

indir = '/directory for CR_num file/'
infncrnum = indir+'aqua_CRnum_map.MODISc6_b42_CR{:02d}.{}-{}_{}x{}.int16dat'.format(ncl,idate.strftime('%Y%m%d'),edate.strftime('%Y%m%d'),nlat,nlon)

### Centroid file info
indir2 = '/directory for centroid/'
inctd = indir2+'MODIS_T+A_b42.cent_km{:02d}_sid01.dpdat'.format(ncl)

### Define object
rm=rmdisp(ncl=ncl)

### Read centroid and get CF
ctd=rm.read_bin_data(inctd,dtp=np.float64).reshape([ncl,nelem])
cf=ctd.sum(axis=1)*100.
print(cf)

### Read CRnum file
crnums=rm.read_bin_data(infncrnum,dtp=np.int16).reshape([nday,nlat,nlon])
rfomap=rm.calc_rfomap(crnums)
print(rfomap.shape)

### Basic figure info
ncol=2; nrow=4  ### parameter to set up figure layout

psize=[ncol*4+1,nrow*1.+1.5]  ## [lx,ly]
lon_limit=[0.,360.1]
lat_limit=[-20.,20.]

rm.set_parameters(nrow,ncol,lons2d,lats2d,lon_limit,lat_limit)

### Draw figure
suptit="K={} Clustering RFO [MODIS Aqua TR 1yr]".format(ncl)
rm.map_display_main(rfomap,cf,page_size_inches=psize,suptit=suptit)

### Save the figure
outdir="./"
fnout = outdir+"MODIS_T+A_b42.cent_km{:02d}_sid01.RFO.png".format(ncl)
rm.savefig(fnout)
