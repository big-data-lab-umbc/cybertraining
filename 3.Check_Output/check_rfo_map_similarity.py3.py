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
cent=rm.read_bin_data(inctd,dtp=np.float64).reshape([ncl,nelem])
cf=cent.sum(axis=1)*100.
print(cf)

### Read CRnum file
crnums=rm.read_bin_data(infncrnum,dtp=np.int16).reshape([nday,nlat,nlon])
rfomap=rm.calc_rfomap(crnums)
print(rfomap.shape)

### Calculate correlations of rfo maps
#rfomap2=rfomap.reshape([ncl,nlat/2,2,nlon/3,3]).sum(axis=4).sum(axis=2)
print("RFO Map Similarity")
corrs,rmsd,pairs=rm.corr_rmsd_mtx_pairs(rfomap.reshape([ncl,-1]))

print("Corr. Coef: Mean={:.3f}, Median={:.3f}".format(corrs.mean(),np.median(corrs)))
print("RMSD: Mean={:.5f}, Median={:.5f}".format(rmsd.mean(),np.median(rmsd)))
idx=np.argsort(corrs)[::-1]
print("Highest Corr, Pair, RMSD")
for i in range(int(ncl/2.5)):
    print("{}: {:.3f}, {}, {:.5f}".format(i+1,corrs[idx[i]],pairs[idx[i]],rmsd[idx[i]]))

### Calculate correlations of centroids
print("Centroid Similarity")
corrs,rmsd,pairs=rm.corr_rmsd_mtx_pairs(cent)

print("Corr. Coef: Mean={:.3f}, Median={:.3f}".format(corrs.mean(),np.median(corrs)))
print("RMSD: Mean={:.5f}, Median={:.5f}".format(rmsd.mean(),np.median(rmsd)))
idx=np.argsort(corrs)[::-1]
print("Highest Corr, Pair, RMSD")
for i in range(int(ncl/2.5)):
    print("{}: {:.3f}, {}, {:.5f}".format(i+1,corrs[idx[i]],pairs[idx[i]],rmsd[idx[i]]))
