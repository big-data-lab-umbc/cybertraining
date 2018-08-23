"""
#
# Read the centroid files, calculating RMSD, and grouping them by the RMSD value
# for clustering results of multi-trials
#
# Daeho Jin, Jun. 30, 2014
#   Revised on 2017.05.08: For MODIS TR data
"""

import numpy as np
import os.path
import sys

###----Funcitons

def bin_file_read2mtx(fname,dtp=np.float32):
    """ Open a binary file, and read data
        fname : file name
        dtp   : data type; np.float32 or np.float64, etc. """

    if not os.path.isfile(fname):
        print("File does not exist:"+fname)
        sys.exit()

    with open(fname,'rb') as fd:
        bin_mat = np.fromfile(file=fd,dtype=dtp)

    return bin_mat

def rmsd(lside, rside): return np.sqrt(((lside - rside) ** 2).mean())
def rmsd_mean(a,b,ncl):
    """
    Calculate RMSD between two sets of centroids [ncl,nelem]
    Since we can not be sure the order of centroid, calculate all possible
    paris and collect the one with minimum rmsd value.
    """
    blist=range(ncl)
    sum=0.
    for i in range(ncl):
        rmsdlist=[rmsd(a[i,:],b[j,:]) for j in blist]
        rmin=min(rmsdlist)
        sum+=rmin
        del blist[rmsdlist.index(rmin)]
#    print "mean",sum/float(len(alist))
    return sum/float(ncl)

####--------------------------------------

if len(sys.argv) < 3:
    sys.exit("Need 2 inputs: km and total number of trials(IDs)")
km    = int(sys.argv[1])       # Number of Clusters
tsid  = int(sys.argv[2])      # Total number of different seeds (1 to tsid)

dir1 = './CTD/'
mdnm = 'Aqua_b42_TR'
#ctdfnm = 'MODIS_{}.cent_km{:02d}_sid{:02d}'.format(mdnm,km,sid)
#fnm=dir1+ctdfnm+'.dpdat'
ncl=km; nelem=42
rmsd_criterion=0.001

print("\n*** Start ***")
print("[mdnm,km,tsid]=",mdnm,km,tsid)

a=np.empty(shape=(tsid,ncl,nelem))
for i in range(tsid):
    ### Read Centroid File
    fnm =  dir1+'MODIS_{}.cent_km{:02d}_sid{:02d}.dpdat'.format(mdnm,ncl,i+1)
#    fnm =  'MODIS_{}.cent_km10_sid24.subCR.sub_k{:02d}_id{:02d}.dpdat'.format(mdnm[mid],ncl,i+1)
    mtx=bin_file_read2mtx(fnm,dtp=np.float64)
    a[i,:] = mtx.reshape([ncl,nelem])

ini_list=range(tsid)
group_idx=0;group_list=[]
point1=0
while (point1<tsid):
#    print "point1=",point1,ini_list
    one_group=[i for i in ini_list if rmsd_mean(a[point1,:,:],a[i,:,:],ncl) <= rmsd_criterion]
#    print "one_group",one_group
    for i in one_group: del ini_list[ini_list.index(i)]

    group_list.append(one_group)
    one_group=[]
    if len(ini_list)<=0 :
        point1=tsid+1
    else:
        point1=min(ini_list)

print(group_list)
nn=len(group_list)
print("\n***"+"-"*15)
print("Total %d groups over %d seeds" % (nn,tsid))
for i in range(nn):
    print("Group%d: %d members" % (i+1,len(group_list[i])))
    print([i+1 for i in group_list[i]])
