from numpy import zeros,empty,linalg
from math import sqrt
import sys

def calc_dist(arr1,arr2):
    # This should be a BLAS operation
    return linalg.norm(arr1-arr2)

def assign_and_get_newsum(indata,ctd,nk):
    nelem = indata.shape[0]
    nrec = indata.shape[1]
    ncl = ctd.shape[1]
    cl = empty(nrec,dtype=int)
    cl.fill(ncl)
    outsum=zeros(shape=(nelem,ncl))

    for ii in range(0,nrec,nk):
        mindd=1e5
        idx=ncl
        for kk in range(ncl):
            # tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
            # Use BLAS here
            tmpdd = calc_dist(indata[:,ii],ctd[:,kk]) ** 2
            if (tmpdd < mindd):
                mindd=tmpdd
                idx=kk
        if (idx == ncl):
            print("Not assigned",idx,ii,indata[:,ii])
            sys.exit()
        cl[ii]=idx

    # !!!--- Sum for New Centroid
    # This should be done separately
    for ii in range(0,nelem):
        # do jj=1,nrec,nk
        for jj in range(0,nrec,nk):
           # outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)

           outsum[ii,cl[jj]]=outsum[ii,cl[jj]]+indata[ii,jj]
    return cl,outsum

def get_wcv_sum(indata,ctd,cl):
    nelem, nrec = indata.shape
    ncl = ctd.shape[1]
    outsum=zeros(nelem)
    for mm in range(nelem):
       for ii in range(nrec):
          outsum[mm,cl[ii]]=outsum[mm,cl[ii]]+(indata[mm,ii]-ctd[mm,cl[ii]])**2
