from numpy import zeros,empty,copy
cimport cython
# from cython.view cimport cvarray
from libc.stdlib cimport malloc, free, pow
import sys
print("Using Cythonized K-means")

# BLAS!
cdef extern from 'cblas.h':
    double ddot 'cblas_ddot'(int n, double* x, int incx, double* y, int incy) nogil
    void daxpby 'cblas_daxpby'(const int n, const double a, const double *x, const int incx, const double b, double *y, const int incy) nogil
    void dcopy 'cblas_dcopy'(const int n, const double *x, const int incx, double *y, const int incy) nogil
    double dnrm2 'cblas_dnrm2'(const int n, const double *x, const int incx) nogil





def calc_dist(arr1, arr2):
    # Get two memory views
    cdef double [:] arr_1 = arr1
    cdef double [:] arr_2 = arr2
    cdef long elem = len(arr1)
    cdef double res = calc_dist_blas(arr_1, arr_2,elem)
    # print(res)
    return res

# This should be a BLAS operation
@cython.boundscheck(False)
cdef double calc_dist_blas(double [:] arr1, double [:] arr2, long elem) nogil:
    # Buffer space
    cdef double *d = <double *>malloc(elem * sizeof(double))
    dcopy(elem, &arr2[0], 1, d, 1)
    daxpby(elem, 1.0, &arr1[0], 1 , -1.0, d, 1)

    cdef double res = dnrm2(elem,d,1)
    free(d)
    # return linalg.norm(arr1-arr2)
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def assign_and_get_newsum(indata,ctd,nk):
    cdef int ii,kk,jj
    # cdef int kk
    # cdef int jj
    cdef int nelem = indata.shape[0]
    cdef int nrec = indata.shape[1]
    cdef int ncl = ctd.shape[1]
    # Predefined variables for efficiency
    cdef int idx = ncl
    cdef double mindd = 1e5
    cdef double tempdd = 1e6

    cl = empty(nrec,dtype=int)
    cl.fill(ncl)
    outsum=zeros(shape=(nelem,ncl))
    
    cdef long clj = 0
    cdef long [:] cl_mview = cl
    cdef double [:,:] outsum_mview = outsum
    cdef double [:,:] indata_mview = indata
    cdef double [:,:] ctd_mview = ctd

    for ii in range(0,nrec,nk):
        mindd=1e5
        idx=ncl
        for kk in range(ncl):
            # tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
            # Use BLAS here
            tmpdd = calc_dist_blas(indata_mview[:,ii],ctd_mview[:,kk],nelem)
            # Simple, Safe, Fast Squaring
            tmpdd *= tmpdd
            if (tmpdd < mindd):
                mindd=tmpdd
                idx=kk
        if (idx == ncl):
            print("Not assigned",idx,ii,indata[:,ii])
            sys.exit()
        cl_mview[ii]=idx

    # !!!--- Sum for New Centroid
    # This should be done separately
    for ii in range(nelem):
        # do jj=1,nrec,nk
        for jj in range(0,nrec,nk):
            # outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)
            clj = cl_mview[jj] 
            # outsum_mview[ii,clj] +=outsum_mview[ii,clj] + indata_mview[ii,jj]
            outsum_mview[ii,clj] += indata_mview[ii,jj]
    return cl,outsum

def get_wcv_sum(indata,ctd,cl):
    nelem, nrec = indata.shape
    ncl = ctd.shape[1]
    outsum=zeros(nelem)
    for mm in range(nelem):
       for ii in range(nrec):
          outsum[mm,cl[ii]]=outsum[mm,cl[ii]]+(indata[mm,ii]-ctd[mm,cl[ii]])**2
