from numpy import zeros,empty,copy
cimport cython
# from cython.view cimport cvarray
from libc.stdlib cimport malloc, free
from libc.math cimport pow, sqrt
import sys
print("Using Cythonized K-means")

# BLAS!
cdef extern from 'cblas.h':
    double ddot 'cblas_ddot'(int n, double* x, int incx, double* y, int incy) nogil
    void daxpby 'cblas_daxpby'(const int n, const double a, const double *x, const int incx, const double b, double *y, const int incy) nogil
    void dcopy 'cblas_dcopy'(const int n, const double *x, const int incx, double *y, const int incy) nogil
    double dnrm2 'cblas_dnrm2'(const int n, const double *x, const int incx) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double calc_dist_simp(double [:] arr1, double  [:] arr2, long l) nogil:
    cdef double d
    cdef long i
    for i in range(l):
        d += (arr1[i] - arr2[i])*(arr1[i] - arr2[i])
    d = sqrt(d)
    return d


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double calc_dist_noslice(double [:,::1] arr1, long ii, double [:,::1] arr2, long kk, long l) nogil:
    cdef double d = 0.0
    cdef double tmp
    cdef long i
    for i in range(l):
        # tmp = arr1[i,ii]-arr2[i,kk]
        tmp = arr1[ii,i]-arr2[kk,i]
        d += tmp*tmp

    d = sqrt(d)
    return d

@cython.boundscheck(False)
@cython.wraparound(False)
def calc_dist(arr1, arr2):
    # Get two memory views
    cdef double [:] arr_1 = arr1
    cdef double [:] arr_2 = arr2
    cdef long elem = len(arr1)
    # cdef double res = calc_dist_blas(arr_1, arr_2,elem)
    cdef double res = calc_dist_simp(arr_1, arr_2, elem)
    # print(res)
    return res


@cython.boundscheck(False)
@cython.wraparound(False)
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
    cdef int nelem = indata.shape[1]
    cdef int nrec = indata.shape[0]
    cdef int ncl = ctd.shape[0]
    # Predefined variables for efficiency
    cdef int idx = ncl
    cdef double mindd = 1e5
    cdef double tempdd = 1e6

    cl = empty(nrec,dtype=int)
    cl.fill(ncl)
    outsum=zeros(shape=(ncl,nelem),order='C')
    
    cdef long clj = 0
    cdef long [::1] cl_mview = cl
    cdef double [:,::1] outsum_mview = outsum
    cdef double [:,::1] indata_mview = indata
    cdef double [:,::1] ctd_mview = ctd
    # print("cl.shape = {}, ctd.shape = {}, indata.shape = {}, outsum.shape = {}".format(
                # cl.shape, ctd.shape, indata.shape, outsum.shape))



    for ii in range(0,nrec,nk):
        mindd=1e5
        idx=ncl
        for kk in range(ncl):
            # tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
            # Use BLAS here
            # tmpdd = calc_dist_simp(indata_mview[:,ii],ctd_mview[:,kk],nelem)
            tmpdd = calc_dist_noslice(indata_mview,ii,ctd_mview,kk,nelem)
            # Simple, Safe, Fast Squaring
            tmpdd *= tmpdd
            if (tmpdd < mindd):
                mindd=tmpdd
                idx=kk
        if (idx == ncl):
            # print("Not assigned",idx,ii,indata[:,ii])
            print("Not assigned",idx,ii,indata[ii,:])
            sys.exit()
        cl_mview[ii]=idx

    # !!!--- Sum for New Centroid
    # This should be done separately
    # for jj in range(0,nrec,nk):
    for jj in range(0,nrec,nk):
        # do jj=1,nrec,nk
        for ii in range(nelem):
            # outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)
            clj = cl_mview[jj] 
            # outsum_mview[ii,clj] +=outsum_mview[ii,clj] + indata_mview[ii,jj]
            # outsum_mview[ii,clj] += indata_mview[ii,jj]
            outsum_mview[clj,ii] += indata_mview[jj,ii]
    return cl,outsum

@cython.boundscheck(False)
@cython.wraparound(False)
def get_wcv_sum(indata,ctd,cl):
    cdef int nelem = indata.shape[1]
    cdef int nrec = indata.shape[0]
    cdef int ncl = ctd.shape[0]
    cdef int ii,mm
    cdef long cli
    outsum=zeros(shape=(ncl,nelem),order='C')
    cdef long [:] cl_mview = cl
    cdef double [:,::1] ctd_mview = ctd
    cdef double [:,::1] outsum_mview = outsum
    cdef double [:,::1] indata_mview = indata
    cdef double tmp
    for ii in range(nrec):
        for mm in range(nelem):
            cli = cl_mview[ii]
            # outsum_mview[mm,cli]=pow(outsum_mview[mm,cli]+(indata_mview[mm,ii]-ctd_mview[mm,cli]), 2)
            tmp = (indata_mview[ii,mm] - ctd_mview[cli,mm])
            outsum_mview[cli,mm] += tmp*tmp
    return outsum

@cython.cdivision(True)
def get_record_spans(long nrec, int rank, int tprocs):
    cdef:
        long l_nrec, rem, startRec, stopRec

    # Calculate the total number of records for each process
    l_nrec = nrec // tprocs
    rem = nrec %  tprocs
    if (rem == 0):
        startRec = l_nrec*rank 
    else:
        if (rank < rem):
            # Pick up an extra record
            l_nrec = l_nrec + 1
            startRec = rank*l_nrec
        else:
            # Accounts for additional records
            startRec = l_nrec*rank + rem
    return startRec, stopRec
