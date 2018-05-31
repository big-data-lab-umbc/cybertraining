from numpy import zeros,empty,copy
cimport cython
cimport openmp
# from cython.view cimport cvarray
from libc.stdlib cimport malloc, free
from libc.math cimport pow, sqrt
# For multithreading
from cython.parallel cimport prange
from cython.parallel cimport parallel
# For quitting if need be
import sys
print("Using Cythonized K-means")

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
    cdef double res = calc_dist_simp(arr_1, arr_2, elem)
    return res

# Used in assign and get newsum
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void calculate_cl(double [:,::1] indata, double [:,::1] ctd, long [::1] cl, int ncl, int nrec, int nk, int nelem) nogil:
    cdef:
        int ii
        int kk
        double mindd = 1.e5
        double tmpdd = 1.e5
        int idx = ncl

        # double mindd
        # double tempdd
    for ii in prange(0, nrec, nk, schedule='static', nogil=True):
        mindd=1.e5
        idx=ncl
        for kk in range(ncl):
            # tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
            tmpdd = calc_dist_noslice(indata,ii,ctd,kk,nelem)
            # Simple, Safe, Fast Squaring
            tmpdd = tmpdd * tmpdd
            if (tmpdd < mindd):
                mindd=tmpdd
                idx=kk
        if (idx == ncl):
            with gil:
                sys.exit()
        cl[ii]=idx
    return

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void calculate_outsum(double [:,::1] indata, long [::1] cl, double [:,::1] outsum, int nrec, int nk, int nelem) nogil:
    cdef:
        int jj
        int ii
        int clj
    # for jj in range(0,nrec,nk):
    # cdef compatible way since the range method is sketchy
    for jj from 0 <= jj < nrec by nk:
        for ii in range(nelem):
            # outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)
            clj = cl[jj] 
            outsum[clj,ii] = outsum[clj,ii] + indata[jj,ii]
    return

@cython.boundscheck(False)
@cython.wraparound(False)
def assign_and_get_newsum(indata,ctd,nk):
    cdef int nelem = indata.shape[1]
    cdef int nrec = indata.shape[0]
    cdef int ncl = ctd.shape[0]
    # Predefined variables for efficiency and multithreading

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


    # 
    # OpenMP is wrapped in this function
    calculate_cl(indata_mview, ctd_mview, cl_mview, ncl, nrec, nk, nelem)


    # !!!--- Sum for New Centroid
    calculate_outsum(indata, cl, outsum, nrec, nk, nelem)
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
    stopRec  = startRec + l_nrec
    return startRec, stopRec
