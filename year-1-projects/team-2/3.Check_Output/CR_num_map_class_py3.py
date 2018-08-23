import os.path
import sys
import numpy as np

class get_crnums(object):
    """
    Assign each input data to one of centroid, and write the assigned Cluster numbers

    Input data: [nday,nlat,nlon,nelem], float32
    Output file: [nday,nlat,nlon], int16

    parameters
    ----------
    nelem: dimension of individual input vector. For cloud histogram, nelem=42
    ncl: # of clusters
    nday, nlat, nlon: dimensions of input data

    Attributes
    ----------


    """

    def __init__(self,nday=1,nlat=180,nlon=360,ncl=5,nelem=42):
        self.nday=nday
        self.nlat=nlat
        self.nlon=nlon
        self.ncl=ncl
        self.nelem=nelem

    def file_existence(self,fname,existence):
        """
        existence=True: if not available, it's error
        existence=False: if exists, it's error
        """

        if os.path.isfile(fname) != existence:
            if existence:
                msg="File does not exist: "+fname
            else:
                msg="File already exists: "+fname
            sys.exit(msg)
        return

    def read_bin_data(self,fname,dtp=np.float32):
        """
        Open a binary file, and read data

        fname : file name
        dtp   : data type; np.float32 or np.float64, etc.
        """

        self.file_existence(fname,True)
        with open(fname,'rb') as fd:
            bin_mat = np.fromfile(file=fd,dtype=dtp)

        return bin_mat



    def crnums_main(self,ctd,indata):
        """
        For each indata, decide the closet centroid, and assign crnum
        """

        crnums=np.full([self.nday,self.nlat,self.nlon],-1,dtype=np.int16)
        for k in range(self.nday):
            tmpdata=indata[k,:,:,:]
            tmpcf=tmpdata.sum(axis=2)
            idx_cf0= tmpcf==0.
            crnums[k,idx_cf0]=0

            idx_cf= np.logical_and(tmpcf>0.,tmpcf<1.000001)
            crnums[k,idx_cf]=self._assign2cr(tmpdata[idx_cf,:],ctd)

            if k%100==0:
                print('nday={}'.format(k+1))

        return crnums

    def _assign2cr(self,tgt,ctd):
        """
        in: tgt[some length of record,nelem], ctd[ncl,nelem]
        out: crs[some length of record], dtype=int16

        Find the index of closest cluster centroid.
        """

        nrec=tgt.shape[0]
        sqdist=np.empty([self.ncl,nrec])
        for k in range(self.ncl):
            sqdist[k,:]=np.power((tgt-ctd[k,:].reshape([1,-1])),2).sum(axis=1)

        crs=sqdist.argmin(axis=0)+1
        return crs.astype(np.int16)

    def savefile(self,fname,outdata):

        with open(fname,'wb') as f:
            outdata.tofile(f)
        return
