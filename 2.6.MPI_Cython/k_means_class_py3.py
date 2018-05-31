import os.path
import sys
import numpy as np
import k_means_cython as km_mod
import datetime
from mpi4py import MPI

class K_means:
    """ 
    K-means clustering importing fortran module

    $ f2py3 -c --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -m k_means_mod k-means_mod.f90

    parameters
    ----------
    domain_size: 2-dim list, [nlat, nlon], Needed for getting initial centroid
    nelem: dimension of individual input vector. For cloud histogram, nelem=42
    nrec: # of features or # of data points of input data
    epsilon: criterion if the clustering is converged

    Attributes
    ----------
    knum: # of clusters or K
    id_: ID for each trial with different initial codition


    """

    def __init__(self, domain_size=[180,360],nelem=42,nrec=-999,epsilon=1.e-6):
        self.domain_size=domain_size
        self.nelem=nelem
        self.nelem=nelem
        self.nrec=nrec
        self.epsilon=epsilon
        # MPI Standards
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.tprocs = self.comm.Get_size()

    def set_knum_id(self,knum,id_):
        self.knum=knum
        self.id_=id_ 
       
    def read_bin_data(self,fname,dtp=np.float32):
        """ 
        Open a binary file, and read data
        
        fname : file name
        dtp   : data type; np.float32 or np.float64, etc. 
        """
        if not os.path.isfile(fname):
            self.print("File does not exist:"+fname)
            sys.exit()

        with open(fname,'rb') as fd:
            bin_mat = np.fromfile(file=fd,dtype=dtp)
        return bin_mat

    def initialize(self,indata, num_threads=1):
        """
        Initialize:
        1. Initialize input data
        2. Set number of threads for OpenMP
        """

        indata=self._initialize_indata(indata)
        return indata

    def _initialize_indata(self,indata):
        """
        Initialize input data
        1. float32 => float64
        2. reshape
        """

        self.nrec=indata.shape[0]/self.nelem
        data = np.reshape(indata,newshape=[int(self.nrec),int(self.nelem)],order='C').astype(float)
        self.print(data.shape)
        self.startRec, self.stopRec = km_mod.get_record_spans(self.nrec, self.rank, self.tprocs)
        self.print("[{:03d}] {}::{}".format(self.rank,self.startRec,self.stopRec))
        return data

    def get_initial_ctd(self, indata, ini_ctd_dist_min=0.125):
        """ 
        Return initial centroid to start the K-means iteration

        After inital assignment, test if each centroid is far enough from others

        Parameters
        ----------
        ini_ctd_dist_min: test criterion for minimum distance between centroids
        ncount_per_day: approximated data ratio, available/(available+missing)
        ntodd: Unit days: odd number

        Input data
        ----------
        indata: original dataset to run K-means clustering, dimension=[nelements,nvectot]

        Output
        ------
        ctd: centroid, dimension(nelements,knum)  
        """

        ny,nx=self.domain_size
        ncount_per_day=nx*ny*0.87   # 0.87: approximated data ratio, available/(available+missing)
        ntodd=17.   # Unit days: odd number

        idx=int(ncount_per_day*(ntodd+(self.id_+13.)/(self.id_+37.)))
        ctd=[]

        ctd.append(indata[idx,:])
        while len(ctd)<self.knum:
            idx+=idx
            if idx>self.nrec:
                idx-=self.nrec
                self.print("idx is over total record")

            tmpctd=indata[int(idx),:]
            if self._test_dist(ctd,tmpctd,ini_ctd_dist_min):
                ctd.append(tmpctd)

        return np.asfarray(ctd)

    def _test_dist(self,group1,elem1,criterion):
        """
        From all elements in group, elem1 should be far enough (>criterion) 
        """

        test=True
        for elem0 in group1:
            distsq=np.sum(np.power(elem0-elem1,2))
            if distsq < criterion**2:
                test=False
                break

        return test

    def K_means_main(self,indata,ctd,iter_max=999):
        """
        Repeat loop until getting converged centroid
        """
        self.print("Start: K={}, ID={}".format(self.knum,self.id_))
        n10=3; nk=2**n10
        self.print("***** nk= {}".format(nk))
        totalTime = 0
        for it in range(iter_max):
            self.print("***** {}".format(it+1))
            startTime = MPI.Wtime()
            ### Assign data to centroid and get new sum (not mean yet)
            cl,outsum=km_mod.assign_and_get_newsum(indata,ctd, self.startRec, self.stopRec, nk)
            maxmove=0.
            cl_count= np.zeros(self.knum)
            for ic in range(self.knum):
                # Pure Python should use this and initialize with ncl
                idx= cl==ic

                # MPI Reduce for the partial sum of cl_count
                val = idx.sum()
                tmp_count = self.comm.allreduce(val, op=MPI.SUM)
                cl_count[ic] = tmp_count
                
                # Merge the partial sums
                fodder = outsum[ic,:]/float(cl_count[ic])
                tmpctd = np.empty(shape=fodder.shape,order='C')

                # MPI Reduce on tmpctd -> since cl==ncl doesn't matter, this works
                self.comm.Allreduce(fodder, tmpctd, op=MPI.SUM)

                # Parallel distance calculation? -> Data structure is too small to matter
                move=km_mod.calc_dist(tmpctd,ctd[ic,:])

                self.print("* {:02d} {}".format(ic+1,move))

                maxmove=max(maxmove,move)
    
                ctd[ic,:]=tmpctd
            # 4 Performance Tests!
            endTime = MPI.Wtime()
            thisTime = endTime-startTime
            totalTime += thisTime

            self.print("** {}".format(datetime.timedelta(seconds=(thisTime))))
            if n10>0 and maxmove < self.epsilon*10.**n10:
                ### Speeding up trick1
                ### : using only part of samples in initial stages
                n10-=1; nk=2**n10
                self.print("***** nk is changed to {}".format(nk))
            elif nk==1 and maxmove < self.epsilon:
                self.print("*** Converged ***",it+1)
                self.print("*** {} ***".format(datetime.timedelta(seconds=(totalTime))))
                break

        if it==iter_max-1:
            self.print("!!!*** Not Converged ***!!!")
            self.print("** Knum= {}, ID= {}, WCV= N/A".format(self.knum,self.id_))
        else:
            startTime = MPI.Wtime()
            # Wcv is a partial sum of only one set of records
            wcvsum  = km_mod.get_wcv_sum(indata,ctd,cl,self.startRec,self.stopRec)

            endTime = MPI.Wtime()
            thisTime = endTime - startTime
            # Merge the partial sum into a smaller sum
            fodder = wcvsum.sum(axis=1)
            wcv = np.empty(shape=fodder.shape)
            self.comm.Allreduce(fodder, wcv, op=MPI.SUM)
            wcv=wcv/cl_count
            
            cf=ctd.sum(axis=1)
            self.print("** Knum= {}, ID= {}, Total WCV= {}, LowestCF WCV={}, WCV Time= {}".format(
                self.knum,self.id_,wcv.sum(),wcv[np.argsort(cf)[0]],datetime.timedelta(seconds=thisTime)))

        return ctd

    def write_centroid(self,fnamehead,ctd,ftype='b'):
        """ 
        Sorting the centroid and then write to a file
        
        ftype='b': binary
        ftype='t': text
        
        """
        if self.rank == 0:
            ctd=ctd.T  #[knum,nelem]
            ctd=self._sort_centroid(ctd)
            self.print('Sorted_CF: ',ctd.sum(axis=1))

            fname=fnamehead+'.cent_k{:02d}_id{:02d}_{}x{}'.format(self.knum,self.id_,self.knum,self.nelem)
            if ftype=='b':
                with open(fname+'.float64_dat','wb') as fd:
                    ctd.tofile(fd)
            elif ftype=='t':
                np.savetxt(fname+'.txt',ctd,fmt='%.8f',delimiter=' ')

        return

    def _sort_centroid(self,nctd):
        """
        Sort the centroid
        
        Thick and high first, thin high second, and thin low last.
        The lowest CF one (less than 50%) is at the end.

        Input: centriod, dimension=[knum,nelem]
        Output: sorted centroid
        """
        ctd = nctd.T
        cf=ctd.sum(axis=1)
        idx= cf<0.5
        ctd2=ctd[~idx,:].reshape([-1,7,3,2]).sum(axis=3)
        ctd2[:,0,:]=ctd2[:,0:3,:].sum(axis=1)
        ctd2[:,1,:]=ctd2[:,3:5,:].sum(axis=1)
        ctd2[:,2,:]=ctd2[:,5:7,:].sum(axis=1)
        ctd2=ctd2[:,0:3,:].reshape([-1,9])

        wt=np.arange(1,10,1).reshape([3,3])[::-1,:].reshape(-1)
        wcf=np.average(ctd2,weights=wt,axis=1)
        ctd0=ctd[~idx,:][np.argsort(wcf)[::-1],:]

        if idx.sum()>0:
            xx=np.argsort(cf[idx])[::-1]
            ctd2=ctd[idx,:].reshape([-1,self.nelem])[xx,:]
            ctd0=np.concatenate((ctd0,ctd2))
        return ctd0
            
    def print(self, *args):
        """(str, [any]) -> None
        
        A lazy way to restrict MPI printing. Checks it's process rank and only prints
        if self.rank == 0."""
        if self.rank == 0:
            print(*args)
            sys.stdout.flush()

    def Allprint(self,*args):
        print("[{:03d}] ".format(self.rank),*args)
        sys.stdout.flush()
