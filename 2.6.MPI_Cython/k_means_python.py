from numpy import zeros,empty
from math import sqrt
def calc_sqdist(arr1,arr2,sz):
  # !!! Calculate squared sum of distance
   # integer :: nn,ii
    # real(8), dimension(nn),intent(in) :: aa,bb
    # sz = arr1.shape[0]
    calc_sqdist=0.0
    # do ii=1,sz
    # print("ARR1 SHAPE SQDIST = {}\n\tSIZE SQDIST = {}".format(arr1.shape,sz))
    # print("ARR2 SHAPE SQDIST = {}\n\tSIZE SQDIST = {}".format(arr2.shape,sz))
    for ii in range(sz):
       calc_sqdist=calc_sqdist+(arr1[ii]-arr2[ii])**2
    # enddo


    return calc_sqdist

def calc_dist(arr1,arr2):
    # integer :: nn,ii
    # real(8),intent(in) :: aa(nn),bb(nn) 
# !    real(8) :: calc_sqdist
    # !F2PY INTENT(HIDE) :: nn
    # !F2PY INTENT(OUT) :: calc_dist
    sz = arr1.shape[0]
    calc_dist=sqrt(calc_sqdist(arr1,arr2,sz))

    return calc_dist

# def assign_and_get_newsum(indata,ctd,nk,cl,outsum,ncl,nelem,nrec):
def assign_and_get_newsum(indata,ctd,nk):
  # !!! Calculate sum of data by clusters
  # !!! Need to get new centroid
    # integer :: nelem,nrec,nk,ncl
    # integer :: ii,jj,kk,idx
    # integer, intent(out) :: cl(nrec)
    # real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    # real(8), intent(out) :: outsum(nelem,ncl)
    # real(8) :: mindd,tmpdd
    # !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    # !F2PY INTENT(OUT) :: cl,outsum
    # !F2PY INTENT(IN) :: nk
    # print("indata.shape = {}\n\tctd.shape = {}".format(indata.shape, ctd.shape))
    nelem = indata.shape[0]
    nrec = indata.shape[1]
    ncl = ctd.shape[1]
    cl = empty(nrec,dtype=int)
    cl.fill(-1)
    # nrec=indata.shape[0]/self.nelem
    # print("NELEM = {}".format(nelem))
    outsum=zeros(shape=(nelem,ncl))

    # !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(indata,ctd,cl,outsum,ncl,nk,nrec,nelem)

    # !!!--- Assigning Clusters
    # !$OMP DO
    # do ii=1,nrec,nk
    for ii in range(1,nrec,nk):
       # idx=-1 because python indexes at 0
       mindd=10000.0;idx=-1
       # do kk=1,ncl
       for kk in range(1,ncl):
          # tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
          # tmpdd=calc_sqdist(indata[:][ii],ctd[:][kk],nelem)
          tmpdd=calc_sqdist(indata[:,ii],ctd[:,kk],nelem)
          if (tmpdd < mindd):
             mindd=tmpdd; idx=kk
          # endif
       # enddo
       if (idx ==-1):
          print("Not assigned",idx,ii,indata[:,ii])
          # stop
          break
       # endif
       cl[ii]=idx
    # enddo
    # !$OMP END DO
    # !$OMP BARRIER
    # !!!--- Sum for New Centroid

    # !$OMP DO
    # do ii=1,nelem
    for ii in range(1,nelem):
       # do jj=1,nrec,nk
       for jj in range(1,nrec,nk):
          # outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)

          outsum[ii,cl[jj]]=outsum[ii,cl[jj]]+indata[ii,jj]
       # enddo
    # enddo
    # !$OMP END DO
    # !print*, outsum,cluster(1000:1100)
    # !$OMP END PARALLEL
    return cl,outsum

  # end subroutine assign_and_get_newsum

def get_wcv_sum(indata,ctd,cl,outsum,ncl,nelem,nrec):
  # !!! Calculate sum of data by clusters
  # !!! Need to get new centroid
    # integer :: nelem,nrec,ncl
    # integer :: mm,ii
    # integer, intent(in) :: cl(nrec)
    # real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    # real(8), intent(out) :: outsum(nelem,ncl)
    # !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    # !F2PY INTENT(OUT) :: outsum

    outsum=zeros(nelem)
    # !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(outsum,cl,indata,ctd,nrec,nelem)
    # do mm=1,nelem
    for mm in range(1,nelem):
       # do ii=1,nrec
       for ii in range(1,nrec):
          # outsum(mm,cl(ii))=outsum(mm,cl(ii))+(indata(mm,ii)-ctd(mm,cl(ii)))**2
          outsum[mm,cl[ii]]=outsum[mm,cl[ii]]+(indata[mm,ii]-ctd[mm,cl[ii]])**2
       # enddo
    # enddo
    # !$OMP END PARALLEL DO
  # end subroutine get_wcv_sum


# end module kmeans_omp_module
