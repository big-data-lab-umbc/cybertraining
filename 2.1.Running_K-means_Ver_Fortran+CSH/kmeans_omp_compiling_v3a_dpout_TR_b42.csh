#!/bin/tcsh

###
### This is a script to run K-means clustering.
###
### Input data: cloud 2D joint histogram data. No missings at all.
###             dimension= [42 vars(nelements), nvectot] with real(4)
###             ## From Python: [nvectot,nvar] with float32
### Output data: centroid result
###             dimension= [42 vars(nelements), ncluster(knum)] with double precision (real(8))
###             ## For Python: [ncluster,nvar] with float
###
### Note:
### 1. Initializing: Starting data points are randomly selected while distances between each starting data points are kept large enough.
### 2. Speed boosting: For initial a few rounds, only part of data are used in order to get plausible centroid faster.
### 3. The loop will stop when the centroid movement is less than epslion value.
### 4. Now the epsilon is quite small, a.k.a, tight version
### 5. openmp is used.

set knums = ( 4 4 ) # Number of clusters (k number)
set sids  = ( 1 1 ) # IDs for each k numbers (number of trials)

set dset = MODIS
set mdid = $dset"_Aqua_b42_TR"

# Directory of input data file
set dir0 = ""
set nvectot = 3445612
set infn = $dir0/aqua_d3_c6_tvp_pcl.noMissing.20050101-20051231_{$nvectot}x42.float32.dat

# Directory for output centroid file
set outdir = "./CTD" #$dir1
set outctd = $mdid

set fn = $mdid.2005_dpout.$knums[1]ompx

set threads_num = 6   # Number of threads for openmp

cat > $fn.f90 << EOF
program kmeans_clustering
  implicit none
  integer, parameter :: knst=$knums[1],knen=$knums[2],sidst=$sids[1],siden=$sids[2]
  integer, parameter ::nelements=42,niter=999, nvectot=$nvectot

  real(8), parameter :: EPSILON=0.000001d0 !0.00001d0
  real :: nmapgd,ntm
  integer :: ncl,id
  real(8) :: tmpdd,tmpcent(nelements),mindd ! Initial mindd should be large enough
  real(8) :: wcv,calc_dist,calc_sqdist,calc_wcv
  integer :: iloop,nk,nrec
  integer :: ii,jj,kk,mm,idx,ict,io1

  integer, allocatable :: cluster(:)
  real, allocatable :: rdvar(:,:)
  real(8), allocatable :: centroid(:,:),move(:),vsum(:,:)
  character*150 :: ctdfnm
  character*2 :: char_k,char_sid
  character*1 :: char_k1
  logical :: lgc,test_dist

!!!
!!! Allocate array and loading data
!!!
  allocate(rdvar(nelements,nvectot))
  nrec=nvectot/4*nelements   ! Dividing reading size of input file if the input file is too large
  print*,nrec,nvectot,nelements,nvectot*nelements
  open(unit=11 &
      ,file='$infn' &
      ,iostat=io1,status='old',access='direct',recl=nrec*4)
  if ( io1.ne.0 ) then
     print*,"Error: data file is not read.",io1
     print*,"$infn"
     stop
  endif

  kk=1
  do jj=1,4
     print*,jj,kk
     read(11,rec=jj) rdvar(:,kk:kk+nvectot/4-1)
     kk=kk+nvectot/4
  enddo

  print*,"File reading is done. ",kk

!!!
!!! Loop: K-number first, SID second, and iloop inner-most
!!!
  do ncl=knst,knen
     do id=sidst,siden

        !!! Initializing
        print*,"Start: K=",ncl,"SID=",id
        allocate(centroid(nelements,ncl),move(ncl))
        allocate(cluster(nvectot),vsum(nelements,ncl))

        nmapgd=360*30*0.87  !!! Estimated average grid points per day
        ntm=17.      !!! Unit days: odd number

        idx=int(nmapgd*(ntm+(id+13.)/real(id+37.)))
        centroid(:,1)=dble(rdvar(:,idx))
        do ii=2,ncl
           lgc=.true.
           do while (lgc)
              idx=idx+idx
              if (idx.gt.nvectot) then
                 idx=idx-nvectot
                 print*,"idx is over total record"
              endif
              tmpcent(:)=dble(rdvar(:,idx))
              !!! Test if distance is large enough
              lgc=test_dist(nelements,tmpcent,ii-1,centroid(:,1:ii-1))
          enddo
          centroid(:,ii)=tmpcent(:)
        enddo
        print*,"Initial assignment of centroid is done."


        print*,"Loop starts"
        iloop=1
        do while (iloop.le.niter)
           print*,"*****",iloop

           !!! Using smaller samples for intial roudns, every nk'th data
           if (iloop.le.10) then
              nk=8
           elseif (iloop.le.20) then
              nk=4
           elseif (iloop.le.32) then
              nk=2
           else
              nk=1
           endif

           vsum=0.d0
!\$OMP PARALLEL DEFAULT(PRIVATE) SHARED(rdvar,centroid,cluster,vsum,ncl,nk)

!!!--- Assigning Clusters
!\$OMP DO
           do ii=1,nvectot,nk
              mindd=10000.d0;idx=0
              do kk=1,ncl
                 tmpdd=calc_sqdist(nelements,dble(rdvar(1:nelements,ii)),centroid(1:nelements,kk))
                 if (tmpdd.lt.mindd) then
                    mindd=tmpdd; idx=kk
                 endif
              enddo
              if (idx.eq.0) then
                 print*,"Not assigned",idx,ii,rdvar(:,ii)
                 stop
              endif
              cluster(ii)=idx
           enddo
!\$OMP END DO

!!!--- New Centroid
!\$OMP DO
           do mm=1,nelements
              do ii=1,nvectot,nk
                 vsum(mm,cluster(ii))=vsum(mm,cluster(ii))+dble(rdvar(mm,ii))
              enddo
           enddo
!\$OMP END DO
           !print*, vsum,cluster(1000:1100)
!\$OMP END PARALLEL

           move=0.d0
           do kk=1,ncl
              ict=COUNT(cluster==kk)
              vsum(:,kk)=vsum(:,kk)/dble(ict)
              move(kk)=calc_dist(nelements,centroid(:,kk),vsum(:,kk))
              print*, "*", kk,real(move(kk))
           enddo
           centroid(:,:)=vsum(:,:)

!!! check how much cluster sizes change from one loop to the next
           if ((MAXVAL(move).le.epsilon).and.(nk.eq.1)) then
              print*, "*** Converged",iloop
              exit
           endif
           iloop=iloop+1
        enddo

        if (iloop.gt.niter) then
           print*,"!!!*** Not converged ***!!!"
           print*,"* KN = ",ncl," SID = ",id," WCV = N/A"
        else
           wcv=calc_wcv(ncl,id,nelements,nvectot,cluster,centroid,rdvar)
           print*,"K=",ncl," SID=",id," WCV= ",wcv
        endif

!!! write out results
        call Look_up_num(2,ncl,char_k)
        call Look_up_num(2,id,char_sid)
        OPEN(unit=21 &
            ,file="$outdir/$outctd.cent_km"//char_k//"_sid"//char_sid//".dpdat" &
            ,status='unknown',access='direct',recl=nelements*ncl*8)
        write(21,rec=1) centroid
        CLOSE(21)

        deallocate(centroid,move,vsum,cluster)
     enddo
  enddo
  deallocate(rdvar)
end program kmeans_clustering

logical function test_dist(nn,tmp,mm,ctd)
!!! Test if the distance is larger or smaller than std
!!! Return True if distance is smaller than std
  implicit none
  real, parameter :: std=0.125
  integer :: nn,mm,kk
  real*8 :: tmp(nn),ctd(nn,mm),calc_dist

  test_dist=.false.
  do kk=1,mm
     if (calc_dist(nn,tmp,ctd(:,kk)).lt.std) then
        test_dist=.true.
        exit
     endif
  enddo
  return
end function test_dist

double precision function calc_sqdist(nn,aa,bb)
!!! Calculate squared sum of distance
  implicit none
  integer :: nn,ii
  real(8) :: aa(nn),bb(nn),dd

  dd=0.d0
  do ii=1,nn
     dd=dd+(aa(ii)-bb(ii))**2
  enddo
  calc_sqdist=dd
  return
end function calc_sqdist

double precision function calc_dist(nn,aa,bb)
  implicit none
  integer :: nn,ii
  real(8) :: aa(nn),bb(nn),dd,calc_sqdist

  calc_dist=sqrt(calc_sqdist(nn,aa,bb))
  return
end function calc_dist

double precision function calc_wcv(knum,sid,nelem,nrec,cl,ctd,rdvec)
!!! Calculate within-cluster-variation
  implicit none
  integer :: knum,sid,nelem,nrec
  integer :: tmprec,ii,i2,jj,kk,mm,ct,ict(knum)
  integer :: cl(nrec),newclnum(knum)
  real :: rdvec(nelem,nrec)
  real(8) :: ctd(nelem,knum)
  real(8) :: vecsum(nelem,knum),ww(knum)


  vecsum=0.d0

!\$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(vecsum,cl,rdvec,ctd,nrec,nelem)
  do mm=1,nelem
     do ii=1,nrec
        vecsum(mm,cl(ii))=vecsum(mm,cl(ii))+(dble(rdvec(mm,ii))-ctd(mm,cl(ii)))**2
     enddo
  enddo
!\$OMP END PARALLEL DO

  !!! Sort the centroid results
  newclnum=0
  call sort_centroid(nelem,knum,ctd,newclnum)

  !!! Within Cluster Variation
  calc_wcv=0.d0; ww=0.d0
  do kk=1,knum
     do mm=1,nelem
        ww(kk)=ww(kk)+vecsum(mm,kk)
     enddo
     calc_wcv=calc_wcv+ww(kk)
  enddo

  ct=0
  do kk=1,knum
     ii=kk
     ict(kk)=COUNT(cl==ii)
     ct=ct+ict(kk)
     print*, "*", kk, ict(kk), ww(ii)
  enddo

  print*,"count[last]+wcv/pop:",knum,sid,ict(knum),calc_wcv/dble(ct)

  return
end function calc_wcv

subroutine Look_up_num(nn,intnum,chars)
  implicit none
  integer :: nn,intnum
  character :: chars(nn)
  character :: num(0:9)
  integer :: ii,inum

  data num /"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"/

  inum=intnum
  do ii=nn,1,-1
     chars(ii)=num(mod(inum,10))
     inum=int(inum/10.0)
  enddo
end subroutine Look_up_num

subroutine sort_centroid(nn,nk,ctd,inum)
!!!
!!! Sorting centroid by cloud top height and thickness
!!!
!!! Obs data
!!! Assume that nelements=42=6*7
!!! Order: 17,27,37,...,67,16,...,61
!!!
  implicit none
  integer :: nn,nk
  real(8)  :: ctd(nn,nk),tmp(nn,nk),centsum(nk)
  real    :: wgt(nn),wctd(nk),wcsorted(nk)
  integer :: inum(nk),idx,irec
  integer :: ii,jj,kk,mm,iost1,nvectors
  logical :: mask(nn)

  !!! Weight
  do kk=1,nn
     ii=7-mod(kk-1,6)
     jj=7-int((kk-1)/6.0)
     wgt(kk) = ii*5+jj*50
  enddo

  !!! Apply weight
  wctd=0.d0; centsum=0.d0
  do jj=1,nk
     do ii=1,nn
        centsum(jj)=centsum(jj)+ctd(ii,jj)
     enddo
     mask=.true.
     wctd(jj)=wgt(MAXLOC(ctd(:,jj),1,mask))
     mask(MAXLOC(ctd(:,jj),1,mask))=.false.
     wctd(jj)=wctd(jj)+wgt(MAXLOC(ctd(:,jj),1,mask))*0.7
!     print*,jj,wctd(jj)
  enddo

!!! Find minimum cent_sum, and send it to the last (lowest weighted value)
    irec=MINLOC(centsum,1)
    if (centsum(irec).lt.0.5d0) wctd(irec)=wctd(irec)/100.0

!!! Sorting high to low
  wcsorted(1)=wctd(1)
  inum(1)=1;idx=1
  do jj=2,nk
     kk=1
     do while (kk.le.idx)
        if (wctd(jj).gt.wcsorted(kk)) exit
        kk=kk+1
     enddo

     if (kk.le.idx) then
        wcsorted(kk+1:idx+1)=wcsorted(kk:idx)
        inum(kk+1:idx+1)=inum(kk:idx)
     endif
     wcsorted(kk)=wctd(jj)
     inum(kk)=jj

     idx=idx+1
  enddo

  print*,"Sort is completed"

   do jj=1,nk
     tmp(:,jj)=ctd(:,inum(jj))
  enddo
  ctd=tmp
  return
end subroutine sort_centroid

EOF

#ifort -o $fn $fn.f90 -mp1 -openmp -check bounds || exit
gfortran -o $fn $fn.f90 -fopenmp -fbounds-check || exit
#gfortran -o $fn $fn.f90 -fbounds-check || exit

#rm $fn.f90
setenv OMP_NUM_THREADS $threads_num
./$fn || exit 0
rm $fn

exit 0
