module kmeans_omp_module
  !$ use omp_lib
  implicit none

contains

  function set_num_threads(nn) result(nt)
    integer :: nn,nt
    !F2PY INTENT(IN) :: nn
    !F2PY INTENT(OUT) :: nt
    
    !$ nt = omp_get_max_threads()
    print*,"Max # of Threads=",nt
    !$ call omp_set_num_threads(nn)
    !$ nt = omp_get_max_threads()
    print*, "Set Max # of Threads=",nt
  end function set_num_threads

  real(8) function calc_sqdist(aa,bb,nn)
  !!! Calculate squared sum of distance
    integer :: nn,ii
    real(8), dimension(nn),intent(in) :: aa,bb
    !F2PY INTENT(HIDE) :: nn
    !F2PY INTENT(OUT) :: calc_sqdist

    calc_sqdist=0.d0
    do ii=1,nn
       calc_sqdist=calc_sqdist+(aa(ii)-bb(ii))**2
    enddo

    return
  end function calc_sqdist

  real(8) function calc_dist(aa,bb,nn) 
    integer :: nn,ii
    real(8),intent(in) :: aa(nn),bb(nn) 
!    real(8) :: calc_sqdist
    !F2PY INTENT(HIDE) :: nn
    !F2PY INTENT(OUT) :: calc_dist

    calc_dist=sqrt(calc_sqdist(aa,bb,nn))

    return
  end function calc_dist

  subroutine assign_and_get_newsum(indata,ctd,nk,cl,outsum,ncl,nelem,nrec)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
    integer :: nelem,nrec,nk,ncl
    integer :: ii,jj,kk,idx
    integer, intent(out) :: cl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    real(8) :: mindd,tmpdd
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: cl,outsum
    !F2PY INTENT(IN) :: nk

    outsum=0.d0
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(indata,ctd,cl,outsum,ncl,nk,nrec,nelem)

    !!!--- Assigning Clusters
    !$OMP DO
    do ii=1,nrec,nk
       mindd=10000.d0;idx=0
       do kk=1,ncl
          tmpdd=calc_sqdist(indata(:,ii),ctd(:,kk),nelem)
          if (tmpdd.lt.mindd) then
             mindd=tmpdd; idx=kk
          endif
       enddo
       if (idx.eq.0) then
          print*,"Not assigned",idx,ii,indata(:,ii)
          stop
       endif
       cl(ii)=idx
    enddo
    !$OMP END DO
    !$OMP BARRIER
    !!!--- Sum for New Centroid

    !$OMP DO
    do ii=1,nelem
       do jj=1,nrec,nk
          outsum(ii,cl(jj))=outsum(ii,cl(jj))+indata(ii,jj)
       enddo
    enddo
    !$OMP END DO
    !print*, outsum,cluster(1000:1100)
    !$OMP END PARALLEL

  end subroutine assign_and_get_newsum

  subroutine get_wcv_sum(indata,ctd,cl,outsum,ncl,nelem,nrec)
  !!! Calculate sum of data by clusters
  !!! Need to get new centroid
    integer :: nelem,nrec,ncl
    integer :: mm,ii
    integer, intent(in) :: cl(nrec)
    real(8), intent(in) :: indata(nelem,nrec),ctd(nelem,ncl)
    real(8), intent(out) :: outsum(nelem,ncl)
    !F2PY INTENT(HIDE) :: ncl,nelem,nrec
    !F2PY INTENT(OUT) :: outsum

    outsum=0.d0
    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(outsum,cl,indata,ctd,nrec,nelem)
    do mm=1,nelem
       do ii=1,nrec
          outsum(mm,cl(ii))=outsum(mm,cl(ii))+(indata(mm,ii)-ctd(mm,cl(ii)))**2
       enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine get_wcv_sum


end module kmeans_omp_module
