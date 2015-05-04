#include <misc.h>
#include <params.h>


subroutine a_c_switching( fu, fv, t2, beglat, endlat)

!!----------------------------------------------------------------------------------
!!  accumulate su,sv,st and update q    (wanhui 2003.10.28-29)
!!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plond, plat, plev, i1, numbnd, beglatex, endlatex
   use comfm1,       only: su, sv, st, q
   use prognostics,  only: qminus

#if ( defined SPMD )
   use mpishorthand, only: mpicom
#endif

   implicit none

#if (defined SPMD)
#include <mpif.h>
#include <commpi.h>

   real(r8) :: fvtmp(plond,plev)
   integer  :: isend, irecv
   integer  :: istatus(mpi_status_size)
#endif
   
   integer  :: beglat, endlat
   real(r8) :: fu (plond,plev,beglat:endlat)
   real(r8) :: fv (plond,plev,beglat:endlat)
   real(r8) :: t2 (plond,plev,beglat:endlat)

   real(r8) :: qsums,qsumn,t2sums,t2sumn

   integer  :: begj, endj
   integer  :: i,    jdyn, jcam, k
   
!!-----------------------------------------------------------------------------

     begj = beglat
     endj = endlat

#if (defined SPMD)

      call mpi_isend( fv(1,1,endlat), plond*plev, mpi_double_precision, &
                                         itop, 1, mpicom, isend, ierr) 
      call mpi_irecv( fvtmp(1,1),     plond*plev, mpi_double_precision, &
                                         ibot, 1, mpicom, irecv, ierr)
#endif

!!(wh 2004.02.25)

#if ( ! defined SPMD )
    do k=1,plev
       qsumn=0.0
       qsums=0.0
       t2sumn=0.0
       t2sums=0.0
       do i=1,plond-2
          qsumn=qsumn+qminus(i1+i-1,k,1,endlat)
          qsums=qsums+qminus(i1+i-1,k,1,beglat)
          t2sumn=t2sumn+t2(i,k,endlat)
          t2sums=t2sums+t2(i,k,beglat)
       enddo
       qsumn=qsumn/dble(plond-2)
       qsums=qsums/dble(plond-2)
       t2sumn=t2sumn/dble(plond-2)
       t2sums=t2sums/dble(plond-2)
       do i=1,plond
          qminus(i,k,1,endlat)=qsumn
          qminus(i,k,1,beglat)=qsums
          t2(i,k,endlat)=t2sumn
          t2(i,k,beglat)=t2sums
       enddo
    enddo
#else

    if (myrank.eq.0) then
     jcam = beglat
     do k=1,plev
       qsums=0.0
       t2sums=0.0
       do i=1,plond-2
          qsums=qsums+qminus(i1+i-1,k,1,jcam)
          t2sums=t2sums+t2(i,k,jcam)
       enddo
       qsums=qsums/dble(plond-2)
       t2sums=t2sums/dble(plond-2)
       do i=1,plond
          qminus(i,k,1,jcam)=qsums
          t2(i,k,jcam)=t2sums
       enddo
     enddo
    endif
!--
   if (myrank.eq.nprocs-1) then
    jcam = endlat
    do k=1,plev
       qsumn=0.0
       t2sumn=0.0
       do i=1,plond-2
          qsumn=qsumn+qminus(i1+i-1,k,1,jcam)
          t2sumn=t2sumn+t2(i,k,jcam)
       enddo
       qsumn=qsumn/dble(plond-2)
       t2sumn=t2sumn/dble(plond-2)
       do i=1,plond
          qminus(i,k,1,jcam)=qsumn
          t2(i,k,jcam)=t2sumn
       enddo
     enddo
    endif
!--
#endif
!!(wh)


     do k=1,plev

!--su--

        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=2,plon
            su(i,jdyn,k) = su(i,jdyn,k) + 0.5*( fu(i-1,k,jcam)+fu(i,k,jcam) )
         enddo
            su(1,jdyn,k) = su(1,jdyn,k) + 0.5*( fu(plon,k,jcam)+fu(1,k,jcam) )

            su(plon+1,jdyn,k)= su(1,jdyn,k)
            su(plon+2,jdyn,k)= su(2,jdyn,k)
        enddo

!--st--

        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=1,plon
            st(i,jdyn,k)= st(i,jdyn,k) + t2(i,k,jcam)
         enddo
            st(plon+1,jdyn,k)= st(1,jdyn,k)
            st(plon+2,jdyn,k)= st(2,jdyn,k)
        enddo
!--q--
 
        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=1,plon
            q(i,jdyn,k) = qminus(i1+i-1,k,1,jcam)
         enddo
            q(plon+1,jdyn,k) = q(1,jdyn,k)
            q(plon+2,jdyn,k) = q(2,jdyn,k)
        enddo

!------------

     enddo

!--sv--

#if (!defined SPMD)

     do k=1,plev
        do jcam=2,plat
		  
          jdyn = plat + 1 - jcam

         do i=1,plon
            sv(i,jdyn,k)= sv(i,jdyn,k)-0.5*( fv(i,k,jcam)+fv(i,k,jcam-1) )
         enddo
            sv(plon+1,jdyn,k) = sv(1,jdyn,k)
            sv(plon+2,jdyn,k) = sv(2,jdyn,k)
        enddo

         do i=1,plond
            sv(i,plat,k) = 0.0
         enddo
     enddo
#else

    begj = beglat + 1

    do k=1,plev
      do jcam=begj,endj

         jdyn = plat + 1 - jcam

         do i=1,plon
            sv(i,jdyn,k)= sv(i,jdyn,k)-0.5*( fv(i,k,jcam)+fv(i,k,jcam-1) )
         enddo
            sv(plon+1,jdyn,k) = sv(1,jdyn,k)
            sv(plon+2,jdyn,k) = sv(2,jdyn,k)
      enddo
    enddo
!!-------------------------------  for the lat next to the south boundary

    call mpi_wait(isend,istatus,ierr)
    call mpi_wait(irecv,istatus,ierr)

    jdyn = plat+1-beglat                   !! (jdyn = endlatexdyn-1)

    if (myrank.eq.0) then                  !! (jdyn = the south pole)

       do k=1,plev
        do i=1,plond
           sv(i,jdyn,k) = 0.0
        enddo
       enddo

    else                                   !! (jdyn = endlatexdyn-1)
       do k=1,plev
        do i= 1,plon
           sv(i,jdyn,k) = sv(i,jdyn,k)-0.5*( fv(i,k,beglat)+fvtmp(i,k) )
        enddo
           sv(plon+1,jdyn,k) = sv(1,jdyn,k)
           sv(plon+2,jdyn,k) = sv(2,jdyn,k)
		    
       enddo
    endif
!!------------------------

#endif

     return
end subroutine a_c_switching
