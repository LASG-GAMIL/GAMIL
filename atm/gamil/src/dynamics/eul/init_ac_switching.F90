#include <misc.h>
#include <params.h>

!!(wanhui 2003.10.28-29)
!!(wanhui 2003.11.24)
!!-----------------------


subroutine init_ac_switching( pmtop )

!!---------------------------------------------------------------------------
!!
!! Purpose : 
!!    switch from A-grid to C-grid
!!
!!---------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond,plat,plev,i1,numbnd,beglatexdyn,endlatexdyn,beglatex,endlatex
   use prognostics
   use comfm1

#if ( defined SPMD )
   use mpishorthand, only: mpicom
#endif
   
   implicit none 

#if (defined SPMD)
#include <mpif.h>
#include <commpi.h>

   real(r8) :: workst(plond*plev),workrb(plond*plev)
   integer  :: isend, irecv
   integer  :: istatus(mpi_status_size)
   integer  :: ii
#endif
   
   real(r8) :: pmtop
   integer  :: begj, endj
   integer  :: i,    jdyn, jcam, k
   real(r8) :: sumn,sums
    
!!--------------------------------------------------------------------------

     begj = beglatexdyn + numbnd
     endj = endlatexdyn - numbnd
!!(wh 2004.02.17)

#if ( ! defined SPMD )
    do k=1,plev
       sumn=0.0
       sums=0.0
       do i=2,plond-1
          sumn=sumn+t3(i,k,1,n3m2)
          sums=sums+t3(i,k,plat,n3m2)
       enddo
       sumn=sumn/dble(plond-2)
       sums=sums/dble(plond-2)
       do i=1,plond
          t3(i,k,1,n3m2)=sumn
          t3(i,k,plat,n3m2)=sums
       enddo
    enddo

#else
   if (myrank.eq.0 ) then
    jcam = beglat
    do k=1,plev
       sums=0.0
       do i=2,plond-1
          sums=sums+t3(i,k,jcam,n3m2)
       enddo
       sums=sums/dble(plond-2)
       do i=1,plond
          t3(i,k,jcam,n3m2)=sums
       enddo
    enddo
   endif
!--
   if (myrank.eq.nprocs-1) then
    jcam = endlat
    do k=1,plev
       sumn=0.0
       do i=2,plond-1
          sumn=sumn+t3(i,k,jcam,n3m2)
       enddo
       sumn=sumn/dble(plond-2)
       do i=1,plond
          t3(i,k,jcam,n3m2)=sumn
       enddo
    enddo
   endif
!--
#endif
!!(wh)


#if (defined SPMD)

      ii = 1
      do k=1,plev
       do i=1,plond
          workst(ii) = v3(i,k,endlatex-1,n3m2)
          ii = ii+1
       enddo
      enddo
      call mpi_isend(workst, plond*plev, mpi_double_precision, itop,1, mpicom, isend, ierr) 
      call mpi_irecv(workrb, plond*plev, mpi_double_precision, ibot,1, mpicom, irecv, ierr)

#endif


!!--- u ----

     do k = 1,plev
      do jdyn = begj,endj
         jcam = plat + 1 - jdyn

         do i=2,plond-2
            u (i,jdyn,k) = 0.5*( u3(i1+i-2, k,jcam,n3m2)+u3(i1+i-1,k,jcam,n3m2) )
         enddo
            u (1,jdyn,k) = 0.5*( u3(plond-1,k,jcam,n3m2)+u3(i1,k,jcam,n3m2) )

            u (plond-1,jdyn,k) = u(1,jdyn,k)
            u (plond,  jdyn,k) = u(2,jdyn,k)
      enddo
     enddo

!!-- t,q ----

     do k = 1,plev
      do jdyn = begj,endj
         jcam = plat + 1 - jdyn

         do i=1,plond-2
            t (i,jdyn,k) = t3 (i1+i-1,k,  jcam,n3m2)
            q (i,jdyn,k) = q3 (i1+i-1,k,1,jcam,n3m2)
         enddo
            t (plond-1,jdyn,k) = t (1,jdyn,k)
            t (plond  ,jdyn,k) = t (2,jdyn,k)
            q (plond-1,jdyn,k) = q (1,jdyn,k)
            q (plond  ,jdyn,k) = q (2,jdyn,k)
           
      enddo
     enddo

!!-- pes,ghs ----

      do jdyn = begj,endj
         jcam = plat + 1 - jdyn

         do i=1,plond-2
            pes (i,jdyn) = ps  (i,jcam,n3m2)*0.01d0 - pmtop
            ghs (i,jdyn) = phis(i,jcam)
         enddo
            pes(plond-1,jdyn) = pes (1,jdyn)     
            pes(plond  ,jdyn) = pes (2,jdyn)     
            ghs(plond-1,jdyn) = ghs (1,jdyn)     
            ghs(plond  ,jdyn) = ghs (2,jdyn)     
      enddo

!!-- v ----

#if (defined SPMD)

      call mpi_wait( isend,istatus,ierr )
      call mpi_wait( irecv,istatus,ierr )

      ii = 1
      do k=1,plev
       do i=1,plond
          v3(i,k,beglatex,n3m2)=workrb(ii)
          ii = ii+1
       enddo
      enddo

#endif

#if (!defined SPMD)

     do k=1,plev
      do jdyn=begj,endj-1
         jcam = plat+1-jdyn
         
         do i=1,plond-2
            v(i,jdyn,k) = -0.5*( v3(i1+i-1,k,jcam,n3m2)+v3(i1+i-1,k,jcam-1,n3m2) )
         enddo
            v(plond-1,jdyn,k) = v(1,jdyn,k)
            v(plond,  jdyn,k) = v(2,jdyn,k)
      enddo

         do i=1,plond
            v(i,endj,k) = 0.0
         enddo
     enddo
          
#else

   if ( myrank == 0 ) then

     do k=1,plev
      do jdyn=begj,endj-1
         jcam = plat+1-jdyn
         
         do i=1,plond-2
            v(i,jdyn,k) = -0.5*( v3(i1+i-1,k,jcam,n3m2)+v3(i1+i-1,k,jcam-1,n3m2) )
         enddo
            v(plond-1,jdyn,k) = v(1,jdyn,k)
            v(plond,  jdyn,k) = v(2,jdyn,k)
      enddo

         do i=1,plond
            v(i,endj,k) = 0.0
         enddo
     enddo

   else
          
     do k=1,plev
      do jdyn=begj,endj
         jcam = plat+1-jdyn
         
         do i=1,plond-2
            v(i,jdyn,k) = -0.5*( v3(i1+i-1,k,jcam,n3m2)+v3(i1+i-1,k,jcam-1,n3m2) )
         enddo
            v(plond-1,jdyn,k) = v(1,jdyn,k)
            v(plond,  jdyn,k) = v(2,jdyn,k)
      enddo
     enddo
   
   endif

#endif      

     return
end subroutine init_ac_switching
