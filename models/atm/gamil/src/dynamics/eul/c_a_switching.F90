#include <misc.h>
#include <params.h>


subroutine c_a_switching( pmtop )

!---------------------------------------------------------------------
! prepare data for physics    (wanhui 2003.10.28-29
!---------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon,plond,plat,plev,i1,numbnd,beglatexdyn,endlatexdyn,beglat,endlat
   use comfm1,       only: u,  v,  t,  q,  wpa,  pes
   use prognostics,  only: u3, v3, t3, q3, omga, ps, n3

#if ( defined SPMD )
   use mpishorthand, only: mpicom
#endif

   implicit none

#if (defined SPMD)
#include <mpif.h>
#include <commpi.h>

   real(r8) :: worksb(plond*plev),workrt(plond*plev)
   integer  :: isend, irecv
   integer  :: istatus(mpi_status_size)
   integer  :: ii
#endif

   real(r8) :: pmtop
   integer  :: begj, endj
   integer  :: i,    jdyn, jcam, k

!---------------------------------------------------------------------
!!     write(*,*) '!!!! c_a_switching----new'
     begj = beglatexdyn + numbnd
     endj = endlatexdyn - numbnd

#if (defined SPMD)

      ii = 1
      do k=1,plev
       do i=1,plond
          worksb(ii) = v(i,endj,k)
          ii = ii+1
       enddo
      enddo
      call mpi_isend(worksb, plond*plev, mpi_double_precision, ibot,1, mpicom, isend, ierr) 
      call mpi_irecv(workrt, plond*plev, mpi_double_precision, itop,1, mpicom, irecv, ierr)

#endif

!!--- 3-d vars----

     do k=1,plev
       do jdyn = begj,endj

          jcam = plat + 1 - jdyn

         do i=1,plond-2
            u3 (i1+i-1, k,  jcam,n3) =  0.5*( u(i,jdyn,k)+u(i+1,jdyn,k) )
            t3 (i1+i-1, k,  jcam,n3) =  t (i,jdyn,k)
            q3 (i1+i-1, k,1,jcam,n3) =  q (i,jdyn,k)
           omga(i,      k,  jcam)    = wpa(i,jdyn,k)*100.0
         enddo
            u3 (1,      k,  jcam,n3) = u3 (plond-1,k,  jcam,n3)
            t3 (1,      k,  jcam,n3) = t3 (plond-1,k,  jcam,n3)
            q3 (1,      k,1,jcam,n3) = q3 (plond-1,k,1,jcam,n3)
           omga(plond-1,k,  jcam)   = omga(1,      k,  jcam)

            u3 (plond,  k,  jcam,n3) = u3 (2,      k,  jcam,n3)
            t3 (plond,  k,  jcam,n3) = t3 (2,      k,  jcam,n3)
            q3 (plond,  k,1,jcam,n3) = q3 (2,      k,1,jcam,n3)
           omga(plond,  k,  jcam)   = omga(2,      k,  jcam)
       enddo
     enddo

!--- 2-d vars----

       do jdyn = begj,endj

          jcam = plat + 1 - jdyn

         do i=1,plond
            ps(i,jcam,n3) = (pes(i,jdyn)+ pmtop)*100.0
         enddo
       enddo

!---------------------------------------------------------------------
#if ( ! defined SPMD )

     do k=1,plev

       do jdyn=2,plat-1
          jcam=plat+1-jdyn
         do i=1,plond-2
            v3 (i1+i-1,k,jcam,n3) = -0.5*( v(i,jdyn-1,k)+v(i,jdyn,k) )
         enddo
            v3 (1    ,k,jcam,n3) = v3 (plond-1,k,jcam,n3)
            v3 (plond,k,jcam,n3) = v3 (2      ,k,jcam,n3)
       enddo

         do i=1,plond
            v3 (i,k,1   ,n3) = 0.0
            v3 (i,k,plat,n3) = 0.0
         enddo

     enddo
#else

    if (myrank.eq.nprocs-1) then
       do k=1,plev
        do i=1,plond
          v3 (i,k,endlat,n3) = 0.0
        enddo
       enddo
       begj = beglatexdyn + numbnd +1
    elseif (myrank.eq.0 ) then
       do k=1,plev
        do i=1,plond
          v3 (i,k,beglat,n3) = 0.0
        enddo
       enddo
       endj = endlatexdyn - numbnd -1
    endif

!!------------------------------------unpack the received data
	     
    call mpi_wait(isend,istatus,ierr)
    call mpi_wait(irecv,istatus,ierr)

    ii = 1 
    if (myrank.ne.nprocs-1) then
       do k=1,plev
        do i=1,plond
           v(i,beglatexdyn,k) = workrt(ii)
           ii = ii+1
        enddo
       enddo
    endif
!!------------------------------------------------

      do k=1,plev
       do jdyn=begj, endj
          jcam= plat+1-jdyn
          do i=1,plond-2
              v3 (i1+i-1,k,jcam,n3) = -0.5*( v(i,jdyn-1,k)+v(i,jdyn,k) )
          enddo
              v3 (1    ,k,jcam,n3) = v3 (plond-1,k,jcam,n3)
              v3 (plond,k,jcam,n3) = v3 (2      ,k,jcam,n3)
       enddo
      enddo

#endif
		    

     return
end subroutine c_a_switching
