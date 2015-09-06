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
   use mpi_gamil

   implicit none

   real(r8) :: pmtop
   integer  :: begj, endj
   integer  :: i,    jdyn, jcam, k

!---------------------------------------------------------------------
!!     write(*,*) '!!!! c_a_switching----new'
     begj = beglatexdyn + numbnd
     endj = endlatexdyn - numbnd

!!--- 3-d vars----

    call gamil_arrays_comm(COMM_TO_LEFT,1,u(:,beglatexdyn,1)) 

!$OMP PARALLEL DO PRIVATE (k, jdyn, jcam, i)
     do k=1,plev
       do jdyn = begj,endj

          jcam = plat + 1 - jdyn

         do i=beglonex,iend2
            u3 (i, k,  jcam,n3) =  0.5*( u(i,jdyn,k)+u(i+1,jdyn,k) )
            t3 (i, k,  jcam,n3) =  t (i,jdyn,k)
            q3 (i, k,1,jcam,n3) =  q (i,jdyn,k)
           omga(i,      k,  jcam)    = wpa(i,jdyn,k)*100.0
         enddo
       enddo
     enddo

!--- 2-d vars----

!$OMP PARALLEL DO PRIVATE (jdyn, jcam, i)
       do jdyn = begj,endj

          jcam = plat + 1 - jdyn

         do i=beglonex,endlonex
            ps(i,jcam,n3) = (pes(i,jdyn)+ pmtop)*100.0
         enddo
       enddo

    call gamil_arrays_comm(COMM_TO_BOT,1,v(:,beglatexdyn,1)) 

!$OMP PARALLEL DO PRIVATE (k, jdyn, jcam, i)
     do k=1,plev
       do jdyn=jbeg0,jend0
         jcam=plat+1-jdyn
         if (jdyn .eq. 1 .or. jdyn .eq. plat) then
           do i=beglonex,endlonex
             v3 (i,k,jcam,n3) = 0.0
           enddo
         else
           do i=beglonex,iend2
             v3 (i,k,jcam,n3) = -0.5*( v(i,jdyn-1,k)+v(i,jdyn,k) )
           enddo
         endif
       enddo
     enddo

     return
end subroutine c_a_switching
