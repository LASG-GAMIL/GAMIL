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
   use mpi_gamil

  
   implicit none 

   real(r8) :: pmtop
   integer  :: i,    jdyn, jcam, k
   real(r8) :: u3_v3_tmp(beglonex:endlonex,beglatexdyn:endlatexdyn,plev)
   real(r8) :: tmp_3d(plond,plat,plev*4)
    
!!--------------------------------------------------------------------------


   call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,u3_v3_tmp(:,beglatexdyn,1))
   call register_comm_array(beglonex,endlonex,beglatex,endlatex,1,plev,1,1,u3(:,1,beglatex,n3m2),is_phys_array=.true.)
   call gamil_arrays_comm(COMM_ROTATE_LEFT,2,u3(:,1,beglatex,n3m2))

     do jdyn=jbeg0, jend0
       if (jdyn .eq. 1 .or. jdyn .eq. plat) then
         jcam = plat-jdyn+1
         call gamil_average_pole_data_phys(jdyn,beglonex,endlonex,plev,t3(beglonex,1,jcam,n3m2))
       endif
     enddo


   do k=1, plev
     do jdyn = jbeg0,jend0
       jcam = plat + 1 - jdyn
       do i=beglonex,endlonex
         u3_v3_tmp(i,jdyn,k) = u3(i,k,jcam,n3m2)
       enddo
     enddo
   enddo

      call gamil_arrays_comm(COMM_TO_RIGHT,1,u3_v3_tmp(:,beglatexdyn,1)) 

     do k = 1,plev
      do jdyn = jbeg0,jend0
         jcam = plat + 1 - jdyn
         do i=ibeg1,iend1
            u (i,jdyn,k) = 0.5*( u3_v3_tmp(i-1,jdyn,k)+u3_v3_tmp(i,jdyn,k) )
         enddo

!         u (1,jdyn,k) = u(plon+1,jdyn,k)
!         u (plond,  jdyn,k) = u(2,jdyn,k)
      enddo
     enddo

     call gamil_arrays_comm(COMM_TO_LEFT,1,u(:,beglatexdyn,1))
     call gamil_arrays_comm(COMM_TO_RIGHT,1,u(:,beglatexdyn,1))

!!-- t,q ----

     do k = 1,plev
      do jdyn = jbeg0,jend0
         jcam = plat + 1 - jdyn

         do i=beglonex,endlonex
            t (i,jdyn,k) = t3 (i,k,  jcam,n3m2)
            q (i,jdyn,k) = q3 (i,k,1,jcam,n3m2)
         enddo
      enddo
     enddo

     call gamil_arrays_comm(COMM_ROTATE_LEFT,2,t(:,beglatexdyn,1),q(:,beglatexdyn,1))

!!-- pes,ghs ----

      do jdyn = jbeg0,jend0
         jcam = plat + 1 - jdyn

         do i=beglonex,endlonex
            pes (i,jdyn) = ps  (i,jcam,n3m2)*0.01d0 - pmtop
            ghs (i,jdyn) = phis(i,jcam)
         enddo
      enddo

     call gamil_arrays_comm(COMM_ROTATE_LEFT,2,pes(:,beglatexdyn),ghs(:,beglatexdyn))

!!-- v ----

      do k=1, plev
      do jdyn = jbeg0,jend0
         jcam = plat + 1 - jdyn
         do i=beglonex,endlonex
           u3_v3_tmp(i,jdyn,k) = v3(i,k,jcam,n3m2)
         enddo
      enddo
      enddo
      call gamil_arrays_comm(COMM_TO_TOP,1,u3_v3_tmp(:,beglatexdyn,1)) 

      do k=1,plev
        do jdyn=jbeg0,jend0
          if (jdyn .eq. plat) then
            do i=beglonex,endlonex
              v(i,jdyn,k) = 0.0
            enddo
          else
            jcam = plat+1-jdyn
            do i=beglonex,endlonex
              v(i,jdyn,k) = -0.5*( u3_v3_tmp(i,jdyn,k)+u3_v3_tmp(i,jdyn+1,k) )
            enddo
          endif
        enddo
      enddo

     call gamil_arrays_comm(COMM_ROTATE_LEFT,2,v(:,beglatexdyn,1))
     call remove_comm_array(u3(:,1,beglatex,n3m2))
     call remove_comm_array(u3_v3_tmp(:,beglatexdyn,1))

     return
end subroutine init_ac_switching
