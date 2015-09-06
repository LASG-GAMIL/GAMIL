#include <misc.h>
#include <params.h>


subroutine a_c_switching( fu, fv, t2, beglat, endlat)

!!----------------------------------------------------------------------------------
!!  accumulate su,sv,st and update q    (wanhui 2003.10.28-29)
!!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plond, plat, plev, numbnd, beglatex, endlatex, beglatexdyn, endlatexdyn
   use comfm1,       only: su, sv, st, q
   use prognostics,  only: qminus
   use mpi_gamil

   implicit none
   
   integer  :: beglat, endlat
   real(r8) :: fu (beglonex:endlonex,plev,beglat:endlat)
   real(r8) :: fv (beglonex:endlonex,plev,beglat:endlat)
   real(r8) :: t2 (beglonex:endlonex,plev,beglat:endlat)
   real(r8) :: fu_fv_tmp(beglonex:endlonex,beglatexdyn:endlatexdyn,plev)

   integer  :: begj, endj
   integer  :: i,    jdyn, jcam, k, j
   
!!-----------------------------------------------------------------------------

   call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,fu_fv_tmp(:,beglatexdyn,1))

     begj = beglat
     endj = endlat

     do jdyn=jbeg0, jend0
       if (jdyn .eq. 1 .or. jdyn .eq. plat) then
         jcam = plat-jdyn+1
         call gamil_average_pole_data_phys(jdyn,beglonex,endlonex,plev,t2(beglonex,1,jcam),qminus(beglonex,1,1,jcam))
      endif
     enddo


!$OMP PARALLEL DO PRIVATE (k, jcam, jdyn, i)
    do k=1, plev
    do jdyn = jbeg0,jend0
      jcam = plat + 1 - jdyn
      do i=beglonex,endlonex
        fu_fv_tmp(i,jdyn,k) = fu(i,k,jcam)
      enddo
    enddo
    enddo

    call gamil_arrays_comm(COMM_ROTATE_LEFT,2,fu_fv_tmp(:,beglatexdyn,1))
    call gamil_arrays_comm(COMM_TO_RIGHT,1,fu_fv_tmp(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (k, jcam, jdyn, i)
     do k=1,plev

!--su--

        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=ibeg1,iend1
            su(i,jdyn,k) = su(i,jdyn,k) + 0.5*( fu_fv_tmp(i-1,jdyn,k)+fu_fv_tmp(i,jdyn,k) )
         enddo
        enddo

!--st--

        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=beglonex,iend2
            st(i,jdyn,k)= st(i,jdyn,k) + t2(i,k,jcam)
         enddo
        enddo
!--q--
 
        do jcam = begj,endj

           jdyn = plat + 1 - jcam

         do i=beglonex,iend2
            q(i,jdyn,k) = qminus(i,k,1,jcam)
         enddo
        enddo

!------------

     enddo



!--sv--

!$OMP PARALLEL DO PRIVATE (k, jcam, jdyn, i)
    do k=1, plev
    do jdyn = jbeg0,jend0
      jcam = plat + 1 - jdyn
      do i=beglonex,endlonex
        fu_fv_tmp(i,jdyn,k) = fv(i,k,jcam)
      enddo
    enddo
    enddo
    call gamil_arrays_comm(COMM_TO_TOP,1,fu_fv_tmp(:,beglatexdyn,1))

!$OMP PARALLEL DO PRIVATE (k, jcam, jdyn, i)
     do k=1,plev
       do jdyn=jbeg0,jend0
         if (jdyn .eq. plat) then
           do i=beglonex,endlonex
             sv(i,jdyn,k) = 0.0
           enddo
         else
           jcam = plat + 1 - jdyn
           do i=beglonex,iend2
             sv(i,jdyn,k)= sv(i,jdyn,k)-0.5*( fu_fv_tmp(i,jdyn,k)+fu_fv_tmp(i,jdyn+1,k) )
           enddo
         endif
       enddo
     enddo

   call gamil_arrays_comm(COMM_TO_RIGHT,1,su(:,beglatexdyn,1)) 
   call gamil_arrays_comm(COMM_ROTATE_LEFT,2,su(:,beglatexdyn,1),st(:,beglatexdyn,1), &
                           q(:,beglatexdyn,1),sv(:,beglatexdyn,1))

   call remove_comm_array(fu_fv_tmp(:,beglatexdyn,1))


     return
end subroutine a_c_switching
