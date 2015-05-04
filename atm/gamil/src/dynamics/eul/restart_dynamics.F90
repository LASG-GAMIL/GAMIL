#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use prognostics
   use ppgrid, only: pcols, pver
!! use comslt
   use binary_io
   use comfm1

   implicit none

#include <pdyn.h> 
#include <comfm2.h> 

CONTAINS

   subroutine write_restart_dynamics (nrg)

!!#include <comqfl.h>

      integer, intent(in) :: nrg     ! unit number

      integer begj    ! starting latitude
      integer ioerr   ! error status

      begj = beglatex+numbnd
!
! prognostics of cam2
!
      call wrtout_r8(nrg, phis, plond)
      call wrtout_r8(nrg, omga, plndlv)

      call wrtout_r8(nrg, u3(1,1,begj,n3m2), plndlv)
      call wrtout_r8(nrg, v3(1,1,begj,n3m2), plndlv)

      call wrtout_r8(nrg, t3(1,1,begj,n3m2), plndlv)
      call wrtout_r8(nrg, q3(1,1,1,begj,n3m2), plndlv*(pcnst+pnats))
      call wrtout_r8(nrg, t31(1,1,begj), plndlv)
      call wrtout_r8(nrg, q31(1,1,begj), plndlv)

      call wrtout_r8(nrg, ps(1,beglat,n3m2), plond)
!
! 'prognostics' of fm2003 : u,v,t,q,wpa,pes,ghs (from module comfm1)
!
!!(wh 2003.10.18)

      call wrtout_r8_dyn(nrg, u, plond, numlatsex, plev)
      call wrtout_r8_dyn(nrg, v, plond, numlatsex, plev)
      call wrtout_r8_dyn(nrg, t, plond, numlatsex, plev)
      call wrtout_r8_dyn(nrg, q, plond, numlatsex, plev)

      call wrtout_r8_dyn (nrg,ws  ,plond,numlatsex, plevp)
      call wrtout_r8_dyn (nrg,wpa ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,ghi ,plond,numlatsex, plevp)
      call wrtout_r8_dyn (nrg,pes ,plond,numlatsex, 1)
      call wrtout_r8_dyn (nrg,ghs ,plond,numlatsex, 1)
!
!  tendencies of fm2003 : su,sv,st  (sq will always be 0.0)
!
      call wrtout_r8_dyn (nrg,su  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,sv  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,st  ,plond,numlatsex, plev)
!
!  variables used in 'dynamics' (from module comfm2)
!
      call wrtout_r8_dyn (nrg,du  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dv  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dtt ,plond,numlatsex, plev)

      call wrtout_r8_dyn (nrg,du0 ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dv0 ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dtt0,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dps0,plond,numlatsex, 1)

      call wrtout_r8_dyn (nrg,du1 ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dv1 ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dtt1,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dps1,plond,numlatsex, 1)

      call wrtout_r8_dyn (nrg,uu  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,vv  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,tt  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,p   ,plond,numlatsex, 1)
      call wrtout_r8_dyn (nrg,ply ,plond,numlatsex, plev)
       
      call wrtout_r8_dyn (nrg,up  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,vp  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,ttp ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,pps ,plond,numlatsex, 1)
      call wrtout_r8_dyn (nrg,dps ,plond,numlatsex, 1)
       
      call wrtout_r8_dyn (nrg,uk  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,vk  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,ttk ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,psk ,plond,numlatsex, 1)

      call wrtout_r8_dyn (nrg,tb  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,cb  ,plond,numlatsex, plev)
      call wrtout_r8_dyn (nrg,dcb ,plond,numlatsex, plev)

      call wrtout_r8_dyn (nrg,hps ,plond,numlatsex, 1)
      call wrtout_r8_dyn (nrg,c0  ,plond,numlatsex, 1)
 
      call wrtout_r8_dyn (nrg,cb0 ,plond,numlatsex, plev)
      call wrtout_int_dyn(nrg,nigw,    1,numlatsex, 1)
!       
!  some other parameters
!
      if (masterproc) then
         write(nrg, iostat=ioerr) itime ,dlt1 ,dlt2

         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if

      return
   end subroutine write_restart_dynamics

!#######################################################################

   subroutine read_restart_dynamics (nrg)

#if ( defined SPMD )
      use mpishorthand
#endif

!!#include <comqfl.h>
!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: begj    ! starting latitude
      integer :: ioerr   ! error status

      call initialize_prognostics
      call initialize_comfm1           !!(wh 2003.10.23)

      begj = beglatex + numbnd

!
! prognostics of cam2
!
      call readin_r8 (nrg,phis  ,plond )
      call readin_r8 (nrg,omga  ,plndlv)

      call readin_r8 (nrg,u3(1,1,begj,n3m2)  ,plndlv)
      call readin_r8 (nrg,v3(1,1,begj,n3m2)  ,plndlv)

      call readin_r8(nrg, t3(1,1,begj,n3m2), plndlv)
      call readin_r8(nrg, q3(1,1,1,begj,n3m2), plndlv*(pcnst+pnats))
      call readin_r8(nrg, t31(1,1,begj), plndlv)
      call readin_r8(nrg, q31(1,1,begj), plndlv)

      call readin_r8(nrg, ps(1,beglat,n3m2), plond)
!
! 'prognostics' of fm2003 : u,v,t,q,wpa,pes,ghs  ( in module comfm1 )
!
      call readin_r8_dyn (nrg,u   ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,v   ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,t   ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,q   ,plond,numlatsex, plev)

      call readin_r8_dyn (nrg,ws  ,plond,numlatsex, plevp)
      call readin_r8_dyn (nrg,wpa ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,ghi ,plond,numlatsex, plevp)
      call readin_r8_dyn (nrg,pes ,plond,numlatsex, 1)
      call readin_r8_dyn (nrg,ghs ,plond,numlatsex, 1)
!
!  tendencies of fm2003 : su,sv,st  (sq will always be 0.0)
!
      call readin_r8_dyn (nrg,su  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,sv  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,st  ,plond,numlatsex, plev)
!
!  variables used in 'dynamics' (from module comfm2)
!
      call readin_r8_dyn (nrg,du  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dv  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dtt ,plond,numlatsex, plev)

      call readin_r8_dyn (nrg,du0 ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dv0 ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dtt0,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dps0,plond,numlatsex, 1)

      call readin_r8_dyn (nrg,du1 ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dv1 ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dtt1,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dps1,plond,numlatsex, 1)

      call readin_r8_dyn (nrg,uu  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,vv  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,tt  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,p   ,plond,numlatsex, 1)
      call readin_r8_dyn (nrg,ply ,plond,numlatsex, plev)
       
      call readin_r8_dyn (nrg,up  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,vp  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,ttp ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,pps ,plond,numlatsex, 1)
      call readin_r8_dyn (nrg,dps ,plond,numlatsex, 1)
       
      call readin_r8_dyn (nrg,uk  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,vk  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,ttk ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,psk ,plond,numlatsex, 1)

      call readin_r8_dyn (nrg,tb  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,cb  ,plond,numlatsex, plev)
      call readin_r8_dyn (nrg,dcb ,plond,numlatsex, plev)

      call readin_r8_dyn (nrg,hps ,plond,numlatsex, 1)
      call readin_r8_dyn (nrg,c0  ,plond,numlatsex, 1)

      call readin_r8_dyn (nrg,cb0 ,plond,numlatsex, plev)
      call readin_int_dyn (nrg,nigw,   1,numlatsex, 1)
!       
!  some other parameters
!
      if (masterproc) then
         read(nrg, iostat=ioerr) itime ,dlt1 ,dlt2

         if (ioerr /= 0 ) then
            write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if


#if ( defined SPMD )
!!   call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
!!   call mpibcast (fixmas,1         ,mpir8  ,0,mpicom)
!!   call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
!!   call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
!!   call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
!!   call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)

     call mpibcast (dlt1  ,1         ,mpir8  ,0,mpicom)      
     call mpibcast (dlt2  ,1         ,mpir8  ,0,mpicom)
     call mpibcast (itime ,1         ,mpiint ,0,mpicom)
#endif

      return

   end subroutine read_restart_dynamics

end module restart_dynamics
