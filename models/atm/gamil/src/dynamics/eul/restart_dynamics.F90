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
      use mpi_gamil

!!#include <comqfl.h>

      integer, intent(in) :: nrg     ! unit number

      integer begj    ! starting latitude
      integer ioerr   ! error status
      real(r8) :: local_3d_array(beglonex:endlonex,beglatexdyn:endlatexdyn,plev+2)
      integer :: i,jdyn,jcam,k,m

      begj = beglatex+numbnd
!
! prognostics of cam2
!

      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,1) = phis(i,jcam)
        enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array)
 
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = omga(i,k,jcam)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)
 
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = u3(i,k,jcam,n3m2)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)

      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = v3(i,k,jcam,n3m2)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)

      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = t3(i,k,jcam,n3m2)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)

      do m=1,pcnst+pnats
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = q3(i,k,m,jcam,n3m2)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)
      enddo

      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = t31(i,k,jcam)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)

      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,k) = q31(i,k,jcam)
        enddo
      enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array)

      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          local_3d_array(i,jdyn,1) = ps(i,jcam,n3m2)
        enddo
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array)

! 'prognostics' of fm2003 : u,v,t,q,wpa,pes,ghs (from module comfm1)
!
!!(wh 2003.10.18)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,u)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,v)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,t)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,q)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,ws)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,wpa)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,ghi)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,pes)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,ghs)

!
!  tendencies of fm2003 : su,sv,st  (sq will always be 0.0)
!
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,su)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,sv)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,st)
!
!  variables used in 'dynamics' (from module comfm2)
!
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du0)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv0)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt0)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps0)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du1)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv1)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt1)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps1)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,uu)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,vv)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,tt)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,p)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,ply)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,uk)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,vk)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,ttk)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,psk)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,tb)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,cb)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dcb)

      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,hps)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,c0)
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,cb0)
      do jdyn=beglatexdyn, endlatexdyn
         local_3d_array(:,jdyn,1) = nigw(jdyn)
      enddo
      call wrtout_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array)
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
      use mpi_gamil

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
      real(r8) :: local_3d_array(beglonex:endlonex,beglatexdyn:endlatexdyn,plev+2)
      integer :: i,jdyn,jcam,k,m

      call initialize_prognostics
      call initialize_comfm1           !!(wh 2003.10.23)

      begj = beglatex + numbnd

!
! prognostics of cam2
!
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array,.false.)
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          phis(i,jcam) = local_3d_array(i,jdyn,1)
        enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          omga(i,k,jcam) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          u3(i,k,jcam,n3m2) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          v3(i,k,jcam,n3m2) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          t3(i,k,jcam,n3m2) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      do m=1,pcnst+pnats
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          q3(i,k,m,jcam,n3m2) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          t31(i,k,jcam) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,local_3d_array,.false.)
      do k=1, plev
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          q31(i,k,jcam) = local_3d_array(i,jdyn,k) 
        enddo
      enddo
      enddo

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array,.false.)
      do jdyn = jbeg0,jend0
        jcam = plat + 1 - jdyn
        do i=beglonex,endlonex
          ps(i,jcam,n3m2) = local_3d_array(i,jdyn,1)
        enddo
      enddo
!
! 'prognostics' of fm2003 : u,v,t,q,wpa,pes,ghs  ( in module comfm1 )
!
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,u,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,v,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,t,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,q,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,ws,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,wpa,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,ghi,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,pes,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,ghs,.true.)
!
!  tendencies of fm2003 : su,sv,st  (sq will always be 0.0)
!
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,su,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,sv,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,st,.true.)
!
!  variables used in 'dynamics' (from module comfm2)
!
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du0,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv0,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt0,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps0,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,du1,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dv1,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dtt1,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,dps1,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,uu,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,vv,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,tt,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,p,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,ply,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,uk,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,vk,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,ttk,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,psk,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,tb,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,cb,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,dcb,.true.)

      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,hps,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,c0,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,cb0,.true.)
      call readin_3D_array_dyn(nrg,beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,local_3d_array,.true.)
      do jdyn=jbeg0, jend0
         nigw(jdyn) = local_3d_array(beglonex+1,jdyn,1)
      enddo
!       
!  some other parameters
!
      itime = 0
      dlt1 = 0.0
      dlt2 = 0.0
      if (masterproc) then
         read(nrg, iostat=ioerr) itime ,dlt1 ,dlt2

         if (ioerr /= 0 ) then
            write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if


#if ( defined SPMD )
     call mpibcast (dlt1  ,1         ,mpir8  ,0,mpicom)      
     call mpibcast (dlt2  ,1         ,mpir8  ,0,mpicom)
     call mpibcast (itime ,1         ,mpiint ,0,mpicom)
#endif

      return

   end subroutine read_restart_dynamics

end module restart_dynamics
