#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use prognostics
   use ppgrid, only: pcols, pver
   use comslt
   use binary_io
#ifdef SPMD
   use pspect, only: psp
#endif
   use comspe, only: lnpstar

   implicit none

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: begj    ! starting latitude
      integer :: ioerr   ! error status
      integer :: num     ! number of values
!
      call wrtout_r8 (nrg,div(1,1,beglat,n3  ) , plndlv)
      call wrtout_r8 (nrg,div(1,1,beglat,n3m1) , plndlv)

      call wrtout_r8 (nrg,urhs,plndlv)
      call wrtout_r8 (nrg,vrhs,plndlv)
      call wrtout_r8 (nrg,trhs,plndlv)
      call wrtout_r8 (nrg,prhs,plndlv)

      call wrtout_r8 (nrg,dpsl  ,plond )
      call wrtout_r8 (nrg,dpslm1,plond )
      call wrtout_r8 (nrg,dpsm  ,plond )
      call wrtout_r8 (nrg,dpsmm1,plond )
      call wrtout_r8 (nrg,ql    ,plndlv)
      call wrtout_r8 (nrg,qm    ,plndlv)
      call wrtout_r8 (nrg,dps   ,plond )
      call wrtout_r8 (nrg,phis  ,plond )
      call wrtout_r8 (nrg,phisl ,plond )
      call wrtout_r8 (nrg,phism ,plond )
      call wrtout_r8 (nrg,tlm1  ,plndlv)
      call wrtout_r8 (nrg,tl    ,plndlv)
      call wrtout_r8 (nrg,tmm1  ,plndlv)
      call wrtout_r8 (nrg,tm    ,plndlv)
      call wrtout_r8 (nrg,omga  ,plndlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      begj = beglatex + numbnd

      call wrtout_r8 (nrg,u3(1,1,begj,n3  )  ,plndlv)
      call wrtout_r8 (nrg,v3(1,1,begj,n3  )  ,plndlv)
      call wrtout_r8 (nrg,t3(1,1,begj,n3  )  ,plndlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3  )  ,plond)

      call wrtout_r8 (nrg,u3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,v3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,t3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m1)  ,plond)
      
      call wrtout_r8 (nrg,q3(1,1,1,begj,n3  ),plndlv*(pcnst+pnats))
      call wrtout_r8 (nrg,q3(1,1,1,begj,n3m1),plndlv*(pcnst+pnats))

      num = plndlv
      call wrtout_r8 (nrg,parrsld(1,1,begj),num)

      num = plond*plevp
      call wrtout_r8 (nrg,etadot(1,1,begj,n3m1),num)

      num = plndlv*pcnst
      call wrtout_r8  (nrg,qfcst(1,1,1,beglat),num)
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0, hw1, hw2, hw3, alpha, lnpstar
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if

      return
   end subroutine write_restart_dynamics

   subroutine read_restart_dynamics (nrg)

#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: begj    ! starting latitude
      integer :: ioerr   ! error status
      integer :: num     ! number of values
!
      call initialize_prognostics
      call readin_r8 (nrg,div(1,1,beglat,n3  ) , plndlv)
      call readin_r8 (nrg,div(1,1,beglat,n3m1) , plndlv)

      call readin_r8 (nrg,urhs,plndlv)
      call readin_r8 (nrg,vrhs,plndlv)
      call readin_r8 (nrg,trhs,plndlv)
      call readin_r8 (nrg,prhs,plndlv)

      call readin_r8 (nrg,dpsl  ,plond )
      call readin_r8 (nrg,dpslm1,plond )
      call readin_r8 (nrg,dpsm  ,plond )
      call readin_r8 (nrg,dpsmm1,plond )
      call readin_r8 (nrg,ql    ,plndlv)
      call readin_r8 (nrg,qm    ,plndlv)
      call readin_r8 (nrg,dps   ,plond )
      call readin_r8 (nrg,phis  ,plond )
      call readin_r8 (nrg,phisl ,plond )
      call readin_r8 (nrg,phism ,plond )
      call readin_r8 (nrg,tlm1  ,plndlv)
      call readin_r8 (nrg,tl    ,plndlv)
      call readin_r8 (nrg,tmm1  ,plndlv)
      call readin_r8 (nrg,tm    ,plndlv)
      call readin_r8 (nrg,omga  ,plndlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      begj = beglatex + numbnd

      call readin_r8 (nrg,u3(1,1,begj,n3  )  ,plndlv)
      call readin_r8 (nrg,v3(1,1,begj,n3  )  ,plndlv)
      call readin_r8 (nrg,t3(1,1,begj,n3  )  ,plndlv)
      call readin_r8 (nrg,ps(1,beglat,n3  )  ,plond)

      call readin_r8 (nrg,u3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,v3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,t3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,ps(1,beglat,n3m1)  ,plond)
      
      call readin_r8 (nrg,q3(1,1,1,begj,n3  ),plndlv*(pcnst+pnats))
      call readin_r8 (nrg,q3(1,1,1,begj,n3m1),plndlv*(pcnst+pnats))

      num = plndlv
      call readin_r8 (nrg,parrsld(1,1,begj),num)

      num = plond*plevp
      call readin_r8 (nrg,etadot(1,1,begj,n3m1),num)

!
! Write slt forcast of moisture and constituents
!
      call initialize_comslt
      num = plndlv*pcnst
      call readin_r8 (nrg,qfcst(1,1,1,beglat),num)
!
! Read global integrals
!
      if (masterproc) then
         read(nrg, iostat=ioerr) tmass0, hw1, hw2, hw3, alpha, lnpstar
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if

#if ( defined SPMD )
      call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
      call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
      call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)
      call mpibcast (lnpstar, psp     ,mpir8  ,0,mpicom)
#endif

      return
   end subroutine read_restart_dynamics

end module restart_dynamics
