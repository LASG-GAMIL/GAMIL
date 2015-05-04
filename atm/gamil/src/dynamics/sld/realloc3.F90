#include <misc.h>
#include <params.h>

subroutine realloc3(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
                    grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
                    grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for Gaussian quadrature.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: realloc3.F90,v 1.4.22.1 2002/06/15 13:48:28 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

#ifdef SPMD

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use mpishorthand
  use spmd_dyn

  implicit none

!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgid = 3000
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: grlps1(     2*pmmax,plat/2) ! ----------------------------
  real(r8), intent(in)   :: grt1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grq1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grz1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grd1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grfu1 (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grfv1 (plev,2*pmmax,plat/2) ! |  definitions: these variables
  real(r8), intent(in)   :: grlps2(     2*pmmax,plat/2) ! |  are declared here for data
  real(r8), intent(in)   :: grt2  (plev,2*pmmax,plat/2) ! |  scoping
  real(r8), intent(in)   :: grq2  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grz2  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grd2  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grfu2 (plev,2*pmmax,plat/2) ! |
  real(r8), intent(in)   :: grfv2 (plev,2*pmmax,plat/2) ! ----------------------------
!
!---------------------------Local workspace-----------------------------
!
  integer procid
  integer num
!
!-----------------------------------------------------------------------
!
  procid = npes - iam - 1

  if (cut(2,iam).le.plat/2) then
!
! Send gr..1 to, receive gr..2 from, my neighbor in the northern hemisphere
!
     num = 2*pmmax*numlats
     call mpisendrecv (grlps1(1,begirow),num,mpir8,procid,msgid, &
                       grlps2(1,begirow),num,mpir8,procid,msgid,mpicom)

     num = plev*2*pmmax*numlats
     call mpisendrecv (grt1 (1,1,begirow),num,mpir8,procid,msgid, &
                       grt2 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grq1 (1,1,begirow),num,mpir8,procid,msgid, &
                       grq2 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grz1 (1,1,begirow),num,mpir8,procid,msgid, &
                       grz2 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grd1 (1,1,begirow),num,mpir8,procid,msgid, &
                       grd2 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grfu1(1,1,begirow),num,mpir8,procid,msgid, &
                       grfu2(1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grfv1(1,1,begirow),num,mpir8,procid,msgid, &
                       grfv2(1,1,begirow),num,mpir8,procid,msgid,mpicom)

  else

!
! Send gr..2 to, receive gr..1 from, my neighbor in the southern hemisphere
!
     num = 2*pmmax*numlats
     call mpisendrecv (grlps2(1,begirow),num,mpir8,procid,msgid, &
                       grlps1(1,begirow),num,mpir8,procid,msgid,mpicom)

     num = plev*2*pmmax*numlats
     call mpisendrecv (grt2 (1,1,begirow),num,mpir8,procid,msgid, &
                       grt1 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grq2 (1,1,begirow),num,mpir8,procid,msgid, &
                       grq1 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grz2 (1,1,begirow),num,mpir8,procid,msgid, &
                       grz1 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grd2 (1,1,begirow),num,mpir8,procid,msgid, &
                       grd1 (1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grfu2(1,1,begirow),num,mpir8,procid,msgid, &
                       grfu1(1,1,begirow),num,mpir8,procid,msgid,mpicom)
     call mpisendrecv (grfv2(1,1,begirow),num,mpir8,procid,msgid, &
                       grfv1(1,1,begirow),num,mpir8,procid,msgid,mpicom)
  end if
#endif
  return
end subroutine realloc3


