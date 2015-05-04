#include <misc.h>
#include <params.h>
subroutine realloc7 (vmax2d  ,vmax2dt ,vcour   )
#ifdef SPMD
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for energy and log stats
!
! Author:  J. Rosinski
!
!----------------------------------------------------------------------------
!
! $Id: realloc7.F90,v 1.5.8.1 2002/06/15 13:48:28 erik Exp $
! $Author: erik $
!
!----------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plat, plev, iam, numlats, beglat
  use mpishorthand
  use spmd_dyn

  implicit none

#include <comsta.h>

!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype=1000 ! message passing id
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: vmax2d (plev,plat) ! Max. wind at each level, latitude
  real(r8), intent(in)   :: vmax2dt(plev,plat) ! Max. truncated wind at each lvl,lat
  real(r8), intent(in)   :: vcour  (plev,plat) ! Maximum Courant number in slice
!
!---------------------------Local workspace-----------------------------
!
  integer procid
  integer bufpos
  integer procj
  integer beglat_p,numlats_p
  integer num
  integer i
!
!-----------------------------------------------------------------------
!
  do procj=1,ceil2(npes)-1
     procid = pair(npes,procj,iam)
     if (procid.ge.0) then
        bufpos = 0

        call mpipack (beglat ,1,mpiint,buf1,bsiz,bufpos,mpicom)
        call mpipack (numlats,1,mpiint,buf1,bsiz,bufpos,mpicom)

        call mpipack(psurf(beglat),numlats,mpir8,buf1,bsiz,bufpos,mpicom)
        call mpipack(stq(beglat)  ,numlats,mpir8,buf1,bsiz,bufpos,mpicom)
        call mpipack(rmst(beglat) ,numlats,mpir8,buf1,bsiz,bufpos,mpicom)
        call mpipack(rmsd(beglat) ,numlats,mpir8,buf1,bsiz,bufpos,mpicom)

        num = numlats*plev
        call mpipack(vmax2d(1,beglat) ,num,mpir8,buf1,bsiz,bufpos,mpicom)
        call mpipack(vmax2dt(1,beglat),num,mpir8,buf1,bsiz,bufpos,mpicom)
        call mpipack(vcour(1,beglat)  ,num,mpir8,buf1,bsiz,bufpos,mpicom)

        call mpisendrecv(buf1,bufpos,mpipk,procid,msgtype, &
                         buf2,bsiz  ,mpipk,procid,msgtype,mpicom)

        bufpos = 0

        call mpiunpack (buf2,bsiz,bufpos,beglat_p ,1,mpiint,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,numlats_p,1,mpiint,mpicom)

        call mpiunpack (buf2,bsiz,bufpos,psurf(beglat_p),numlats_p,mpir8,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,stq(beglat_p)  ,numlats_p,mpir8,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,rmst(beglat_p) ,numlats_p,mpir8,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,rmsd(beglat_p) ,numlats_p,mpir8,mpicom)

        num = numlats_p*plev
        call mpiunpack (buf2,bsiz,bufpos,vmax2d(1,beglat_p) ,num,mpir8,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,vmax2dt(1,beglat_p),num,mpir8,mpicom)
        call mpiunpack (buf2,bsiz,bufpos,vcour(1,beglat_p)  ,num,mpir8,mpicom)
     endif
  enddo
#endif
  return
end subroutine realloc7
