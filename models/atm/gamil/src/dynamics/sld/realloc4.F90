#include <misc.h> 
#include <params.h> 

subroutine realloc4(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
                    grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
                    grz2    ,grd2    ,grfu2   ,grfv2   )
#ifdef SPMD
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for the Fourier coefficients
!
! Author:  J. Rosinski
!
!----------------------------------------------------------------------------
!
! $Id: realloc4.F90,v 1.4.22.1 2002/06/15 13:48:28 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use spmd_dyn
  use mpishorthand

  implicit none

#include <comsta.h>

!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype  = 1000
!
!----------------------------- Arguments ----------------------------
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
  integer procid,length,j,mstrt,j_p
  integer length_p,mstrt_p
  integer bpos
  integer procj
  integer begirow_p,endirow_p
  integer num
!
!-----------------------------------------------------------------------
!     
!     Send gr..1 "m" values to processor(s) owning that wavenumber
!     
  length_p = 2*numm(iam)
  mstrt_p = 2*begm(iam)-1

  do procj=1,ceil2(npes)-1
     procid = pair(npes,procj,iam)
     length = 2*numm(procid)
     mstrt = 2*begm(procid)-1

     if (length > 0 .and. length_p > 0 .and. procid >= 0) then
        bpos = 0
        call mpipack (begirow,1,mpiint,buf1,bsiz,bpos,mpicom)
        call mpipack (endirow,1,mpiint,buf1,bsiz,bpos,mpicom)

        do j=begirow,endirow
           call mpipack (grlps1(mstrt,j),length,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grlps2(mstrt,j),length,mpir8,buf1,bsiz,bpos,mpicom)

           num = length*plev
           call mpipack (grt1 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grq1 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grz1 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grd1 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grfu1(1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grfv1(1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grt2 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grq2 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grz2 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grd2 (1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grfu2(1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
           call mpipack (grfv2(1,mstrt,j),num,mpir8,buf1,bsiz,bpos,mpicom)
        end do

        call mpisendrecv(buf1,bpos,mpi_packed,procid,msgtype, &
                         buf2,bsiz,mpi_packed,procid,msgtype,mpicom)

        bpos = 0
        call mpiunpack (buf2,bsiz,bpos,begirow_p,1,mpiint,mpicom)
        call mpiunpack (buf2,bsiz,bpos,endirow_p,1,mpiint,mpicom)
        do j_p=begirow_p,endirow_p
           call mpiunpack (buf2,bsiz,bpos,grlps1(mstrt_p,j_p),length_p,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grlps2(mstrt_p,j_p),length_p,mpir8,mpicom)

           num = length_p*plev
           call mpiunpack (buf2,bsiz,bpos,grt1 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grq1 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grz1 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grd1 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grfu1(1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grfv1(1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grt2 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grq2 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grz2 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grd2 (1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grfu2(1,mstrt_p,j_p),num,mpir8,mpicom)
           call mpiunpack (buf2,bsiz,bpos,grfv2(1,mstrt_p,j_p),num,mpir8,mpicom)
        end do
     end if
!JR        call barrier (mpicom)
  end do

#endif
  return
end subroutine realloc4
