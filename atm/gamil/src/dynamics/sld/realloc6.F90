#include <misc.h>
#include <params.h>

subroutine realloc6

#ifdef SPMD
!-----------------------------------------------------------------------
!
! Purpose:
! Reallocation routine for spectral prognostics.
!
! Author:  J. Rosinski
!
!----------------------------------------------------------------------------
!
! $Id: realloc6.F90,v 1.3.22.1 2002/06/15 13:48:28 erik Exp $
! $Author: erik $
!
!----------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use spmd_dyn
  use mpishorthand

  implicit none

!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtype=1001 ! message passing id
!
!---------------------------Local workspace-----------------------------
!
  integer m,k
  integer mstrt,length
  integer i,mask,procid
  integer length_p,mstrt_p
  integer bpos
  integer mb,me
!
!----------------------------------------------------------------------------
!
! Spectral processors ship their "m" values
!
  if(iam.le.npessp-1) then
     mb = begm(iam)
     me = endm(iam)
     mstrt = 2*nstart(mb)+1
     length = 2*nstart(me)-2*nstart(mb)+2*nlen(me)
  else
     mstrt = 2*psp
     length = 0
  endif
  mask = 1

  do while (mask.lt.ceil2(npes))
     procid = pair(npes,iam,mask)
     if (procid.ge.0) then
        bpos = 0
        call mpipack (length,1,mpiint,buf1,bsiz,bpos,mpicom)
        call mpipack (mstrt,1,mpiint,buf1,bsiz,bpos,mpicom)
        if (length.gt.0) then
           call mpipack (alps(mstrt),length,mpir8,buf1,bsiz,bpos,mpicom)
           do k=1,plev
              call mpipack (t (mstrt,k),length,mpir8,buf1,bsiz,bpos,mpicom)
              call mpipack (q (mstrt,k),length,mpir8,buf1,bsiz,bpos,mpicom)
              call mpipack (d (mstrt,k),length,mpir8,buf1,bsiz,bpos,mpicom)
              call mpipack (vz(mstrt,k),length,mpir8,buf1,bsiz,bpos,mpicom)
           enddo
        endif

        call mpisendrecv(buf1,bpos,mpipk,procid,msgtype, &
                         buf2,bsiz,mpipk,procid,msgtype, &
                         mpicom)

        bpos = 0
        call mpiunpack (buf2,bsiz,bpos,length_p,1,mpiint,mpicom)
        call mpiunpack (buf2,bsiz,bpos,mstrt_p,1,mpiint,mpicom)
        if (length_p.gt.0) then
           call mpiunpack (buf2,bsiz,bpos,alps(mstrt_p),length_p,mpir8,mpicom)
           do k=1,plev
              call mpiunpack (buf2,bsiz,bpos,t (mstrt_p,k),length_p,mpir8,mpicom)
              call mpiunpack (buf2,bsiz,bpos,q (mstrt_p,k),length_p,mpir8,mpicom)
              call mpiunpack (buf2,bsiz,bpos,d (mstrt_p,k),length_p,mpir8,mpicom)
              call mpiunpack (buf2,bsiz,bpos,vz(mstrt_p,k),length_p,mpir8,mpicom)
           enddo
        endif
     endif

     if (npes.eq.ceil2(npes)) then
        mstrt = min(mstrt,mstrt_p)
        length = length+length_p
        mask = mask*2
     else
        mask = mask+1
     endif
!JR         call barrier(mpicom)
  end do
#endif
  return
end subroutine realloc6
