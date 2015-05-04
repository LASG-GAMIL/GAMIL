#include <misc.h>
#include <params.h>
subroutine trunc
!-----------------------------------------------------------------------
!
! Purpose:
! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: trunc.F90,v 1.3.6.1 2002/06/15 13:48:34 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe

  implicit none

!---------------------------Local variables-----------------------------
!
  integer n              ! Loop index over diagonals
  integer ik2            ! K+2
  integer m              ! loop index
!
!-----------------------------------------------------------------------
!
! trunc first evaluates truncation parameters for a general pentagonal 
! truncation for which the following parameter relationships are true
!
! 0 .le. |m| .le. ptrm
!
! |m| .le. n .le. |m|+ptrn for |m| .le. ptrk-ptrn
!
! |m| .le. n .le. ptrk     for (ptrk-ptrn) .le. |m| .le. ptrm
!
! Most commonly utilized truncations include:
!  1: triangular  truncation for which ptrk=ptrm=ptrn
!  2: rhomboidal  truncation for which ptrk=ptrm+ptrn
!  3: trapezoidal truncation for which ptrn=ptrk .gt. ptrm
!
! Simple sanity check
! It is necessary that ptrm .ge. ptrk-ptrn .ge. 0
!
  if (ptrm.lt.(ptrk-ptrn)) then
     write(6,*)'TRUNC: Error in truncation parameters'
     write(6,*)'       ntrm.lt.(ptrk-ptrn)'
     call endrun
  end if
  if (ptrk.lt.ptrn) then
     write(6,*)'TRUNC: Error in truncation parameters'
     write(6,*)'       ptrk.lt.ptrn'
     call endrun
  end if
!
! Evaluate pointers and displacement info based on truncation params
!
! The following ifdef logic seems to have something do with SPMD 
! implementation, although it's not clear how this info is used
! Dave, can you check this with JR?
!
  ncoefi(1) = 1
  ik2 = ptrk + 2
  do n=1,pmax
     ncoefi(n+1) = ncoefi(n) + min0(pmmax,ik2-n)
     nalp(n) = ncoefi(n) - 1
     nco2(n) = ncoefi(n)*2
     nm(n) = ncoefi(n+1) - ncoefi(n)
  end do
  nstart(1) = 0
  nlen(1) = ptrn + 1
  do m=2,pmmax
     nstart(m) = nstart(m-1) + nlen(m-1)
     nlen(m) = min0(ptrn+1,ptrk+2-m)
  end do
!      write(6,*)'Starting index  length'
!      do m=1,ptrm+1
!         write(6,'(1x,i14,i8)')nstart(m),nlen(m)
!      end do
!
! Define break-even point for vector lengths in GRCALC.  Don't implement
! for non-PVM machines
!
  ncutoff = pmax
#if ( defined PVP )
  if (2*nm(1).lt.plev) ncutoff = 0
  do n=2,pmax,2
     if (2*nm(n).lt.plev) then
        ncutoff = n
        goto 100
     end if
  end do
100 continue
  write(6,*)'TRUNC: n cutoff for GRCALC vectorization = ',ncutoff
#endif
!
  return
end subroutine trunc

