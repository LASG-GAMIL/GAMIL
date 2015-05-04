#include <misc.h>
#include <params.h>

subroutine qneg3a (subnam  ,lchnk     ,q       ,mfirst  ,mlast   , &
                   ncol    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, and return information
! to allow warning message to be printed.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: E. Kluzek
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use constituents, only: ppcnst, qmin

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam     ! name of calling routine
   integer, intent(in) :: lchnk            ! chunk index
   integer, intent(in) :: mfirst           ! Starting constituent index
   integer, intent(in) :: mlast            ! Last constituent index
   integer, intent(in) :: ncol
!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(pcols,pver,ppcnst) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer indx(pcols)      ! array of indices of points < qmin
   integer nval             ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst               ! biggest violator
!
!-----------------------------------------------------------------------
!
! Loop over constituents
!
   do m=mfirst, mlast
      nvals = 0
      found = .false.
      worst = 1.e35
!
! Test all field values for being less than minimum value.
! Trace offenders and identify worst one.
!
      do k=1,pver
         nval = 0
         do i=1,ncol
            if (q(i,k,m) < qmin(m)) then
               nval = nval + 1
               indx(nval) = i
            end if
         end do
         if (nval > 0) then
            found = .true.
            nvals = nvals + nval
            do ii=1,nval
               i = indx(ii)
               if (q(i,k,m) < worst) then
                  worst = q(i,k,m)
                  kw = k
                  iw = i
               end if
            end do
         end if
      end do
      if (found) then
         !write(6,9000)subnam,m,get_lat_p(lchnk,iw),nvals,worst,get_lon_p(lchnk,iw),kw  !--sxj--
      end if
   end do
!
   return
9000 format(' QNEG3A from ',a,':m=',i3,' lat=',i3, &
            ' Min. mixing ratio violated at ',i4,' points.  ', &
            ' Worst =',e8.1,' at i,k=',i4,i3)
end subroutine qneg3a

