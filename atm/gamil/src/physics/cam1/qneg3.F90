#include <misc.h>
#include <params.h>
subroutine qneg3 (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst  , &
                  qmin    ,q       )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst       ! number of constituents

   real(r8), intent(in) :: qmin(lconst)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(ncold,lver,lconst) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer indx(ncol )      ! array of indices of points < qmin
   integer nval             ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst           ! biggest violator
!
!-----------------------------------------------------------------------
!
   do m=1,lconst
#ifdef HADVTEST
!jr Disable this routine for purposes of advection test
      return
#endif
      nvals = 0
      found = .false.
      worst = 1.e35
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
      do k=1,lver
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
               q(i,k,m) = qmin(m)
            end do
         end if
      end do
      if (found) then
         !write(6,9000)subnam,m,idx,nvals,qmin(m),worst,iw,kw !--sxj--
      end if
   end do
!
   return
9000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i3, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
end subroutine qneg3

