#include <misc.h>
#include <params.h>
subroutine virtem(ncol    ,ncold   ,nver    ,t       ,q       ,zvir    , &
                  tv      )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the virtual temperature.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol  ! number of atmospheric columns
   integer, intent(in) :: ncold ! atmospheric column index dimension
   integer, intent(in) :: nver  ! number of vertical levels in a column

   real(r8) t(ncold,nver)       ! temperature
   real(r8) q(ncold,nver)       ! specific humidity
   real(r8) zvir                ! virtual temperature constant
!
! Output arguments
!
   real(r8), intent(out) :: tv(ncold,nver)      ! virtual temperature
!
!---------------------------Local storage-------------------------------
!
   integer i,k              ! column and level indexes
!
   do k=1,nver
      do i=1,ncol
         tv(i,k) = t(i,k)*(1.0 + zvir*q(i,k))
      end do
   end do
!
   return
end subroutine virtem

