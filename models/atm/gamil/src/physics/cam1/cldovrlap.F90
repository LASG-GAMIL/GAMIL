#include <misc.h>
#include <params.h>
subroutine cldovrlap(lchnk   ,ncol    ,pint    ,cld     ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Partitions each column into regions with clouds in neighboring layers.
! This information is used to implement maximum overlap in these regions
! with random overlap between them.
! On output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! Method: 

! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pint(pcols,pverp)   ! Interface pressure
   real(r8), intent(in) :: cld(pcols,pver)     ! Fractional cloud cover
!
! Output arguments
!
   real(r8), intent(out) :: pmxrgn(pcols,pverp)! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer nmxrgn(pcols)                    ! Number of maximally overlapped regions
!
!---------------------------Local variables-----------------------------
!
   integer i                    ! Longitude index
   integer k                    ! Level index
   integer n                    ! Max-overlap region counter

   real(r8) pnm(pcols,pverp)    ! Interface pressure

   logical cld_found            ! Flag for detection of cloud
   logical cld_layer(pver)      ! Flag for cloud in layer
!
!------------------------------------------------------------------------
!

   do i = 1, ncol
      cld_found = .false.
      cld_layer(:) = cld(i,:) > 0.0_r8
      pmxrgn(i,:) = 0.0
      pnm(i,:)=pint(i,:)*10.
      n = 1
      do k = 1, pver
         if (cld_layer(k) .and.  .not. cld_found) then
            cld_found = .true.
         else if ( .not. cld_layer(k) .and. cld_found) then
            cld_found = .false.
            if (count(cld_layer(k:pver)) == 0) then
               exit
            endif
            pmxrgn(i,n) = pnm(i,k)
            n = n + 1
         endif
      end do
      pmxrgn(i,n) = pnm(i,pverp)
      nmxrgn(i) = n
   end do

   return
end subroutine cldovrlap

