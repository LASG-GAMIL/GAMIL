#include <misc.h>
#include <params.h>
subroutine cldsav(lchnk   ,ncol    , &
                  cld     ,pmid    ,cldtot  ,cldlow  ,cldmed  , &
                  cldhgh  ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
! Pressure ranges for the 3 cloud levels are specified.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none
!------------------------------Parameters-------------------------------
   real(r8) plowmax             ! Max prs for low cloud cover range
   real(r8) plowmin             ! Min prs for low cloud cover range
   real(r8) pmedmax             ! Max prs for mid cloud cover range
   real(r8) pmedmin             ! Min prs for mid cloud cover range
   real(r8) phghmax             ! Max prs for hgh cloud cover range
   real(r8) phghmin             ! Min prs for hgh cloud cover range
!
   parameter (plowmax = 120000.,plowmin = 70000., &
              pmedmax =  70000.,pmedmin = 40000., &
              phghmax =  40000.,phghmin =  5000.)

   real(r8) ptypmin(4)
   real(r8) ptypmax(4)

   data ptypmin /phghmin, plowmin, pmedmin, phghmin/
   data ptypmax /plowmax, plowmax, pmedmax, phghmax/
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: cld(pcols,pver)     ! Cloud fraction
   real(r8), intent(in) :: pmid(pcols,pver)    ! Level pressures
   real(r8), intent(in) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc

   integer, intent(in) :: nmxrgn(pcols)        ! Number of maximally overlapped regions
!
! Output arguments
!
   real(r8), intent(out) :: cldtot(pcols)       ! Total random overlap cloud cover
   real(r8), intent(out) :: cldlow(pcols)       ! Low random overlap cloud cover
   real(r8), intent(out) :: cldmed(pcols)       ! Middle random overlap cloud cover
   real(r8), intent(out) :: cldhgh(pcols)       ! High random overlap cloud cover

!
!---------------------------Local workspace-----------------------------
!
   integer i,k                  ! Longitude,level indices
   integer irgn                 ! Max-overlap region index
   integer ityp                 ! Type counter
   real(r8) clrsky              ! Max-random clear sky fraction
   real(r8) clrskymax           ! Maximum overlap clear sky fraction
   real(r8) cldamt(4)           ! Low, medium, high, and total cld amounts
!
!-----------------------------------------------------------------------
   do i=1,ncol
      do ityp = 1, 4
!
! Initialize region number
!
         irgn = 1
         do while (pmxrgn(i,irgn) < ptypmin(ityp) .and. irgn < nmxrgn(i))
            irgn = irgn + 1
         end do
!
! Compute cloud amount by estimating clear-sky amounts
!
         clrsky = 1.0
         clrskymax = 1.0

         do k = 1, pver
            if (pmid(i,k) >= ptypmin(ityp) .and. pmid(i,k) <= ptypmax(ityp)) then
               if (pmxrgn(i,irgn) < pmid(i,k) .and. irgn < nmxrgn(i)) then
                  irgn = irgn + 1
                  clrsky = clrsky * clrskymax
                  clrskymax = 1.0
               endif
               clrskymax = min(clrskymax,1.0-cld(i,k))
            endif
         end do
         cldamt(ityp) = 1.0 - (clrsky * clrskymax)
      end do
      cldtot(i) = cldamt(1)
      cldlow(i) = cldamt(2)
      cldmed(i) = cldamt(3)
      cldhgh(i) = cldamt(4)
   end do

   return
end subroutine cldsav

