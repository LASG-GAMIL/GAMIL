#include <misc.h>
#include <params.h>

subroutine radozn (lchnk, ncol, pmid, o3vmr)
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone from current time-interpolated values to model levels
! 
! Method: Use pressure values to determine interpolation levels
! 
! Author: Bruce Briegleb
! 
!--------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,     only: get_lat_all_p, get_lon_all_p
   use comozp
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)

   real(r8), intent(out) :: o3vmr(pcols,pver) ! ozone volume mixing ratio
!
! local storage
!
   integer i                   ! longitude index
   integer k, kk, kkstart      ! level indices
   integer kupper(pcols)       ! Level indices for interpolation
   integer kount               ! Counter
   integer lats(pcols)         ! latitude indices
   integer lons(pcols)         ! latitude indices

   real(r8) dpu                ! upper level pressure difference
   real(r8) dpl                ! lower level pressure difference
!
! Initialize latitude indices
!
   call get_lat_all_p(lchnk, ncol, lats)
   call get_lon_all_p(lchnk, ncol, lons)
!
! Initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do

   do k=1,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=1,ncol
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=1,ncol
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,k) = (ozmix(lons(i),kupper(i)  ,lats(i))*dpl + &
                             ozmix(lons(i),kupper(i)+1,lats(i))*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
      do i=1,ncol
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,k) = ozmix(lons(i),1,lats(i))*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,k) = ozmix(lons(i),levsiz,lats(i))
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,k) = (ozmix(lons(i),kupper(i)  ,lats(i))*dpl + &
                          ozmix(lons(i),kupper(i)+1,lats(i))*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
         write(6,*)'RADOZN: Bad ozone data: non-monotonicity suspected'
         call endrun
      end if
35    continue
   end do
!!!!!!!!!!!!!!ljli2010for cmip5 ozone
   do k=10,pver
      do i=1,ncol
         if(o3vmr(i,k) > 1.e-5.and.o3vmr(i,k-1)<=1.e-5)  o3vmr(i,k) = o3vmr(i,k-1)
      end do
   end do


   return
end subroutine radozn

