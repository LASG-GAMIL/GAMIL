#include <misc.h>
#include <params.h>

subroutine oznint
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
! Method: Find next month of ozone data to interpolate.  Linearly interpolate 
!         vertically and horizontally
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use comozp
   use pspect
   use rgrid
   use commap
   use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                           is_perpetual
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

#include <comctl.h>
#include <comlun.h>
#include <forcingtimediff.h>
!
! Local workspace
!
   integer cnt4(4)                ! array of counts for each dimension
   integer strt4(4)               ! array of starting indices
   integer i, k, lat              ! longitude, level, latitude indices
   integer ntmp                   ! temporary
   integer :: yr, mon, day        ! components of a date
   integer :: ncdate              ! current date in integer format [yyyymmdd]
   integer :: ncsec               ! current time of day [seconds]

   real(r8) fact1, fact2          ! time interpolation factors
   real(r8) :: calday             ! current calendar day
   real(r8) caldayloc             ! calendar day (includes yr if no cycling)
   real(r8) tmpozmix(levsiz,plat) ! temporary ozmix array
   real(r8) deltat                ! time (days) between interpolating ozone data

   real(r8), allocatable :: oznbdyp(:,:,:) ! ozone data from next time sample
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
   if (masterproc) then
!
! Use year information only if a multiyear dataset
!
      calday = get_curr_calday()
      if ( is_perpetual() ) then
         call get_perp_date(yr, mon, day, ncsec)
      else
         call get_curr_date(yr, mon, day, ncsec)
      end if
        yr = yr + timediff  !ljli
      ncdate = yr*10000 + mon*100 + day
      if (ozncyc) then
         caldayloc = calday
      else
         caldayloc = calday + yr*365.
      end if

      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1
!
! If model time is past current forward ozone timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for ozncyc = .true. 
! interpolation between December and January (np1 == 1).  Note that 
! np1 is never 1 when ozncyc is .false.
!
      if (caldayloc > cdayozp .and. .not. (np1 == 1 .and. caldayloc > cdayozm)) then
!
! Allocate memory for dynamic local workspace
!
         allocate (oznbdyp(lonsiz,levsiz,latsiz))

         if (ozncyc) then
            np1 = mod(np1,12) + 1
         else
            np1 = np1 + 1
         end if
         if (np1 > timesiz) then
             np1 = np1-12
             yr = yr+1
!            write(6,*)'OZNINT: Attempt to read past end of O3 dataset'
!            call endrun
         end if
         cdayozm = cdayozp
         call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
         if (.not.ozncyc) then
            if (yr <= 2099) yr = date_oz(np1)/10000
            cdayozp = cdayozp + yr*365.
         end if
         if (np1 == 1 .or. caldayloc <= cdayozp) then
            ntmp = nm
            nm = np
            np = ntmp
            strt4(4) = np1
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            write(6,*)'OZNINT: Read ozone for date (yyyymmdd) ', date_oz(np1),' sec ',sec_oz(np1)
!
! Spatial interpolation.  If ozone dataset is only 2-d (i.e. lonsiz = 1) and 
! thus only latitude interpolation is necessary, expand to 3-d after 
! interpolation.
!
            if (lonsiz == 1) then
               call lininterp (oznbdyp, ozlat, levsiz, latsiz, tmpozmix, &
                               latdeg, plat)
               do lat=1,plat
                  do k=1,levsiz
                     do i=1,nlon(lat)
                        ozmixm(i,k,lat,np) = tmpozmix(k,lat)
                     end do
                  end do
               end do
            else
               call bilin (oznbdyp ,ozlon   ,ozlat   ,lonsiz  ,lonsiz  , &
                           levsiz  ,levsiz  ,latsiz  ,ozmixm(1,1,1,np) , londeg  , &
                           latdeg  ,plond   ,nlon    ,levsiz  ,plat    )
            end if
         else
            write(6,*)'OZNINT: Input ozone for date',date_oz(np1),' sec ',sec_oz(np1), &
                 'does not exceed model date',ncdate,' sec ',ncsec,' Stopping.'
            call endrun
         end if
!
! Deallocate dynamic memory for local workspace.
!
         deallocate (oznbdyp)
      end if
!
! Determine time interpolation factor.  Account for December-January 
! interpolation if cycling ozone dataset.  Again note that np1 is never 1 
! when ozncyc is false
!
      if (np1 == 1) then                    ! Dec-Jan interpolation
         deltat = cdayozp + 365. - cdayozm
         if (caldayloc > cdayozp) then      ! We're in December
            fact1 = (cdayozp + 365. - caldayloc)/deltat
            fact2 = (caldayloc - cdayozm)/deltat
         else                                ! We're in January
            fact1 = (cdayozp - caldayloc)/deltat
            fact2 = (caldayloc + 365. - cdayozm)/deltat
         end if
      else
         deltat = cdayozp - cdayozm
         fact1 = (cdayozp - caldayloc)/deltat
         fact2 = (caldayloc - cdayozm)/deltat
      end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
      if (abs(fact1+fact2-1.) > 1.e-6 .or. fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
          fact2 > 1.000001 .or. fact2 < -1.e-6) then

         write(6,*)'OZNINT: Bad fact1 and/or fact2=',fact1,fact2
         call endrun
      end if
!
! Time interpolation.
!
      do lat=1,plat
         do k=1,levsiz
            do i=1,nlon(lat)
               ozmix(i,k,lat) = ozmixm(i,k,lat,nm)*fact1 + ozmixm(i,k,lat,np)*fact2
            end do
         end do
      end do
   end if

#if ( defined SPMD )
   call mpibcast(ozmix, plond*levsiz*plat, mpir8, 0, mpicom)
#endif

   return
end subroutine oznint

