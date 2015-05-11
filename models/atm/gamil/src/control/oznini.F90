#include <misc.h>
#include <params.h>

subroutine oznini
!----------------------------------------------------------------------- 
! 
! Purpose: Do initial read of time-variant ozone boundary dataset, containing
!          ozone mixing ratios as a function of latitude and pressure.  Read two
!          consecutive months between which the current date lies.  Routine
!          RADOZ2 then evaluates the two path length integrals (with and without
!          pressure weighting) from zero to the interfaces between the input
!          levels.  It also stores the contribution to the integral from each
!          layer.
! 
! Method: Call appropriate netcdf wrapper routines and interpolate to model grid
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
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Local workspace
!
   integer dateid                          ! netcdf id for date variable
   integer secid                           ! netcdf id for seconds variable
   integer londimid                        ! netcdf id for longitude dimension
   integer latdimid                        ! netcdf id for latitude dimension
   integer levdimid                        ! netcdf id for level dimension
   integer lonid                           ! netcdf id for longitude variable
   integer latid                           ! netcdf id for latitude variable
   integer levid                           ! netcdf id for level variable
   integer timeid                          ! netcdf id for time variable
   integer dimids(nf_max_var_dims)         ! variable shape
   integer cnt4(4)                         ! array of counts for each dimension
   integer strt4(4)                        ! array of starting indices
   integer i, k, lat, n                    ! longitude, level, latitude, time indices
   integer  :: yr, mon, day                ! components of a date
   integer  :: ncdate                      ! current date in integer format [yyyymmdd]
   integer  :: ncsec                       ! current time of day [seconds]
   real(r8) :: calday                      ! current calendar day
   real(r8) caldayloc                      ! calendar day (includes yr if no cycling)
   real(r8), allocatable :: tmpozmix(:,:)  ! temporary ozmix array
   real(r8), allocatable :: oznbdym(:,:,:) ! ozone data previous time sample
   real(r8), allocatable :: oznbdyp(:,:,:) ! ozone data next time sample

   nm = 1
   np = 2
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
   if (masterproc) then
!
! Use year information only if not cycling ozone dataset
!
      calday = get_curr_calday()
      if ( is_perpetual() ) then
         call get_perp_date(yr, mon, day, ncsec)
      else
         call get_curr_date(yr, mon, day, ncsec)
      end if
        yr  = yr+timediff  ! ljli
      ncdate = yr*10000 + mon*100 + day
      if (ozncyc) then
         caldayloc = calday
      else
         caldayloc = calday + yr*365.
      end if
!
! Get and check dimension info
!
      CALL WRAP_INQ_DIMID( ncid_oz, 'lon', londimid   )
      CALL WRAP_INQ_DIMID( ncid_oz, 'lev', levdimid   )
      CALL WRAP_INQ_DIMID( ncid_oz, 'time', timeid  )
      CALL WRAP_INQ_DIMID( ncid_oz, 'lat', latdimid   )

      CALL WRAP_INQ_DIMLEN( ncid_oz, londimid, lonsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, levdimid, levsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, latdimid, latsiz   )
      CALL WRAP_INQ_DIMLEN( ncid_oz, timeid, timesiz   )

      CALL WRAP_INQ_VARID( ncid_oz, 'date', dateid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'datesec', secid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'OZONE', oznid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lon', lonid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lat', latid   )
      CALL WRAP_INQ_VARID( ncid_oz, 'lev', levid   )

      CALL WRAP_INQ_VARDIMID (ncid_oz, oznid, dimids)
      if (dimids(1) /= londimid .or. dimids(2) /= levdimid .or. dimids(3) /= latdimid) then
         write(6,*)'OZNINI: Data must be ordered lon, lev, lat, time'
         call endrun
      end if
!
! Dynamically allocated memory for module comozp 
!
      allocate (ozlon(lonsiz))
      allocate (ozlat(latsiz))
      allocate (date_oz(timesiz))
      allocate (sec_oz(timesiz))
      allocate (pin(levsiz))
      allocate (ozmixm(plond,levsiz,plat,2))
      allocate (ozmix(plond,levsiz,plat))
!
! Locally dynamic that will be deallocated before "return"
!
      allocate (oznbdym(lonsiz,levsiz,latsiz))
      allocate (oznbdyp(lonsiz,levsiz,latsiz))
      allocate (tmpozmix(levsiz,plat))
!
! Retrieve longitude, latitude and level arrays for interpolation.
!
      CALL WRAP_GET_VAR_REALX (NCID_OZ, lonid,ozlon)
      CALL WRAP_GET_VAR_REALX (NCID_OZ, latid,ozlat)
      CALL WRAP_GET_VAR_REALX (NCID_OZ, levid,pin)
!
! Convert from millibars to pascals
!
      do k=1,levsiz
         pin(k) = pin(k)*100.
      end do
!
! Retrieve entire date and sec variables.
!
      CALL WRAP_GET_VAR_INT (ncid_oz,dateid,date_oz)
      CALL WRAP_GET_VAR_INT (ncid_oz,secid,sec_oz)
      if (ozncyc) then
         if (timesiz < 12) then 
            write(6,*)'OZNINI: When cycling ozone, dataset must have 12 consecutive ', &
                      'months of data starting with Jan'
            write(6,*)'Current dataset has only ',timesiz,' months'
            call endrun
         end if
         do n = 1,12
            if (mod(date_oz(n),10000)/100 /= n) then
               write(6,*)'OZNINI: When cycling ozone, dataset must have 12 consecutive ', &
                         'months of data starting with Jan'
               write(6,*)'Month ',n,' of dataset says date=',date_oz(n)
               call endrun
            end if
         end do
      end if

      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1
!
! Special code for interpolation between December and January
!
      if (ozncyc) then
         n = 12
         np1 = 1
         call bnddyi(date_oz(n  ), sec_oz(n  ), cdayozm)
         call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
         if (caldayloc <= cdayozp .or. caldayloc > cdayozm) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdym)

            strt4(4) = np1
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            goto 10
         end if
      end if
!
! Normal interpolation between consecutive time slices.
!
      do n=1,timesiz-1
         np1 = n + 1
         call bnddyi(date_oz(n  ), sec_oz(n  ), cdayozm)
         call bnddyi(date_oz(np1), sec_oz(np1), cdayozp)
         if (.not.ozncyc) then
            yr = date_oz(n)/10000
            cdayozm = cdayozm + yr*365.
            yr = date_oz(np1)/10000
            cdayozp = cdayozp + yr*365.
         end if
         if (caldayloc > cdayozm .and. caldayloc <= cdayozp) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdym)

            strt4(4) = np1
            call wrap_get_vara_realx (ncid_oz,oznid,strt4,cnt4,oznbdyp)
            goto 10
         end if
      end do
      write(6,*)'OZNINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
      call endrun
10    continue
      write(6,*)'OZNINI: Read ozone data for dates ',date_oz(n), &
                sec_oz(n),' and ',date_oz(np1),sec_oz(np1)
!
! Spatial interpolation.  If ozone dataset is 2-d (i.e. lonsiz = 1) and 
! thus only latitude interpolation is necessary, expand to 3-d after 
! interpolation.
!
      if (lonsiz == 1) then
         call lininterp (oznbdym ,ozlat   ,levsiz  ,latsiz  ,tmpozmix, &
                         latdeg  ,plat    )
         do lat=1,plat
            do k=1,levsiz
               do i=1,nlon(lat)
                  ozmixm(i,k,lat,nm) = tmpozmix(k,lat)
               end do
            end do
         end do

         call lininterp (oznbdyp ,ozlat   ,levsiz  ,latsiz  ,tmpozmix, &
                         latdeg  ,plat)
         do lat=1,plat
            do k=1,levsiz
               do i=1,nlon(lat)
                  ozmixm(i,k,lat,np) = tmpozmix(k,lat)
               end do
            end do
         end do

      else
         write(*,*) 'oznini, oznbdym'
         call bilin (oznbdym, ozlon, ozlat, lonsiz, lonsiz, &
                     levsiz, levsiz, latsiz, ozmixm(1,1,1,nm), londeg, &
                     latdeg, plond, nlon, levsiz, plat)

         write(*,*) 'oznini, oznbdyp'
         call bilin (oznbdyp, ozlon, ozlat, lonsiz, lonsiz, &
                     levsiz, levsiz, latsiz, ozmixm(1,1,1,np), londeg, &
                     latdeg, plond, nlon, levsiz, plat)
      end if
!
! Deallocate dynamic memory for local workspace.  NOT for pointers in common.
!
      deallocate (oznbdym)
      deallocate (oznbdyp)
      deallocate (tmpozmix)
   end if

#if ( defined SPMD )
   call mpibcast (levsiz, 1, mpiint, 0, mpicom)
   if (.not.masterproc) then
      allocate (pin(levsiz))
      allocate (ozmix(plond,levsiz,plat))
   end if
   call mpibcast (pin, levsiz, mpir8, 0, mpicom)
#endif

   return
end subroutine oznini


