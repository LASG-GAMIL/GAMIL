#include <misc.h>
#include <params.h>

subroutine rampnl_ghg( year )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc
   implicit none

#include <ramp.h>
!------------------------- Input args. ---------------------------------
   integer, intent(in) :: year ! Ramped gases fixed at this year
!-----------------------------------------------------------------------
   rampYear_ghg = year
   fixYear_ghg = .false.
   if ( year > 0 ) then
      fixYear_ghg = .true.
      if (masterproc) &
         write(6,*) 'RAMP_GHG: Ramped gases being fixed at year ',rampYear_ghg
   end if
   return
end subroutine rampnl_ghg

!##############################################################################

subroutine ramp_ghg
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes greenhouse gas volume mixing ratios via interpolation of
! yearly input data.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use physconst,    only: mwdry, mwco2
   use pmgrid,       only: masterproc
   use constituents, only: co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr
   use time_manager, only: get_curr_date, get_curr_calday

   implicit none

#include <ramp.h>
#include <crdcon.h>
#include <ramp_ghg_bau.h>
!---------------------------Local variables-----------------------------

   real(r8), parameter :: rmwco2 = mwco2/mwdry   ! ratio of molecular weights of co2 to dry air

   integer yrmodel           ! model year
   integer nyrm              ! year index
   integer nyrp              ! year index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) doymodel             ! model day of year
   real(r8) doydatam             ! day of year for input data yrdata(nyrm)
   real(r8) doydatap             ! day or year for input data yrdata(nyrp)
   real(r8) deltat               ! delta time
   real(r8) fact1, fact2         ! time interpolation factors
   real(r8) cfcscl               ! cfc scale factor for f11

!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   if (ramp_write) then
      if ( masterproc ) &
         write(6,*) ramp_type
      ramp_write = .false.
   endif
!
! determine index into input data
!
   if ( fixYear_ghg ) then
      yrmodel  = rampYear_ghg
   else
      yrmodel  = ncdate/10000
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is before yrdata(1), quit
!
   if (nyrm < 1) then
      if ( masterproc )then
         write(6,*)'RAMP_GHG: data time index is out of bounds'
         write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      end if
      call endrun
   endif
!
! if current date later than yrdata(ntim), just use ntim values
! if want to use just use ntim values - uncomment the following lines
! below and comment the call to endrun and previous write
!
   if (nyrp > ntim) then
      if ( masterproc )then
         write(6,*)'RAMP: error - current date is past the end of ', &
                ' valid sulfate scale factor data'
      end if
      call endrun
!         write(6,*)'RAMP_GHG: using ghg data for ',yrdata(ntim)
!         co2vmr = co2(ntim)*1.e-06
!         ch4vmr = ch4(ntim)*1.e-09
!         n2ovmr = n2o(ntim)*1.e-09
!         f11vmr = f11(ntim)*1.e-12*(1.+cfcscl)
!         f12vmr = f12(ntim)*1.e-12
!         co2mmr = rmwco2 * co2vmr
!         return
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365.    + calday
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat

   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. &
       fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. &
       fact2 < -1.e-6) then
      if ( masterproc ) &
         write(6,*)'RAMP_GHG: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! do time interpolation:
!   co2     in ppmv
!   n2o,ch4 in ppbv
!   f11,f12 in pptv
!
   co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06
   ch4vmr = (ch4(nyrm)*fact1 + ch4(nyrp)*fact2)*1.e-09
   n2ovmr = (n2o(nyrm)*fact1 + n2o(nyrp)*fact2)*1.e-09

   cfcscl = (adj(nyrm)*fact1 + adj(nyrp)*fact2)
   f11vmr = (f11(nyrm)*fact1 + f11(nyrp)*fact2)*1.e-12*(1.+cfcscl)
   f12vmr = (f12(nyrm)*fact1 + f12(nyrp)*fact2)*1.e-12

   co2mmr = rmwco2 * co2vmr
!
! output statements of ramping
!
   if ( masterproc )then
      write(6,'(a,f8.2,6(1pe22.14))') 'calday1 = ',calday,co2vmr/1.e-06,ch4vmr/1.e-09, &
                                      n2ovmr/1.e-09
      write(6,'(a,f8.2,6(1pe22.14))') 'calday2 = ',calday,cfcscl, &
                                      (f11(nyrm)*fact1 + f11(nyrp)*fact2),f12vmr/1.e-12
   end if
   return
end subroutine ramp_ghg
