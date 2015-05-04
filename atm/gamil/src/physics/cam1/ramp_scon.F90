#include <misc.h>
#include <params.h>

subroutine rampnl_scon( year )

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

   integer, intent(in) :: year ! Ramped gases fixed at this year
!-----------------------------------------------------------------------
   rampYear_scon = year
   fixYear_scon = .false.
   if ( year > 0 ) then
      fixYear_scon = .true.
      if (masterproc) &
         write(6,*) 'RAMP_SCON: Ramped gases being fixed at year ',rampYear_scon
   end if
   return
end subroutine rampnl_scon

!##############################################################################

subroutine ramp_scon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes ramping of solar constant
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: masterproc
   use time_manager, only: get_curr_date, get_curr_calday

   implicit none

#include <ramp.h>
#include <comsol.h>
#include <ramp_scon.h>
!---------------------------Local variables-----------------------------

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
!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   if (ramp_write) then
      if (masterproc) &
        write(6,*) ramp_type
      ramp_write = .false.
   endif
!
! determine index into input data
!
   if ( fixYear_scon ) then
      yrmodel  = rampYear_scon
   else
      yrmodel  = ncdate/10000
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is before 1870, quit
!
   if (nyrm < 1) then
      if (masterproc) then
         write(6,*)'RAMP_SCON: data time index is out of bounds'
         write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      end if
      call endrun
   endif
!
! if current date later than ntim, just use ntim values
!
   if (nyrp > ntim) then
      scon = sconst(ntim)
      return
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365.    + calday
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam !365
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat

   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. &
       fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. &
       fact2 < -1.e-6) then
      if (masterproc) &
         write(6,*)'RAMP_SCON: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! do time interpolation:
!
   scon = (sconst(nyrm)*fact1 + sconst(nyrp)*fact2)

   return
end subroutine ramp_scon

