#include <misc.h>
#include <params.h>
subroutine initesttr( q3,nlon )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Initialize test tracers.  The test tracers are:
!
!    1) Radon, init to zero, surface fluxes from WCRP95, 5.5
!       day e-folding decay.
!    2) conserved unit tracer
!    3) ozone-like tracer, init to 1.e-9 above ~100mb, zero
!       elsewhere, re-zero the bottom level at each timestep.
! Note that:
!    o ixtrct   = index of radon advected tracer
!    o ixtrct+1 = index of conserved unit tracer
!    o ixtrct+2 = index of ozone-like tracer
!
!-------------------------Code History----------------------------------
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! Reviewed:
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use tracers, only: pcnst, pnats, ixtrct
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!
! Output arguments:
!
   real(r8), intent(out) :: q3(plond,plev,pcnst+pnats)    ! kg tracer/kg dry air
   integer, intent(in) :: nlon
!
!--------------------------Local Variables------------------------------
!
   integer i, k                       !  loop counters
!
!-----------------------------------------------------------------------
!
!
! Initialize radon tracer to zero.
!
   if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
      do k = 1, plev
         do i = 1, nlon
            q3(i,k,ixtrct) = 0.0
         end do
      end do
   end if
!
! Initialize conserved unit tracer.
!
   if ( trace_test2 .or. trace_test3 ) then
      do k = 1, plev
         do i = 1, nlon
            q3(i,k,ixtrct+1) = 1.0
         end do
      end do
   end if
!
! Initialize strat tracer to 1.e-9 above 100mb
!
   if ( trace_test3 ) then
      do k = 1, plev
         do i = 1, nlon
            q3(i,k,ixtrct+2) = 0.0
         end do
      end do
      do k = 1, 5
         do i = 1, nlon
            q3(i,k,ixtrct+2) = 1.e-9
         end do
      end do
   end if

   return
end subroutine initesttr
!
!#######################################################################
!
subroutine rndecay( lchnk, ncol, rn, deltat, rnsnk)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Radon decay.
!
!-------------------------Code History----------------------------------
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! Reviewed:
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   implicit none
!-------------------------Arguments--------------------------------------
!
! input args
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: rn(pcols,pver)       ! radon mixing ratio (kg/(kg moist air))
   real(r8), intent(in) :: deltat               ! time step
!
! output args
   real(r8), intent(out) :: rnsnk(pcols,pver)    ! conversion rate
!                               !              (kg rn /(s kg moist air))
!
!--------------------------Local Variables------------------------------
!
   integer i                 ! x index
   integer k                 ! z index
!
   real(r8) a                    ! lifetime
   parameter( a = 2.10e-6 )
!
!-----------------------------------------------------------------------
!
!   calculate tendencies using Euler Backward
!
   do k = 1,pver
      do i = 1,ncol
         rnsnk(i,k) = -rn(i,k)*a / (1. + a*deltat)
      end do
   end do
!
!
   return
!
end subroutine rndecay
!
!########################################################################
!
subroutine rnsfwcrp( lchnk, ncol, landfrac, flux )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Set surface fluxes for radon for WCRP95 RN-PB simulation.
!
!  The flux is specified non-zero over land between 60S - 70N, except
!  exclude Greenland.
!
!  Flux strength:
!  60S - 60N:  3.69e-21 kg/m^2/s
!  60N - 70N:  (3.69e-21)/2 kg/m^2/s
!
!  This land source is has been adjusted so that the total radon flux is
!  15 kg/yr for a T42 grid.
!
!-------------------------Code History----------------------------------
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! Reviewed:
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments-------------------------------------
!
! Input arguments:
!
   integer, intent(in) :: lchnk           ! chunk identifier
   integer, intent(in) :: ncol            ! number of atmospheric columns

   real(r8), intent(in) :: landfrac(pcols)! landfraction
!
! Output arguments:
!
   real(r8), intent(out) :: flux(pcols)    ! specified radon flux in kg/m^2/s
!
!--------------------------Local Variables------------------------------
!
   integer i      ! loop counter

   real(r8) rlat(pcols)                  ! current latitudes(radians)
   real(r8) rlon(pcols)                  ! current longitudes(radians)
   real(r8) landflx   ! land flux
   real(r8) landflxn  ! (land flux)/2
   real(r8) rad2deg   ! convert radians to degrees
   real(r8) latdeg    ! latitude in degrees
!
!--------------------------Statement functions--------------------------
!
   logical land
   land(i) = nint(landfrac(i)).gt.0.9999_r8
!
!-----------------------------------------------------------------------
!
!
!--------------------------Parameters-----------------------------------
!
   parameter( rad2deg = 360. / 6.283185308)
!
!------------------------------------------------------------------------
!
!      landflx = 3.69e-21
   landflx = 3.7796e-21   ! rescaled so total flux is 15 kg/yr (T42)
   landflxn = landflx/2.
!
   call get_rlat_all_p(lchnk, ncol, rlat)
   call get_rlon_all_p(lchnk, ncol, rlon)
   do i = 1, ncol
!
      flux(i) = 0.
      latdeg = rlat(i) * rad2deg
      if ( latdeg .ge. -60.  .and.  latdeg .le. 60. ) then    ! 60S - 60N
         if ( land(i) ) flux(i) = landflx
      else if ( latdeg .gt. 60. .and. latdeg .le. 70 ) then   ! 60N - 70N
         if (rlon(i)*rad2deg .le. 300.0) then            ! 0 - 300E excludes Greenland
            if ( land(i) ) flux(i) = landflxn
         end if
      end if
!
   end do
!
   return
!
end subroutine rnsfwcrp
