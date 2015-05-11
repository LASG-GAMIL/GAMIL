#include <misc.h>
#include <params.h>

subroutine preset

!! (wh 2003.04.28)
!! (wh 2004.04.14)
!-----------------------------------------------------------------------
!
! Purpose: Preset namelist ATMEXP input variables and initialize some other variables
!
! Method: Hardwire the values
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use infnan,       only: inf
   use history,      only: fincl, fexcl, fhstpr, fwrtpr
   use pmgrid
   use tracers,      only: nusr_adv, nusr_nad
   use constituents, only: ch4vmr, n2ovmr, f11vmr, f12vmr, co2vmr
!! use pspect
   use rgrid
   use shr_orb_mod
   use dycore
   use comhd, only: dfs0          !! (wh 2004.04.14)
#if ( ! defined COUP_SOM ) && ( ! defined COUP_CSM )
   use ice_dh, only: prognostic_icesnow,reset_csim_iceprops
#endif
   use time_manager, only: timemgr_preset
!-----------------------------------------------------------------------
   implicit none
!------------------------------Commons----------------------------------
#include <comadj.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
!!#include <comhd.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comtfc.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
#include <perturb.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Preset character history variables here because module initialization of character arrays
! does not work on all machines
!
   fincl(:,:)  = ' '
   fexcl(:,:)  = ' '
   fhstpr(:,:) = ' '
   fwrtpr(:,:) = ' '
!
! Flags
!
   nlend       = .false.       ! end of run flag
   nlres       = .false.       ! continuation run flag
   nlhst       = .false.       ! regen run or branch run flag
   lbrnch      = .false.       ! branch run flag
   adiabatic   = .false.       ! no physics
   ideal_phys  = .false.       ! "idealized" model configuration
   aqua_planet = .false.       ! global oceans/analytical SST's
   print_step_cost = .false.   ! print per timestep cost info
!
! Tracers
!
   trace_gas   = .false.          ! greenhouse gas code not implemented
   trace_test1 = .false.          ! test tracer 1 code not implemented
   trace_test2 = .false.          ! test tracer 2 code not implemented
   trace_test3 = .false.          ! test tracer 3 code not implemented
   readtrace   = .true.           ! initialize from initial conditions file
   nusr_adv    = 0                ! zero user advected tracers
   nusr_nad    = 0                ! zero user non-advected tracers
!
! Ice flags
!
#if ( ! defined COUP_SOM ) && ( ! defined COUP_CSM )
   prognostic_icesnow = .true.    ! snow falls on ice by default but
                                  ! it is limited to 0.5 meter.
   reset_csim_iceprops= .false.   ! use initial condition info unless
                                  ! need to reset ice properties in csim
#endif
!
   trace_gas   = .false.          ! greenhouse gas code not implemented
!
! Default run type is initialization
!
   nsrest = 0
!
! Default value for writing restart files
!
   nrefrq = 1      ! normal run, dispose with full history file
!
! Cycling flags for input boundary data files
!
   sstcyc = .true.
   icecyc = .true.
   ozncyc = .true.
!
! Model time defaults
!
   call timemgr_preset()
!
! Frequency in iterations of absorptivity/emissivity calc (negative
! values in model hours)
!
   iradae = -12
!
! Frequency of annual cycle sst update
!
   itsst  =  1
!
! Default frequency of shortwave and longwave radiation computations:
! once per hour (negative value implies model hours)
!
   iradsw = -1
   iradlw = -1
!
! Numerical scheme default values
!
   eps    = 0.06
!
! The Courant limiter is only used in eul dycore.
!
!!   kmxhdc = 5
   nlvdry = 3
!
!! Horizontal diffusion is used in both eul and sld dycores, but for sld
!! the parameters are zero in the production version.
!! Horizontal diffusion is not used in the fv dycore.
!!
!!   if (dycore_is ('EUL')) then
!!
!!     Set dif2 and dif4 for T42 resolution otherwise require setting on namelist
!!
!!      dif2   = 2.5e5   ! Leave dif2 to this value
!!      if ( (plon == 128) .and. (plat == 64) )then
!!        dif4   = 1.e16
!!      else
!!        dif4   = inf
!!      end if
!!   else
!!      dif2   = 0.
!!      dif4   = 0.
!!   end if
!!
    if (dycore_is ('EUL')) dfs0 = 0.02           !!(wh 2004.04.14)
!
!! No divergence damping
!!
!!   divdampn = 0.
!
! Precipitation threshold for PRECCINT, PRECLINT, PRECCFRQ, and PRECLFRQ output fields
! (mm/hr)
!
   precc_thresh = 0.1
   precl_thresh = 0.05
!
! Initialize volume mixing ratios
!
   ch4vmr = 1.714e-6
   n2ovmr = 0.311e-6
   f11vmr = 0.280e-9
   f12vmr = 0.503e-9
   co2vmr = 3.550e-4
!
! Orbital parameters.
! NOTE: if iyear_AD is set to SHR_ORB_UNDEF_INT after namelist input
! then namelist values of obliq,eccen,and mvelp are used otherwise
! obliq,eccen and mvelp are calculated based on iyear_AD
!
   iyear_ad = shr_orb_undef_int
   obliq    = shr_orb_undef_real
   eccen    = shr_orb_undef_real
   mvelp    = shr_orb_undef_real
!
! Solar constant
!
   scon       = 1.367e6
!
! Visible optical depth
!
   tauvis = .14

#ifdef COUP_CSM
!
! Communications with the flux coupler
!
   flxave = .true.
#endif
!
! rgrid: set default to full grid
!
   nlon(:) = plon
!
! Unit numbers: set to invalid
!
   nsds     = -1
   nrg      = -1
   nrg2     = -1
   ncid_ini = -1
   ncid_oz  = -1
   ncid_sst = -1
   ncid_trc = -1
   luhrest  = -1
!
! /perturb/
!
  pertlim = 0.0

   return
end subroutine preset
