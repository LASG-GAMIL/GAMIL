#include <misc.h>
#include <preproc.h>

#if (defined OFFLINE)

PROGRAM program_off

!----------------------------------------------------------------------- 
! 
! Purpose: 
! This program is the driver for CLM2.0 to work as "off-line" code 
! to mimic coupling to an atmospheric model.
!
! Method: 
! This code can be used to run CLM2 in offline mode, where the 
! appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
!
! o The land surface model may use a different grid than the input 
!   atmospheric data. The atmospheric data is then interpolated to the 
!   land model grid inside the atmospheric driver module [atmdrvMod.F90].
!
! o To map from the atmospheric grid to the land grid, the atmospheric
!   datasets must provide latitudes and longitudes (degrees) for each grid
!
! o The zenith angle calculation is for the NEXT time step rather 
!   than the current time step. The calendar day must therefore be
!   the NEXT time step and is for Greenwich time.
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: program_off.F90,v 1.6.2.6 2002/06/21 05:44:39 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_orb_mod  , only : SHR_ORB_UNDEF_REAL, shr_orb_params
  use clm_varctl   , only : irad, nsrest   
  use initializeMod, only : initialize
  use atmdrvMod    , only : atmdrv
#if (defined SPMD)
  use spmdMod      , only : masterproc, iam, spmd_init
  use mpishorthand , only : mpicom
#else
  use spmdMod      , only : masterproc, iam
#endif
  use time_manager , only : is_last_step, advance_timestep, get_nstep 
  implicit none

#include <gpt.inc>

! -------------------- local variables ---------------------------
  logical doalb     !true if surface albedo calculation time step
  integer nstep     !time step 
  integer ier       !error code

! Earth's orbital characteristics

  integer iyear_AD  !Year (AD) to simulate above earth's orbital parameters for
  real(r8) eccen    !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) obliq    !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) mvelp    !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)

! Orbital information after call to routine shr_orbit_params

  real(r8) obliqr   !Earth's obliquity in radians
  real(r8) lambm0   !Mean longitude (radians) of perihelion at the vernal equinox 
  real(r8) mvelpp   !Earth's moving vernal equinox longitude
  logical log_print !true=> print diagnostics  
! -----------------------------------------------------------------

!
! Initialize timing library.  2nd arg 0 means disable, 1 means enable
!
  call t_setoptionf (usrsys, 1)
  call t_initializef ()

#if (defined SPMD)

!
! Initialize intra-MPI communication stuff 
!
  call spmd_init
#endif

! 
! Initialize orbital parameters.
! variables obliq, eccen and mvelp determined based on value of 
! iyear_AD
! 

   if (masterproc) then
     log_print = .true.
   else
     log_print = .false.
   end if
   iyear_AD = 1950
   obliq    = SHR_ORB_UNDEF_REAL
   eccen    = SHR_ORB_UNDEF_REAL
   mvelp    = SHR_ORB_UNDEF_REAL
   call shr_orb_params (iyear_AD, eccen, obliq, mvelp, obliqr, &
                        lambm0, mvelpp, log_print)
!
! Initialize CLM2 
!
  call initialize (eccen, obliqr, lambm0, mvelpp)   

  call t_startf('total')
!
! Begin time stepping loop
!
  do 
!
! Obtain atmospheric state and fluxes.
!
     nstep = get_nstep()
     call atmdrv (nstep)    
!
! Determine if albedo calculation is to be done (only done when
! the next time step is a radiation time step). For example: 
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!
     doalb = (irad==1 .or. (mod(nstep,irad)==0 .and. nstep+1/=1))
!
! Call CLM2 driver (only sets land points)
!
     call driver (doalb, eccen, obliqr, lambm0, mvelpp)
!
! Determine if time to stop
!
     if (is_last_step()) exit
!
! Increment time step
!
     call advance_timestep()

  end do
  call t_stopf('total')

!
! Exit gracefully
!

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
#if (defined SPMD) 
  call mpi_barrier (mpicom, ier)
  call mpi_finalize(ier)
#endif

  stop
end PROGRAM program_off

#else

!The following is only here since empty file won't compile
subroutine program_off_stub
  write(6,*) 'PROGRAM_OFF: this routine should not be called'
  return
end subroutine program_off_stub

#endif
