#include <misc.h>
#include <preproc.h>

#ifdef COUP_CSM

PROGRAM program_csm

!-----------------------------------------------------------------------
!
! Purpose:
! driver for CLM2.0 as the land component of CCSM
!
! Method:
! This program is the driver for CLM to work as the land component of
! CCSM.  The flux coupler will provide all the appropriate atmospheric
! forcing for the land model to run.
!
! o Currently, the land surface model grid is must be on the same model
!   grid as the atmospheric model when running in flux coupled mode.
!
! o CLM returns to the CCSM flux coupler surface
!   fluxes, temperatures, and albedos for only the land points on the
!   [lsmlon x lsmlat] grid
!
! o CLM uses its own grid dimensions (lsmlon and lsmlat)
!   and its own surface type data set. Since it processes only land points,
!   this data set must correctly define what points are land and what
!   are not land.
!
! o The zenith angle calculation is performed for the NEXT model
!   time step rather than the current time step.
!
! o The land model calculates its own net solar radiation. The net solar
!   radiation should equal that calculated by the atmospheric component.
!   If not, there is a problem.
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: program_csm.F90,v 1.10.8.5 2002/06/21 05:44:39 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar         !parameters
  use clm_varctl         !run control variables
  use shr_orb_mod        !orbital parameters and routines
  use shr_msg_mod        !csm message passing routines and variables
#if (defined SPMD)
  use spmdMod      , only : masterproc, iam, spmd_init
#else
  use spmdMod      , only : masterproc, iam
#endif
  use initializeMod, only : initialize
  use clm_csmMod   , only : csmstop_now
  use time_manager , only : advance_timestep, get_nstep
  implicit none
#include <gpt.inc>

! ----------------local variables ---------------------------------
  integer :: i,j        !loop indices
  integer :: nstep      !time step
  logical :: doalb      !true if surface albedo calculation time step

! Earth's orbital characteristics

  integer  :: iyear_AD  !Year (AD) to simulate above earth's orbital parameters for
  real(r8) :: eccen     !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) :: obliq     !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) :: mvelp     !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)

! Orbital information after call to routine shr_orbit_params

  real(r8) :: obliqr    !Earth's obliquity in radians
  real(r8) :: lambm0    !Mean longitude (radians) of perihelion at the vernal equinox
  real(r8) :: mvelpp    !Earth's moving vernal equinox longitude
  logical  :: log_print !true=> print diagnostics
! -----------------------------------------------------------------

!
! Initialize timing library.  2nd arg 0 means disable, 1 means enable
!
  call t_setoptionf (usrsys, 0)
  call t_initializef ()
!
! Determine input/output units
!
  call shr_msg_stdio ('lnd')
!
! Initialize MPI communication groups for flux coupler
!
  call shr_msg_init ('lnd')
  call shr_msg_groups ('lnd')
!
! Initialize intra-MPI communication stuff or
! set masterproc if not running in SPMD mode
!
#if (defined SPMD)
  call spmd_init()
#endif
!
! Initialize land model - initialize communication with flux coupler
!
  call initialize (eccen, obliqr, lambm0, mvelpp)

  call t_startf('lnd_timeloop')
!
! begin time stepping loop
!
  do
!
! doalb is true when the next time step is a radiation time step
! For example:
!      nstep dorad doalb
!        1     F     F
!        2     F     T
!        3     T     F
!
     nstep = get_nstep()
     doalb = ((irad==1 .and. nstep+1/=1) .or. (mod(nstep,irad)==0 .and. nstep+1/=1))
!
! Call land surface model driver
!
     call driver (doalb, eccen, obliqr, lambm0 ,mvelpp)
!
! determine if time to stop
!
     if (csmstop_now) exit
!
! increment time step
!
     call advance_timestep()

  end do

  call t_stopf ('lnd_timeloop')

!
! Exit gracefully
!
  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
  call shr_msg_finalize

  stop
end PROGRAM program_csm

#else

!The following is only here since empty file won't compile
subroutine program_csm_stub
  write(6,*) 'PROGRAM_CSM: this routine should not be called'
  stop 99
end subroutine program_csm_stub

#endif





