#include <misc.h>
#include <params.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  inital --- Define initial conditions for first run of case
!
! !INTERFACE:
subroutine inital

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use inidat
   use buffer
   use radae, only: initialize_radbuffer
   use pspect
   use prognostics
   use comsrf
   use dynamics_vars, only : dynamics_init
   use phys_grid, only: phys_grid_init
   use time_manager, only: timemgr_init, get_step_size
   use filenames, only: ncdata
#ifdef COUP_CSM
   use ccsm_msg, only: initialize_ccsm_msg
#endif
   use ioFileMod

!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comqfl.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'

! !DESCRIPTION:
!
!   Define initial conditions for first run of case
!
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.03.01      Acker         Modifications
!   96.04.01      Boville       Reviewed
!   01.06.17      Sawyer        Added call to dynamics_init
!   01.07.12      Sawyer        Added arguments to dynamics_init
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer k                  ! indices
   character(len=256) locfn   ! local filename
   real(r8) :: dtime          ! timestep size
!
!-----------------------------------------------------------------------
!
! Obtain initial dataset
!
   if (masterproc) then
      call getfil(ncdata, locfn)
      call wrap_open(locfn, NF_NOWRITE, ncid_ini)
   end if
!
! Check for consistent settings on initial dataset
!
   call readinitial(ncid_ini)

! Initialize time manager.

   call timemgr_init()

   dtime = get_step_size()
   call dynamics_init( dtime, iord, jord, nsplit )
!
! Initialize prognostics variables
!
   call initialize_prognostics
!
! Initialize commons
!
   call initcom
!
! Define physics data structures
!
   call phys_grid_init

#ifdef COUP_CSM
!
! Initialize ccsm arrays (must be done after phys_grid_init where
! begchunk and endchunk are defined
!
   call initialize_ccsm_msg
#endif
!
! Initialize buffer, comsrf, and radbuffer variables
! (which must occur after the call to phys_grid_init)
!
   call initialize_buffer
   call initialize_comsrf
   call initialize_radbuffer
!
! Read in initial data
!
   call read_inidat

   return
!EOC
end subroutine inital
!-----------------------------------------------------------------------
