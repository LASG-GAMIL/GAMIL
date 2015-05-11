#include <misc.h>
#include <params.h>

logical function rstwr ()
!-----------------------------------------------------------------------
!
! Purpose: Determine whether it is time to write restart files
!
! Method:
!    If the model is run through flux coupler:
!     o At end of current day after receiving message from coupler.
!    If the model is not run through flux coupler:
!     o If end of run
!     o If monthly average primary and are on month boundary
!     o If not a monthly average primary and writing primary tape
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc
   use history, only:  nhtfrq, hstwr, nfils, mfilt
   use time_manager, only: get_nstep, is_end_curr_day
#ifdef COUP_CSM
   use ccsm_msg, only: csmstop, csmrstrt
#endif

   implicit none

#include <comctl.h>
!-----------------------------------------------------------------------

   rstwr = .false.

#ifdef COUP_CSM

   if (csmrstrt) then
      if (is_end_curr_day()) rstwr = .true.
   end if

#else

   if (nlend) then
      rstwr = .true.
      if (masterproc) then
         write (6,*) 'RSTWR: returning true due to nlend'
      end if
   else if (nhtfrq(1) == 0) then
      if (hstwr(1) .and. get_nstep() > 1) rstwr=.true.
   else
      if (hstwr(1) .and. nfils(1) >= mfilt(1) .and. get_nstep() > 1) then
         rstwr = .true.
         if (masterproc) then
            write (6,*) 'RSTWR: returning true due to full volume'
         end if
      end if
   end if

#endif

end function rstwr
