#include <misc.h>
#include <params.h>

subroutine inital

!-----------------------------------------------------------------------
!
! Purpose:
! Define initial conditions for first run of case
!
! Method:
!
! Author:
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          B. Boville, April 1996
!
! $Id: inital.F90,v 1.13.2.6 2002/08/28 15:35:32 erik Exp $
! $Author: erik $
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, i1, plond, masterproc
   use inidat,       only: read_inidat
   use prognostics,  only: ps, u3, v3, t3, q3, div, ptimelevels, initialize_prognostics
   use buffer,       only: initialize_buffer
   use comslt,       only: initialize_comslt
   use comsrf,       only: initialize_comsrf
   use ioFileMod,    only: getfil
   use radae,        only: initialize_radbuffer
   use phys_grid,    only: phys_grid_init
   use time_manager, only: timemgr_init
   use filenames,    only: ncdata
#ifdef COUP_CSM
   use ccsm_msg,     only: initialize_ccsm_msg
#endif
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!---------------------------Local variables-----------------------------
!
   integer k,n              ! indices
   character(len=256) locfn ! local filename
!
!-----------------------------------------------------------------------
!
!
! Obtain initial dataset
!
   if (masterproc) then
      call getfil (ncdata, locfn)
      call wrap_open (locfn, NF_NOWRITE, ncid_ini)
   end if
!
! Check for consistent settings on initial dataset
!
   call readinitial (ncid_ini)

! Initialize time manager.

   call timemgr_init()
!
! Initialize comslt and prognostics variables
!
   call initialize_prognostics
   call initialize_comslt
!
! Set commons
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
!
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!
   do n=2,ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(i1:i1+plon-1,:,:,:,n) = q3(i1:i1+plon-1,:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
   end do

   return
end subroutine inital
