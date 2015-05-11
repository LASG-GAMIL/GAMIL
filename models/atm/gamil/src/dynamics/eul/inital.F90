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
    ! $Id: inital.F90,v 1.15.2.6 2002/08/28 15:35:14 erik Exp $
    ! $Author: erik $
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid,       only: plon, plond, masterproc, &
        i1 ! added by ljli 2005/08/12
    use inidat,       only: read_inidat
    use prognostics,  only: ps, u3, v3, t3, q3, ptimelevels, initialize_prognostics, &
        q31, t31, q32, t32 ! added by LI Lijuan
    use buffer,       only: initialize_buffer
    use comsrf,       only: initialize_comsrf
    use ioFileMod,    only: getfil
    use radae,        only: initialize_radbuffer
    use phys_grid,    only: phys_grid_init
    use time_manager, only: timemgr_init
    use filenames,    only: ncdata
    use comfm1,       only: initialize_comfm1 ! WAN Hui 2003/10/23
#ifdef COUP_CSM
    use ccsm_msg, only: initialize_ccsm_msg
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
    integer n                ! index
    character(len=256) locfn ! local filename

    !-----------------------------------------------------------------------
    !
    ! Obtain initial dataset
    !
    if (masterproc) then
        call getfil (ncdata, locfn)
        call wrap_open (locfn, NF_NOWRITE, ncid_ini)
    end if

    call initialize_prognostics

    call initialize_comfm1      !!(wh 2003.10.23)
    !
    ! Check for consistent settings on initial dataset
    !
    call readinitial(ncid_ini)

    ! Initialize time manager.

    call timemgr_init

    !
    ! Initialize comslt and prognostics variables
    !
!!!   call initialize_prognostics
    !!   call initialize_comslt
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
    q31(:,:,:)       = q3(:,:,1,:,1)
    q32(:,:,:)       = q31(:,:,:)
    t31(:,:,:)       = t3(:,:,:,1)
    t32(:,:,:)       = t31(:,:,:)
    do n = 2, ptimelevels
        ps(:,:,n)     = ps(:,:,1)
        u3(:,:,:,n)   = u3(:,:,:,1)
        v3(:,:,:,n)   = v3(:,:,:,1)
        t3(:,:,:,n)   = t3(:,:,:,1)
        q3(:,:,:,:,n) = q3(:,:,:,:,1)
    end do

    return
end subroutine inital
