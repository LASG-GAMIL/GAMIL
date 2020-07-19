#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
!
!   Initialize external models and/or boundary dataset information
!
! Method:
!
! Author:
!
!   CCM Core Group
!
!-----------------------------------------------------------------------

subroutine initext
!!(wh 2003.12.27)

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid
    use ppgrid,         only: begchunk, endchunk
    use phys_grid,      only: get_ncols_p, get_rlat_all_p, get_rlon_all_p,get_lat_all_p, get_lon_all_p
    use comsrf
    use rgrid
    use shr_orb_mod
    use ioFileMod
    use so4bnd
    use so4bnd_IPCC ! added by WAN Hui
    use commap
#ifndef COUP_CSM
    use ice_constants,  only: Tffresh
#endif
    use filenames,      only: bndtvo, bndtvs
    use physconst,      only: stebol
    use time_manager,   only: is_first_step, is_perpetual, &
                              get_curr_calday, get_curr_date, get_perp_date
#ifdef SPMD
    use mpishorthand
#endif

#ifdef COUP_CSM
    use ccsm_msg,       only: ccsmini
#else
    use atm_lndMod,     only: atmlnd_ini
#ifndef COUP_SOM
    use sst_data,       only: sstini, sstint, sstan, sst
    use ice_data,       only: iceini, iceint
    use atm_lndMod
#endif
#endif

    implicit none

#include <comlun.h>
#include <comctl.h>
#include <comsol.h>
    include 'netcdf.inc'

#ifndef COUP_CSM
    integer  i                  ! indices
    integer  ncol               ! number of columns in current chunk
    real(r8) coszrs(pcols)      ! Cosine solar zenith angle
    real(r8) clat1(pcols)       ! Current latitude(radians)
    real(r8) clon1(pcols)       ! Current longitude(radians)
    integer  sghid              ! NetCDF sgh field id
    logical  oro_hires          ! true => ORO came from high res topo file
    logical  log_print          ! Flag to print out log information or not
    integer  ret                ! NetCDF returned status
    integer  attlen             ! NetCDF attribute length
    character(256) text         ! NetCDF attribute
#endif
    character(256) locfn        ! netcdf local filename to open
    character(4) ncnam(5)
    integer  yr, mon, day, tod  ! components of a date
    real(r8) calday             ! current calendar day
    integer  lchnk
    integer  lats(pcols)
    integer  lons(pcols)
    real(r8) tssav(pcols,begchunk:endchunk) ! cam surface temperatures

    calday = get_curr_calday()
    !
    !----------------------------------------------------------------------
    ! 1. Obtain datasets
    !----------------------------------------------------------------------
    !
    ! Obtain time-variant ozone and sst datatsets and do initial read of
    ! ozone dataset
    !
    if (.not. ideal_phys) then
        if (masterproc) then
            call getfil(bndtvo, locfn)
            call wrap_open(locfn, 0, ncid_oz)
            write(6, "('Notice: initext: ')", advance="no")
            write(6, "('wrap_open returns ncid ', I8)", advance="no") ncid_oz
            write(6, "(' for file ', A)") trim(locfn)
        end if
#ifndef COUP_CSM
        if (.not. aqua_planet .and. .not. adiabatic) then
            if (masterproc) then
                call getfil(bndtvs, locfn)
                call wrap_open(locfn, 0, ncid_sst)
                write(6, "('Notice: initext: ')", advance="no")
                write(6, "('wrap_open returns ncid ', I8)", advance="no") ncid_sst
                write(6, "(' for file ', A)") trim(locfn)
            end if
        end if
#endif
        call oznini
    end if
    !
    !----------------------------------------------------------------------
    ! 2. Obtain sulfate aerosol datasets
    !----------------------------------------------------------------------
    !
    if (doRamp_so4) then
        call sulfini
    end if

    if (doIPCC_so4) then
        call sulfini_IPCC ! added by WAN Hui
    end if

#ifndef COUP_CSM
    !
    !----------------------------------------------------------------------
    ! 3. Determine if SGH field came from hi-res dataset
    !----------------------------------------------------------------------
    !
    if (is_first_step()) then
        if (masterproc) then
            call wrap_inq_varid(ncid_ini, 'SGH', sghid)
            ret = nf_inq_attlen(ncid_ini, sghid, 'from_hires', attlen)
            if (ret == nf_noerr .and. attlen > 256) then
                write(6, "('Error: initext: attribute length of ""from_hires"" is too long')")
                call endrun
            end if
            ret = nf_get_att_text(ncid_ini, sghid, 'from_hires', text)
            if (ret == nf_noerr .and. text(1:4) == 'true') then
                oro_hires = .true.
                write(6, "('Notice: initext: attribute ""from_hires"" is true')")
                write(6, "('Notice: initext: ""tssub"" will be used to guess sea ice')")
            else
                oro_hires = .false.
                write(6, "('Notice: initext: attribute ""from_hires"" is either false or not present')")
                write(6, "('Notice: initext: where sea ice exists, its initial temperature will be just below freezing')")
            end if
        end if
#if ( defined SPMD )
        call mpibcast(oro_hires, 1, mpilog, 0, mpicom)
#endif
    end if
    !
    !----------------------------------------------------------------------
    ! 4. Setup the characteristics of the orbit (Based on the namelist parameters)
    !----------------------------------------------------------------------
    !
    if (masterproc) then
        log_print = .true.
    else
        log_print = .false.
    end if
    call shr_orb_params(iyear_AD, eccen, obliq , mvelp, obliqr, lambm0, mvelpp, log_print)
    !
    !----------------------------------------------------------------------
    ! 5. Initialize land model
    !----------------------------------------------------------------------
    !
    ! This involves initializing land albedos, surface temperature, lwup and snowh.
    ! NOTE: On restart, lwup, ts, albedos and snowh, come from the atm restart data.
    !
    if (is_first_step()) then
        call srfflx_state_reset(srfflx_state2d)
    end if
    if (.not. adiabatic .and. .not. ideal_phys .and. .not. aqua_planet) then
        call atmlnd_ini(srfflx_parm2d)
    end if
    !
    ! Save off ts here because it is needed below for calculating lwup.  The
    ! updated state values are summed according to the fractional area of the
    ! underlying surface and only represent an actual grid box value after the
    ! last surface process has been called. TS is a special case as it is
    ! calculated from lwup and overwritten each time update surface
    ! fluxes is called.  The intermediate ts values returned from the update
    ! routine are wrong until ts is calculated from the full gridbox value of lwup
    ! lwup is only complete after the last surface process is called or if we
    ! are calculating a grid point that is all land.  We will save the
    ! intermediate value returned from land for our calculations below.
    !
    do lchnk = begchunk, endchunk
        tssav(:,lchnk) = srfflx_parm2d(lchnk)%ts(:)
    end do
    ! LIU Li: only call update_srf_fluxes at initial run
    if (is_first_step()) then
        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, landfrac)
    end if

    !----------------------------------------------------------------------
    ! 6. Initialize ocean and ice model
    !----------------------------------------------------------------------
#ifdef COUP_SOM
    !
    ! Slab ocean model: set initial surf temps for initial run. Read in 2 time slices of
    ! mixed layer depths and q fluxes from boundary dataset whether initial or restart
    !
    call somini(oro_hires)

#else
    !
    ! Data ocean model: Initialize ocean/sea-ice surface datasets and determine initial sea surface
    ! temperature
    !
    if (.not. adiabatic .and. .not. ideal_phys) then
        call sstini
        call iceini
        call sstint
        call iceint
    else
        icefrac(:pcols,begchunk:endchunk) = 0.0
        call update_srf_fractions
    end if
    !
    !----------------------------------------------------------------------
    ! 7. Initialize surface and sub-surface temperatures, set new sea ice
    !    concentrations and compute longwave up over non-land
    !----------------------------------------------------------------------
    !
    if (is_first_step()) then
        do lchnk = begchunk, endchunk
            if (.not. adiabatic .and. .not. ideal_phys) then
                ncol = get_ncols_p(lchnk)
                do i = 1, ncol
                    srfflx_state2d(lchnk)%ts(i) = &
                        landfrac(i,lchnk)*tssav(i,lchnk) + &
                        icefrac(i,lchnk)*tsice(i,lchnk) + &
                        ocnfrac(i,lchnk)*(sst(i,lchnk)+Tffresh)
                    if (landfrac(i,lchnk).ne.1.) then
                        srfflx_state2d(lchnk)%lwup(i) = &
                            stebol*(srfflx_state2d(lchnk)%ts(i)**4)
                    end if
                end do
            end if
        end do
    end if
#endif
    !
    !----------------------------------------------------------------------
    ! 8. Initialize non-land albedos at NSTEP = 0.  At NSTEP = 1 and
    !    beyond, albedos will be computed for the *next* timestep to
    !    accomodate coupling with a single interface.
    !----------------------------------------------------------------------
    !
    if (is_first_step()) then
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            call get_rlat_all_p(lchnk, ncol, clat1)
            call get_rlon_all_p(lchnk, ncol, clon1)
            call zenith(calday, clat1, clon1, coszrs, ncol)
            call albocean(lchnk, ncol, coszrs, &
                          srfflx_parm2d(lchnk)%asdir, srfflx_parm2d(lchnk)%aldir, &
                          srfflx_parm2d(lchnk)%asdif, srfflx_parm2d(lchnk)%aldif)
        end do

        call update_srf_fluxes(srfflx_state2d, srfflx_parm2d, ocnfrac)

        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            call get_lat_all_p(lchnk, ncol, lats)
            call get_lon_all_p(lchnk, ncol, lons)
            call get_rlat_all_p(lchnk, ncol, clat1)
            call get_rlon_all_p(lchnk, ncol, clon1)
            call zenith (calday, clat1, clon1, coszrs, ncol)
            call albice(lchnk, ncol, tsice(1,lchnk), snowhice(1,lchnk), coszrs, &
                        srfflx_parm2d(lchnk)%asdir, srfflx_parm2d(lchnk)%aldir, &
                        srfflx_parm2d(lchnk)%asdif, srfflx_parm2d(lchnk)%aldif)
            !
            ! fill in ice albedoes for therm ice model
            !
            asdirice(:ncol,lchnk)= srfflx_parm2d(lchnk)%asdir(:ncol)
            aldirice(:ncol,lchnk)= srfflx_parm2d(lchnk)%aldir(:ncol)
            asdifice(:ncol,lchnk)= srfflx_parm2d(lchnk)%asdif(:ncol)
            aldifice(:ncol,lchnk)= srfflx_parm2d(lchnk)%aldif(:ncol)
        end do
        call update_srf_fluxes(srfflx_state2d,srfflx_parm2d,icefrac)
    end if
#endif

#ifdef COUP_CSM
    !
    ! 1. Initial communications with coupler
    !
    call ccsmini
#endif

    return
end subroutine initext
