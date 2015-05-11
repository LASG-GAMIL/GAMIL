#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose: Module to handle all of the message passing to/from
!          the CCSM coupler for coupled simulations.
!
! Author: Erik Kluzek
!
! Adapted from: The "ccm_csm*.F90" series of Mariana Vertenstein
! Extensively modified by Rob Jacob and Tony Craig, 3/2002 to
!   interface with cpl6.
!
! CVS Id: $Id: ccsm_msg.F90,v 1.11.2.17.36.2 2006/01/26 22:32:03 mvr Exp $
!
!-----------------------------------------------------------------------

module ccsm_msg

#ifdef COUP_CSM
    use shr_kind_mod,   only: r8 => shr_kind_r8                                     ! atmospheric model precision
    use pmgrid,         only: plat, plon, beglat, endlat, plond, masterproc, iam    ! model grid
    use ppgrid,         only: pcols, pver, begchunk, endchunk                       ! physics grid
    use phys_grid,      only: get_ncols_p, nlcols, ngcols, &                        ! physics parallel decomposition
                              scatter_field_to_chunk, gather_chunk_to_field, &
                              buff_to_chunk, chunk_to_buff, get_chunk_owner_p, &
                              read_chunk_from_field, write_field_from_chunk
    use shr_sys_mod,    only: shr_sys_flush, shr_sys_irtc                           ! standardized system subroutines
    use shr_kind_mod,   only: shr_kind_in                                           ! defines CCSM real & integer kinds
    use shr_const_mod,  only: shr_const_spval
    use cpl_contract_mod
    use cpl_interface_mod
    use cpl_fields_mod
!--------------------------------------------------------------------------
! NOTE:  if nlon is not the same as nlon_p in phys_grid, this module will
!            need to be modified  -RLJ
!--------------------------------------------------------------------------
    use rgrid,          only: nlon                                                  ! reduced grid
!
#if ( defined SPMD )
    use mpishorthand,   only: mpicom, mpir8, mpiint, mpilog                         ! MPI interface
#endif
    use history,        only: outfld                                                ! history output

    implicit none

    public ccsmini                  ! Initialization
    public ccsmsnd                  ! Send information to coupler
    public ccsmrcv                  ! Receive information from coupler
    public ccsmfin                  ! Finalization, shut down model
    public ccsmave                  ! Average CCSM data (when flxave set)
    public read_restart_ccsm        ! Read the CCSM restart information
    public write_restart_ccsm       ! Write the CCSM restart information
    public initialize_ccsm_msg      ! Initialize ccsm_msg data
    !
    ! When to send/receive messages to coupler and when to make restart and stop
    !
    logical, public :: dorecv       ! receive data from coupler this step
    logical, public :: dosend       ! send data to coupler this step
    logical, public :: csmstop      ! stop signal received from coupler
    logical, public :: csmrstrt     ! restart write signal received from coupler
    !
    ! Surface data important for CCSM only
    !
    real(r8), public, allocatable :: rho(:,:)      ! surface air density
    real(r8), public, allocatable :: netsw(:,:)    ! net shortwave
    real(r8), public, allocatable :: psl(:,:)      ! sea-level pressure

    private
    !
    ! Buffer information
    !
    integer(shr_kind_in), private :: ibuff(cpl_fields_ibuf_total) ! integer buffer to/from cpl
    real(r8), private             :: rbuff(cpl_fields_rbuf_total) ! floating point buffer to/from cpl

    real(r8) :: spval = shr_const_spval          ! Special value for real msg data
    !
    ! Timing information
    !
    logical, private :: csm_timing            ! turn timing of CCSM messages on
    integer, private :: irtc_w                ! rtc ticks when waiting for msg
    integer, private :: irtc_r                ! rtc ticks when msg recved
    integer, private :: irtc_s                ! rtc ticks when msg sent
    !
    ! send/recv buffers
    !
    integer(shr_kind_in), private, parameter :: nsnd=cpl_fields_a2c_total   ! number of send variables
    integer(shr_kind_in), private, parameter :: nrcv=cpl_fields_c2a_total   ! number of recv variables
    real(r8), allocatable, private :: rbuffs(:,:)                           ! real buffer for send array
    real(r8), allocatable, private :: rbuffr(:,:)                           ! real buffer for recv array

    real(r8), private :: arget(plon,plat,nrcv)          ! recv array on 1 processor
    real(r8), private :: recv2d(plon,nrcv,plat)         ! permutation of recv array on 1 proc

    real(r8), allocatable :: recv2d_chunk(:,:,:)        ! chunked recv array
    real(r8), allocatable :: send2d_chunk(:,:,:)        ! chunked send array
    type(cpl_contract), save :: contractS
    type(cpl_contract), save :: contractR
    !
    ! flux accumulator
    !
    integer, private :: countfa      ! counter for flux accumulators
    !
    ! Surface data that needs to be averaged
    !
    real(r8), allocatable:: precca(:,:)   ! average convective precipitation
    real(r8), allocatable:: precla(:,:)   ! average large-scale precipation
    real(r8), allocatable:: precsca(:,:)  ! average convective snow-fall
    real(r8), allocatable:: precsla(:,:)  ! average large-scale snow-fall
    real(r8), allocatable:: rainconv(:,:) ! convective rainfall
    real(r8), allocatable:: rainlrsc(:,:) ! large-scale rainfall
    real(r8), allocatable:: snowconv(:,:) ! convective snowfall
    real(r8), allocatable:: snowlrsc(:,:) ! larse-scale snowfall
    real(r8), allocatable:: prc_err(:,:)  ! error in precipitation sent to coupler

contains

    !-----------------------------------------------------------------------
    !
    ! Purpose: Initialize ccsm coupler communications
    !
    ! Method:
    !
    ! Author: Mariana Vertenstein
    !
    !-----------------------------------------------------------------------

    subroutine ccsmini

        use comsrf,         only: srfflx_state2d
        use physconst,      only: stebol
        use constituents,   only: pcnst, pnats
        use time_manager,   only: is_first_step

#include <comctl.h>

        integer i, m, lchnk, n  ! indices
        integer ncols           ! number of columns
        integer ierr            ! allocation error signal
        integer sizebuf         ! size of buffer for sending grid data to coupler

        call t_startf('ccsm_initialization')
        !
        ! set the CCSM stop and restart to false
        !
        csmstop  = .false.
        csmrstrt = .false.
        !
        ! allocate chunked send/receive buffers
        !
        if (.not. allocated(recv2d_chunk))then
            allocate(recv2d_chunk(pcols,nrcv,begchunk:endchunk), stat=ierr)
            if (ierr /= 0) then
                call endrun('ccsmini: recv2d_chunk allocation error')
            end if
        end if
        allocate(send2d_chunk(pcols,nsnd,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            call endrun('ccsmini: send2d_chunk allocation error')
        end if
        !
        ! allocate long/lat send/receive buffers
        !
        if (.not. allocated(rbuffr)) then
            allocate(rbuffr(nlcols,nrcv), stat=ierr)
            if (ierr /= 0) then
                call endrun('ccsmini: rbuffr allocation error')
            end if
        end if
        allocate(rbuffs(nlcols,nsnd), stat=ierr)
        if (ierr /= 0) then
            call endrun ('ccsmini: rbuffs allocation error')
        end if
        !
        ! for now set all tracer fluxes to zero
        !
        do m = 2, pcnst+pnats
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i=1,ncols
                    srfflx_state2d(lchnk)%cflx(i,m) = 0.
                end do
            end do
        end do
        !
        ! require the short and longwave radiation frequencies to match, since these
        ! fluxes will be sent as instantaneous fluxes to the coupler, valid over the
        ! next interval.
        !
        if (masterproc) then
            if (flxave) then
                if (iradsw == iradlw) then
                    write(6, *) 'ccsmini: coupling will take place every ',iradsw, ' steps'
                else
                    write(6, *) 'ccsmini: iradsw != iradlw ', iradsw, iradlw
                    call endrun("ccsmini: bad irad")
                end if
            else
                write(6, *) 'ccsmini: coupling will take place every time step'
            end if
            call shr_sys_flush(6)
        end if
        !
        ! send grid and initial ibuf to flux coupler
        !
        if (masterproc) then
            write(6, *) 'ccsmini: send grid to coupler'
            call shr_sys_flush(6)
        end if
        call ccsm_msg_sendgrid
        !
        ! receive orbital parameters
        !
        if (masterproc) then
            write(6, *) 'ccsmini: get orbital parameters from coupler'
            call shr_sys_flush(6)
        end if
        call ccsm_msg_getorb
        !
        ! For initial run only:
        !
        if (is_first_step()) then
            !
            ! initial run only: get albedos and ice fraction
            !
            if (masterproc) then
                write(6, *) 'ccsmini: initial run, get albedos from coupler'
                call shr_sys_flush(6)
            end if
            call ccsm_msg_getalb
            !
            ! Initial run only: determine longwave up flux from the surface temperature.
            !
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    srfflx_state2d(lchnk)%lwup(i) = stebol*(srfflx_state2d(lchnk)%ts(i)**4)
                end do
            end do
        end if

        if (masterproc) then
            write(6, *) 'ccsmini: CCSM initialization complete!'
            call shr_sys_flush(6)
        end if

        call t_stopf('ccsm_initialization')
        call t_startf('ccsm_rcvtosnd')
        call t_startf('ccsm_runtotal')

        return
    end subroutine ccsmini

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Get the message array from the csm driver and extract the data
    !
    ! Method:
    !
    ! Author: Byron Boville
    !        modified for cpl6 by Robert Jacob
    !
    !-----------------------------------------------------------------------

    subroutine ccsmrcv

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                          snowhice, snowhland, verify_fractions

#include <comctl.h>

        integer i, lat, n, lchnk    ! indices
        integer ncols               ! Number of columns
        integer len                 ! temporary variable length

        real(r8) newlandfrac        ! land fraction computed as residual of icefrac + ocnfrac
        real(r8) delta              ! land fraction difference across timesteps
                                    ! code needs to be rewritten to guarantee that this is zero

        call t_startf('ccsm_rcv')
        call t_stopf('ccsm_sndtorcv')
        !
        ! get data from flux coupler.
        !
        if (csm_timing) irtc_w = shr_sys_irtc()
        call cpl_interface_contractRecv(cpl_fields_cplname, contractR, ibuff, rbuffr)
        if (masterproc) then
            if (csm_timing) then
                irtc_r = shr_sys_irtc()
                write(6, "('[mp timing]  irtc = ',i20,' ',a)") irtc_w, 'd->a waiting'
                write(6, "('[mp timing]  irtc = ',i20,' ',a)") irtc_r, 'd->a received'
            end if
        end if
        !
        ! rearrange data from lon-lat buffer into chunk structure
        !
        call buff_to_chunk(nrcv, plon, rbuffr, recv2d_chunk)
        !
        ! split buffer into component arrays. Change signs as required.
        ! Note that coupler has convention that fluxes are positive downward.
        !

        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do i = 1, ncols
                srfflx_state2d(lchnk)%wsx(i)    = -recv2d_chunk(i,cpl_fields_c2a_taux ,lchnk)  ! wind stress, zonal
                srfflx_state2d(lchnk)%wsy(i)    = -recv2d_chunk(i,cpl_fields_c2a_tauy ,lchnk)  ! wind stress, meridional
                srfflx_state2d(lchnk)%lhf(i)    = -recv2d_chunk(i,cpl_fields_c2a_lat  ,lchnk)  ! latent          heat flux
                srfflx_state2d(lchnk)%shf(i)    = -recv2d_chunk(i,cpl_fields_c2a_sen  ,lchnk)  ! sensible        heat flux
                srfflx_state2d(lchnk)%lwup(i)   = -recv2d_chunk(i,cpl_fields_c2a_lwup ,lchnk)  ! upward longwave heat flux
                srfflx_state2d(lchnk)%cflx(i,1) = -recv2d_chunk(i,cpl_fields_c2a_evap ,lchnk)  !  evaporation    water flux
                srfflx_state2d(lchnk)%asdir(i)  =  recv2d_chunk(i,cpl_fields_c2a_avsdr,lchnk)  ! albedo, visible, direct
                srfflx_state2d(lchnk)%aldir(i)  =  recv2d_chunk(i,cpl_fields_c2a_anidr,lchnk)  ! albedo, near-ir, direct
                srfflx_state2d(lchnk)%asdif(i)  =  recv2d_chunk(i,cpl_fields_c2a_avsdf,lchnk)  ! albedo, visible, diffuse
                srfflx_state2d(lchnk)%aldif(i)  =  recv2d_chunk(i,cpl_fields_c2a_anidf,lchnk)  ! albedo, near-ir, diffuse
                srfflx_state2d(lchnk)%ts(i)     =  recv2d_chunk(i,cpl_fields_c2a_t    ,lchnk)  ! surface temperature
                srfflx_state2d(lchnk)%sst(i)    =  recv2d_chunk(i,cpl_fields_c2a_sst  ,lchnk)  ! sea surface temperature
                snowhland(i,lchnk)              =  recv2d_chunk(i,cpl_fields_c2a_snowh,lchnk)  ! surface snow depth
                icefrac(i,lchnk)                =  recv2d_chunk(i,cpl_fields_c2a_ifrac,lchnk)  ! surface ice fraction
                ocnfrac(i,lchnk)                =  recv2d_chunk(i,cpl_fields_c2a_ofrac,lchnk)  ! surface ocn fraction
                !
                ! Bit-for-bit stuff.  Much of the code between "bit-for-bit" comments was originally added
                ! to fix an inexact restart problem in CSM.  Checking on change in landfrac between timesteps
                ! is new.
                !
                ! Verify that landfrac is within roundoff of its previous value.  A better solution is to
                ! communicate landfrac only on startup, and never recompute it.  Also, icefrac and ocnfrac
                ! should NEVER be changed by CAM in coupled mode.  Once this is done, code between the
                ! "bit-for-bit" comments below and similarly in ccsm_msg_getalb (below) can be excised.
                !
                newlandfrac = 1.0-icefrac(i,lchnk)-ocnfrac(i,lchnk)
                delta = newlandfrac-landfrac(i,lchnk)
                if (abs(delta) > 10.*epsilon(1.0_r8)) then
                    write(6, *) 'ccsmrcv: new landfrac differs beyond roundoff from its previous value'
                    write(6, *) 'i, lchnk, oldlandfrac, newlandfrac =', i, lchnk, landfrac(i,lchnk), newlandfrac
                    call endrun
                end if

                landfrac(i,lchnk) = newlandfrac
                if (icefrac(i,lchnk)+landfrac(i,lchnk) > 1.0) then
                    icefrac(i,lchnk) = 1.0-landfrac(i,lchnk)
                end if

                ocnfrac(i,lchnk) = 1.0-landfrac(i,lchnk)-icefrac(i,lchnk)
                !
                ! End of bit-for-bit stuff.
                !
                srfflx_state2d(lchnk)%tref(i) = recv2d_chunk(i,cpl_fields_c2a_tref,lchnk)  ! 2m reference temperature
                srfflx_state2d(lchnk)%qref(i) = recv2d_chunk(i,cpl_fields_c2a_qref,lchnk)  ! 2m reference specific humidity
            end do
            !
            ! Ensure that fractions are valid
            !
            call verify_fractions(lchnk, ncols)
        end do
        !
        ! Set snowh over ice to zero since flux coupler only returns snowh over land
        !
        snowhice(:,:) = 0.0
        !
        ! Determine if stop at end of day
        !
        if (csmstop == .false. .and. ibuff(cpl_fields_ibuf_stopeod) /= 0) then
            csmstop = .true.
            if (masterproc) write(6, *) 'ccsmrcv: received stop at end of day signal from flux coupler'
        end if
        !
        ! Determine if write restart at end of day
        !
        if (csmrstrt == .false. .and. ibuff(cpl_fields_ibuf_resteod) /= 0) then
            csmrstrt = .true.
            if (masterproc) write(6, *) 'ccsmrcv: received write restart at end of day signal from flux coupler'
        else if (ibuff(cpl_fields_ibuf_resteod) == 0)then
            csmrstrt = .false.
        end if

        call t_stopf('ccsm_rcv')
        call t_startf('ccsm_rcvtosnd')

        return
    end subroutine ccsmrcv

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Send the message array to the csm driver.
    !
    ! Method:
    ! On steps where the data is to be sent to the coupler, fill the
    ! message passing array with instantaneous atmospheric states,
    ! instantaneous downward radiative fluxes, averaged precipitation,
    ! instantaneous surface states and averaged surface fluxes.
    ! Condense the data into one array. The coupler has the convention that
    ! fluxes are positive downward. Note that precc and precl precipitation rates
    ! in units of m/sec. They are turned into fluxes by multiplying by 1000 kg/m^3.
    !
    ! Author: Byron Boville
    !         modified for cpl6 by Robert Jacob
    !
    !-----------------------------------------------------------------------

    subroutine ccsmsnd

        use comsrf, only: surface_state2d, srfflx_state2d

#include <comctl.h>

        integer i, lchnk, n, lat    ! indices
        integer ncols               ! Number of columns
        integer len                 ! temporary length variable

        call t_startf('ccsm_snd')
        call t_stopf('ccsm_rcvtosnd')
        !
        ! Divide total precipitation and snowfall into rain and snowfall
        !
        if (flxave) then
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    rainconv(i,lchnk) = ((precca(i,lchnk)-precsca(i,lchnk)))*1000.
                    rainlrsc(i,lchnk) = ((precla(i,lchnk)-precsla(i,lchnk)))*1000.
                    snowconv(i,lchnk) = precsca(i,lchnk)*1000.
                    snowlrsc(i,lchnk) = precsla(i,lchnk)*1000.
                end do
            end do
        else
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    rainconv(i,lchnk) = ((surface_state2d(lchnk)%precc(i)-surface_state2d(lchnk)%precsc(i)))*1000.
                    rainlrsc(i,lchnk) = ((surface_state2d(lchnk)%precl(i)-surface_state2d(lchnk)%precsl(i)))*1000.
                    snowconv(i,lchnk) = surface_state2d(lchnk)%precsc(i)*1000.
                    snowlrsc(i,lchnk) = surface_state2d(lchnk)%precsl(i)*1000.
                end do
            end do
        end if
        !
        ! If averaging flux over several timesteps, ensure rain and snow do not
        ! exist simultaneously to satisfy a limitation in LSM.
        ! This code removed 3/2003, TC
        ! CLM2.2 now accepts both rain and snow on same grid point, prc_err == 0.
        !
        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            prc_err(1:ncols,lchnk)  = 0.
        end do
        !
        ! Copy from component arrays into one chunk array.
        ! Note that coupler has convention that fluxes are positive downward.
        !
        do lchnk = begchunk, endchunk
            ncols = get_ncols_p(lchnk)
            do i = 1, ncols
                send2d_chunk(i,cpl_fields_a2c_z,lchnk)     = surface_state2d(lchnk)%zbot(i) ! Atmospheric state variable m
                send2d_chunk(i,cpl_fields_a2c_u,lchnk)     = surface_state2d(lchnk)%ubot(i) ! Atmospheric state variable m/s
                send2d_chunk(i,cpl_fields_a2c_v,lchnk)     = surface_state2d(lchnk)%vbot(i) ! Atmospheric state variable m/s
                send2d_chunk(i,cpl_fields_a2c_tbot,lchnk)  = surface_state2d(lchnk)%tbot(i) ! Atmospheric state variable K
                send2d_chunk(i,cpl_fields_a2c_ptem,lchnk)  = surface_state2d(lchnk)%thbot(i)! Atmospheric state variable K
                send2d_chunk(i,cpl_fields_a2c_pbot,lchnk)  = surface_state2d(lchnk)%pbot(i) ! Atmospheric state variable Pa
                send2d_chunk(i,cpl_fields_a2c_pslv,lchnk)  = psl(i,lchnk)                   ! Atmospheric state variable Pa
                send2d_chunk(i,cpl_fields_a2c_shum,lchnk)  = surface_state2d(lchnk)%qbot(i) ! Atmospheric state variable kg/kg
                send2d_chunk(i,cpl_fields_a2c_dens,lchnk)  = rho(i,lchnk)                   ! Atmospheric state variable kg/m^3
                send2d_chunk(i,cpl_fields_a2c_swnet,lchnk) = netsw(i,lchnk)                 ! Atmospheric flux W/m^2
                send2d_chunk(i,cpl_fields_a2c_lwdn,lchnk)  = surface_state2d(lchnk)%flwds(i)! Atmospheric flux W/m^2
                send2d_chunk(i,cpl_fields_a2c_rainc,lchnk) = rainconv(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,cpl_fields_a2c_rainl,lchnk) = rainlrsc(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,cpl_fields_a2c_snowc,lchnk) = snowconv(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,cpl_fields_a2c_snowl,lchnk) = snowlrsc(i,lchnk)              ! Atmospheric flux kg/s/m^2
                send2d_chunk(i,cpl_fields_a2c_swndr,lchnk) = surface_state2d(lchnk)%soll(i) ! Atmospheric flux W/m^2
                send2d_chunk(i,cpl_fields_a2c_swvdr,lchnk) = surface_state2d(lchnk)%sols(i) ! Atmospheric flux W/m^2
                send2d_chunk(i,cpl_fields_a2c_swndf,lchnk) = surface_state2d(lchnk)%solld(i)! Atmospheric flux W/m^2
                send2d_chunk(i,cpl_fields_a2c_swvdf,lchnk) = surface_state2d(lchnk)%solsd(i)! Atmospheric flux W/m^2
            end do
        end do
        !
        ! Rearrange data from chunk structure into couplers lon-lat buffer
        !
        call chunk_to_buff(nsnd, plon, send2d_chunk, rbuffs)
        !
        ! Output to history file the snow and rain actually sent to coupler as well as the
        ! error between what is sent and what is reported on history file in PRECT/PRECS
        !
        do lchnk = begchunk, endchunk
            call outfld('CPLRAINC', rainconv(1,lchnk), pcols, lchnk)
            call outfld('CPLRAINL', rainlrsc(1,lchnk), pcols, lchnk)
            call outfld('CPLSNOWC', snowconv(1,lchnk), pcols, lchnk)
            call outfld('CPLSNOWL', snowlrsc(1,lchnk), pcols, lchnk)
            call outfld('CPLPRCER', prc_err (1,lchnk), pcols, lchnk)
        end do
        !
        ! Send buffer to coupler
        !
        call msgsnd

        call t_stopf('ccsm_snd')
        call t_startf('ccsm_sndtorcv')

        return
    end subroutine ccsmsnd

    !-----------------------------------------------------------------------
    !
    ! Purpose: Send and receive final msgs at end of run.
    !
    ! Method:
    ! The coupler currently expects a final when nlend is true -
    ! this data is only written out the coupler restart file
    ! restart file and is not used upon restart by the coupler
    ! for the cam component. The coupler also sends a final msg.
    ! This data is put into a dummy array
    !
    ! Author: Mariana Vertenstein
    !         modified for cpl6 by Robert Jacob
    !
    !-----------------------------------------------------------------------

    subroutine ccsmfin

        use time_manager, only: get_nstep, get_prev_date

#include <comctl.h>

        integer nstepcsm            ! time steps sent to flux coupler
        integer cdatecsm, cseccsm   ! date, sec at beginning of current time step
        integer yr, mon, day        ! year, month, day components of cdatecsm

        !
        ! Determine final date to send to coupler
        !
        nstepcsm = get_nstep()-1
        call get_prev_date(yr, mon, day, cseccsm)
        cdatecsm = yr*10000+mon*100+day
        !
        ! Initialize ibuff
        !
        ibuff(:)  = 0
        ibuff(cpl_fields_ibuf_cdate) = cdatecsm   ! model date (yyyymmdd)
        ibuff(cpl_fields_ibuf_sec)   = cseccsm    ! elapsed seconds in current day
        !ibuff(6)  = nstepcsm   ! ending model time step
        !
        ! Send final message
        !
        rbuffs(:,:) = spval
        call cpl_interface_contractSend(cpl_fields_cplname, contractS, ibuff, rbuffs)
        !
        ! Receive final message
        !
        call cpl_interface_contractRecv(cpl_fields_cplname, contractR, ibuff, rbuffr)

        call t_stopf('ccsm_runtotal')
        call t_stopf('ccsm_rcvtosnd')

        return
    end subroutine ccsmfin

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Write COUP_CSM specific variables to restart dataset
    !
    ! Author:
    !         modified for cpl6 by Robert Jacob
    !
    !-----------------------------------------------------------------------

    subroutine write_restart_ccsm

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, snowhland

#include <comlun.h>
#include <comctl.h>

        integer ioerr
        integer n, lat, i
        integer ncols, lchnk
        !
        ! Write out flux averaging flag
        !
        if (masterproc) then
            write(nrg, iostat=ioerr) flxave
            if (ioerr /= 0) then
                write(6, "('Error: ccsm_msg: write_restart_ccsm: ')", advance="no")
                write(6, "('I/O error ', I5, ' on unit ', I5)") ioerr, nrg
                call endrun
            end if
        end if
        !
        ! If flux averaging is enabled write out necessary info
        !
        if (flxave) then
            ! Update recv2d_chunk to the current state, not the state last
            ! received.  CAM very infrequently changes some of the coupling fields
            ! but when it does, exact restart testing in CCSM fails.
            ! Added June, 2003 (TC,MV)

            !new stuff
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    recv2d_chunk(i,cpl_fields_c2a_taux,lchnk)  =  -srfflx_state2d(lchnk)%wsx(i)
                    recv2d_chunk(i,cpl_fields_c2a_tauy,lchnk)  =  -srfflx_state2d(lchnk)%wsy(i)
                    recv2d_chunk(i,cpl_fields_c2a_lat,lchnk)   =  -srfflx_state2d(lchnk)%lhf(i)
                    recv2d_chunk(i,cpl_fields_c2a_sen,lchnk)   =  -srfflx_state2d(lchnk)%shf(i)
                    recv2d_chunk(i,cpl_fields_c2a_lwup,lchnk)  =  -srfflx_state2d(lchnk)%lwup(i)
                    recv2d_chunk(i,cpl_fields_c2a_evap,lchnk)  =  -srfflx_state2d(lchnk)%cflx(i,1)
                    recv2d_chunk(i,cpl_fields_c2a_avsdr,lchnk) =  srfflx_state2d(lchnk)%asdir(i)
                    recv2d_chunk(i,cpl_fields_c2a_anidr,lchnk) =  srfflx_state2d(lchnk)%aldir(i)
                    recv2d_chunk(i,cpl_fields_c2a_avsdf,lchnk) =  srfflx_state2d(lchnk)%asdif(i)
                    recv2d_chunk(i,cpl_fields_c2a_anidf,lchnk) =  srfflx_state2d(lchnk)%aldif(i)
                    recv2d_chunk(i,cpl_fields_c2a_t,lchnk)     =  srfflx_state2d(lchnk)%ts(i)
                    recv2d_chunk(i,cpl_fields_c2a_sst,lchnk)   =  srfflx_state2d(lchnk)%sst(i)
                    recv2d_chunk(i,cpl_fields_c2a_snowh,lchnk) =  snowhland(i,lchnk)
                    recv2d_chunk(i,cpl_fields_c2a_ifrac,lchnk) =  icefrac(i,lchnk)
                    recv2d_chunk(i,cpl_fields_c2a_ofrac,lchnk) =  ocnfrac(i,lchnk)
                    recv2d_chunk(i,cpl_fields_c2a_tref,lchnk)  =  srfflx_state2d(lchnk)%tref(i)
                    recv2d_chunk(i,cpl_fields_c2a_qref,lchnk)  =  srfflx_state2d(lchnk)%qref(i)
                end do
            end do
            call gather_chunk_to_field(1, nrcv, 1, plon, recv2d_chunk, recv2d)
            if (masterproc) then
                do n = 1, nrcv
                    do lat = 1, plat
                        do i = 1, plon
                            arget(i,lat,n) = recv2d(i,n,lat)
                        end do
                    end do
                end do

                write(nrg, iostat=ioerr) dosend, countfa, arget
                if (ioerr /= 0 ) then
                    write(6, "('Error: ccsm_msg: write_restart_ccsm: ')", advance="no")
                    write(6, "('I/O error ', I5, ' on unit ', I5)") ioerr, nrg
                    call endrun
                end if
            end if
            call write_field_from_chunk(nrg, 1, 1, 1, precca)
            call write_field_from_chunk(nrg, 1, 1, 1, precla)
            call write_field_from_chunk(nrg, 1, 1, 1, precsca)
            call write_field_from_chunk(nrg, 1, 1, 1, precsla)
        end if

        return
    end subroutine write_restart_ccsm

    !-----------------------------------------------------------------------
    ! Read in COUP_CSM specific variables and determine surface state
    ! variables and fluxes from arget. NOTE: are not assured that will
    ! do a recv upon restart if flux averaging is used so must split
    ! these variables off from arget since most are not individually
    ! written out to the restart file.
    !-----------------------------------------------------------------------

    subroutine read_restart_ccsm

        use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                          snowhland, snowhice

#include <comlun.h>
#include <comctl.h>

        integer i, n, lat, lchnk     ! indices
        integer ncols                ! Number of columns
        integer ioerr                ! error return
        integer len                  ! length for spmd
        logical flxave_res           ! flux averaging flag from restart file
        integer ierr                 ! error flag
        !
        ! Read in flux averaging flag
        !
        if (masterproc) then
            read(nrg, iostat=ioerr) flxave_res
            if (ioerr /= 0 ) then
                write(6, "('Error: ccsm_msg: read_restart_ccsm: ')", advance="no")
                write(6, "('I/O error ', I5, ' on unit ', I5)") ioerr, nrg
                call endrun
            end if
            if ((flxave_res .and. .not. flxave) .or. &
                (.not. flxave_res .and. flxave)) then
                write(6, "('Error: ccsm_msg: read_restart_ccsm: ')", advance="no")
                write(6, "('namelist flxave (', L1, ') must equal restart flxave (', L1, ')')") flxave, flxave_res
                call endrun
            end if
        end if
        !
        ! If flux averaging is enabled read in necessary info
        !
        if (flxave) then
            if (masterproc) then
                read(nrg, iostat=ioerr) dosend, countfa, arget
                if (ioerr /= 0 ) then
                    write(6, "('Error: ccsm_msg: read_restart_ccsm: ')", advance="no")
                    write(6, "('I/O error ', I5, ' on unit ', I5)") ioerr, nrg
                    call endrun
                end if
            end if
            call read_chunk_from_field(nrg, 1, 1, 1, precca)
            call read_chunk_from_field(nrg, 1, 1, 1, precla)
            call read_chunk_from_field(nrg, 1, 1, 1, precsca)
            call read_chunk_from_field(nrg, 1, 1, 1, precsla)
#ifdef SPMD
            call mpibcast(dosend , 1, mpilog, 0, mpicom)
            call mpibcast(countfa, 1, mpiint, 0, mpicom)
#endif
            if (masterproc) then
                do n = 1, nrcv
                    do lat = 1, plat
                        do i = 1, plon
                            recv2d(i,n,lat) = arget(i,lat,n)
                        end do
                    end do
                end do
            end if
            if (.not. allocated(recv2d_chunk)) then
                allocate(recv2d_chunk(pcols,nrcv,begchunk:endchunk), stat=ierr)
                if (ierr /= 0) then
                    write(6, "('Error: ccsm_msg: read_restart_ccsm: ')", advance="no")
                    write(6, "('I/O error ', I5, ' on unit ', I5)") ioerr, nrg
                    call endrun
                end if
            end if

            call scatter_field_to_chunk(1, nrcv, 1, plon, recv2d, recv2d_chunk)

            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                do i = 1, ncols
                    srfflx_state2d(lchnk)%wsx(i)    = -recv2d_chunk(i, cpl_fields_c2a_taux,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%wsy(i)    = -recv2d_chunk(i, cpl_fields_c2a_tauy,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%lhf(i)    = -recv2d_chunk(i, cpl_fields_c2a_lat ,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%shf(i)    = -recv2d_chunk(i, cpl_fields_c2a_sen ,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%lwup(i)   = -recv2d_chunk(i, cpl_fields_c2a_lwup,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%cflx(i,1) = -recv2d_chunk(i, cpl_fields_c2a_evap,lchnk)   ! Atmosphere-surface flux
                    srfflx_state2d(lchnk)%asdir(i)  =  recv2d_chunk(i, cpl_fields_c2a_avsdr,lchnk)  ! Surface state variable
                    srfflx_state2d(lchnk)%aldir(i)  =  recv2d_chunk(i, cpl_fields_c2a_anidr,lchnk)  ! Surface state variable
                    srfflx_state2d(lchnk)%asdif(i)  =  recv2d_chunk(i, cpl_fields_c2a_avsdf,lchnk)  ! Surface state variable
                    srfflx_state2d(lchnk)%aldif(i)  =  recv2d_chunk(i, cpl_fields_c2a_anidf,lchnk)  ! Surface state variable
                    srfflx_state2d(lchnk)%ts(i)     =  recv2d_chunk(i,cpl_fields_c2a_t,lchnk)       ! Surface state variable
                    srfflx_state2d(lchnk)%sst(i)    =  recv2d_chunk(i,cpl_fields_c2a_sst,lchnk)     ! Surface state variable
                    snowhland(i,lchnk)              =  recv2d_chunk(i,cpl_fields_c2a_snowh,lchnk)   ! Surface state variable
                    icefrac(i,lchnk)                =  recv2d_chunk(i,cpl_fields_c2a_ifrac,lchnk)   ! Surface type fraction
                    ocnfrac(i,lchnk)                =  recv2d_chunk(i,cpl_fields_c2a_ofrac,lchnk)   ! Surface type fraction
                    srfflx_state2d(lchnk)%tref(i)   =  recv2d_chunk(i,cpl_fields_c2a_tref,lchnk)    ! Surface state variable
                    !
                    ! Do not update ocnfrac or icefrac (e.g. by calling update_ocnice), as these are read from
                    ! the restart file.  Physics buffer at restart already had the correct fractions.  Don't
                    ! even update landfrac!
                    !
                    !landfrac(i,lchnk) = 1.0_r8 - icefrac(i,lchnk) - ocnfrac(i,lchnk)
                end do
            end do
        end if

        snowhice(:,:) = 0.0

        return
    end subroutine read_restart_ccsm

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Initialize data so that if data used before set the program will die.
    !
    ! Method:
    !
    ! Author: Mariana Vertenstein
    !
    !-----------------------------------------------------------------------

    subroutine initialize_ccsm_msg

        use infnan

        integer ierr

        allocate(rho(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('rho allocation error')")
            call endrun
        end if
        allocate(netsw(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('netsw allocation error')")
            call endrun
        end if
        allocate(psl(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('psl allocation error')")
            call endrun
        end if
        allocate(precca(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('precca allocation error')")
            call endrun
        end if
        allocate(precla(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('precla allocation error')")
            call endrun
        end if
        allocate(precsca(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('precsca allocation error')")
            call endrun
        end if
        allocate(precsla(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('precsla allocation error')")
            call endrun
        end if
        allocate(rainconv(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('rainconv allocation error')")
            call endrun
        end if
        allocate(rainlrsc(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('rainlrsc allocation error')")
            call endrun
        end if
        allocate(snowconv(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('snowconv allocation error')")
            call endrun
        end if
        allocate(snowlrsc(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('snowlrsc allocation error')")
            call endrun
        end if
        allocate(prc_err(pcols,begchunk:endchunk), stat=ierr)
        if (ierr /= 0) then
            write(6, "('Error: ccsm_msg: initialize_ccsm_msg: ')", advance="no")
            write(6, "('prc_err allocation error')")
            call endrun
        end if
        !
        ! Initialize to NaN or Inf
        !
        rho     (:,:) = inf
        netsw   (:,:) = inf
        psl     (:,:) = inf
        precca  (:,:) = inf
        precla  (:,:) = inf
        precsca (:,:) = inf
        precsla (:,:) = inf
        snowconv(:,:) = inf
        snowlrsc(:,:) = inf
        rainconv(:,:) = inf
        rainlrsc(:,:) = inf
        prc_err (:,:) = inf

    end subroutine initialize_ccsm_msg

!===============================================================================
! The following subroutines private to this module!
!===============================================================================

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Send message to flux coupler
    !
    ! Method:
    !
    ! Author: Mariana Vertenstein
    !         modified for cpl6 by Robert Jacob
    !
    !-----------------------------------------------------------------------

    subroutine msgsnd

        use time_manager, only: get_nstep, get_step_size, get_curr_date, &
                                get_prev_date

#include <comctl.h>

        integer n                      ! count indices
        integer nstep                  ! current time step
        integer nstepcsm               ! time step sent to flux coupler
        integer yr, mon, day           ! year, month, day components of cdatecsm
        integer cdatecsm, cseccsm      ! current date,sec
        integer msgpday                ! number of send/recv msgs per day
        logical nextsw                 ! set to true for next sw calculation
        real(r8) dtime                 ! timestep size
        real(r8) albshift              ! albedo calculation time shift

        nstep = get_nstep()
        dtime = get_step_size()
        !
        ! Determine time step sent to flux coupler and corresponding date.
        !
        if (nstep==0) then
            nstepcsm = nstep
            call get_curr_date(yr, mon, day, cseccsm)
            cdatecsm = yr*10000+mon*100+day
        else
            nstepcsm = nstep - 1
            call get_prev_date(yr, mon, day, cseccsm)
            cdatecsm = yr*10000+mon*100+day
        end if
        !
        ! Determine albedo calculation time shift, which is the time interval
        ! from nstepcsm until the next short wave calculation.
        !
        if (nstep /= 0) then
            if (flxave) then
                albshift = nint((nstep+iradsw-nstepcsm)*dtime)
            else
                nextsw = .false.
                n = 1
                do while (.not. nextsw)
                    nextsw = (mod((nstep+n-1),iradsw)==0)
                    if (nextsw) albshift = nint((nstep+n-nstepcsm)*dtime)
                    n = n+1
                end do
            end if
        else
            albshift = nint(iradsw*dtime)+dtime
        end if
        !
        ! Determine number of send/recv msgs per day
        !
        if (flxave) then
            msgpday = nint(86400./dtime)/iradsw
        else
            msgpday = nint(86400./dtime)
        end if
        !
        ! Determine ibuff array
        !
        ibuff(:)  = 0
        ibuff(cpl_fields_ibuf_cdate) = cdatecsm             ! model date (yyyymmdd)
        ibuff(cpl_fields_ibuf_sec)   = cseccsm              ! elapsed seconds in current day
        !ibuff(6)  = nstepcsm                               ! model time step
        ibuff(cpl_fields_ibuf_gisize)  = plon               ! number of model longitudes
        ibuff(cpl_fields_ibuf_gjsize)  = plat               ! number of model latitudes
        ibuff(cpl_fields_ibuf_ncpl)    = msgpday            ! number of send/recv msgs per day
        ibuff(cpl_fields_ibuf_ashift ) = albshift           ! albedo calculation time shift
        ibuff(cpl_fields_ibuf_userest) = 1                  ! if ccsm_msg_getalb was not called, this may be used
                                                            ! to send the first message.  So tell coupler
                                                            ! again to not use its restart data.
        !
        ! Send data to coupler. rbuffs was set above in ccsmsnd
        !
        call cpl_interface_contractSend(cpl_fields_cplname, contractS, ibuff, rbuffs)

        if(masterproc) then
            if (csm_timing) then
                irtc_s = shr_sys_irtc()
                write(6, "('[mp timing]  irtc = ',i20,' ',a)") irtc_s, 'a->d sending'
            end if
            end if

        return
    end subroutine msgsnd

    subroutine ccsmave(iradsw, nstep, dosw)

!-----------------------------------------------------------------------
!
! Purpose:
! Average the input fluxes to lsm between solar radiation times.
!
! Method:
! Currently, the only flux requiring averaging is the precipitation,
! since the radiative fluxes are constant over the averaging interval.
!
! Author: Byron Boville
!
!-----------------------------------------------------------------------

    use comsrf, only: surface_state2d

!------------------------------Arguments--------------------------------
    integer, intent(in) :: iradsw  ! solar radiation interval
    integer, intent(in) ::  nstep  ! time step number
    logical, intent(in) ::  dosw   ! time to compute averages (solar radiation time)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
    integer i,lchnk  ! indices
    integer ncols    ! Number of columns
    real(r8) rcount  ! reciprocal of count
!-----------------------------------------------------------------------
!
! If iradsw == 1, then no averaging is required
!
    if (iradsw == 1) return
!
! Set the counter and normalizing factor
!
    if (nstep == 0) countfa = 0
    countfa = countfa + 1
    if (dosw) then
       rcount = 1./countfa
    end if

    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       if (countfa == 1) then
          do i = 1, ncols
             precca(i,lchnk)  = surface_state2d(lchnk)%precc(i)
             precla(i,lchnk)  = surface_state2d(lchnk)%precl(i)
             precsca(i,lchnk) = surface_state2d(lchnk)%precsc(i)
             precsla(i,lchnk) = surface_state2d(lchnk)%precsl(i)
          end do
!
! Final call of averaging interval, complete averaging and copy data back
!
       else if (dosw) then
          do i = 1, ncols
             precca(i,lchnk)  = rcount*(precca(i,lchnk) + surface_state2d(lchnk)%precc(i))
             precla(i,lchnk)  = rcount*(precla(i,lchnk) + surface_state2d(lchnk)%precl(i))
             precsca(i,lchnk) = rcount*(precsca(i,lchnk) + surface_state2d(lchnk)%precsc(i))
             precsla(i,lchnk) = rcount*(precsla(i,lchnk) + surface_state2d(lchnk)%precsl(i))
          end do
!
! Intermediate call, add data to accumulators
!
       else
          do i = 1, ncols
             precca(i,lchnk)  = precca(i,lchnk) + surface_state2d(lchnk)%precc(i)
             precla(i,lchnk)  = precla(i,lchnk) + surface_state2d(lchnk)%precl(i)
             precsca(i,lchnk) = precsca(i,lchnk) + surface_state2d(lchnk)%precsc(i)
             precsla(i,lchnk) = precsla(i,lchnk) + surface_state2d(lchnk)%precsl(i)
          end do
       end if
    end do
!
! Reset the counter if the average was just computed
!
    if (dosw) then
       countfa = 0
    end if

    return
  end subroutine ccsmave

!===============================================================================

  subroutine ccsm_msg_getorb

!-----------------------------------------------------------------------
!
! Purpose: Get orbital values from flux coupler
!
! Method:
!
! Author: Erik Kluzek
!         modified for cpl6 by Robert Jacob
!
!-----------------------------------------------------------------------

     use physconst, only:

#include <comctl.h>
#include <comsol.h>

!--------------------------Local Variables------------------------------
    integer cplcdate           ! current date from coupler
    integer cplcsec            ! elapsed sec on current date
    integer info_time          ! T => turn on msg-passing timing
    integer ierr               ! Return error
!
!-----------------------------------------------------------------------
!
!
! Receive first ibuff and rbuff message from coupler. ibuff is currently only used
! to determine if output csm timing will occur.
!
       ibuff(:) = 0
    rbuff(:) = 0.
    call cpl_interface_ibufRecv(cpl_fields_cplname,ibuff,rbuff)

! unload the integer parts
    ierr      = ibuff(cpl_fields_ibuf_rcode)    ! error code
    cplcdate  = ibuff(cpl_fields_ibuf_cdate)    ! current date from coupler
    cplcsec   = ibuff(cpl_fields_ibuf_sec)      ! elapsed sec on current date
    info_time = ibuff(cpl_fields_ibuf_infotim)  ! T => turn on msg-passing timing
    write(6,*)'(CCSM_MSG_GET_ORB): recd d->a initial ibuf '
       call shr_sys_flush(6)


! unload the real parts
    spval  = rbuff(cpl_fields_rbuf_spval)      !Special flag value for data
    eccen  = rbuff(cpl_fields_rbuf_eccen)      !Earth's eccentricity of orbit
    obliqr = rbuff(cpl_fields_rbuf_obliqr)     !Earth's Obliquity radians
    lambm0 = rbuff(cpl_fields_rbuf_lambm0)     !longitude of perihelion at v-equinox
    mvelpp = rbuff(cpl_fields_rbuf_mvelpp)     !Earth's Moving vernal equinox of orbit + pi
!
! Check that data sent is good data and not the special value
!
       call ccsm_compat_check_spval(spval, eccen ,'Eccentricity' )
       call ccsm_compat_check_spval(spval, obliqr,'Obliquity' )
       call ccsm_compat_check_spval(spval, lambm0,'long of perh.' )
       call ccsm_compat_check_spval(spval, mvelpp,'Moving lon of perh')
       write(6,*)'(CCSM_MSG_GET_ORB): eccen:  ', eccen
       write(6,*)'(CCSM_MSG_GET_ORB): obliqr: ', obliqr
       write(6,*)'(CCSM_MSG_GET_ORB): lambm0: ', lambm0
       write(6,*)'(CCSM_MSG_GET_ORB): mvelpp: ', mvelpp
    write(6,*)'(CCSM_MSG_GET_ORB): recd d->a initial real buf'
          call shr_sys_flush(6)
!
! Determine if will output csm timing info.
!
       if (info_time == 0) then
          csm_timing = .false.
       else
          csm_timing = .true.
       end if

    return
  end subroutine ccsm_msg_getorb

!===============================================================================

  subroutine ccsm_msg_sendgrid

!-----------------------------------------------------------------------
!
! Purpose:
! Send grid to flux coupler
!
! Method:
!:
! Author: Mariana Vertenstein
!        modified for cpl6 by Robert Jacob
!
!-----------------------------------------------------------------------

    use infnan
    use commap, only: latdeg, londeg
    use dycore, only: dycore_is
    use time_manager, only: get_nstep, get_step_size

#include <comctl.h>

!--------------------------Local Variables------------------------------
    integer lat, lon, i, j, n     ! loop indices
    integer nstep                  ! current time step
    integer msgpday               ! number of send/recv msgs per day
    integer sizebuf               ! size of buffer for sending grid data to coupler
    integer startpoint            ! starting value for grid numbering scheme
    integer(SHR_KIND_IN) ::  mask(plon,plat)       ! Mask of valid data
    real(r8),allocatable :: sbuf(:,:)  ! array for holding grid data to be sent to coupler
    real(r8) dtime                ! timestep size [s]
    real(r8) area(plon,plat)      ! Area in radians squared for each grid point
    real(r8) clondeg(plon,plat)   ! Longitude grid
    real(r8) clatdeg(plon,plat)   ! latitude grid as 2 dimensional array
    real(r8) ns_vert(4,plon,plat) ! latitude grid vertices
    real(r8) ew_vert(4,plon,plat) ! longitude grid vertices
    real(r8) del_theta            ! difference in latitude at a grid point
    real(r8) del_phi              ! difference in longitude at a grid point
    real(r8) pie                  ! mathmatical constant 3.1415...
    real(r8) degtorad             ! convert degrees to radians
!-----------------------------------------------------------------------


       nstep = get_nstep()
       dtime = get_step_size()
!
! Determine number of send/recv msgs per day
!
       if (flxave) then
          msgpday = nint(86400./dtime)/iradsw
       else
          msgpday = nint(86400./dtime)
       end if
       write(6,*)'(CCSM_MSG_SENDGRID): there are ',msgpday,' send/recv msgs per day'
       call shr_sys_flush(6)
!
! Determine ibuff sent to coupler
!

       ibuff(:) = 0
       ibuff(cpl_fields_ibuf_lsize ) = nlcols             ! local number of columns (gridpoints)
       ibuff(cpl_fields_ibuf_gsize ) = ngcols             ! global number of columns (gridpoints)
       ibuff(cpl_fields_ibuf_lisize) = maxval(nlon(beglat:endlat))  ! send the maximum local longitudes
       ibuff(cpl_fields_ibuf_ljsize) = endlat-beglat+1    ! number of local latitudes
       ibuff(cpl_fields_ibuf_gisize) = plon               ! number of model longitudes
       ibuff(cpl_fields_ibuf_gjsize) = plat               ! number of model latitudes
       ibuff(cpl_fields_ibuf_ncpl)   = msgpday            ! number of send/recv msgs per day
       ibuff(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
       ibuff(cpl_fields_ibuf_dead   ) = 0                 ! not a dead model

! Constants
!
       pie       = acos(-1.)
       degtorad  = pie / 180.0
!
! Mask for which cells are active and inactive and 2D latitude grid
!
       mask(:,:)    = 0        ! Initialize mask so that cells are inactive
       clatdeg(:,:) = spval
       clondeg(:,:) = spval
       do lat = 1, plat
         mask(1:nlon(lat),lat)    = 1     ! Active cells
         clatdeg(1:nlon(lat),lat) = latdeg(lat) ! Put latitude in 2D array
         clondeg(1:nlon(lat),lat) = londeg(1:nlon(lat),lat)
       end do
!
! Send vertices of each grid point
! Verticies are ordered as follows:
! 1=lower left, 2 = upper left, 3 = upper right, 4 = lower right
!
       ns_vert(:,:,:) = spval
       ew_vert(:,:,:) = spval
!
! Longitude vertices
!
       do lat = 1, plat
         ew_vert(1,1,lat)             = (londeg(1,lat) - 360.0 + londeg(nlon(lat),lat))*0.5
         ew_vert(1,2:nlon(lat),lat)   = (londeg(1:nlon(lat)-1,lat) + &
                                         londeg(2:nlon(lat),lat))*0.5
         ew_vert(2,:nlon(lat),lat)    = ew_vert(1,:nlon(lat),lat)  ! Copy lowleft corner to upleft
         ew_vert(3,:nlon(lat)-1,lat)  = ew_vert(1,2:nlon(lat),lat)
         ew_vert(3,nlon(lat),lat)     = (londeg(nlon(lat),lat) + (360.0 + londeg(1,lat)))*0.5
         ew_vert(4,:nlon(lat),lat)    = ew_vert(3,:nlon(lat),lat)  ! Copy lowright corner to upright
       end do
!
! Latitude
!
       if ( dycore_is('LR') )then
         ns_vert(1,:nlon(1),1)         = -90.0 + (latdeg(1) - latdeg(2))*0.5
         ns_vert(2,:nlon(plat),plat)   =  90.0 + (latdeg(plat) - latdeg(plat-1))*0.5
       else
         ns_vert(1,:nlon(1),1)         = -90.0
         ns_vert(2,:nlon(plat),plat)   =  90.0
       end if
       ns_vert(4,:nlon(1),1)       = ns_vert(1,nlon(1),1)        ! Copy lower left to lower right
       ns_vert(3,:nlon(plat),plat) = ns_vert(2,nlon(plat),plat)  ! Copy up left to up right
       do lat = 2, plat
         ns_vert(1,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat-1) )*0.5
         ns_vert(4,:nlon(lat),lat) = ns_vert(1,:nlon(lat),lat)
       end do
       do lat = 1, plat-1
         ns_vert(2,:nlon(lat),lat) = (latdeg(lat) + latdeg(lat+1) )*0.5
         ns_vert(3,:nlon(lat),lat) = ns_vert(2,:nlon(lat),lat)
       end do
!
! Get area of grid cells (as radians squared)
!
       area(:,:) = 0.0
       do lat = 1, plat
         do lon = 1, nlon(lat)
           del_phi = sin( ns_vert(2,lon,lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat) = del_theta*del_phi
         end do
       end do
!
! If grid has a pole point (as in Lin-Rood dynamics
!
      if ( dycore_is('LR') )then
         lat = 1
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi = -sin( latdeg(lat)*degtorad ) + sin( ns_vert(2,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
         lat = plat
!         mask(2:nlon(lat),lat) = 0   ! Only active one point on pole
         do lon = 1, nlon(lat)
           del_phi =  sin( latdeg(lat)*degtorad ) - sin( ns_vert(1,lon,lat)*degtorad )
           del_theta = ( ew_vert(4,lon,lat) - ew_vert(1,lon,lat) )*degtorad
           area(lon,lat)  = del_theta*del_phi
         end do
       end if
       if ( abs(sum(area) - 4.0*pie) > 1.e-12 )then
         write (6,*) 'CCSM_MSG_SENDGRID: sum of areas on globe does not = 4*pi'
         write (6,*) ' sum of areas = ', sum(area)
         call endrun
       end if

! find size of sbuf
       sizebuf=0
       do j=beglat,endlat
	 sizebuf=sizebuf+nlon(j)
       enddo

! NOTE:  Numbering scheme is: West to East and South to North
! starting at south pole.  Should be the same as what's used
! in SCRIP

       allocate(sbuf(sizebuf,cpl_fields_grid_total))

! load in the lat, lon, area, mask, and compute gridpoint numbers for
! points on this processor
       n=0
       startpoint = 0
       do j=1,plat
          do i=1,nlon(j)
            if(get_chunk_owner_p(i,j) .eq. iam) then
	      n=n+1
	      sbuf(n,cpl_fields_grid_lon) = clondeg(i,j)
	      sbuf(n,cpl_fields_grid_lat) = clatdeg(i,j)
	      sbuf(n,cpl_fields_grid_area) = area(i,j)
	      sbuf(n,cpl_fields_grid_mask) = float(mask(i,j))
	      sbuf(n,cpl_fields_grid_index) = startpoint + i
            end if
          enddo
	  startpoint = startpoint + nlon(j)
       enddo
!
! Send ibuff and local grid information to flux coupler.  Initialize cpl6 contracts
!

       call cpl_interface_contractInit(contractS,cpl_fields_atmname, &
	   cpl_fields_cplname,cpl_fields_a2c_fields,ibuff,sbuf)
       call cpl_interface_contractInit(contractR,cpl_fields_atmname, &
	   cpl_fields_cplname,cpl_fields_c2a_fields,ibuff,sbuf)

       write(6,*)'(CCSM_MSG_SENDGRID): sent a->d startup'
       call shr_sys_flush(6)

       deallocate(sbuf)


    return
  end subroutine ccsm_msg_sendgrid

!===============================================================================

  subroutine ccsm_msg_getalb

!-----------------------------------------------------------------------
!
! Purpose:
! Send first time of albedo calculation (along with dummy data) to
! coupler and get albedos along with snow and ocn/ice fractions back
!
! Method:
!
! Author: Mariana Vertenstein
!         modified for cpl6 by Robert Jacob
!
!-----------------------------------------------------------------------

    use comsrf, only: srfflx_state2d, icefrac, ocnfrac, landfrac, &
                      snowhice, snowhland, verify_fractions
    use time_manager, only: get_start_date

#include <comctl.h>

!--------------------------Local Variables------------------------------
    integer i,m,n,lat,lchnk             ! indices
    integer ncols                       ! Number of columns
    integer yr, mon, day                ! year, month, day components of cdatecsm
    integer cdatecsm,cseccsm            ! current date,sec
    integer msgpday                     ! number of send/recv msgs per day
!-----------------------------------------------------------------------
!
!
! Send first time of albedo calculation (along with dummy data) to the flux coupler.
!
       call get_start_date(yr, mon, day, cseccsm)
       cdatecsm = yr*10000 + mon*100 + day

       ibuff(:)  = 0
    ibuff(cpl_fields_ibuf_cdate)   = cdatecsm   ! model date (yyyymmdd)
    ibuff(cpl_fields_ibuf_sec)     = cseccsm    ! elapsed seconds in current day
!    ibuff()  = 0          ! current time step
    ibuff(cpl_fields_ibuf_gisize)  = plon       ! number of model longitudes
    ibuff(cpl_fields_ibuf_gjsize)  = plat       ! number of model latitudes
    ibuff(cpl_fields_ibuf_ncpl)    = msgpday    ! number of send/recv msgs per day
    ibuff(cpl_fields_ibuf_ashift)  = 0          ! albedo calculation time shift
    ibuff(cpl_fields_ibuf_xalbic) = 1           ! tell coupler to do extra albedo calculation
						!   on startup (this routine is already doing that)
    ibuff(cpl_fields_ibuf_userest) = 1          ! use own restart info, not coupler's

! send a dummy message to coupler which is expecting an initial data message.
! Coupler will proceed to the albedo init portion.  After it sends the albedo
! to us, it will wait for the true initial message from physpkg.
    rbuffs(:,:) = spval
    call cpl_interface_contractSend(cpl_fields_cplname,contractS,ibuff,rbuffs)
       if (csm_timing) irtc_s = shr_sys_irtc()
!
! Receive merged surface state from flux coupler.
!
       ibuff(:) = 0
       if (csm_timing) irtc_w = shr_sys_irtc()
    call cpl_interface_contractRecv(cpl_fields_cplname,contractR,ibuff,rbuffr)

    if(masterproc) then
       if (csm_timing) then
          irtc_r = shr_sys_irtc()
          write(6,9099) irtc_s,'a->d sending'
          write(6,9099) irtc_w,'d->a waiting'
          write(6,9099) irtc_r,'d->a received'
       end if
      write(6,*) '(CCSM_MSG_GETALB) recd d->a surface state'
    end if
       call shr_sys_flush(6)

!
! Extract the surface state variables and surface type fractions.
! NOTE: coupler sends a full buffer but only has good data in
! the variables needed below (albedos, t, snow, ifrac, ofrac)
!

! Copy data from coupler buffer into local chunk structure
    call buff_to_chunk(nrcv,plon,rbuffr,recv2d_chunk)

! Copy data from chunks to chunks
    do lchnk=begchunk,endchunk
       ncols = get_ncols_p(lchnk)
       do i=1,ncols
          srfflx_state2d(lchnk)%asdir(i)  = recv2d_chunk(i,cpl_fields_c2a_avsdr,lchnk)  ! Surface state variable
          srfflx_state2d(lchnk)%aldir(i)  = recv2d_chunk(i,cpl_fields_c2a_anidr ,lchnk) ! Surface state variable
          srfflx_state2d(lchnk)%asdif(i)  = recv2d_chunk(i,cpl_fields_c2a_avsdf ,lchnk) ! Surface state variable
          srfflx_state2d(lchnk)%aldif(i)  = recv2d_chunk(i,cpl_fields_c2a_anidf,lchnk)  ! Surface state variable
          srfflx_state2d(lchnk)%ts(i)     = recv2d_chunk(i,cpl_fields_c2a_t,lchnk)      ! Surface state variable
          srfflx_state2d(lchnk)%sst(i)    = recv2d_chunk(i,cpl_fields_c2a_sst,lchnk)    ! Surface state variable
          snowhland(i,lchnk)              = recv2d_chunk(i,cpl_fields_c2a_snowh,lchnk)  ! Surface state variable
          icefrac(i,lchnk)                = recv2d_chunk(i,cpl_fields_c2a_ifrac,lchnk) ! Surface type fraction
          ocnfrac(i,lchnk)                = recv2d_chunk(i,cpl_fields_c2a_ofrac,lchnk) ! Surface type fraction
          landfrac(i,lchnk)               = 1.0 - icefrac(i,lchnk) - ocnfrac(i,lchnk)
!
! Start of bit-for-bit stuff.
! Code between the "bit-for-bit" comments is only there to allow identical
! results with ongoing IPCC coupled runs.  See earlier "bit-for-bit" comments
! in ccsmrcv (above).
!
          if (icefrac(i,lchnk)+landfrac(i,lchnk) > 1.0) then
             icefrac(i,lchnk) = 1.0 - landfrac(i,lchnk)
          end if

          ocnfrac(i,lchnk) = 1.0 - landfrac(i,lchnk) - icefrac(i,lchnk)
!
! end of bit-for-bit stuff.
!
       end do
!
! Ensure that fractions are valid
!
       call verify_fractions (lchnk, ncols)
    end do
!
! Set snowh over ice to zero since flux coupler only returns snowh over land
!
    snowhice(:,:) = 0.0

9099 format('[mp timing]  irtc = ',i20,' ',a)

        return
    end subroutine ccsm_msg_getalb

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Checks that the message recieved from the coupler is compatable
    ! with the type of message that I expect to recieve.
    !
    ! Method:
    ! If the minor version numbers differ I print a warning message.  If the major
    ! numbers differ I abort since that means that the change is drastic enough that
    ! I can't run with the differences.
    !
    ! Author: Erik Kluzek
    !
    !-----------------------arguments---------------------------------------

    subroutine ccsm_msg_compat( cpl_maj_vers, cpl_min_vers, expect_maj_vers, expect_min_vers )

        integer, intent(in) :: cpl_maj_vers                 ! major version from coupler initial ibuff array
        integer, intent(in) :: cpl_min_vers                 ! minor version from coupler initial ibuff array
        integer(SHR_KIND_IN), intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
        integer(SHR_KIND_IN), intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting

        write(6,*)'(CCSM_MSG_COMPAT): This is revision: $Revision: 1.11.2.17.36.2 $'
        write(6,*)'              Tag: $Name: ccsm3_0_1_beta14 $'
        write(6,*)'              of the message compatability interface:'
        if (cpl_min_vers /= expect_min_vers) then
            write(6,*)'WARNING(cpl_compat):: Minor version of coupler messages different than expected: '
            write(6,*)'The version of the coupler being used is: ', cpl_min_vers
            write(6,*)'The version I expect is:                  ', expect_min_vers
        end if
        if ( cpl_maj_vers /= expect_maj_vers )then
            write(6,*) 'ERROR(cpl_compat):: Major version of coupler messages different than expected: '
            write(6,*) 'The version of the coupler being used is: ', cpl_maj_vers
            write(6,*) 'The version I expect is:                  ', expect_maj_vers
            call endrun
        end if

        return
    end subroutine ccsm_msg_compat

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Check that the given piece of real data sent from the coupler is valid data
    ! and not the couplers special data flag.  This ensures that the data
    ! you expect is actually being sent by the coupler.
    !
    ! Method:
    !
    ! Author: Erik Kluzek
    !
    !-----------------------------------------------------------------------

    subroutine ccsm_compat_check_spval( spval, data, string )

        real(r8) , intent(in) ::  spval, data
        character, intent(in) ::  string*(*)

        if ( spval == data )then
            write(6, *)'ERROR::( lsm_compat_check_spval) msg incompatibility'
            write(6, *)'ERROR:: I expect to recieve the data type: ', string
            write(6, *)'from CPL, but all I got was the special data flag'
            write(6, *)'coupler must not be sending this data, you are'
            write(6, *)'running with an incompatable version of the coupler'
            call endrun ('CCSM_COMPAT_CHECK_SPVAL')
        end if

        return
    end subroutine ccsm_compat_check_spval

#endif

end module ccsm_msg
