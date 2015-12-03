#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
! Loop over time, calling driving routines for physics, dynamics,
! transport
!
! Method:
!
! Author:
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Restructured:      J. Truesdale, May 1999
!
!-----------------------------------------------------------------------

subroutine stepon

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use history,        only: wshist, wrapup
    use pmgrid
    ! DONG Li: What is "rgrid" for?
    use rgrid
    use prognostics
    use comfm1
    use buffer
    use restart,        only: write_restart
#ifdef COUP_CSM
    use ccsm_msg,       only: csmstop, ccsmfin
#endif

    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state, physics_tend
    use dp_coupling,    only: d_p_coupling, p_d_coupling
    use commap
    use physconst,      only: gravit
    use time_manager,   only: advance_timestep, get_step_size, get_nstep, &
                              is_first_step, is_first_restart_step, &
                              is_last_step, is_end_curr_day, get_curr_calday, &
                              dtdy ! added by WANG Hui

    implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>

    type(physics_state), allocatable :: phys_state(:)
    type(physics_state), allocatable :: phys_state0(:) ! added by WAN Hui, according to P Liu 2003)
    type(physics_tend ), allocatable :: phys_tend(:)

    real(r8), allocatable :: t2(:,:,:) ! temp tendency
    real(r8), allocatable :: fu(:,:,:) ! u wind tendency
    real(r8), allocatable :: fv(:,:,:) ! v wind tendency

    real(r8) rpmid(plond,plev)
    real(r8) pdel(plond,plev)
    real(r8) pint(plond,plevp)
    real(r8) pmid(plond,plev)
    real(r8) dtime               ! timestep size  (physics package)
    real(r8) ztodt               ! twice time step unless nstep=0
    real(r8) wcstart, wcend      ! wallclock timestamp at start, end of timestep
    real(r8) usrstart, usrend    ! user timestamp at start, end of timestep
    real(r8) sysstart, sysend    ! sys timestamp at start, end of timestep

    real(r8) calday              ! current calendar day
    integer nseq

    integer i, k, lat, j, begj   ! longitude,level,latitude indices
    !
    ! Externals
    !
    logical, external :: rstwr  ! whether or not to write restart files
    !
    !-----------------------------------------------------------------------
    call t_startf('stepon_startup'); if(masterproc) write(6,*) '+++++ stepon_startup +++++'
    dtime = get_step_size();         if(masterproc) write(6,*) 'dtime = ', dtime
                                     if(masterproc) write(6,*) 'dtdy  = ', dtdy
    nseq  = dtime/(dtdy-0.01);       if(masterproc) write(6,*) 'nseq  = ', nseq

    pmtop = pmtop*0.01d0

    ! WAN Hui 2003.07.08)

    if (is_first_step()) then
        if(masterproc) write(6, "('Notice: stepon: first step start')")
        itime = 0
        !
        ! Calculate vertical motion field
        !
        omga(:,:,:) = 0.0
        if(masterproc) write(6, "('Notice: stepon: set omga to zero')")
        call init_ac_switching(pmtop) ! added by WAN Hui 2003.10.28
    else if (is_first_restart_step()) then
        !sq(:,:,:) = 0.0
        ! DONG Li: clean this out
        if(masterproc) write(6, "('Notice: stepon: sq set to 0.0 in comfm1 (Please clarify this)')")
    else
        if(masterproc) write(*, "('Error: neither first_step nor first_restart_step')")
        call endrun
    end if

    allocate(phys_state(begchunk:endchunk))
    allocate(phys_state0(begchunk:endchunk)) ! added by WAN Hui
    allocate(phys_tend(begchunk:endchunk))
    allocate(t2(beglonex:endlonex,plev,beglat:endlat))
    allocate(fu(beglonex:endlonex,plev,beglat:endlat))
    allocate(fv(beglonex:endlonex,plev,beglat:endlat))
    !
    ! Beginning of basic time step loop
    !
    call t_stopf ('stepon_startup')

    ! Begin time loop.

    do

        call t_startf('stepon_st')
        if (masterproc .and. print_step_cost) then
            call t_stampf(wcstart, usrstart, sysstart)
        end if

        ! DONG Li: clarify this
        !ztodt = 2.0*dtime
        ztodt = dtime

        calday = get_curr_calday()

        if (masterproc) then
            write(6, *)
            write(6, *) 'date:', calday
        end if

        !----------------------------------------------------------
        ! PHYSPKG  Call the Physics package
        !----------------------------------------------------------
        if (masterproc) write(6, *) '------physpkg------'

        begj = beglatex+numbnd

        call t_stopf('stepon_st')
        call t_startf('d_p_coupling')
        call d_p_coupling(ps(beglonex,beglat,n3m2), t3(beglonex,1,begj,n3m2), u3(beglonex,1,begj,n3m2), &
                          v3(beglonex,1,begj,n3m2), q3(beglonex,1,1,begj,n3m2), &
                          q31(beglonex,1,begj), t31(beglonex,1,begj), q32(beglonex,1,begj), t32(beglonex,1,begj),&
                          omga, phis, phys_state)
        call t_stopf('d_p_coupling')
        !!(wh, according to P Liu 2003)
        if (is_first_restart_step()) then
            call t_startf('d_p_coupling')

            call d_p_coupling(ps(beglonex,beglat,n3m2), t3(beglonex,1,begj,n3m2), u3(beglonex,1,begj,n3m2), &
                              v3(beglonex,1,begj,n3m2), q3(beglonex,1,1,begj,n3m2), &
                              q31(beglonex,1,begj), t31(beglonex,1,begj), q32(beglonex,1,begj), t32(beglonex,1,begj),&!(ljli)
                              omga, phis, phys_state0)
            call t_stopf('d_p_coupling')
        else
            !for Tiedtke scheme
            call t_startf('d_p_coupling')
            call d_p_coupling(ps(beglonex,beglat,n3), t3(beglonex,1,begj,n3), u3(beglonex,1,begj,n3), &
                              v3(beglonex,1,begj,n3), q3(beglonex,1,1,begj,n3), &
                              q31(beglonex,1,begj), t31(beglonex,1,begj), q32(beglonex,1,begj), t32(beglonex,1,begj),&!(ljli)
                              omga, phis, phys_state0)
            call t_stopf('d_p_coupling')
        end if
        !!(wh, according to P Liu 2003)
        call t_startf('phys_driver')

        if (masterproc) then
            write(6, *) 'ideal_phys =', ideal_phys
            write(6, *) 'adiabatic  =', adiabatic
        end if

        if (ideal_phys) then
            call phys_idealized(phys_state, phys_tend, ztodt, sigl)
        else if (adiabatic) then
            call phys_adiabatic(phys_state, phys_tend)
        else
            call physpkg( &
                phys_state, phys_state0, w, ztodt, phys_tend,     &
                cld(1,1,begchunk,n3m2),   cld(1,1,begchunk,n3),   &
                tcwat(1,1,begchunk,n3m2), tcwat(1,1,begchunk,n3), &
                qcwat(1,1,begchunk,n3m2), qcwat(1,1,begchunk,n3), &
                lcwat(1,1,begchunk,n3m2), lcwat(1,1,begchunk,n3))
        end if
        call t_stopf('phys_driver')

        call t_startf('p_d_coupling')
        call p_d_coupling(phys_state, phys_tend, t2, fu, fv, &
            qminus(beglonex,1,1,begj), q3(beglonex,1,1,begj,n3), q31(beglonex,1,begj), t31(beglonex,1,begj))
        call t_stopf('p_d_coupling')

        !----------------------------------------------------------
        ! DYNPKG Call the Dynamics Package
        !----------------------------------------------------------
        if (masterproc) write(6, *) '------dynpkg------'

        call t_startf('dynpkg')

        ! accumulate su, sv, st and update q

        call a_c_switching(fu, fv, t2, beglat, endlat)   !!(wh 2003.10.28)

        call dynpkg(dtdy, nseq)        !!(wh 2003.10.23)

        ! prepare data for physics

        call c_a_switching(pmtop)            !!(wh 2003.10.28)

        call t_stopf('dynpkg')
        !
        ! Shift time pointers
        !
        call shift_time_indices

        call t_startf('stepon_st')
        if (is_first_restart_step()) then
            call print_memusage
        end if

        ! Set end of run flag.

#ifndef COUP_CSM
        if (is_last_step()) nlend = .true.
#else
        if (csmstop) then
            if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
            if (is_end_curr_day()) nlend = .true.
        end if
#endif
        !
        !----------------------------------------------------------
        ! History and restart logic: Write and/or dispose history tapes if required
        !----------------------------------------------------------
        !
        call t_startf ('wshist')
        call wshist ()
        call t_stopf ('wshist')
        !
        ! Write restart file
        !
        if (rstwr() .and. nrefrq /= 0) then
            call t_startf ('write_restart')
            call write_restart
            call t_stopf ('write_restart')
        end if
        !
        ! Dispose necessary files
        !
        call t_startf ('wrapup')
        call wrapup
        call t_stopf ('wrapup')

        if (masterproc .and. print_step_cost) then
            call t_stampf (wcend, usrend, sysend)
            write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
        end if
        !
        ! Advance timestep before returning to top of loop
        !
        call advance_timestep()
        call t_stopf('stepon_st')
        !
        ! Check for end of run
        !
        if (nlend) then
            deallocate(phys_state)
            deallocate(phys_state0)   !!(wh)
            deallocate(phys_tend)
            deallocate(t2)
            deallocate(fu)
            deallocate(fv)
#ifdef COUP_CSM
            call ccsmfin
#endif
            return
        end if

    end do  ! End of timestep loop

end subroutine stepon
