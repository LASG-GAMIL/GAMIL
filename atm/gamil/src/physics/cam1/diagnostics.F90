#include <misc.h>
#include <params.h>

module diagnostics

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,       only: pcols, pver, pvermx
    use history,      only: outfld
    use constituents, only: pcnst, pnats, cnst_name

    implicit none

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

contains

    !-----------------------------------------------------------------------
    !
    ! Purpose: record dynamics variables on physics grid
    !
    !-----------------------------------------------------------------------

    subroutine diag_dynvar(lchnk, ncol, state)

        use physics_types, only: physics_state
        use physconst,     only: gravit, rga, rair
        use wv_saturation, only: aqsat
#ifdef COUP_CSM
        use ccsm_msg,      only: psl   ! Store sea-level pressure for CCSM
#endif

        integer, intent(in) :: lchnk            ! chunk identifier
        integer, intent(in) :: ncol             ! longitude dimension
        type(physics_state), intent(inout) :: state

        real(r8) ftem(pcols,pver) ! temporary workspace
        real(r8) psl_tmp(pcols)   ! Sea Level Pressure
        real(r8) z3(pcols,pver)   ! geo-potential height
        real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
        real(r8) tem2(pcols,pver) ! temporary workspace

        integer k, m              ! index

        call outfld('T       ', state%t , pcols, lchnk)
        call outfld('PS      ', state%ps, pcols, lchnk)
        call outfld('U       ', state%u , pcols, lchnk)
        call outfld('V       ', state%v , pcols, lchnk)
        do m = 1, pcnst+pnats
            call outfld(cnst_name(m), state%q(1,1,m), pcols, lchnk)
        end do
        call outfld('PHIS    ', state%phis, pcols, lchnk)
        !
        ! Add height of surface to midpoint height above surface
        !
        do k = 1, pver
            z3(:ncol,k) = state%zm(:ncol,k)+state%phis(:ncol)*rga
        end do
        call outfld('Z3      ', z3, pcols, lchnk)

        call outfld('PMID    ', state%pmid, pcols, lchnk) ! added by DONG Li
        !
        ! Output Z3 on 500mb, 300, 50 and 700 mb surface
        !
        call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
        call outfld('Z700    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
        call outfld('Z500    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid,   500._r8, z3, p_surf)
        call outfld('Z050    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
        call outfld('Z300    ', p_surf, pcols, lchnk)
        !
        ! Quadratic height fiels Z3*Z3
        !
        ftem(:ncol,:) = z3(:ncol,:)*z3(:ncol,:)
        call outfld('ZZ      ', ftem, pcols, lchnk)

        ftem(:ncol,:) = z3(:ncol,:)*state%v(:ncol,:)*gravit
        call outfld('VZ      ', ftem, pcols, lchnk)
        !
        ! Meridional advection fields
        !
        ftem(:ncol,:) = state%v(:ncol,:)*state%t(:ncol,:)
        call outfld('VT      ', ftem, pcols, lchnk)

        ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)
        call outfld('VQ      ', ftem, pcols, lchnk)

        ftem(:ncol,:) = state%v(:ncol,:)**2
        call outfld('VV      ', ftem, pcols, lchnk)

        ftem(:ncol,:) = state%v(:ncol,:) * state%u(:ncol,:)
        call outfld('VU      ', ftem, pcols, lchnk)
        !
        ! zonal advection
        !
        ftem(:ncol,:) = state%u(:ncol,:)**2
        call outfld('UU      ', ftem, pcols, lchnk)
        !
        ! Wind speed
        !
        ftem(:ncol,:) = sqrt( state%u(:ncol,:)**2 + state%v(:ncol,:)**2)
        call outfld('WSPEED  ', ftem, pcols, lchnk)
        !
        ! Vertical velocity and advection
        !
        call outfld('OMEGA   ', state%omega, pcols, lchnk)
        ftem(:ncol,:) = state%omega(:ncol,:)*state%t(:ncol,:)
        call outfld('OMEGAT  ', ftem, pcols, lchnk)
        ftem(:ncol,:) = state%omega(:ncol,:)*state%u(:ncol,:)
        call outfld('OMEGAU  ', ftem, pcols, lchnk)
        !
        ! Output omega at 850 and 600 mb pressure levels
        !
        call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
        call outfld('OMEGA850', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 60000._r8, state%omega, p_surf)
        call outfld('OMEGA600', p_surf, pcols, lchnk)
        !
        ! Mass of q, by layer and vertically integrated
        !
        ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
        call outfld('MQ      ', ftem, pcols, lchnk)

        do k = 2, pver
            ftem(:ncol,1) = ftem(:ncol,1)+ftem(:ncol,k)
        end do
        call outfld('TMQ     ', ftem, pcols, lchnk)
        !
        ! Relative humidity
        !
        call aqsat(state%t, state%pmid, tem2, ftem, pcols, ncol, pver, 1, pver)
        ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100.
        call outfld('RELHUM  ', ftem, pcols, lchnk)
        !
        ! Sea level pressure
        !
        call cpslec(ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair)
        call outfld('PSL     ', psl_tmp, pcols, lchnk)
#ifdef COUP_CSM
        psl(:ncol,lchnk) = psl_tmp(:ncol)
#endif
        !
        ! Output T,q,u,v fields on pressure surfaces
        !
        call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
        call outfld('T850    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
        call outfld('T300    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
        call outfld('Q850    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
        call outfld('Q200    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
        call outfld('U850    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
        call outfld('U200    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
        call outfld('V850    ', p_surf, pcols, lchnk)
        call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
        call outfld('V200    ', p_surf, pcols, lchnk)

        ftem(:ncol,:) = state%t(:ncol,:)*state%t(:ncol,:)
        call outfld('TT      ', ftem, pcols, lchnk)

        call outfld('UBOT    ', state%u(1,pver),   pcols, lchnk)
        call outfld('VBOT    ', state%v(1,pver),   pcols, lchnk)
        call outfld('QBOT    ', state%q(1,pver,1), pcols, lchnk)
        call outfld('TBOT    ', state%t(1,pver),   pcols, lchnk)
        call outfld('ZBOT    ', state%zm(1,pver),  pcols, lchnk)

        return
    end subroutine diag_dynvar

    subroutine diag_surf(lchnk, ncol,     shflx,     lhflx,   cflx,    &
                         tref,  trefmxav, trefmnav,                    &
#if (defined COUP_CSM)
                         qref,  rhref,    ps,                          &
#endif
                         taux,  tauy, icefrac, ocnfrac, landfrac,      &
                         tssub, tsnam, ts, sicthk, snowhland, snowhice)

!-----------------------------------------------------------------------
!
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------
#if (defined COUP_CSM)
        use wv_saturation, only: calc_qs
        use time_manager,  only: is_first_step, is_end_curr_day
#endif
!-----------------------------------------------------------------------
!
! Input arguments
!
        integer,  intent(in) :: lchnk               ! chunk identifier
        integer,  intent(in) :: ncol                ! longitude dimension

        real(r8), intent(in) :: shflx(pcols)        ! sensible heat flux (w/m^2)
        real(r8), intent(in) :: lhflx(pcols)        ! latent heat flux (w/m^2)
        real(r8), intent(in) :: cflx(pcols)         ! surface water flux (kg/m^2/s)
        real(r8), intent(in) :: tref(pcols)         ! 2m surface air temperature (not skin temp)
        real(r8), intent(inout) :: trefmnav(pcols)  ! daily minimum tref
        real(r8), intent(inout) :: trefmxav(pcols)  ! daily maximum tref
#if (defined COUP_CSM)
        real(r8), intent(in) :: qref(pcols)         ! 2m surface specific humidity    ! added by DONG Li
        real(r8), intent(inout) :: rhref(pcols)     ! 2m surface relative humidity    ! for FGOALS2.0
        real(r8), intent(in) :: ps(pcols)
        real(r8) qs
#endif
        real(r8), intent(in) :: taux(pcols)         ! x surface stress (zonal) (N/m2)
        real(r8), intent(in) :: tauy(pcols)         ! y surface stress (meridional) (N/m2)
        real(r8), intent(in) :: icefrac(pcols)      ! ice fraction
        real(r8), intent(in) :: ocnfrac(pcols)      ! ocean fraction
        real(r8), intent(in) :: landfrac(pcols)     ! land fraction
        real(r8), intent(in) :: tssub(pcols,pvermx) ! sub-surface soil temperatures
        character(8), intent(in) :: tsnam(pvermx)
        real(r8), intent(in) :: ts(pcols)           ! surface temperature
        real(r8), intent(in) :: sicthk(pcols)       ! sea-ice thickness
        real(r8), intent(in) :: snowhland(pcols)    ! equivalent liquid water snow depth
        real(r8), intent(in) :: snowhice(pcols)     ! equivalent liquid water snow depth

        integer i, k

        call outfld('SHFLX   ', shflx, pcols, lchnk)
        call outfld('LHFLX   ', lhflx, pcols, lchnk)
        call outfld('QFLX    ', cflx,  pcols, lchnk)

        call outfld('TAUX    ', taux,  pcols, lchnk)
        call outfld('TAUY    ', tauy,  pcols, lchnk)
        call outfld('TREFHT  ', tref,  pcols, lchnk)
        call outfld('TREFHTMX', tref,  pcols, lchnk)
        call outfld('TREFHTMN', tref,  pcols, lchnk)
#if (defined COUP_CSM)
        ! added by DONG Li
        call outfld('QREFHT  ', qref,  pcols, lchnk)
        do i = 1, ncol
            call calc_qs(tref(i), ps(i), qs)
            rhref(i) = qref(i)/qs
        end do
        call outfld('RHREFHT ', rhref, pcols, lchnk)
#endif

        call outfld('LANDFRAC', landfrac, pcols, lchnk)
        call outfld('ICEFRAC ', icefrac,  pcols, lchnk)
        call outfld('OCNFRAC ', ocnfrac,  pcols, lchnk)

#if (defined COUP_CSM)
        !
        ! Compute daily minimum and maximum of TREF
        !
        if (is_first_step()) then  ! initialize trefmxav and trefmnav
            trefmxav(:ncol) = -1.0e36
            trefmnav(:ncol) =  1.0e36
        end if
        do i = 1,ncol
            trefmxav(i) = max(tref(i), trefmxav(i))
            trefmnav(i) = min(tref(i), trefmnav(i))
        end do
        if (is_end_curr_day()) then
            ! DONG Li: "TREFMXAV" and "TREFMNAV" are only output in coupled mode?
            call outfld('TREFMXAV', trefmxav, pcols, lchnk)
            call outfld('TREFMNAV', trefmnav, pcols, lchnk)
            trefmxav(:ncol) = -1.0e36
            trefmnav(:ncol) =  1.0e36
        end if
#else
        do k = 1, pvermx
            call outfld(tsnam(k), tssub(1,k), pcols, lchnk)
        end do
        call outfld('SICTHK  ', sicthk, pcols, lchnk)
#endif

        call outfld('TS      ',ts,      pcols,   lchnk     )
        call outfld('TSMN    ',ts,      pcols,   lchnk     )
        call outfld('TSMX    ',ts,      pcols,   lchnk     )
        call outfld('SNOWHLND',snowhland,   pcols,   lchnk     )
        call outfld('SNOWHICE',snowhice ,   pcols,   lchnk     )

        return
    end subroutine diag_surf

end module diagnostics
