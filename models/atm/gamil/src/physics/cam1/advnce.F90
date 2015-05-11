#include <misc.h>
#include <params.h>

subroutine advnce(state)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Advance time information
    !
    ! Method:
    !
    ! Author: CCM1, CMS Contact: J. Truesdale
    !
    !-----------------------------------------------------------------------
    !! (wh 2003.12.27)

    use shr_kind_mod,        only: r8 => shr_kind_r8
    use chemistry,           only: chem_time_interp
    use so4bnd
    use so4bnd_IPCC                                     ! added by WAN Hui
    use ramp_so4_mod
    use time_manager,        only: get_nstep

    use prescribed_aerosols, only:aerint                ! added by SHI Xiangjun 2009-0311
    use physics_types,       only: physics_state        ! added by SHI Xiangjun 2009-0311
    use ppgrid,              only: begchunk, endchunk   ! added by SHI Xiangjun 2009-0311

#if (! defined COUP_CSM) && (! defined COUP_SOM)
    use ice_data,            only: iceint
    use sst_data,            only: sstint
#endif

    ! *** added by DONG Li *** !
#if (defined TRIAL_RUN_FORCING)
    use ForcingScenario
    use SolarForcing, only: SolarForcing_scenario => scenario, SolarForcing_UpdateSolarConstant
    use GHGForcing,   only: GHGForcing_scenario => scenario, GHGForcing_UpdateGHGvmr
#endif
    ! ************************ !

    implicit none

#include <comctl.h>
#include <comlun.h>
#include<RK_or_MG.h>
    !
    !-----------------------------------------------------------------------
    !
    ! Local workspace
    type(physics_state), intent(in), dimension(begchunk:endchunk) :: state !!! sxj-2009-0311
    !
    integer :: nstep             ! current timestep number
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    !
    ! Determine whether it is time for a shortwave or longwave radiation
    ! calculation
    !
    dosw = nstep.eq.0 .or. iradsw.eq.1 .or. (mod(nstep-1,iradsw).eq.0 .and. nstep.ne.1)
    dolw = nstep.eq.0 .or. iradlw.eq.1 .or. (mod(nstep-1,iradlw).eq.0 .and. nstep.ne.1)
    !
    ! Determine whether it is time for an absorptivity/emissivity calculation
    !
    doabsems = nstep.eq.0 .or. iradae.eq.1 .or. (mod(nstep-1,iradae).eq.0 .and. nstep.ne.1)
    aeres = (mod(nstep,iradae).ne.0)
    !
    ! Update ozone data on shortwave or longwave time step.
    ! Note that the ozone data is not needed on a longwave time step unless the
    ! absorptivities are being updated ("doabsems").
    !
    if (dosw .or. dolw) then
        call oznint
    end if
    if (RK_or_MG=='MG')then  !sxj----
        call aerint(state)
    endif
    !write(*,*) "advnce.F90-line67"
    !call endrun

    !
    ! Time interpolate sst data
    !
#ifndef COUP_CSM
#ifdef COUP_SOM
    call t_startf ('somint')
    call somint ()
    call t_stopf ('somint')
#else
    call sstint ()
    call iceint ()
#endif
#endif
    !
    ! Ramping ghg if appropriate
    !
    if (doRamp_ghg ) call ramp_ghg
    if (doCmip5_ghg ) call cmip5_ghg
    !
    ! Ramp solar constant if appropraite
    !
#if (defined TRIAL_RUN_FORCING)
    if (SolarForcing_scenario == RampedScenario) call SolarForcing_UpdateSolarConstant
    if (GHGForcing_scenario == RampedScenario) call GHGForcing_UpdateGHGvmr
#else
    if (doRamp_scon) call ramp_scon
    if (doCmip5_scon) call cmip5_scon
#endif
    !
    ! Time interpolate for chemistry, if appropriate
    !
    if (trace_gas) call chem_time_interp
    !
    ! sulfate aerosols
    !
    if ( doRamp_so4 ) then
        call ramp_so4
        call sulfint
    endif

    if ( doIPCC_so4 ) call sulfint_IPCC

end subroutine advnce
