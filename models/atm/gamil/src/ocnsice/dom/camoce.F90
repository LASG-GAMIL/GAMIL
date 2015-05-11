#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
! CAM ocean and sea ice surface fluxes.
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------
!
! $Id: camoce.F90,v 1.1.2.3 2002/06/15 13:48:55 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

subroutine camoce(srf_state,srfflx)

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use pspect
    use comsrf,       only: surface_state,srfflx_parm,ocnfrac
    use phys_grid,    only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use sst_data,     only: sstan
    use time_manager, only: get_nstep, get_step_size, get_curr_calday

    implicit none

    type(surface_state), intent(inout) :: srf_state(begchunk:endchunk)
    type(srfflx_parm),   intent(inout) :: srfflx(begchunk:endchunk)

#include <comctl.h>

    integer nstep               ! current timestep number
    integer dtime               ! timestep size [seconds]
    real(r8) cdaynext           ! calendar day for next timestep
    real(r8) clat(pcols)        ! current latitudes(radians)
    real(r8) clon(pcols)        ! current longitudes(radians)
    real(r8) cosznext(pcols)    ! cosine solar zenith angle next timestep
    integer ncol                ! number of columns in chunk
    integer lchnk               ! chunk index

    !
    ! Calendar day for next time step
    !
    nstep = get_nstep()
    dtime = get_step_size()
    cdaynext = get_curr_calday(offset=dtime)

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, CLAT, CLON, COSZNEXT)

    do lchnk = begchunk, endchunk
        ncol = get_ncols_p(lchnk)
        !
        ! Update ts and for open ocean
        !
        if (anncyc .and. mod(nstep, itsst) == 0) then
            call sstan(lchnk, ncol, ocnfrac(1,lchnk), srfflx(lchnk)%ts)
        end if
        !
        ! Ocean surface fluxes and temperatures
        !
        call srfoce(lchnk, ncol, &
                    ocnfrac(1,lchnk)     , srf_state(lchnk)%ubot,  &
                    srf_state(lchnk)%vbot, srf_state(lchnk)%tbot,  &
                    srf_state(lchnk)%qbot, srf_state(lchnk)%thbot, &
                    srf_state(lchnk)%zbot, srf_state(lchnk)%pbot,  &
                    srfflx(lchnk)%cflx, &
                    srfflx(lchnk)%wsx   , srfflx(lchnk)%wsy ,&
                    srfflx(lchnk)%ts    , srfflx(lchnk)%shf ,&
                    srfflx(lchnk)%lhf   , srfflx(lchnk)%lwup,&
                    srfflx(lchnk)%tref)
        !
        ! Albedos for next time step
        !
        call get_rlat_all_p(lchnk, ncol, clat)
        call get_rlon_all_p(lchnk, ncol, clon)
        call zenith(cdaynext, clat, clon, cosznext, ncol)

        call albocean(lchnk, ncol, cosznext, &
                      srfflx(lchnk)%asdir, srfflx(lchnk)%aldir, &
                      srfflx(lchnk)%asdif, srfflx(lchnk)%aldif)
    end do

    return
end subroutine camoce
