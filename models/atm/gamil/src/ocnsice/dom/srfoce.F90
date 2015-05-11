#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
! Compute ocean to atmosphere sensible/latent heat fluxes and
! stress components
!
! Method:
! Follows the same basic parameterizations as for ocean surfaces
!
! Author: CCM1
!
!-----------------------------------------------------------------------
!
! $Id: srfoce.F90,v 1.5.6.2 2002/06/15 13:48:55 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

subroutine srfoce(lchnk   ,ncol    , &
                  ocnfrac ,ubot    ,vbot    ,tbot    ,qbot    , &
                  thbot   ,zbot    ,pmidm1  ,qflx    , &
                  taux    ,tauy    ,ts      ,shflx   ,lhflx   , &
                  lwup    ,tref    )

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use comsrf,       only: plevmx
    use tracers,      only: pcnst, pnats
    use history,      only: outfld ! added by DONG Li for outputting

    implicit none
#include <comtsc.h>

    !------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                      ! chunk identifier
    integer, intent(in) :: ncol                       ! number of atmospheric columns

    real(r8), intent(in) :: ocnfrac(pcols)            ! ocean fraction
    real(r8), intent(in) :: ubot(pcols)               ! Bottom level u wind
    real(r8), intent(in) :: vbot(pcols)               ! Bottom level v wind
    real(r8), intent(in) :: tbot(pcols)               ! Bottom level temperature
    real(r8), intent(in) :: qbot(pcols)               ! Bottom level specific humidity
    real(r8), intent(in) :: thbot(pcols)              ! Bottom level potential temperature
    real(r8), intent(in) :: zbot(pcols)               ! Bottom level height above surface
    real(r8), intent(in) :: pmidm1(pcols)             ! Bottom level pressure
    real(r8), intent(inout) :: qflx(pcols,pcnst+pnats)! Constituent flux (kg/m2/s)
    real(r8), intent(inout) :: taux(pcols)            ! X surface stress (N/m2)
    real(r8), intent(inout) :: tauy(pcols)            ! Y surface stress (N/m2)
    real(r8), intent(inout) :: ts(pcols)              ! surface temperature (K)
    real(r8), intent(inout) :: shflx(pcols)           ! Surface sensible heat flux (J/m2/s)
    real(r8), intent(inout) :: lhflx(pcols)           ! Surface latent   heat flux (J/m2/s)
    real(r8), intent(inout) :: lwup(pcols)            ! surface longwave up flux (W/m2)
    real(r8), intent(inout) :: tref(pcols)            ! 2m reference temperature

    !---------------------------Local variables-----------------------------
    integer i, ii             ! column indices
    integer m                 ! constituent index
    integer indx(pcols)       ! column index array (land)
    integer npts              ! Number of land points
    real(r8) ltheat(pcols)    ! Latent heat for given sfc conditions
    real(r8) wind10m(pcols)   ! added by DONG Li for outputting "u10n"
    !-----------------------------------------------------------------------
    !
    ! Set up index array of ocean surfaces
    !
    npts = 0
    do i = 1, ncol
        if (ocnfrac(i) > 0.) then
            npts = npts + 1
            indx(npts) = i
        end if
    end do
    if (npts == 0) return
    !
    ! Determine latent heat
    !
    do ii = 1, npts
        i = indx(ii)
        ltheat(i) = latvap
    end do
    !
    ! Compute surface fluxes, derviatives, and exchange coefficiants
    !
    call flxoce(indx    ,npts    ,pmidm1  ,ubot    ,vbot      , &
                tbot    ,qbot    ,thbot   ,zbot    ,ts        , &
                ltheat  ,shflx   ,lhflx   ,taux    ,tauy      , &
                lwup    ,tref    ,wind10m )

    call outfld('WIND10M ', wind10m, pcols, lchnk) ! added by DONG Li
    !
    ! Evaluate contituent fluxes
    !
    do ii = 1, npts
        i = indx(ii)
        qflx(i,1) = lhflx(i)/ltheat(i)
    end do
    !
    ! Set non-water constituent fluxes to zero
    !
    do m = 2, pcnst+pnats
        do ii = 1, npts
            i = indx(ii)
            qflx(i,m) = 0.
        end do
    end do

    return
end subroutine srfoce

