#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose:
!   Compute ocean to atmosphere surface fluxes of sensible, latent heat
!   and stress components:
!
! Method:
! Assume:
!   1) Neutral 10m drag coeff:
!         cdn = .0027/U10N + .000142 + .0000764 U10N
!   2) Neutral 10m stanton number:
!         ctn = .0327 sqrt(cdn), unstable
!         ctn = .0180 sqrt(cdn), stable
!   3) Neutral 10m dalton number:
!         cen = .0346 sqrt(cdn)
!   4) The saturation humidity of air at T(K):
!         qsat(T)  (kg/m^3)
! Note:
!   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   2) wind speeds should all be above a minimum speed (umin)
!
! Author: Bill Large/M.Vertenstein, Sep. 1995
!
! See P147~148 in CAM 3.0 Description for reference
!
!-----------------------------------------------------------------------

subroutine flxoce(indx    ,npts    ,pmidm1  ,ubot    ,vbot    , &
                  tbot    ,qbot    ,thbot   ,zbot    ,ts      , &
                  ltheat  ,shf     ,lhf     ,taux    ,tauy    , &
                  lwup    ,tref    ,wind10m)

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use physconst,    only: rair, cpair, cpvir, zvir, gravit, stebol

    implicit none

#include <parpbl.h>

    !------------------------------Arguments--------------------------------

    integer , intent(in) :: indx(pcols)    ! column index array (land)
    integer , intent(in) :: npts           ! Number of land points
    real(r8), intent(in) :: pmidm1(pcols)  ! Bottom level pressure
    real(r8), intent(in) :: ubot(pcols)    ! Bottom level u wind
    real(r8), intent(in) :: vbot(pcols)    ! Bottom level v wind
    real(r8), intent(in) :: tbot(pcols)    ! Bottom level temperature
    real(r8), intent(in) :: qbot(pcols)    ! Bottom level specific humidity
    real(r8), intent(in) :: thbot(pcols)   ! Bottom level potential temperature
    real(r8), intent(in) :: zbot(pcols)    ! Bottom level height above surface
    real(r8), intent(in) :: ts(pcols)      ! Surface temperature
    real(r8), intent(in) :: ltheat(pcols)  ! Latent heat for given srf conditions

    real(r8), intent(inout) :: shf(pcols)  ! Initial sensible heat flux (W/m2)
    real(r8), intent(inout) :: lhf(pcols)  ! Initial latent heat flux (W/m2)
    real(r8), intent(inout) :: taux(pcols) ! X surface stress (N/m2)
    real(r8), intent(inout) :: tauy(pcols) ! Y surface stress (N/m2)
    real(r8), intent(inout) :: lwup(pcols) ! Longwave up flux at surface (W/m2)
    real(r8), intent(inout) :: tref(pcols) ! 2m reference temperature
    real(r8), intent(out) :: wind10m(pcols)! added by DONG Li for outputting "u10n"
    !-----------------------------------------------------------------------

    !---------------------------Local variables-----------------------------
    integer i,ii            ! column indices
    real(r8) ssq            ! Surface saturation specific humidity
    real(r8) ustar          ! ustar
    real(r8) tstar          ! tstar
    real(r8) qstar          ! qstar
    real(r8) u10n           ! neutral 10 m wind speed over ocean
    real(r8) vmag           ! Surface wind magnitude
    real(r8) thvbot         ! Bottom lev virtual potential temp
    real(r8) delt           ! potential T difference (K)
    real(r8) delq           ! specific humidity difference (kg/kg)
    real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
    real(r8) rhn            ! sqrt of neutral exchange coeff (heat)
    real(r8) ren            ! sqrt of neutral exchange coeff (tracers)
    real(r8) rd             ! sqrt of exchange coeff (momentum)
    real(r8) rh             ! sqrt of exchange coeff (heat)
    real(r8) re             ! sqrt of exchange coeff (tracers)
    real(r8) hol            ! Ref hgt (10m) / monin-obukhov length
    real(r8) xsq            ! Temporary variable
    real(r8) xqq            ! Temporary variable
    real(r8) alz            ! ln(zbot/z10)
    real(r8) cp             ! Specific heat of moist air
    real(r8) tau            ! Reference height stress
    real(r8) psimh          ! Stability function at reference level (momentum)
    real(r8) psixh          ! Stability function at reference level (heat & tracers)
    real(r8) stable         ! Stability factor
    real(r8) rbot(pcols)    ! Density at bottom model level
    real(r8) bn             ! exchange coef funct for interpolation
    real(r8) bh             ! exchange coef funct for interpolation
    real(r8) fac            ! interpolation factor
    real(r8) ln0            ! log factor for interpolation
    real(r8) ln3            ! log factor for interpolation
    real(r8), parameter :: ztref = 2.0 ! reference height for air temperature
    !-----------------------------------------------------------------------


    !--------------------------Statement functions--------------------------
    real(r8) psimhu         ! Unstable part of psimh
    real(r8) psixhu         ! Unstable part of psixh
    real(r8) qsat           ! Saturation specific humidty of air
    real(r8) cdn            ! Neutral drag coeff at bottom model level
    real(r8) xd             ! Dummy argument
    real(r8) Tk             ! Temperature (K)
    real(r8) Umps           ! Wind velocity (m/sec)

    qsat(Tk)   = 640380./exp(5107.4/Tk)
    cdn(Umps)  = 0.0027/Umps+.000142+.0000764*Umps                    ! see (4.447)
    psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.)-2.*atan(xd)+1.571 ! see (4.442)
    psixhu(xd) = 2.*log((1.+xd*xd)/2.)                                ! see (4.443)
    !-----------------------------------------------------------------------
    !
    ! Loop over ocean points
    !
    do ii = 1, npts
        i = indx(ii)
        !
        !---------------------------------------------------------------
        ! Set up necessary variables
        !---------------------------------------------------------------
        !
        rbot(i) = pmidm1(i)/(rair*tbot(i))
        vmag    = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
        thvbot  = thbot(i) * (1.0 + zvir*qbot(i))
        ssq     = 0.98*qsat(ts(i))/rbot(i)
        delt    = thbot(i)-ts(i)
        delq    = qbot(i)-ssq
        alz     = log(zbot(i)/zref)
        cp      = cpair*(1.+cpvir*ssq)
        !
        !---------------------------------------------------------------
        ! First iteration to converge on Z/L and hence the fluxes
        !---------------------------------------------------------------
        !
        ! Initial guess for roots of neutral exchange coefficients,
        ! assume z/L=0. and u10n is approximated by vmag.
        ! Stable if (thbot > ts ).
        !
        stable = 0.5+sign(0.5_r8, delt)
        rdn = sqrt(cdn(vmag)) ! DONG Li: vmag is not at 10m in initial guess
        rhn = (1.0-stable)*0.0327+stable*0.018
        ren = 0.0346
        !
        ! Initial guess of ustar,tstar and qstar
        !
        ustar = rdn*vmag
        tstar = rhn*delt ! DONG Li: This is different from (4.446)
        qstar = ren*delq !
        !
        ! Compute stability and evaluate all stability functions
        ! Stable if (thbot > ts or hol > 0 )
        !
        hol = xkar*gravit*zbot(i)*(tstar/thvbot+qstar/(1./zvir+qbot(i)))/ustar**2
        hol = sign(min(abs(hol), 10.0_r8), hol)
        stable = 0.5+sign(0.5_r8 , hol)
        xsq   = max(sqrt(abs(1.0-16.0*hol)), 1.0_r8)     !
        xqq   = sqrt(xsq)                                ! see (4.444)
        psimh = -5.0*hol*stable+(1.0-stable)*psimhu(xqq) ! see (4.441) and (4.442)
        psixh = -5.0*hol*stable+(1.0-stable)*psixhu(xqq) ! see (4.441) and (4.443)
        !
        ! Shift 10m neutral wind speed using old rdn coefficient
        !
        rd   = rdn/(1.+rdn/xkar*(alz-psimh))
        u10n = vmag*rd/rdn
        !
        ! Update the neutral transfer coefficients at 10m and neutral stability
        !
        rdn = sqrt(cdn(u10n))
        rhn = (1.0-stable)*0.0327+stable*0.018
        ren = 0.0346
        !
        ! Shift all coeffs to measurement height and stability
        !
        rd = rdn/(1.+rdn/xkar*(alz-psimh))
        rh = rhn/(1.+rhn/xkar*(alz-psixh))
        re = ren/(1.+ren/xkar*(alz-psixh))
        !
        ! Update ustar, tstar, qstar using updated, shifted coeffs
        !
        ustar = rd*vmag
        tstar = rh*delt
        qstar = re*delq
        !
        !---------------------------------------------------------------
        ! Second iteration to converge on Z/L and hence the fluxes
        !---------------------------------------------------------------
        !
        ! Recompute stability & evaluate all stability functions
        ! Stable if (thbot > ts or hol > 0 )
        !
        hol = xkar*gravit*zbot(i)*(tstar/thvbot+qstar/(1./zvir+qbot(i)))/ustar**2
        hol = sign(min(abs(hol), 10._r8), hol)
        stable = 0.5+sign(0.5_r8, hol)
        xsq = max(sqrt(abs(1.-16.*hol)), 1._r8)
        xqq = sqrt(xsq)
        psimh = -5.*hol*stable+(1.-stable)*psimhu(xqq)
        psixh = -5.*hol*stable+(1.-stable)*psixhu(xqq)
        !
        ! Shift 10m neutral wind speed using old rdn coefficient
        !
        rd   = rdn/(1.+rdn/xkar*(alz-psimh))
        u10n = vmag*rd/rdn
        wind10m(i) = u10n ! added by DONG Li for outputting "u10n"
        !
        ! Update the neutral transfer coefficients at 10m and neutral stability
        !
        rdn = sqrt(cdn(u10n))
        rhn = (1.0-stable)*0.0327+stable*0.018
        ren = 0.0346
        !
        ! Shift all coeffs to measurement height and stability
        !
        rd = rdn/(1.0+rdn/xkar*(alz-psimh))
        rh = rhn/(1.0+rhn/xkar*(alz-psixh))
        re = ren/(1.0+ren/xkar*(alz-psixh))
        !
        !---------------------------------------------------------------
        ! Compute the fluxes
        !---------------------------------------------------------------
        !
        ! Update ustar, tstar, qstar using updated, shifted coeffs
        !
        ustar = rd*vmag
        tstar = rh*delt
        qstar = re*delq
        !
        ! Compute surface stress components
        !
        tau     =  rbot(i)*ustar*ustar
        taux(i) = -tau*ubot(i)/vmag
        tauy(i) = -tau*vbot(i)/vmag
        !
        ! Compute heat flux components at current surface temperature
        ! (Define positive latent and sensible heat as upwards into the atm)
        !
        shf(i) = -cp*tau*tstar/ustar
        lhf(i) = -ltheat(i)*tau*qstar/ustar
        lwup(i) = stebol*ts(i)**4

        !---------------------------------------------------------------
        ! Following Geleyn(1988), interpolate ts to fixed height zref
        !---------------------------------------------------------------

        ! Compute function of exchange coefficients. Assume that
        ! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore
        ! 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh

        bn = xkar/rdn
        bh = xkar/rh

        ! Interpolation factor for stable and unstable cases

        ln0 = log(1.0+(ztref/zbot(i))*(exp(bn)-1.0))
        ln3 = log(1.0+(ztref/zbot(i))*(exp(bn-bh)-1.0))
        fac = (ln0-ztref/zbot(i)*(bn-bh))/bh*stable+(ln0-ln3)/bh*(1.0-stable)
        fac = min(max(fac, 0.0_r8), 1.0_r8)

        ! Actual interpolation

        tref(i) = ts(i)+(tbot(i)-ts(i))*fac

    end do

    return
end subroutine flxoce


