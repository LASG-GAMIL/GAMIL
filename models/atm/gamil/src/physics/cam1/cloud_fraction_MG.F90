module cloud_fraction_MG

    use shr_kind_mod,   only: r8 => shr_kind_r8

    implicit none

    private

    save

    public cldfrc_MG   !  Computation of cloud fraction

    real(r8), parameter :: rhminl = 0.905_r8         ! minimum rh for low stable clouds (.915_r8)
    real(r8), parameter :: rhminh = 0.78_r8          ! minimum rh for high stable clouds (.80_r8)
    real(r8), parameter :: sh1    = 0.04_r8
    real(r8), parameter :: sh2    = 500.0_r8         ! parameters for shallow convection cloud fraction
    real(r8), parameter :: dp1    = 0.10_r8
    real(r8), parameter :: dp2    = 500.0_r8         ! parameters for deep convection cloud fraction
    real(r8), parameter :: premit = 750.0e2_r8       ! top pressure bound for mid level cloud
    real(r8), parameter :: pnot   = 1.0e5_r8         ! reference pressure
    real(r8), parameter :: lapse  = 6.5e-3_r8        ! U.S. Standard Atmsophere lapse rate
    real(r8), parameter :: premib_uw  = 750.0e2_r8   ! bottom pressure bound of middle cloud for UW
    real(r8), parameter :: premib_cam = 750.0e2_r8   ! bottom pressure bound of middle cloud for CAM
    real(r8), parameter :: premib = premib_cam       ! bottom pressure bound of middle cloud
    real(r8), parameter :: pretop = 1.0e2_r8         ! pressure bounding high cloud

    ! DONG Li: "inversion_cld_off" has not been initialized!!!
    logical inversion_cld_off                        ! Turns off stratification-based cld frc

contains

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Compute cloud fraction
    !
    !
    ! Method:
    ! This calculate cloud fraction using a relative humidity threshold
    ! The threshold depends upon pressure, and upon the presence or absence
    ! of convection as defined by a reasonably large vertical mass flux
    ! entering that layer from below.
    !
    ! Author: Many. Last modified by Jim McCaa
    !
    !-----------------------------------------------------------------------

    subroutine cldfrc_MG(state, cmfmc, cmfmc2, zdu,  &
                         landfrac, ocnfrac,          & ! DONG Li: These variables come from "comsrf", consider to use that module directly?
                         ts, sst, ps, snowh,         & !
                         cloud, clc, concld, cldst,  &
                         relhum, rhu00, dindex)


        use physics_types,  only: physics_state
        use tracers,        only: ixmoist
        use ppgrid
        use physconst,      only: cappa, gravit, rair, tmelt
        use cldconst
        use wv_saturation,  only: aqsat
        use phys_grid,      only: get_rlat_all_p, get_rlon_all_p

#include<misc.h>

        type(physics_state), intent(in) :: state
        real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
        real(r8), intent(in) :: cmfmc2(pcols,pver)    ! shallow convective mass flux--m sub c
        real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection
        real(r8), intent(in) :: landfrac(pcols)       ! Land fraction
        real(r8), intent(in) :: ocnfrac(pcols)        ! Ocean fraction
        real(r8), intent(in) :: ts(pcols)             ! surface temperature
        real(r8), intent(in) :: sst(pcols)            ! sea surface temperature
        real(r8), intent(in) :: ps(pcols)             ! surface pressure
        real(r8), intent(in) :: snowh(pcols)          ! snow depth (liquid water equivalent)
        integer, intent(in) :: dindex                 !

        real(r8), intent(out) :: cloud(pcols,pver)    ! cloud fraction
        real(r8), intent(out) :: clc(pcols)           ! vertical integrated convective cloud fraction
        real(r8), intent(out) :: concld(pcols,pver)   ! convective cloud fraction
        real(r8), intent(out) :: cldst(pcols,pver)    ! stratus cloud fraction
        real(r8), intent(out) :: relhum(pcols,pver)   ! relative humidity
        real(r8), intent(out) :: rhu00(pcols,pver)    ! relative humidity threshold for cloud

        real(r8) dthdpmn(pcols)        ! most stable lapse rate below 750 mb
        real(r8) dthdp                 ! lapse rate (intermediate variable)
        real(r8) es(pcols,pver)        ! saturation vapor pressure
        real(r8) qs(pcols,pver)        ! saturation specific humidity
        real(r8) rhwght                ! weighting function for rhlim transition
        real(r8) rh(pcols,pver)        ! relative humidity
        real(r8) rhdif                 ! intermediate scratch variable
        real(r8) strat                 ! intermediate scratch variable
        real(r8) theta(pcols,pver)     ! potential temperature
        real(r8) rhlim                 ! local rel. humidity threshold estimate
        real(r8) coef1                 ! coefficient to convert mass flux to mb/d
        real(r8) clrsky(pcols)         ! temporary used in random overlap calc
        real(r8) rpdeli(pcols,pver-1)  ! 1./(pmid(k+1)-pmid(k))
        real(r8) rhpert                !the specified perturbation to rh
        real(r8) deepcu                ! deep convection cloud fraction
        real(r8) shallowcu             ! shallow convection cloud fraction
        real(r8) sst_K(pcols)          ! added by SHI Xiangjun

        logical cldbnd(pcols)          ! region below high cloud boundary

        integer i, ierror, k           ! column, level indices
        integer kp1
        integer kdthdp(pcols)          ! index of interface of max inversion
        integer numkcld                ! number of levels in which to allow clouds

        real(r8) thetas(pcols)         ! ocean surface potential temperature
        real(r8) clat(pcols)           ! current latitudes(radians)
        real(r8) clon(pcols)           ! current longitudes(radians)
        !
        ! Statement functions
        !
        logical land

        land(i) = nint(landfrac(i)) == 1

        !==================================================================================
        ! PHILOSOPHY OF PRESENT IMPLEMENTATION
        !
        ! There are three co-existing cloud types: convective, inversion related low-level
        ! stratocumulus, and layered cloud (based on relative humidity).  Layered and
        ! stratocumulus clouds do not compete with convective cloud for which one creates
        ! the most cloud.  They contribute collectively to the total grid-box average cloud
        ! amount.  This is reflected in the way in which the total cloud amount is evaluated
        ! (a sum as opposed to a logical "or" operation)
        !
        !==================================================================================

        ! set defaults for rhu00
        rhu00  = 2.0_r8
        ! define rh perturbation in order to estimate rhdfda
        rhpert = 0.01_r8
        !
        ! Evaluate potential temperature and relative humidity
        !
        call aqsat(state%t, state%pmid, es, qs, pcols, state%ncol, pver, 1, pver)

        do k = 1, pver
            do i = 1, state%ncol
                theta(i,k)  = state%t(i,k)*(pnot/state%pmid(i,k))**cappa
                rh(i,k)     = state%q(i,k,ixmoist)/qs(i,k)*(1.0_r8+real(dindex,r8)*rhpert)
                !
                !  record relhum, rh itself will later be modified related with concld
                !
                relhum(i,k) = rh(i,k)
                cloud(i,k)  = 0.0_r8
                cldst(i,k)  = 0.0_r8
                concld(i,k) = 0.0_r8
            end do
        end do
        !
        ! Initialize other temporary variables
        !
        ierror = 0
        do i = 1, state%ncol
!ljli20121023 this line is for cmip test           sst_K(i) = sst(i)  !! sxj
             sst_K(i) = sst(i)+273.15    !if amip test, 273.15 is required. 
            ! Adjust thetas(i) in the presence of non-zero ocean heights.
            ! This reduces the temperature for positive heights according to a standard lapse rate.
            if (ocnfrac(i) > 0.01_r8) then
                thetas(i) = (sst_K(i)-lapse*state%phis(i)/gravit)*(pnot/ps(i))**cappa
                !thetas(i)  = ts(i)*(pnot/ps(i))**cappa ! SHI Xiangjun
                if (sst_K(i) < 260.0_r8) ierror = i
            end if
            clc(i) = 0.0_r8
        end do
        coef1 = gravit*864.0_r8    ! conversion to millibars/day

        !if (ierror > 0) then
            !write(6, "('Warning: cloud_fraction_MG: cold sst (', F5.1, ') encountered ')", advance="no") sst_K(ierror)
            !write(6, "('at chunk ', I5, ' column ', I5, ' with fraction ', F)") lchnk, ierror, ocnfrac(ierror)
        !end if

        ! DONG Li: Can we use "pdel" here?
        do k = 1, pver-1
            do i = 1, state%ncol
                rpdeli(i,k) = 1.0_r8/(state%pmid(i,k+1)-state%pmid(i,k))
            end do
        end do

        ! DONG Li: "concld" has been zeroed already!!!
        !do k = 1, pver
        !    do i = 1, state%ncol
        !        concld(i,k) = 0.0_r8
        !    end do
        !end do

        ! ---------------------------------------------------------------------------
        ! 1. Calculate convective cloud fraction
        ! ---------------------------------------------------------------------------
        !
        ! Estimate of local convective cloud cover based on convective mass flux
        ! Modify local large-scale relative humidity to account for presence of
        ! convective cloud when evaluating relative humidity based layered cloud amount
        !
        ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
        ! shallow and deep convective cloudiness are evaluated separately (since processes
        ! are evaluated separately) and summed
        !
#ifndef PERGRO
        do k = 1, pver-1
            do i = 1, state%ncol
                shallowcu   = max(0.0_r8, min(sh1*log(1.0_r8+sh2*              cmfmc2(i,k+1) ), 0.30_r8))
                deepcu      = max(0.0_r8, min(dp1*log(1.0_r8+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))), 0.60_r8))
                concld(i,k) = min(shallowcu+deepcu, 0.80_r8)
                rh(i,k)     = (rh(i,k)-concld(i,k))/(1.0_r8 - concld(i,k))
            end do
        end do
#endif

        ! added by SHI Xiangjun
        do i = 1, state%ncol
            clrsky(i) = 1.0
        end do
        do k = pver, 1, -1
            do i = 1, state%ncol
                clrsky(i) = clrsky(i)*(1.0-concld(i,k))
            end do
        end do
        do i = 1, state%ncol
            clc(i) = 1.0-clrsky(i)
        end do

        ! ---------------------------------------------------------------------------
        ! 2. Calculate layered cloud fraction
        ! ---------------------------------------------------------------------------
        !
        ! Begin the evaluation of layered cloud amount based on (modified) RH
        !
        numkcld = pver
        do k = 2, numkcld
            kp1 = min(k+1, pver)
            do i = 1, state%ncol
                cldbnd(i) = state%pmid(i,k) >= pretop
                if (state%pmid(i,k) >= premib) then
                    !==============================================================
                    ! This is the low cloud (below premib) block
                    !==============================================================
                    ! enhance low cloud activation over land with no snow cover
                    if (land(i) .and. (snowh(i) <= 0.000001_r8)) then
                        rhlim = rhminl-0.10_r8
                    else
                        rhlim = rhminl
                    end if
                    rhdif = (rh(i,k)-rhlim)/(1.0_r8-rhlim)
                    cloud(i,k) = min(0.999_r8, max(rhdif, 0.0_r8)**2)
                    ! SJV: decrease cloud amount if very low water vapor content
                    ! (thus very cold): "freeze dry"
                    cloud(i,k) = cloud(i,k)*max(0.15_r8, min(1.0_r8, state%q(i,k,ixmoist)/0.0030_r8))
                else if (state%pmid(i,k) < premit) then
                    !==============================================================
                    ! This is the high cloud (above premit) block
                    !==============================================================
                    !
                    rhlim = rhminh
                    rhdif = (rh(i,k)-rhlim)/(1.0_r8-rhlim)
                    cloud(i,k) = min(0.999_r8, max(rhdif, 0.0_r8)**2)
                else
                    !==============================================================
                    ! This is the middle cloud block
                    !==============================================================
                    !
                    ! linear rh threshold transition between thresholds for low & high cloud
                    !
                    rhwght = (premib-max(state%pmid(i,k), premit))/(premib-premit)
                    if (land(i) .and. (snowh(i) <= 0.000001_r8)) then
                        rhlim = rhminh*rhwght+(rhminl-0.10_r8)*(1.0_r8-rhwght)
                    else
                        rhlim = rhminh*rhwght+rhminl*(1.0_r8-rhwght)
                    end if
                    rhdif = (rh(i,k)-rhlim)/(1.0_r8-rhlim)
                    cloud(i,k) = min(0.999_r8, max(rhdif, 0.0_r8)**2)
                end if
                !==================================================================================
                ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
                !==================================================================================
                !
                ! save rhlim to rhu00, it handles well by itself for low/high cloud
                !
                rhu00(i,k) = rhlim
            end do
        end do

        ! ---------------------------------------------------------------------------
        ! 3. Calculate marine stratus fraction
        ! ---------------------------------------------------------------------------
        !
        ! Marine stratus should be a special case of layered cloud.
        ! Cloud currently contains layered cloud determined by rh criteria.
        ! Take the maximum of the diagnosed layered cloud or stratocumulus
        !
        ! Some observations about the following section of code (missed in earlier look)
        ! K700 is set as a constant based on hybrid coordinate: it does not depend on
        ! local pressure; there is no pressure ramp => looks level dependent and
        ! discontinuous in space (i.e., stratus will end suddenly with no transition)
        !
        ! It appears that stratus is evaluated according to Klein and Hartmann; however,
        ! The actual stratus amount (cldst) appears to depend directly on the rh below
        ! the strongest part of the low level inversion.
        !
        if (.not. inversion_cld_off) then
            ! DONG Li: check out what is the value of "inversion_cld_off"
            !write(6, "('Debug: cloud_fraction_MG: inversion_cld_off is false')")
            !
            ! Find most stable level below 750 mb for evaluating stratus regimes
            !
            do i = 1, state%ncol
                ! Nothing triggers unless a stability greater than this minimum threshold is found
                dthdpmn(i) = -0.125_r8
                kdthdp(i) = 0
            end do
            do k = 2, pver
                do i = 1, state%ncol
                    if (state%pmid(i,k) >= premib .and. ocnfrac(i) > 0.01_r8) then
                        ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
                        dthdp = 100.0_r8*(theta(i,k)-theta(i,k-1))*rpdeli(i,k-1)
                        if (dthdp < dthdpmn(i)) then
                            dthdpmn(i) = dthdp
                            kdthdp(i) = k     ! index of interface of max inversion
                        end if
                    end if
                end do
            end do
            !
            ! Also check between the bottom layer and the surface
            ! Only perform this check if the criteria were not met above
            !
            do i = 1, state%ncol
                if (kdthdp(i) == 0 .and. ocnfrac(i) > 0.01_r8) then
                    dthdp = 100.0_r8*(thetas(i)-theta(i,pver))/(ps(i)-state%pmid(i,pver))
                    if (dthdp < dthdpmn(i)) then
                        dthdpmn(i) = dthdp
                        kdthdp(i) = pver
                    end if
                end if
            end do

            do i = 1, state%ncol
                if (kdthdp(i) /= 0) then
                    k = kdthdp(i)
                    kp1 = min(k+1, pver)
                    ! Note: "strat" will be zero unless ocnfrac > 0.01
                    strat = min(1.0_r8, max(0.0_r8, 2.0*ocnfrac(i)*((theta(i,k700)-thetas(i))*0.057_r8-0.5573_r8)))
                    !
                    ! assign the stratus to the layer just below max inversion
                    ! the relative humidity changes so rapidly across the inversion
                    ! that it is not safe to just look immediately below the inversion
                    ! so limit the stratus cloud by rh in both layers below the inversion
                    !
                    cldst(i,k) = min(strat, max(rh(i,k), rh(i,kp1)))
                end if
            end do
        end if

        do k = 1, pver
            do i = 1, state%ncol
                !
                ! which is greater; standard layered cloud amount or stratocumulus diagnosis
                !
                cloud(i,k) = max(cloud(i,k), cldst(i,k))
                !
                ! add in the contributions of convective cloud (determined separately and accounted
                ! for by modifications to the large-scale relative humidity.
                !
                cloud(i,k) = min(cloud(i,k)+concld(i,k), 1.0_r8)
            end do
        end do

        return
    end subroutine cldfrc_MG

end module cloud_fraction_MG
