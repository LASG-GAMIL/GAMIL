#include <misc.h>
#include <params.h>

subroutine grmult (ztodt   ,rcoslat ,coslat  ,coriol  ,div     , &
                   q3      ,t3      ,u3      ,v3      ,u3sld   , &
                   v3sld   ,ed1     ,ql      ,qm      ,tl      , &
                   tm      ,phis    ,phisl   ,phism   ,dpsl    , &
                   dpsm    ,omga    ,pmid    ,pdel    ,ps      , &
                   logps   ,rpmid   ,etamid  ,fu      ,fv      , &
                   t2      ,t3m1    ,u3m1    ,v3m1    ,etadot  , &
                   dpslon  ,dpslat  ,urhs    ,vrhs    ,trhs    , &
                   prhs    ,etadotm1,lnpssld ,prhssld ,tarrsld , &
                   parrsld ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Compute non-linear dynamics terms in grid point space (in preparation
! for SLD interpolation)
!
!---------------------------Code history--------------------------------
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: grmult.F90,v 1.6.2.2 2002/06/15 13:48:24 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comslt
  use commap
  use physconst, only: rair, zvir, cappa, cpvir
  use time_manager, only: is_first_step

  implicit none

#include <comhyb.h>
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                 ! delta-t
  real(r8), intent(in)   :: rcoslat               ! 1/(cos lat)
  real(r8), intent(in)   :: coslat                ! cos(lat)
  real(r8), intent(in)   :: coriol                ! Coriolis parameter
  real(r8), intent(in)   :: div     (plond,plev)  ! divergence
  real(r8), intent(in)   :: q3      (plond,plev)  ! specific humidity
  real(r8), intent(in)   :: t3      (plond,plev)  ! temperature
  real(r8), intent(in)   :: u3      (plond,plev)  ! zonal wind
  real(r8), intent(in)   :: v3      (plond,plev)  ! meridional wind
  real(r8), intent(out)  :: u3sld   (plond,plev)  ! u-wind (time n) (used for advection)
  real(r8), intent(out)  :: v3sld   (plond,plev)  ! v-wind (time n) (used for advection)
  real(r8), intent(in)   :: ed1     (plond,plev)  ! (1/ps)etadot(dp/deta) (time n-1)
  real(r8), intent(in)   :: ql      (plond,plev)  ! longitudinal derivative of q
  real(r8), intent(in)   :: qm      (plond,plev)  ! latitudinal  derivative of q
  real(r8), intent(in)   :: tl      (plond,plev)  ! zonal derivative of T
  real(r8), intent(in)   :: tm      (plond,plev)  ! meridional derivative of T
  real(r8), intent(in)   :: phis    (plond)       ! Phi at surface
  real(r8), intent(in)   :: phisl   (plond)       ! longitudinal derivative of phis
  real(r8), intent(in)   :: phism   (plond)       ! latitudinal  devivative of phis
  real(r8), intent(in)   :: dpsl    (plond)       ! longitudinal component of grad ln(ps)
  real(r8), intent(in)   :: dpsm    (plond)       ! latitudinal  component of grad ln(ps)
  real(r8), intent(in)   :: omga    (plond,plev)  ! vertical pressure velocity
  real(r8), intent(in)   :: pmid    (plond,plev)  ! pressure at full levels
  real(r8), intent(in)   :: pdel    (plond,plev)  ! layer thicknesses (pressure)
  real(r8), intent(in)   :: ps      (plond)       ! Surface pressure (n)
  real(r8), intent(in)   :: logps   (plond)       ! log(ps)
  real(r8), intent(in)   :: rpmid   (plond,plev)  ! 1./pmid
  real(r8), intent(in)   :: etamid  (plev)        ! midpoint values of eta (a+b)
  real(r8), intent(inout):: fu      (plond,plev)  ! nonlinear term - u momentum eqn
  real(r8), intent(inout):: fv      (plond,plev)  ! nonlinear term - v momentum eqn
  real(r8), intent(inout):: t2      (plond,plev)  ! nonlinear term - temperature
  real(r8), intent(out)  :: t3m1    (plond,plev)  ! T at time n-1
  real(r8), intent(inout):: u3m1    (plond,plev)  ! U at previous time step
  real(r8), intent(inout):: v3m1    (plond,plev)  ! V at previous time step
  real(r8), intent(out)  :: etadot  (plond,plevp) ! vertical velocity in eta coordinates
  real(r8), intent(out)  :: dpslon  (plond,plev)  ! Pressure gradient term
  real(r8), intent(out)  :: dpslat  (plond,plev)  ! Pressure gradient term
  real(r8), intent(inout):: urhs    (plond,plev)  ! RHS of U  eqn valid for mid-point
  real(r8), intent(inout):: vrhs    (plond,plev)  ! RHS of V  eqn valid for mid-point
  real(r8), intent(inout):: trhs    (plond,plev)  ! RHS of T  eqn valid for mid-point
  real(r8), intent(inout):: prhs    (plond,plev)  ! RHS of Ps eqn valid for mid-point
  real(r8), intent(inout):: etadotm1(plond,plevp) ! etadot for time n-1
  real(r8), intent(out)  :: lnpssld (plond,plev)  ! RHS Ps term for SLD
  real(r8), intent(out)  :: prhssld (plond,plev)  ! RHS Ps term for SLD
  real(r8), intent(out)  :: tarrsld (plond,plev)  ! T  at arr. pt. (SLD)
  real(r8), intent(out)  :: parrsld (plond,plev)  ! Ps at arr. pt. (SLD)
  integer , intent(in)   :: nlon                  ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) tmp1                  ! temporary workspace
  real(r8) tmp2                  ! temporary workspace
  real(r8) tmp                   ! temporary workspace
  real(r8) tmpk                  ! workspace
  real(r8) tmpkp1                ! workspace
  real(r8) u3l     (plond,plev)  ! u-wind used locally only
  real(r8) v3l     (plond,plev)  ! v-wind used locally only
  real(r8) tv      (plond,plev)  ! virtual temperature
  real(r8) ddpk    (plond)       ! partial sum of div*delta p
  real(r8) ddpn    (plond)       ! complete sum of div*delta p
  real(r8) vkdp    (plond,plev)  ! V dot grad(ln(ps))
  real(r8) vpdsk   (plond)       ! partial sum  V dot grad(ln(ps)) delta b
  real(r8) vpdsn   (plond)       ! complete sum V dot grad(ln(ps)) delta b
  real(r8) lpsstar (plond)       ! Reference ln(Ps) (used to define a new 
!                                ! perturbation Ps)
  real(r8) lpsstarl(plond)       ! long. grad of reference ln(Ps)
  real(r8) lpsstarm(plond)       ! lat.  grad of reference ln(Ps)
  real(r8) rtv     (plond,plev)  ! rair*(tv+t0)
  real(r8) dt                    ! time step
  real(r8) hsl     (plond,plev)  ! zonal      deriv of hydrostatic term
  real(r8) hsm     (plond,plev)  ! meridional deriv of hydrostatic term
  real(r8) pspsl   (plond)       ! Ps*d(lnPs)/d(long.)
  real(r8) pspsm   (plond)       ! Ps*d(lnPs)/d(lat. )
  real(r8) ed1p    (plond,plevp) ! (1/ps)etadot(dp/deta) (time n-1)
  real(r8) rtvl                  ! zonal      derivative of R*Tv
  real(r8) rtvm                  ! meridional derivative of R*Tv
  real(r8) abp0                  ! constant for grad(H(n)) matrix
!
! Arrays which hold results from the RHS of the prognostic equations.
! Results will be interpolated to trajectory departure points in the routine
! SCANSLT.
!
  real(r8) tsld0a  (plond,plev)  ! RHS of T  eqn valid for mid-point
  real(r8) usldm   (plond,plev)  ! RHS of U  eqn valid for departure pt (n)
  real(r8) vsldm   (plond,plev)  ! RHS of V  eqn valid for departure pt (n)
  real(r8) tsldm   (plond,plev)  ! RHS of T  eqn valid for departure pt (n)
  real(r8) psldm   (plond,plev)  ! RHS of Ps eqn valid for departure pt (n)
  real(r8) urhsl   (plond,plev)  ! RHS of U  eqn valid for mid-point (n+1/2)
  real(r8) vrhsl   (plond,plev)  ! RHS of V  eqn valid for mid-point (n+1/2)
  real(r8) trhsl   (plond,plev)  ! RHS of T  eqn valid for mid-point (n+1/2)
  real(r8) prhsl   (plond,plev)  ! RHS of Ps eqn valid for mid-point (n+1/2)

  real(r8) onemeps               ! 1 - epssld (SLD decentering coefficient)
  real(r8) onepeps               ! 1 + epssld (SLD decentering coefficient)
  real(r8) detai(plevp)          ! interval between interfaces
  real(r8) facm1                 ! interpolation factor for time n-1
  real(r8) facm2                 ! interpolation factor for time n-2

  integer i,l,k                  ! longitude, level indices
  integer npr                    ! index
!
!-----------------------------------------------------------------------
!
  onemeps = 1. - epssld
  onepeps = 1. + epssld
  facm1   =  3./2.
  facm2   = -1./2.
!
  do k = 1,plev
     detai  (k) = (hyai(k+1) + hybi(k+1)) - (hyai(k) + hybi(k))
  end do
!
! Compute U/V for time n + 1/2
! The first formula will be used later for trajectory calculation
! The second will be used locally in the Hortal Temperature correction
!
  do k = 1,plev
     do i = 1,nlon
        u3sld(i,k) = 2.   *u3(i,k) -       u3m1(i,k)
        v3sld(i,k) = 2.   *v3(i,k) -       v3m1(i,k)
        u3l  (i,k) = facm1*u3(i,k) + facm2*u3m1(i,k)
        v3l  (i,k) = facm1*v3(i,k) + facm2*v3m1(i,k)
     end do
  end do
!
! Zero auxiliary fields
!
  tmp = 1./(rair*t0(plev))
  do i=1,nlon
     ddpk    (i) = 0.0
     ddpn    (i) = 0.0
     vpdsk   (i) = 0.0
     vpdsn   (i) = 0.0
     lpsstar (i) = -phis (i)*tmp
     lpsstarl(i) = -phisl(i)*tmp
     lpsstarm(i) = -phism(i)*tmp
     pspsl   (i) = ps(i)*dpsl(i)
     pspsm   (i) = ps(i)*dpsm(i)
     etadot(i,1) = 0.0
     ed1p  (i,1) = 0.0
     etadot(i,plevp) = 0.0
  end do
!
! Virtual temperature
!
  call virtem(nlon, plond, plev, t3      ,q3      ,zvir    ,tv)
!
! calculate some auxiliary quantities
!
  do k=1,plev
     do i=1,nlon
        ed1p(i,k+1) = ed1(i,k)
        rtv(i,k) = rair*tv(i,k)
!
! sum(plev)(div(k)*dp(k))
!
        ddpn(i) = ddpn(i) + div(i,k)*pdel(i,k)
     end do
  end do
!
! sum(plev)(v(k)*grad(lnps)*db(k))
!
  do k=nprlev,plev
     do i=1,nlon
        vkdp(i,k)= rcoslat*(u3(i,k)*pspsl(i) + v3(i,k)*pspsm(i))
        vpdsn(i) = vpdsn(i) + vkdp(i,k)*hybd(k)
     end do
  end do
!
! Compute etadot (top and bottom = 0.)
!
  do k = 1,plev-1
!
! Compute etadot(dp/deta)(k+1/2) and sum(k)(div(j)*dp(j))
!
     do i=1,nlon
        ddpk(i) = ddpk(i) + div(i,k)*pdel(i,k)
#ifdef HADVTEST
!
!jr Set etadot to zero for horizontal advection test
!
        etadot(i,k+1) = 0.
#else
        etadot(i,k+1) = -ddpk(i)
#endif
     end do
!
! sum(k)(v(j)*grad(ps)*db(j))
!
     if (k.ge.nprlev) then
        do i=1,nlon
           vpdsk(i) = vpdsk(i) + vkdp(i,k)*hybd(k)
#ifdef HADVTEST
!
!jr Set etadot to zero for horizontal advection test
!
           etadot(i,k+1) = 0.
#else
           etadot(i,k+1) = etadot(i,k+1) - vpdsk(i) + hybi(k+1)*(ddpn(i)+vpdsn(i))
#endif
        end do
     end if
!
! Convert eta-dot(dp/deta) to eta-dot
!
     tmp = etamid(k+1) - etamid(k)
     do i = 1,nlon
        etadot(i,k+1) = etadot(i,k+1)*tmp/(pmid(i,k+1) - pmid(i,k))
     end do
  end do
!
! Zonal and meridional derivatives of the hydrostatic term in the
! momentum equations.
!
  do k = 1,plev
     do i = 1,nlon
        tmp1     = (1. + zvir*q3(i,k))
        tmp2     = t3(i,k)*zvir
        rtvl     = rair*( tmp1*tl(i,k) + tmp2*ql(i,k) )
        rtvm     = rair*( tmp1*tm(i,k) + tmp2*qm(i,k) )
!
        tmp      = rpmid(i,k)*pdel(i,k)
        hsl(i,k) = -rtvl*tmp
        hsm(i,k) = -rtvm*tmp
     end do
     if(k .ge. nprlev) then
        abp0 = ps0*(hyam(k)*(hybi(k+1) - hybi(k)) - hybm(k)*(hyai(k+1) - hyai(k)))
        do i = 1,nlon
           tmp      = rtv(i,k)*abp0/(pmid(i,k)*pmid(i,k))
           hsl(i,k) = hsl(i,k) - pspsl(i)*tmp
           hsm(i,k) = hsm(i,k) - pspsm(i)*tmp
        end do
     endif
  end do
!
! Calculate RHS of all prognostic eqns
!
  do k = 1,plev
     tmp  = cappa*t0(k)*hypi(plevp)/hypm(k)
     tmp1 = hypd(k)/hypi(plevp)
     do i = 1,nlon
!
! Surface Pressure eqn
!
        prhsl(i,k) =  div(i,k)*(tmp1 - pdel(i,k)/ps(i))
        psldm(i,k) = -div(i,k)*tmp1 - ed1p(i,k+1) + ed1p(i,k)
!
! Temperature eqn
!
        trhsl (i,k) = cappa*tv(i,k)/(1. + cpvir*q3(i,k))*omga(i,k)*rpmid(i,k)
        trhsl (i,k) = trhsl(i,k) - omga(i,k)/ps(i)*tmp
        tsldm (i,k) = 0.5*(ed1p(i,k+1) + ed1p(i,k))*tmp
!
! ... horizontal advection portion of Hortal Temperature correction
!
        trhsl (i,k) = trhsl(i,k) - &
                      rcoslat*( u3(i,k)*lpsstarl(i) + v3(i,k)*lpsstarm(i) )*hortalc(k)
!
! ... Ritchie damping term for Temperature eqn.
!
        tsld0a(i,k) = rcoslat*( u3l(i,k)*lpsstarl(i) + v3l(i,k)*lpsstarm(i) )
!
! U/V eqns (includes only the diagonal portion of the hydrostatic term)
!
        urhsl(i,k) = 0.5*hsl(i,k) + href(k,k)*tl(i,k) + bps(k)*dpsl(i)
        vrhsl(i,k) = 0.5*hsm(i,k) + href(k,k)*tm(i,k) + bps(k)*dpsm(i)
        usldm(i,k) = -phisl(i) + v3(i,k)*coriol*coslat -href(k,k)*tl(i,k) - bps(k)*dpsl(i)
        vsldm(i,k) = -phism(i) - u3(i,k)*coriol*coslat -href(k,k)*tm(i,k) - bps(k)*dpsm(i)
     end do
!
! Add pressure gradient terms to momentum tendencies
!
     if (k.ge.nprlev) then
        do i=1,nlon
           tmp        = rtv(i,k)*hybm(k)*rpmid(i,k)
           dpslon(i,k) = rcoslat*tmp*pspsl(i)
           dpslat(i,k) = rcoslat*tmp*pspsm(i)
           urhsl(i,k) = urhsl(i,k) - dpslon(i,k)*coslat
           vrhsl(i,k) = vrhsl(i,k) - dpslat(i,k)*coslat
        end do
     else
        do i = 1,nlon
           dpslon(i,k) = 0.
           dpslat(i,k) = 0.
        end do
     end if
  end do
!
! Interior levels of the hydrostatic term
!
  do k=1,plev-1
     do l=k+1,plev
        do i=1,nlon
           urhsl(i,k) = urhsl(i,k) + href(l,k)*tl(i,l) + hsl(i,l)
           vrhsl(i,k) = vrhsl(i,k) + href(l,k)*tm(i,l) + hsm(i,l)
!
           usldm(i,k) = usldm(i,k) - href(l,k)*tl(i,l)
           vsldm(i,k) = vsldm(i,k) - href(l,k)*tm(i,l)
        end do
     end do
  end do

  if(is_first_step()) then
     do k = 1,plevp
        do i = 1,nlon
           etadotm1(i,k) = etadot(i,k)
        end do
     end do
  endif
!
! The modified etadotm1 will be used later for trajectory calculation in SCANSLT
!
  do k = 1,plevp
     do i = 1,nlon
        etadotm1(i,k) = 2.*etadot(i,k) - etadotm1(i,k)
     end do
  end do
!
! Compute vertical advection portion of Hortal Temperature correction
!
  npr = nprlev - 1
  if(npr .lt. 1) npr = 1
  do k = npr,plev-1
     tmpk   = 0.5*hdel(k-1)/detai(k)
     tmpkp1 = 0.5*hdel(k  )/detai(k)
     do i = 1,nlon
        trhsl(i,k) = trhsl(i,k) - (etadot(i,k  )*tmpk + etadot(i,k+1)*tmpkp1)*lpsstar(i)
     end do
  end do
!
! ... bottom level
!
  tmpk = 0.5*hdel(plev-1)/detai(plev)
  do i = 1,nlon
     trhsl(i,plev) = trhsl(i,plev) - etadot(i,plev)*tmpk*lpsstar(i)
  end do
!
! Compute (by extrapolation) RHS terms for time n + 1/2
!
  if(is_first_step()) then
     do k = 1,plev
        do i = 1,nlon
           urhs  (i,k) = urhsl (i,k)
           vrhs  (i,k) = vrhsl (i,k)
           trhs  (i,k) = trhsl (i,k)
           prhs  (i,k) = prhsl (i,k)
        end do
     end do
  endif
!
  do k = 1,plev
     do i = 1,nlon
        tmp         =       urhsl (i,k)
        urhsl (i,k) = facm1*urhsl (i,k) + facm2*urhs  (i,k)
        urhs  (i,k) = tmp
        tmp         =       vrhsl (i,k)
        vrhsl (i,k) = facm1*vrhsl (i,k) + facm2*vrhs  (i,k)
        vrhs  (i,k) = tmp
        tmp         =       trhsl (i,k)
        trhsl (i,k) = facm1*trhsl (i,k) + facm2*trhs  (i,k)
        trhs  (i,k) = tmp
        tmp         =       prhsl (i,k)
        prhsl (i,k) = facm1*prhsl (i,k) + facm2*prhs  (i,k)
        prhs  (i,k) = tmp
     end do
  end do
!
! Prepare SLD arrays for interpolation by SCANSLT
!
! Append appropriate coefficients (including decentering epsilon values -- See
! Notes 13,14,15)
!
  dt   = 0.5*ztodt
  tmp1 = 0.5*ztodt/coslat
  do k = 1,plev
     tmp = cappa*t0(k)*hypi(plevp)/hypm(k)
     do i = 1,nlon
        urhsl (i,k) = tmp1*urhsl (i,k)
        vrhsl (i,k) = tmp1*vrhsl (i,k)
        trhsl (i,k) = dt  *trhsl (i,k)
        tsld0a(i,k) = dt  *tsld0a(i,k)
        prhsl (i,k) = dt  *prhsl (i,k)
        t2    (i,k) = ztodt*t2   (i,k)
        fu    (i,k) = ztodt*fu   (i,k)
        fv    (i,k) = ztodt*fv   (i,k)
!
! Combine terms.
! (Time split the midpoint results between arrival and departure
!  points)
!
        tmp2          = lpsstar(i)*hortalc(k)
        u3m1    (i,k) = (u3  (i,k)*coslat + onemeps*usldm(i,k)*dt)/coslat + &
                                            onemeps*urhsl(i,k) + &
                                            onemeps*fu   (i,k)*0.5
        v3m1    (i,k) = (v3  (i,k)*coslat + onemeps*vsldm(i,k)*dt)/coslat + &
                                            onemeps*vrhsl(i,k) + &
                                            onemeps*fv   (i,k)*0.5
#ifdef HADVTEST
        t3m1    (i,k) = t3  (i,k)
#else
        t3m1    (i,k) = t3  (i,k) +         onemeps*tsldm(i,k)*dt + &
                                            onemeps*trhsl(i,k) - tmp2 + &
                                            onemeps*t2(i,k)*0.5
#endif
        prhssld (i,k) = (psldm(i,k)*dt + prhsl(i,k))*onemeps
        tarrsld (i,k) = (trhsl(i,k) + &
                         hybm(k)*tmp*tsld0a(i,k))*onepeps + tmp2 + onepeps*t2(i,k)*0.5
        parrsld (i,k) = (prhsl(i,k) - hybd(k)*tsld0a(i,k))*onepeps
        lnpssld (i,k) = logps(i) - lpsstar(i) - tsld0a(i,k)*onemeps
        fu      (i,k) = urhsl(i,k)*coslat*onepeps + fu(i,k)*coslat*onepeps*0.5
        fv      (i,k) = vrhsl(i,k)*coslat*onepeps + fv(i,k)*coslat*onepeps*0.5
     end do
  end do
#ifdef VADVTEST
  do k=2,plev
     do i=1,nlon
        etadot(i,k) = -0.5/(12.*86400.)
     end do
  end do
#endif
!
  return
end subroutine grmult
