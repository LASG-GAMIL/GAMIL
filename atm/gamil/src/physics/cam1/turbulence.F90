#include <misc.h>
#include <params.h>
module turbulence
!---------------------------------------------------------------------------------
! Module to compute mixing coefficients associated with turbulence in the 
! planetary boundary layer and elsewhere.
!
! calling sequence:
!
!    trbinti       initializes time independent coefficients
!     .
!     .
!    vdiff_eddy
!      trbintr        interface for vertical diffusion and pbl scheme and write output
!          trbintd       initializes time dependent variables
!          pblintd       initializes time dependent variables that depend pbl depth
!          vd_k_freeatm  computes mixing coefficients for free atmosphere
!          vd_k_pbl      computes mixing coefficients for pbl
!
!---------------------------Code history--------------------------------
! Standardized:      J. Rosinski, June 1992
! Reviewed:          P. Rasch, B. Boville, August 1992
! Reviewed:          P. Rasch, April 1996
! Reviewed:          B. Boville, April 1996
! rewritten:         B. Boville, May 2000
! rewritten:         B. Stevens, August 2000
!---------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only    : pver, pverp, pcols
  use pmgrid, only    : masterproc
  use tracers, only   : pcnst, pnats
  use history, only: outfld

  implicit none
!
! PBL limits
!
  real(r8), parameter :: ustar_min = 0.01        ! min permitted value of ustar
  real(r8), parameter :: pblmaxp   = 4.e4        ! pbl max depth in pressure units
  real(r8), parameter :: zkmin     = 0.01        ! Minimum kneutral*f(ri)
!
! PBL Parameters
!
  real(r8), parameter :: onet  = 1./3. ! 1/3 power in wind gradient expression
  real(r8), parameter :: betam = 15.0  ! Constant in wind gradient expression
  real(r8), parameter :: betas =  5.0  ! Constant in surface layer gradient expression
  real(r8), parameter :: betah = 15.0  ! Constant in temperature gradient expression 
  real(r8), parameter :: fakn  =  7.2  ! Constant in turbulent prandtl number
  real(r8), parameter :: fak   =  8.5  ! Constant in surface temperature excess         
  real(r8), parameter :: ricr  =  0.3  ! Critical richardson number
  real(r8), parameter :: sffrac=  0.1  ! Surface layer fraction of boundary layer
  real(r8), parameter :: binm  = betam*sffrac       ! betam * sffrac
  real(r8), parameter :: binh  = betah*sffrac       ! betah * sffrac
!
! Pbl constants set using values from other parts of code
!
  real(r8), save :: cpair      ! Specific heat of dry air
  real(r8), save :: rair       ! Gas const for dry air
  real(r8), save :: zvir       ! rh2o/rair - 1
  real(r8), save :: g          ! Gravitational acceleration
  real(r8), save :: ml2(pverp) ! Mixing lengths squared
  real(r8), save :: vk         ! Von Karman's constant
  real(r8), save :: ccon       ! fak * sffrac * vk

  integer, save  :: npbl       ! Maximum number of levels in pbl from surface
  integer, save  :: ntop_turb  ! Top level to which turbulent vertical diffusion is applied.
  integer, save  :: nbot_turb  ! Bottom level to which turbulent vertical diff is applied.

CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
 subroutine trbinti(gravx,      cpairx, rairx, zvirx, ntop_eddy,  &
                    nbot_eddy,  hypm,   vkx )
!----------------------------------------------------------------------- 
! 
! Purpose:  
! Initialize time independent variables of turbulence/pbl package.
! 
! Author: B. Boville, B. Stevens (August 2000)
! 
!-----------------------------------------------------------------------
    implicit none
   !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: gravx     ! acceleration of gravity
    real(r8), intent(in) :: cpairx    ! specific heat of dry air
    real(r8), intent(in) :: rairx     ! gas constant for dry air
    real(r8), intent(in) :: zvirx     ! rh2o/rair - 1
    real(r8), intent(in) :: hypm(pver)! reference pressures at midpoints
    real(r8), intent(in) :: vkx       ! Von Karman's constant
    integer, intent(in)  :: ntop_eddy ! Top level to which eddy vert diff is applied.
    integer, intent(in)  :: nbot_eddy ! Bottom level to which eddy vert diff is applied.
   !---------------------------Local workspace-----------------------------
    integer :: k                     ! vertical loop index
   !-----------------------------------------------------------------------
   !
   ! Basic constants
   !
    cpair = cpairx
    rair  = rairx
    g     = gravx
    zvir  = zvirx
    vk    = vkx
    ccon  = fak*sffrac*vk
    ntop_turb = ntop_eddy
    nbot_turb = nbot_eddy
   !
   ! Set the square of the mixing lengths.
   !
    ml2(ntop_turb) = 0.
    do k = ntop_turb+1, nbot_turb
       ml2(k) = 30.0**2
    end do
    ml2(nbot_turb+1) = 0.
   !
   ! Limit pbl height to regions below 400 mb
   ! npbl = max number of levels (from bottom) in pbl
   !
    npbl = 0
    do k=nbot_turb,ntop_turb,-1
       if (hypm(k) >= pblmaxp) then
          npbl = npbl + 1
       end if
    end do
    npbl = max(npbl,1)

    if (masterproc) then
       write(6,*)'TRBINTI: PBL height will be limited to bottom ',npbl, &
            ' model levels. Top is ',hypm(pverp-npbl),' pascals'
    end if

    return
 end subroutine trbinti
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
 subroutine trbintr(lchnk   ,ncol    ,                            &
                    th      ,t       ,q       ,z       ,zi      , &
                    pmid    ,u       ,v       ,taux    ,tauy    , &
                    shflx   ,cflx    ,obklen  ,ustar   ,pblh    , &
                    kvm     ,kvh     ,cgh     ,cgs     ,kqfs    , &
                    tpert   ,qpert   ,ktopbl  ,ktopblmn)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Interface routines for calcualtion and diatnostics of turbulence related
!  coefficients
!
! Author: B. Stevens (rewrite August 2000)
! 
!-----------------------------------------------------------------------
    implicit none
   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
    integer, intent(in) :: lchnk                      ! chunk identifier
    integer, intent(in) :: ncol                       ! number of atmospheric columns

    real(r8), intent(in)  :: th(pcols,pver)           ! potential temperature [K]
    real(r8), intent(in)  :: t(pcols,pver)            ! temperature (used for density)
    real(r8), intent(in)  :: q(pcols,pver,pcnst+pnats)! specific humidity [kg/kg]
    real(r8), intent(in)  :: z(pcols,pver)            ! height above surface [m]
    real(r8), intent(in)  :: zi(pcols,pverp)          ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)            ! zonal velocity
    real(r8), intent(in)  :: v(pcols,pver)            ! meridional velocity
    real(r8), intent(in)  :: taux(pcols)              ! zonal stress [N/m2]
    real(r8), intent(in)  :: tauy(pcols)              ! meridional stress [N/m2]
    real(r8), intent(in)  :: shflx(pcols)             ! sensible heat flux
    real(r8), intent(in)  :: cflx(pcols,pcnst+pnats)  ! constituent flux
    real(r8), intent(in)  :: pmid(pcols,pver)         ! midpoint pressures

   !
   ! Output arguments
   !
    real(r8), intent(out) :: kqfs(pcols,pcnst+pnats)  ! kinematic surf constituent flux (kg/m2/s)
    real(r8), intent(out) :: kvm(pcols,pverp)         ! eddy diffusivity for momentum [m2/s]
    real(r8), intent(out) :: kvh(pcols,pverp)         ! eddy diffusivity for heat [m2/s]
    real(r8), intent(out) :: cgh(pcols,pverp)         ! counter-gradient term for heat [J/kg/m]
    real(r8), intent(out) :: cgs(pcols,pverp)         ! counter-gradient star (cg/flux)
    real(r8), intent(out) :: tpert(pcols)             ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)             ! convective humidity excess
    real(r8), intent(out) :: ustar(pcols)             ! surface friction velocity [m/s]
    real(r8), intent(out) :: obklen(pcols)            ! Obukhov length
    real(r8), intent(out) :: pblh(pcols)              ! boundary-layer height [m]

    integer, intent(out)  :: ktopbl(pcols)            ! index of first midpoint inside pbl
    integer, intent(out)  :: ktopblmn                 ! min value of ktopbl
   !
   !---------------------------Local workspace-----------------------------
   !
    real(r8) :: wstar(pcols)            ! convective velocity scale [m/s]
    real(r8) :: khfs(pcols)             ! kinimatic surface heat flux 
    real(r8) :: kbfs(pcols)             ! surface buoyancy flux 
    real(r8) :: kvf(pcols,pverp)        ! free atmospheric eddy diffsvty [m2/s]
    real(r8) :: s2(pcols,pver)          ! shear squared
    real(r8) :: n2(pcols,pver)          ! brunt vaisaila frequency
    real(r8) :: ri(pcols,pver)          ! richardson number: n2/s2
   !
   ! Initialize time dependent variables that do not depend on pbl height
   !
    call trbintd(lchnk   ,ncol    ,                            &
                 th      ,q       ,z       ,u       ,v       , &
                 t       ,pmid    ,cflx    ,shflx   ,taux    , &
                 tauy    ,ustar   ,obklen  ,kqfs    ,khfs    , &
                 kbfs    ,s2      ,n2      ,ri      )
   !
   ! Initialize time dependent variables that do depend on pbl height
   !
    call  pblintd(lchnk   ,ncol    ,                            &
                  th      ,q       ,z       ,u       ,v       , &
                  ustar   ,obklen  ,kbfs    ,pblh    ,wstar   )
   !
   ! Get free atmosphere exchange coefficients
   !
    call austausch_atm(lchnk   ,ncol    ,ri      ,s2      ,kvf     )
   ! 
   ! Get pbl exchange coefficients
   !
    call austausch_pbl(lchnk   ,ncol    ,                            &
                       z       ,kvf     ,kqfs    ,khfs    ,kbfs    , &
                       obklen  ,ustar   ,wstar   ,pblh    ,kvm     , &
                       kvh     ,cgh     ,cgs     ,tpert   ,qpert   , &
                       ktopbl  ,ktopblmn)
   !
    call outfld ('PBLH    ',pblh ,pcols,lchnk)
    call outfld ('TPERT   ',tpert,pcols,lchnk)
    call outfld ('QPERT   ',qpert,pcols,lchnk)
    call outfld ('USTAR   ',ustar,pcols,lchnk)
    call outfld ('KVH     ',kvh,pcols,lchnk)
    call outfld ('KVM     ',kvm,pcols,lchnk)
    call outfld ('CGS     ',cgs,pcols,lchnk)

    return
 end subroutine trbintr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
 subroutine trbintd(lchnk   ,ncol    ,                            &
                    th      ,q       ,z       ,u       ,v       , &
                    t       ,pmid    ,cflx    ,shflx   ,taux    , &
                    tauy    ,ustar   ,obklen  ,kqfs    ,khfs    , &
                    kbfs    ,s2      ,n2      ,ri      )
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Time dependent initialization
! 
! Method: 
!  Diagnosis of variables that do not depend on mixing assumptions or
!  PBL depth
!
! Author: B. Stevens (extracted from pbldiff, August, 2000)
! 
!-----------------------------------------------------------------------
    implicit none
   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
    integer, intent(in) :: lchnk                     ! chunk identifier
    integer, intent(in) :: ncol                      ! number of atmospheric columns

    real(r8), intent(in)  :: th(pcols,pver)          ! potential temperature [K]
    real(r8), intent(in)  :: q(pcols,pver,pcnst+pnats)     ! specific humidity [kg/kg]
    real(r8), intent(in)  :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)           ! windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)           ! windspeed y-direction [m/s]
    real(r8), intent(in)  :: t(pcols,pver)           ! temperature (used for density)
    real(r8), intent(in)  :: pmid(pcols,pver)        ! midpoint pressures
    real(r8), intent(in)  :: cflx(pcols,pcnst+pnats) ! surface constituent flux (kg/m2/s)
    real(r8), intent(in)  :: shflx(pcols)            ! surface heat flux (W/m2)
    real(r8), intent(in)  :: taux(pcols)             ! surface u stress [N/m2]
    real(r8), intent(in)  :: tauy(pcols)             ! surface v stress [N/m2]

   !
   ! Output arguments
   !
    real(r8), intent(out) :: ustar(pcols)            ! surface friction velocity [m/s]
    real(r8), intent(out) :: obklen(pcols)           ! Obukhov length
    real(r8), intent(out) :: khfs(pcols)             ! surface kinematic heat flux [mK/s]
    real(r8), intent(out) :: kbfs(pcols)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(r8), intent(out) :: kqfs(pcols,pcnst+pnats) ! sfc kinematic constituent flux [m/s]
    real(r8), intent(out) :: s2(pcols,pver)          ! shear squared
    real(r8), intent(out) :: n2(pcols,pver)          ! brunt vaisaila frequency
    real(r8), intent(out) :: ri(pcols,pver)          ! richardson number: n2/s2
   !
   !---------------------------Local workspace-----------------------------
   !
    integer  :: i                        ! longitude index
    integer  :: k                        ! level index
    integer  :: m                        ! constituent index

    real(r8) :: thvsrf(pcols)            ! sfc (bottom) level virtual temperature
    real(r8) :: thv(pcols,pver)          ! bulk Richardson no. from level to ref lev
    real(r8) :: rrho(pcols)              ! 1./bottom level density (temporary)
    real(r8) :: vvk                      ! velocity magnitude squared
    real(r8) :: dvdz2                    ! velocity shear squared
    real(r8) :: dz                       ! delta z between midpoints
   !
   ! Compute ustar, and kinematic surface fluxes from surface energy fluxes
   !
    do i=1,ncol
       rrho(i)   = rair*t(i,pver)/pmid(i,pver)
       ustar(i)  = max(sqrt(sqrt(taux(i)**2 + tauy(i)**2)*rrho(i)),ustar_min)
       khfs(i)   = shflx(i)*rrho(i)/cpair
    end do
    do m=1,pcnst+pnats
       do i=1,ncol
          kqfs(i,m)= cflx(i,m)*rrho(i)
       end do
    end do
   !
   ! Compute Obukhov length virtual temperature flux and various arrays for use later:
   !
    do i=1,ncol
       kbfs(i)      = khfs(i) + 0.61*th(i,pver)*kqfs(i,1)
       thvsrf(i)    = th(i,pver)*(1.0 + 0.61*q(i,pver,1))
       obklen(i)    = -thvsrf(i)*ustar(i)**3/(g*vk*(kbfs(i) + sign(1.e-10_r8,kbfs(i))))
    end do
   !
   ! Compute shear squared (s2), brunt vaisaila frequency (n2) and related richardson
   ! number (ri).  For the n2 calcualtion use the dry theta_v derived from virtem
   !
    call virtem(ncol, pcols, pver, th      ,q(1,1,1),0.61_r8 ,thv)
    do k=ntop_turb,nbot_turb-1
       do i=1,ncol
          dvdz2   = (u(i,k)-u(i,k+1))**2 + (v(i,k)-v(i,k+1))**2
          dvdz2   = max(dvdz2,1.e-36_r8)
          dz      = z(i,k) - z(i,k+1)
          s2(i,k) = dvdz2/(dz**2)
          n2(i,k) = g*2.0*( thv(i,k) - thv(i,k+1))/((thv(i,k) + thv(i,k+1))*dz)
          ri(i,k) = n2(i,k)/s2(i,k)
       end do
    end do

    return
 end subroutine trbintd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
 subroutine pblintd(lchnk   ,ncol    ,                            &
                    th      ,q       ,z       ,u       ,v       , &
                    ustar   ,obklen  ,kbfs    ,pblh    ,wstar   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Diagnose standard PBL variables
! 
! Method: 
! Diagnosis of PBL depth and related variables.  In this case only wstar.
! The PBL depth follows:
!    Holtslag, A.A.M., and B.A. Boville, 1993:
!    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
!    Model. J. Clim., vol. 6., p. 1825--1842.
!
! Updated by Holtslag and Hack to exclude the surface layer from the
! definition of the boundary layer Richardson number. Ri is now defined
! across the outer layer of the pbl (between the top of the surface
! layer and the pbl top) instead of the full pbl (between the surface and
! the pbl top). For simiplicity, the surface layer is assumed to be the
! region below the first model level (otherwise the boundary layer depth
! determination would require iteration).
!
! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
! >>>>>>>>>  (Use ricr = 0.3 in this formulation)
! 
! Author: B. Stevens (extracted from pbldiff, August 2000)
! 
!-----------------------------------------------------------------------
    implicit none
   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
    integer, intent(in) :: lchnk                     ! chunk identifier
    integer, intent(in) :: ncol                      ! number of atmospheric columns

    real(r8), intent(in)  :: th(pcols,pver)          ! potential temperature [K]
    real(r8), intent(in)  :: q(pcols,pver,pcnst+pnats)     ! specific humidity [kg/kg]
    real(r8), intent(in)  :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in)  :: u(pcols,pver)           ! windspeed x-direction [m/s]
    real(r8), intent(in)  :: v(pcols,pver)           ! windspeed y-direction [m/s]
    real(r8), intent(in)  :: ustar(pcols)            ! surface friction velocity [m/s]
    real(r8), intent(in)  :: obklen(pcols)           ! Obukhov length
    real(r8), intent(in)  :: kbfs(pcols)             ! sfc kinematic buoyancy flux [m^2/s^3]
 
   !
   ! Output arguments
   !
    real(r8), intent(out) :: wstar(pcols)            ! convective sclae velocity [m/s]
    real(r8), intent(out) :: pblh(pcols)             ! boundary-layer height [m]
   !
   !---------------------------Local parameters----------------------------
   !
    real(r8), parameter   :: tiny = 1.e-36           ! lower bound for wind magnitude
    real(r8), parameter   :: fac  = 100.             ! ustar parameter in height diagnosis 
   !
   !---------------------------Local workspace-----------------------------
   !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index
  
    real(r8) :: thvref(pcols)           ! reference level virtual temperature
    real(r8) :: phiminv(pcols)          ! inverse phi function for momentum
    real(r8) :: phihinv(pcols)          ! inverse phi function for heat
    real(r8) :: rino(pcols,pver)        ! bulk Richardson no. from level to ref lev
    real(r8) :: tlv(pcols)              ! ref. level pot tmp + tmp excess
    real(r8) :: tkv                     ! model level potential temperature
    real(r8) :: vvk                     ! velocity magnitude squared
  
    logical  :: unstbl(pcols)           ! pts w/unstbl pbl (positive virtual ht flx)
    logical  :: check(pcols)            ! True=>chk if Richardson no.>critcal
   !
   ! Compute Obukhov length virtual temperature flux and various arrays for use later:
   !
    do i=1,ncol
       check(i)     = .true.
       rino(i,pver) = 0.0
       thvref(i)    = th(i,pver)*(1.0 + 0.61*q(i,pver,1))
       pblh(i)      = z(i,pver)
    end do
   !
   !
   ! PBL height calculation:  Scan upward until the Richardson number between
   ! the first level and the current level exceeds the "critical" value.
   !
    do k=pver-1,pver-npbl+1,-1
       do i=1,ncol
          if (check(i)) then
             vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
             vvk = max(vvk,tiny)
             tkv = th(i,k)*(1. + 0.61*q(i,k,1))
             rino(i,k) = g*(tkv - thvref(i))*(z(i,k)-z(i,pver))/(thvref(i)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1)) * &
                          (z(i,k) - z(i,k+1))
                check(i) = .false.
             end if
          end if
       end do
    end do
   !
   ! Estimate an effective surface temperature to account for surface fluctuations
   !
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       unstbl(i) = (kbfs(i) > 0.)
       check(i)  = (kbfs(i) > 0.)
       if (check(i)) then
          phiminv(i)   = (1. - binm*pblh(i)/obklen(i))**onet
          rino(i,pver) = 0.0
          tlv(i)       = thvref(i) + kbfs(i)*fak/( ustar(i)*phiminv(i) )
       end if
    end do
   !
   ! Improve pblh estimate for unstable conditions using the convective temperature excess:
   !
    do k=pver-1,pver-npbl+1,-1
       do i=1,ncol
          if (check(i)) then
             vvk = (u(i,k) - u(i,pver))**2 + (v(i,k) - v(i,pver))**2 + fac*ustar(i)**2
             vvk = max(vvk,tiny)
             tkv = th(i,k)*(1. + 0.61*q(i,k,1))
             rino(i,k) = g*(tkv - tlv(i))*(z(i,k)-z(i,pver))/(thvref(i)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) - rino(i,k+1))* &
                          (z(i,k) - z(i,k+1))
                check(i) = .false.
              end if
         end if
       end do
    end do
   !
   ! PBL height must be greater than some minimum mechanical mixing depth
   ! Several investigators have proposed minimum mechanical mixing depth
   ! relationships as a function of the local friction velocity, u*.  We
   ! make use of a linear relationship of the form h = c u* where c=700.
   ! The scaling arguments that give rise to this relationship most often
   ! represent the coefficient c as some constant over the local coriolis
   ! parameter.  Here we make use of the experimental results of Koracin
   ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
   ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
   ! latitude value for f so that c = 0.07/f = 700.  Also, do not allow 
   ! PBL to exceed some maximum (npbl) number of allowable points
   !
    do i=1,ncol
       if (check(i)) pblh(i) = z(i,pverp-npbl)
       pblh(i) = max(pblh(i),700.0*ustar(i))
       wstar(i) = (max(0._r8,kbfs(i))*g*pblh(i)/thvref(i))**onet
    end do

    return
 end subroutine pblintd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
  subroutine austausch_atm(lchnk   ,ncol    ,ri      ,s2      ,kvf     )
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  Computes exchange coefficients for free turbulent flows. 
! 
! Method: 
!
! The free atmosphere diffusivities are based on standard mixing length
! forms for the neutral diffusivity multiplied by functns of Richardson
! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for
! momentum, potential temperature, and constitutents. 
!
! The stable Richardson num function (Ri>0) is taken from Holtslag and
! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))
! The unstable Richardson number function (Ri<0) is taken from  CCM1.
! f = sqrt(1 - 18*Ri)
! 
! Author: B. Stevens (rewrite, August 2000)
! 
!-----------------------------------------------------------------------
    implicit none
   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
    integer, intent(in) :: lchnk                    ! chunk identifier
    integer, intent(in) :: ncol                     ! number of atmospheric columns

    real(r8), intent(in)  ::  s2(pcols,pver)        ! shear squared
    real(r8), intent(in)  ::  ri(pcols,pver)        ! richardson no
   !
   ! Output arguments
   !
    real(r8), intent(out) :: kvf(pcols,pverp)       ! coefficient for heat and tracers
   !
   !---------------------------Local workspace-----------------------------
   !
    real(r8) :: fofri                  ! f(ri)
    real(r8) :: kvn                    ! neutral Kv

    integer  :: i                      ! longitude index
    integer  :: k                      ! vertical index
   !
   !-----------------------------------------------------------------------
   !
   ! The surface diffusivity is always zero
   !
    kvf(:ncol,pverp) = 0.0
   !
   ! Set the vertical diffusion coefficient above the top diffusion level
   ! Note that nbot_turb != pver is not supported
   !
    kvf(:ncol,1:ntop_turb) = 0.0
   !
   ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 
   !
    do k = ntop_turb, nbot_turb-1
       do i=1,ncol
          if (ri(i,k) < 0.0) then
             fofri = sqrt(max(1. - 18.*ri(i,k),0._r8))
          else 
             fofri = 1.0/(1.0 + 10.0*ri(i,k)*(1.0 + 8.0*ri(i,k)))    
          end if
          kvn = ml2(k)*sqrt(s2(i,k))
          kvf(i,k+1) = max(zkmin,kvn*fofri)
       end do
    end do

    return
 end subroutine austausch_atm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================================
 subroutine austausch_pbl(lchnk   ,ncol    ,                            &
                          z       ,kvf     ,kqfs    ,khfs    ,kbfs    , &
                          obklen  ,ustar   ,wstar   ,pblh    ,kvm     , &
                          kvh     ,cgh     ,cgs     ,tpert   ,qpert   , &
                          ktopbl  ,ktopblmn)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Atmospheric Boundary Layer Computation
! 
! Method: 
! Nonlocal scheme that determines eddy diffusivities based on a
! specified boundary layer height and a turbulent velocity scale;
! also, countergradient effects for heat and moisture, and constituents
! are included, along with temperature and humidity perturbations which
! measure the strength of convective thermals in the lower part of the
! atmospheric boundary layer.
!
! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
! Model. J. Clim., vol. 6., p. 1825--1842.
!
! Updated by Holtslag and Hack to exclude the surface layer from the
! definition of the boundary layer Richardson number. Ri is now defined
! across the outer layer of the pbl (between the top of the surface
! layer and the pbl top) instead of the full pbl (between the surface and
! the pbl top). For simiplicity, the surface layer is assumed to be the
! region below the first model level (otherwise the boundary layer depth
! determination would require iteration).
!
! Author: B. Boville, B. Stevens (rewrite August 2000)
! 
!-----------------------------------------------------------------------
    implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: lchnk                    ! chunk identifier
    integer, intent(in) :: ncol                     ! number of atmospheric columns

    real(r8), intent(in) :: z(pcols,pver)           ! height above surface [m]
    real(r8), intent(in) :: kvf(pcols,pverp)        ! free atmospheric eddy diffsvty [m2/s]
    real(r8), intent(in) :: kqfs(pcols,pcnst+pnats) ! kinematic surf cnstituent flux (kg/m2/s)
    real(r8), intent(in) :: khfs(pcols)             ! kinimatic surface heat flux 
    real(r8), intent(in) :: kbfs(pcols)             ! surface buoyancy flux 
    real(r8), intent(in) :: pblh(pcols)             ! boundary-layer height [m]
    real(r8), intent(in) :: obklen(pcols)           ! Obukhov length
    real(r8), intent(in) :: ustar(pcols)            ! surface friction velocity [m/s]
    real(r8), intent(in) :: wstar(pcols)            ! convective velocity scale [m/s]
!
! Output arguments
!
    real(r8), intent(out) :: kvm(pcols,pverp)       ! eddy diffusivity for momentum [m2/s]
    real(r8), intent(out) :: kvh(pcols,pverp)       ! eddy diffusivity for heat [m2/s]
    real(r8), intent(out) :: cgh(pcols,pverp)       ! counter-gradient term for heat [J/kg/m]
    real(r8), intent(out) :: cgs(pcols,pverp)       ! counter-gradient star (cg/flux)
    real(r8), intent(out) :: tpert(pcols)           ! convective temperature excess
    real(r8), intent(out) :: qpert(pcols)           ! convective humidity excess

    integer, intent(out) :: ktopbl(pcols)           ! index of first midpoint inside pbl
    integer, intent(out) :: ktopblmn                ! min value of ktopbl
!
!---------------------------Local workspace-----------------------------
!
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index

    real(r8) :: phiminv(pcols)          ! inverse phi function for momentum
    real(r8) :: phihinv(pcols)          ! inverse phi function for heat
    real(r8) :: wm(pcols)               ! turbulent velocity scale for momentum
    real(r8) :: zp(pcols)               ! current level height + one level up 
    real(r8) :: fak1(pcols)             ! k*ustar*pblh     
    real(r8) :: fak2(pcols)             ! k*wm*pblh
    real(r8) :: fak3(pcols)             ! fakn*wstar/wm
    real(r8) :: pblk(pcols)             ! level eddy diffusivity for momentum
    real(r8) :: pr(pcols)               ! Prandtl number for eddy diffusivities
    real(r8) :: zl(pcols)               ! zmzp / Obukhov length
    real(r8) :: zh(pcols)               ! zmzp / pblh
    real(r8) :: zzh(pcols)              ! (1-(zmzp/pblh))**2
    real(r8) :: zmzp                    ! level height halfway between zm and zp
    real(r8) :: term                    ! intermediate calculation

    logical  :: unstbl(pcols)           ! pts w/unstbl pbl (positive virtual ht flx)
    logical  :: pblpt(pcols)            ! pts within pbl
!
! Initialize height independent arrays
!
    do i=1,ncol
       unstbl(i) = (kbfs(i) > 0.)
       pblk(i) = 0.0
       fak1(i) = ustar(i)*pblh(i)*vk
       if (unstbl(i)) then
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1. - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          fak3(i)    = fakn*wstar(i)/wm(i)
          tpert(i)   = max(khfs(i)*fak/wm(i),0._r8)
          qpert(i)   = max(kqfs(i,1)*fak/wm(i),0._r8)
       else
          tpert(i)   = max(khfs(i)*fak/ustar(i),0._r8)
          qpert(i)   = max(kqfs(i,1)*fak/ustar(i),0._r8)
       end if
    end do
!
! Initialize output arrays with free atmosphere values
!
    do k=1,pverp
       do i=1,ncol
          kvm(i,k) = kvf(i,k)
          kvh(i,k) = kvf(i,k)
          cgh(i,k) = 0.0
          cgs(i,k) = 0.0
       end do
    end do
!
! Main level loop to compute the diffusivities and counter-gradient terms. These terms are 
! only calculated at points determined to be in the interior of the pbl (pblpt(i)==.true.),
! and then calculations are directed toward regime: stable vs unstable, surface vs outer 
! layer.
!
    do k=pver,pver-npbl+2,-1
       do i=1,ncol
          pblpt(i) = (z(i,k) < pblh(i))
          if (pblpt(i)) then
             ktopbl(i) = k
             zp(i)  = z(i,k-1)
             if (zkmin == 0.0 .and. zp(i) > pblh(i)) zp(i) = pblh(i)
             zmzp    = 0.5*(z(i,k) + zp(i))
             zh(i)   = zmzp/pblh(i)
             zl(i)   = zmzp/obklen(i)
             zzh(i)  = zh(i)*max(0._r8,(1. - zh(i)))**2
             if (unstbl(i)) then
                if (zh(i) < sffrac) then
                   term     = (1. - betam*zl(i))**onet
                   pblk(i)  = fak1(i)*zzh(i)*term
                   pr(i)    = term/sqrt(1. - betah*zl(i))
                else
                   pblk(i)  = fak2(i)*zzh(i)
                   pr(i)    = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
                   cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
                   cgh(i,k) = khfs(i)*cgs(i,k)*cpair
                end if
             else
                if (zl(i) <= 1.) then
                   pblk(i) = fak1(i)*zzh(i)/(1. + betas*zl(i))
                else
                   pblk(i) = fak1(i)*zzh(i)/(betas + zl(i))
                end if
                pr(i)    = 1.
             end if
             kvm(i,k) = max(pblk(i),kvf(i,k))
             kvh(i,k) = max(pblk(i)/pr(i),kvf(i,k))
          end if
       end do
    end do
!
! Check whether last allowed midpoint is within pbl, determine ktopblmn
!
    ktopblmn = pver
    k = pver-npbl+1
    do i = 1, ncol
       if (z(i,k) < pblh(i)) ktopbl(i) = k
       ktopblmn = min(ktopblmn, ktopbl(i))
    end do

    return
 end subroutine austausch_pbl
END MODULE turbulence
