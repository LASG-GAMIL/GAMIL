#include <misc.h>
#include <params.h>

module gw_drag

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an 
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!---------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use physics_types,  only: physics_state, physics_ptend
  use pmgrid, only:         masterproc
  use history, only:        outfld

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public gw_inti                  ! Initialization
  public gw_intr                  ! interface to actual parameterization
!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  integer, parameter :: pgwv = 0  ! number of waves allowed

  integer :: kbotbg, kbotoro      ! interface of gwd source
  integer :: ktopbg, ktoporo      ! top interface of gwd region

  real(r8) :: alpha(0:pver)       ! newtonian cooling coefficients
  real(r8) :: c(-pgwv:pgwv)       ! list of wave phase speeds
  real(r8) :: cpair               ! specific heat of dry air (constant p)
  real(r8) :: dback               ! background diffusivity
  real(r8) :: effkwv              ! effective wavenumber (fcrit2*kwv)
  real(r8) :: effgw               ! tendency efficiency
  real(r8) :: fracldv             ! fraction of stress deposited in low level region
  real(r8) :: g                   ! acceleration of gravity
  real(r8) :: kwv                 ! effective horizontal wave number

  real(r8) :: mxasym              ! max asymmetry between tau(c) and tau(-c)
  real(r8) :: mxrange             ! max range of tau for all c
  real(r8) :: n2min               ! min value of bouyancy frequency
  real(r8) :: fcrit2              ! critical froude number
  real(r8) :: oroko2              ! 1/2 * horizontal wavenumber
  real(r8) :: orohmin             ! min surface displacment height for orographic waves
  real(r8) :: orovmin             ! min wind speed for orographic waves
  real(r8) :: r                   ! gas constant for dry air
  real(r8) :: rog                 ! r / g
  real(r8) :: taubgnd             ! background source strength (/tauscal)
  real(r8) :: taumin              ! minimum (nonzero) stress
  real(r8) :: tauscal             ! scale factor for background stress source
  real(r8) :: tndmax              ! maximum wind tendency
  real(r8) :: umcfac              ! factor to limit tendency to prevent reversing u-c
  real(r8) :: ubmc2mn             ! min (u-c)**2
  real(r8) :: zldvcon             ! constant for determining zldv from tau0

contains

!===============================================================================
  subroutine gw_inti (cpairx, cpwv, gx, rx, hypi)
!-----------------------------------------------------------------------
! Time independent initialization for multiple gravity wave parameterization.
!-----------------------------------------------------------------------
    use history,    only: addfld, add_default, phys_decomp

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: cpairx                ! specific heat of dry air (constant p)
    real(r8), intent(in) :: cpwv                  ! specific heat of water vapor (constant p)
    real(r8), intent(in) :: gx                    ! acceleration of gravity
    real(r8), intent(in) :: rx                    ! gas constant for dry air
    real(r8), intent(in) :: hypi(pver+1)          ! reference interface pressures

!---------------------------Local storage-------------------------------
    integer :: k

!-----------------------------------------------------------------------
! Copy model constants
    cpair  = cpairx
    g      = gx
    r      = rx

! Set MGWD constants
    effgw  = 0.125            ! efficiency of the tendencies
    kwv    = 6.28e-5          ! 100 km wave length
    dback  = 0.05             ! background diffusivity
    fcrit2 = 0.5              ! critical froude number squared
    tauscal= 0.001            ! scale factor for background stress
    taubgnd= 6.4              ! background stress amplitude

    fracldv= 0.0              ! fraction of tau0 diverged in low level region
    zldvcon= 10.              ! constant for determining zldv

! Set phase speeds 
    do k = -pgwv, pgwv
       c(k)   = 10. * k       ! 0, +/- 10, +/- 20, ... m/s
    end do

    if (masterproc) then
       write(6,*) ' '
       write(6,*) 'GW_INTI: pgwv = ', pgwv
       write(6,*) 'GW_INTI: c(l) = ', c
       write(6,*) ' '
    end if

! Set radiative damping times
    do k = 0, pver
       alpha(k) = 1.e-6       ! about 10 days.
    end do

! Min and max values to keep things reasonable
    mxasym = 0.1              ! max factor of 10 from |tau(c)| to |tau(-c)|
    mxrange= 0.001            ! factor of 100 from max to min |tau(c)|
    n2min  = 1.e-8            ! min value of Brunt-Vaisalla freq squared
    orohmin= 10.              ! min surface displacement for orographic wave drag
    orovmin=  2.              ! min wind speed for orographic wave drag
    taumin = 1.e-10           ! min stress considered > 0
    tndmax = 500. / 86400.    ! max permitted tendency (500 m/s/day)
    umcfac = 0.5              ! max permitted reduction in u-c
    ubmc2mn= 0.01             ! min value of (u-c)^2

! Determine other derived constants
    oroko2 = 0.5 * kwv
    effkwv = fcrit2 * kwv
    rog    = r/g

! Determine the bounds of the background and orographic stress regions
    ktopbg  = 0
    kbotoro = pver
    do k = 0, pver
       if (hypi(k+1) .lt. 10000.) kbotbg  = k    ! spectrum source at 100 mb
!!$       if (hypi(k+1) .lt.  3000.) ktoporo = k
    end do
    ktoporo = 0

    if (masterproc) then
       write (6,*) 'KTOPBG  =',ktopbg
       write (6,*) 'KBOTBG  =',kbotbg
       write (6,*) 'KTOPORO =',ktoporo
       write (6,*) 'KBOTORO =',kbotoro
    end if

! Declare history variables for orgraphic term
    call addfld ('TTGWORO ','K/s     ',pver, 'A','T tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('UTGWORO ','m/s2    ',pver, 'A','U tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('VTGWORO ','m/s2    ',pver, 'A','V tendency - orographic gravity wave drag',phys_decomp)
    call addfld ('TAUGWX  ','N/m2    ',1,    'A','Zonal gravity wave surface stress',        phys_decomp)
    call addfld ('TAUGWY  ','N/m2    ',1,    'A','Meridional gravity wave surface stress',   phys_decomp)

!!  (wanhui 2003.06.19)
!
!    call add_default ('UTGWORO ', 1, ' ')
!    call add_default ('VTGWORO ', 1, ' ')
!    call add_default ('TAUGWX  ', 1, ' ')
!    call add_default ('TAUGWY  ', 1, ' ')
!
!!  (2003.06.19)

! Declare history variables for spectrum
    if (pgwv > 0) then
       call addfld ('TTGWSPEC','K/s     ',pver, 'A','T tendency - gravity wave spectrum',       phys_decomp)
       call addfld ('UTGWSPEC','m/s2    ',pver, 'A','U tendency - gravity wave spectrum',       phys_decomp)
       call addfld ('VTGWSPEC','m/s2    ',pver, 'A','V tendency - gravity wave spectrum',       phys_decomp)
       call add_default ('UTGWSPEC', 1, ' ')
       call add_default ('VTGWSPEC', 1, ' ')
    end if

    return
  end  subroutine gw_inti

!===============================================================================

  subroutine gw_intr (state,  sgh,  pblh,   dt, ptend)

!-----------------------------------------------------------------------
! Interface for multiple gravity wave drag parameterization.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    real(r8), intent(in) :: pblh(pcols)           ! planetary boundary layer height
    real(r8), intent(in) :: dt                    ! time step

    type(physics_state), intent(in) :: state      ! physics state structure
    type(physics_ptend), intent(inout):: ptend    ! parameterization tendency structure

!---------------------------Local storage-------------------------------
    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns

    integer :: i,k                                ! loop indexes
    integer :: kldv(pcols)                        ! top interface of low level stress divergence region
    integer :: kldvmn                             ! min value of kldv
    integer :: ksrc(pcols)                        ! index of top interface of source region
    integer :: ksrcmn                             ! min value of ksrc

    real(r8) :: ttgw(pcols,pver)                  ! temperature tendency
    real(r8) :: utgw(pcols,pver)                  ! zonal wind tendency
    real(r8) :: vtgw(pcols,pver)                  ! meridional wind tendency

    real(r8) :: ni(pcols,0:pver)                  ! interface Brunt-Vaisalla frequency
    real(r8) :: nm(pcols,pver)                    ! midpoint Brunt-Vaisalla frequency
    real(r8) :: rdpldv(pcols)                     ! 1/dp across low level divergence region
    real(r8) :: rhoi(pcols,0:pver)                ! interface density
    real(r8) :: tau(pcols,-pgwv:pgwv,0:pver)      ! wave Reynolds stress
    real(r8) :: tau0x(pcols)                      ! c=0 sfc. stress (zonal)
    real(r8) :: tau0y(pcols)                      ! c=0 sfc. stress (meridional)
    real(r8) :: ti(pcols,0:pver)                  ! interface temperature
    real(r8) :: ubi(pcols,0:pver)                 ! projection of wind at interfaces
    real(r8) :: ubm(pcols,pver)                   ! projection of wind at midpoints
    real(r8) :: xv(pcols)                         ! unit vectors of source wind (x)
    real(r8) :: yv(pcols)                         ! unit vectors of source wind (y)

!-----------------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

! Profiles of background state variables
    call gw_prof(lchnk, ncol, &
         state%u   , state%v   , state%t   , state%pmid   , state%pint, &
         rhoi      , ni        , ti        , nm)

!-----------------------------------------------------------------------------
! Non-orographic backgound gravity wave spectrum
!-----------------------------------------------------------------------------
    if (pgwv >0) then

! Determine the wave source for a background spectrum at ~100 mb

       call gw_bgnd (lchnk          , ncol       ,                           &
            state%u    , state%v    , state%t    , state%pmid , state%pint , &
            state%pdel , state%rpdel, state%lnpint,kldv       , kldvmn     , &
            ksrc       , ksrcmn     , rdpldv     , tau        , ubi        , &
            ubm        , xv         , yv         , PGWV       , kbotbg     )

! Solve for the drag profile
       call gw_drag_prof (lchnk     , ncol       ,                           &
            PGWV       , kbotbg     , ktopbg     , state%u    , state%v    , &
            state%t    , state%pint , state%pdel , state%rpdel, state%lnpint,&
            rhoi       , ni         , ti         , nm         , dt         , &
            kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
            tau        , ubi        , ubm        , xv         , yv         , &
            utgw       , vtgw       , tau0x      , tau0y                   )

! Add the momentum tendencies to the output tendency arrays
       do k = 1, pver
          do i = 1, ncol
             ptend%u(i,k) = utgw(i,k)
             ptend%v(i,k) = vtgw(i,k)
          end do
       end do

! Write output fields to history file
       call outfld ('UTGWSPEC', utgw, pcols, lchnk)
       call outfld ('VTGWSPEC', vtgw, pcols, lchnk)

! zero net tendencies if no spectrum computed
    else
       ptend%u = 0.
       ptend%v = 0.
    end if
!-----------------------------------------------------------------------------
! Orographic stationary gravity wave
!-----------------------------------------------------------------------------

! Determine the orographic wave source

    call gw_oro (lchnk, ncol,                                             &
         state%u    , state%v    , state%t    , sgh        , state%pmid , &
         state%pint , state%pdel , state%zm   , nm         , pblh       , &
         kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
         tau        , ubi        , ubm        , xv         , yv         )

! Solve for the drag profile
    call gw_drag_prof (lchnk, ncol,                                  &
         0          , kbotoro    , ktoporo    , state%u    , state%v    , &
         state%t    , state%pint , state%pdel , state%rpdel, state%lnpint,&
         rhoi       , ni         , ti         , nm         , dt         , &
         kldv       , kldvmn     , ksrc       , ksrcmn     , rdpldv     , &
         tau        , ubi        , ubm        , xv         , yv         , &
         utgw       , vtgw       , tau0x      , tau0y                   )

! Add the orographic tendencies to the spectrum tendencies
! Compute the temperature tendency from energy conservation (includes spectrum).
    do k = 1, pver
       do i = 1, ncol
          ptend%u(i,k) = ptend%u(i,k) + utgw(i,k)
          ptend%v(i,k) = ptend%v(i,k) + vtgw(i,k)
          ptend%s(i,k) = -(ptend%u(i,k) * (state%u(i,k) + ptend%u(i,k)*0.5*dt) &
                          +ptend%v(i,k) * (state%v(i,k) + ptend%v(i,k)*0.5*dt))
          ttgw(i,k) = ptend%s(i,k) / cpair
       end do
    end do

! Set flags for nonzero tendencies, q not yet affected by gwd
    ptend%name  = "vertical diffusion"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .TRUE.
    ptend%lu    = .TRUE.
    ptend%lv    = .TRUE.

! Write output fields to history file
    call outfld ('UTGWORO', utgw,  pcols, lchnk)
    call outfld ('VTGWORO', vtgw,  pcols, lchnk)
    call outfld ('TTGWORO', ttgw,  pcols, lchnk)
    call outfld ('TAUGWX',  tau0x, pcols, lchnk)
    call outfld ('TAUGWY',  tau0y, pcols, lchnk)
    call outfld ('SGH    ', sgh,   pcols, lchnk)

    return
  end  subroutine gw_intr

!===============================================================================
  subroutine gw_prof (lchnk, ncol, u, v, t, pm, pi, rhoi, ni, ti, nm)
!-----------------------------------------------------------------------
! Compute profiles of background state quantities for the multiple
! gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures

    real(r8), intent(out) :: rhoi(pcols,0:pver)   ! interface density
    real(r8), intent(out) :: ni(pcols,0:pver)     ! interface Brunt-Vaisalla frequency
    real(r8), intent(out) :: ti(pcols,0:pver)     ! interface temperature
    real(r8), intent(out) :: nm(pcols,pver)       ! midpoint Brunt-Vaisalla frequency

!---------------------------Local storage-------------------------------
    integer :: i,k                                ! loop indexes

    real(r8) :: dtdp
    real(r8) :: n2                                ! Brunt-Vaisalla frequency squared

!-----------------------------------------------------------------------------
! Determine the interface densities and Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------

! The top interface values are calculated assuming an isothermal atmosphere 
! above the top level.
    k = 0
    do i = 1, ncol
       ti(i,k) = t(i,k+1)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = sqrt (g*g / (cpair*ti(i,k)))
    end do

! Interior points use centered differences
    do k = 1, pver-1
       do i = 1, ncol
          ti(i,k) = 0.5 * (t(i,k) + t(i,k+1))
          rhoi(i,k) = pi(i,k) / (r*ti(i,k))
          dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
          n2 = g*g/ti(i,k) * (1./cpair - rhoi(i,k)*dtdp)
          ni(i,k) = sqrt (max (n2min, n2))
       end do
    end do

! Bottom interface uses bottom level temperature, density; next interface
! B-V frequency.
    k = pver
    do i = 1, ncol
       ti(i,k) = t(i,k)
       rhoi(i,k) = pi(i,k) / (r*ti(i,k))
       ni(i,k) = ni(i,k-1)
    end do

!-----------------------------------------------------------------------------
! Determine the midpoint Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------
    do k=1,pver
       do i=1,ncol
          nm(i,k) = 0.5 * (ni(i,k-1) + ni(i,k))
       end do
    end do

    return
  end subroutine gw_prof

!===============================================================================

  subroutine gw_oro (lchnk, ncol,                &
       u, v, t, sgh, pm, pi, dpm, zm, nm, pblh,  &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv)
!-----------------------------------------------------------------------
! Orographic source for multiple gravity wave drag parameterization.
! 
! The stress is returned for a single wave with c=0, over orography.
! For points where the orographic variance is small (including ocean),
! the returned stress is zero. 
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: sgh(pcols)            ! standard deviation of orography
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: zm(pcols,pver)        ! midpoint heights
    real(r8), intent(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency
    real(r8), intent(in) :: pblh(pcols)           ! planetary boundary layer height

    integer, intent(out) :: kldv(pcols)           ! top interface of low level stress div region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real(r8), intent(out) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    real(r8), intent(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    real(r8), intent(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    real(r8), intent(out) :: xv(pcols)            ! unit vectors of source wind (x)
    real(r8), intent(out) :: yv(pcols)            ! unit vectors of source wind (y)

!---------------------------Local storage-------------------------------
    integer :: i,k                                ! loop indexes

    real(r8) :: hdsp(pcols)                       ! surface streamline displacment height (2*sgh)
    real(r8) :: sghmax                            ! max orographic sdv to use
    real(r8) :: tauoro(pcols)                     ! c=0 stress from orography
!    real(r8) :: zldv(pcols)                       ! top of the low level stress divergence region
    real(r8) :: nsrc(pcols)                       ! b-f frequency averaged over source region
    real(r8) :: psrc(pcols)                       ! interface pressure at top of source region
    real(r8) :: rsrc(pcols)                       ! density averaged over source region
    real(r8) :: usrc(pcols)                       ! u wind averaged over source region
    real(r8) :: vsrc(pcols)                       ! v wind averaged over source region

!---------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the apropiate
! values of wind, stability, etc. for determining the wave source are 
! averages over the depth of the atmosphere pentrated by the typical mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
! 
! Also determine the depth of the low level stress divergence region, as
! the max of the boundary layer depth and the source region depth. This
! can be done here if the stress magnitude does not determine the depth,
! otherwise it must be done below.
!---------------------------------------------------------------------------
    k = pver
    do i = 1, ncol
       ksrc(i) = k-1
       kldv(i) = k-1
       psrc(i) = pi(i,k-1)
       rsrc(i) = pm(i,k)/(r*t(i,k)) * dpm(i,k)
       usrc(i) = u(i,k) * dpm(i,k)
       vsrc(i) = v(i,k) * dpm(i,k)
       nsrc(i) = nm(i,k)* dpm(i,k)
       hdsp(i) = 2.0 * sgh(i)
    end do
    do k = pver-1, pver/2, -1
       do i = 1, ncol
          if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
             ksrc(i) = k-1
             kldv(i) = k-1
             psrc(i) = pi(i,k-1)
             rsrc(i) = rsrc(i) + pm(i,k) / (r*t(i,k))* dpm(i,k)
             usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
             vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
             nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
          elseif (pblh(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
             kldv(i) = k-1
          end if
       end do
    end do
    do i = 1, ncol
       rsrc(i) = rsrc(i) / (pi(i,pver) - psrc(i))
       usrc(i) = usrc(i) / (pi(i,pver) - psrc(i))
       vsrc(i) = vsrc(i) / (pi(i,pver) - psrc(i))
       nsrc(i) = nsrc(i) / (pi(i,pver) - psrc(i))

!!     (2003.06.12)

!!     ubi(i,pver) = sqrt (usrc(i)**2 + vsrc(i)**2)
       ubi(i,pver) = max( 1.0, sqrt (usrc(i)**2 + vsrc(i)**2) )

       xv(i) = usrc(i) / ubi(i,pver)
       yv(i) = vsrc(i) / ubi(i,pver)
    end do

! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
       do i = 1, ncol
          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
       end do
    end do

! Compute the interface wind projection by averaging the midpoint winds.
! Use the top level wind at the top interface.
    do i = 1, ncol
!!     ubi(i,0) = ubm(i,1)
       ubi(i,0) = max( 1.0, ubm(i,1) )
    end do
    do k = 1, pver-1
       do i = 1, ncol
!!        ubi(i,k) = 0.5 * (ubm(i,k) + ubm(i,k+1))
          ubi(i,k) = max( 1.0, 0.5 * (ubm(i,k) + ubm(i,k+1)) )
       end do
    end do

! Determine the orographic c=0 source term following McFarlane (1987).
! Set the source top interface index to pver, if the orographic term is zero.
    do i = 1, ncol
       if ((ubi(i,pver) .gt. orovmin) .and. (hdsp(i) .gt. orohmin)) then
          sghmax = fcrit2 * (ubi(i,pver) / nsrc(i))**2
          tauoro(i) = oroko2 * min(hdsp(i)**2, sghmax) * rsrc(i) * nsrc(i) * ubi(i,pver)
       else
          tauoro(i) = 0.
          ksrc(i) = pver
          kldv(i) = pver
       end if
    end do

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).
    do i = 1, ncol
       tau(i,0,pver) = tauoro(i)
    end do
! +
! Find the top interface of the low level stress divergence region according
! to the maximum depth of three criterion.
! 1. source region depth
! 2. planetary boundary layer depth
! 3. 10 * (u_*) / N where u_* is defined from the gravity wave stresss
! = sqrt(tau/rho) using source region values
! -
!   if (kbot .lt. pver) then
!   do i = 1, ncol
!   kldv(i) = kbot
!   end do
!   else
!   do i = 1, ncol
!   zldv(i) = max (pblh(i), sgh(i)
!   zldv(i) = max (zdlv(i),
!   $           zldvcon * sqrt(tau(i,0,k)/rsrc(i)) / nsrc(i))
!   end do
!   kldv(i) = pver-1
!   do k = pver-1, pver/2, -1
!   do i = 1, ncol
!   if (zldv(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
!   kldv(i)  = k-1
!   end do
!   end do
!   end if

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver
    do i = 1, ncol
       ksrcmn = min(ksrcmn, ksrc(i))
       kldvmn = min(kldvmn, kldv(i))
       if (kldv(i) .ne. pver) then
          rdpldv(i) = 1. / (pi(i,kldv(i)) - pi(i,pver))
       end if
    end do
    if (fracldv .le. 0.) kldvmn = pver

    return
  end subroutine gw_oro

!===============================================================================
  subroutine gw_bgnd (lchnk, ncol, &
       u, v, t, pm, pi, dpm, rdpm, piln,     &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv, &
       ngwv, kbot)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
    use phys_grid, only: get_rlat_all_p
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns

    integer, intent(in) :: kbot                   ! index of bottom (source) interface
    integer, intent(in) :: ngwv                   ! number of gravity waves to use

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pm(pcols,pver)        ! midpoint pressures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real(r8), intent(in) :: piln(pcols,0:pver)    ! ln(interface pressures)

    integer, intent(out) :: kldv(pcols)           ! top interface of low level stress divergence region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc(pcols)           ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real(r8), intent(in) :: rdpldv(pcols)        ! 1/dp across low level divergence region
    real(r8), intent(out) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress
    real(r8), intent(out) :: ubi(pcols,0:pver)    ! projection of wind at interfaces
    real(r8), intent(out) :: ubm(pcols,pver)      ! projection of wind at midpoints
    real(r8), intent(out) :: xv(pcols)            ! unit vectors of source wind (x)
    real(r8), intent(out) :: yv(pcols)            ! unit vectors of source wind (y)

!---------------------------Local storage-------------------------------
    integer :: i,k,l                              ! loop indexes

    real(r8) :: rlat(pcols)                       ! latitude in radians for columns
    real(r8) :: tauback(pcols)                    ! background stress at c=0
    real(r8) :: usrc(pcols)                       ! u wind averaged over source region
    real(r8) :: vsrc(pcols)                       ! v wind averaged over source region
    real(r8) :: al0                               ! Used in lat dependence of GW spec. 
    real(r8) :: dlat0                             ! Used in lat dependence of GW spec.
    real(r8) :: flat_gw                           ! The actual lat dependence of GW spec.
    real(r8) :: pi_g                              ! 3.14........

!---------------------------------------------------------------------------
! Determine the source layer wind and unit vectors, then project winds.
!---------------------------------------------------------------------------

! Just use the source level interface values for the source
! wind speed and direction (unit vector).
    k = kbot
    do i = 1, ncol
       ksrc(i) = k
       kldv(i) = k
       usrc(i) = 0.5*(u(i,k+1)+u(i,k))
       vsrc(i) = 0.5*(v(i,k+1)+v(i,k))
       ubi(i,kbot) = sqrt (usrc(i)**2 + vsrc(i)**2)
       xv(i) = usrc(i) / ubi(i,k)
       yv(i) = vsrc(i) / ubi(i,k)
    end do

! Project the local wind at midpoints onto the source wind.
    do k = 1, kbot
       do i = 1, ncol
          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
       end do
    end do

! Compute the interface wind projection by averaging the midpoint winds.
! Use the top level wind at the top interface.
    do i = 1, ncol
       ubi(i,0) = ubm(i,1)
    end do
    do k = 1, kbot-1
       do i = 1, ncol
          ubi(i,k) = 0.5 * (ubm(i,k) + ubm(i,k+1))
       end do
    end do

!-----------------------------------------------------------------------
! Gravity wave sources
!-----------------------------------------------------------------------

! Determine the background stress at c=0
    do i=1,ncol
      tauback(i) = taubgnd * tauscal
    enddo

!	Include dependence on latitude:
! 	The lat function was obtained by RR Garcia as 
!	currently used in his 2D model
!       [Added by F. Sassi on May 30, 1997]

    pi_g=4.*atan(1.)
    al0=40.*pi_g/180.
    dlat0=10.*pi_g/180.
!
    call get_rlat_all_p(lchnk, ncol, rlat)
    do i=1,ncol
      flat_gw= 0.5*(1.+tanh((rlat(i)-al0)/dlat0)) + 0.5*(1.+tanh(-(rlat(i)+al0)/dlat0)) 
      flat_gw=max(flat_gw,0.2_r8)
      tauback(i)=tauback(i)*flat_gw
    enddo

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).

    do l = 1, ngwv
       do i = 1, ncol
          tau(i, l,kbot) = tauback(i) * exp(-(c(l)/30.)**2)
          tau(i,-l,kbot) = tau(i, l,kbot)
       end do
    end do
    do i = 1, ncol
       tau(i,0,kbot) = tauback(i)
    end do

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver

    return
  end  subroutine gw_bgnd

!===============================================================================
  subroutine gw_drag_prof (lchnk, ncol,                           &
       ngwv, kbot, ktop, u, v, t,                                 &
       pi, dpm, rdpm, piln, rhoi, ni, ti, nm, dt,                 &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv, &
       ut, vt, tau0x, tau0y)
!-----------------------------------------------------------------------
! Solve for the drag profile from the multiple gravity wave drag
! parameterization.
! 1. scan up from the wave source to determine the stress profile
! 2. scan down the stress profile to determine the tendencies
!     => apply bounds to the tendency
!          a. from wkb solution
!          b. from computational stability constraints
!     => adjust stress on interface below to reflect actual bounded tendency
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: kbot                   ! index of bottom (source) interface
    integer, intent(in) :: ktop                   ! index of top interface of gwd region
    integer, intent(in) :: ngwv                   ! number of gravity waves to use
    integer, intent(in) :: kldv(pcols)            ! top interface of low level stress  divergence region
    integer, intent(in) :: kldvmn                 ! min value of kldv
    integer, intent(in) :: ksrc(pcols)            ! index of top interface of source region
    integer, intent(in) :: ksrcmn                 ! min value of ksrc

    real(r8), intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real(r8), intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real(r8), intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real(r8), intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real(r8), intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real(r8), intent(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real(r8), intent(in) :: piln(pcols,0:pver)    ! ln(interface pressures)
    real(r8), intent(in) :: rhoi(pcols,0:pver)    ! interface density
    real(r8), intent(in) :: ni(pcols,0:pver)      ! interface Brunt-Vaisalla frequency
    real(r8), intent(in) :: ti(pcols,0:pver)      ! interface temperature
    real(r8), intent(in) :: nm(pcols,pver)        ! midpoint Brunt-Vaisalla frequency
    real(r8), intent(in) :: dt                    ! time step
    real(r8), intent(in) :: rdpldv(pcols)         ! 1/dp across low level divergence region
    real(r8), intent(in) :: ubi(pcols,0:pver)     ! projection of wind at interfaces
    real(r8), intent(in) :: ubm(pcols,pver)       ! projection of wind at midpoints
    real(r8), intent(in) :: xv(pcols)             ! unit vectors of source wind (x)
    real(r8), intent(in) :: yv(pcols)             ! unit vectors of source wind (y)

    real(r8), intent(inout) :: tau(pcols,-pgwv:pgwv,0:pver)! wave Reynolds stress

    real(r8), intent(out) :: ut(pcols,pver)       ! zonal wind tendency
    real(r8), intent(out) :: vt(pcols,pver)       ! meridional wind tendency
    real(r8), intent(out) :: tau0x(pcols)         ! c=0 sfc. stress (zonal)
    real(r8), intent(out) :: tau0y(pcols)         ! c=0 sfc. stress (meridional)

!---------------------------Local storage-------------------------------
    integer :: i,k,l                              ! loop indexes

    real(r8) :: d(pcols)                          ! "total" diffusivity 
    real(r8) :: dsat(pcols,-pgwv:pgwv)            ! saturation diffusivity
    real(r8) :: dscal                             ! fraction of dsat to use
    real(r8) :: mi                                ! imaginary part of vertical wavenumber
    real(r8) :: taudmp                            ! stress after damping
    real(r8) :: taumax(pcols)                     ! max(tau) for any l
    real(r8) :: tausat(pcols,-pgwv:pgwv)          ! saturation stress
    real(r8) :: ubmc(pcols,-pgwv:pgwv)            ! (ub-c)
    real(r8) :: ubmc2                             ! (ub-c)**2
    real(r8) :: ubt(pcols,pver)                   ! ubar tendency
    real(r8) :: ubtl                              ! ubar tendency from wave l
    real(r8) :: ubtlsat                           ! saturation tendency

! Initialize gravity wave drag tendencies to zero

    do k=1,pver
       do i=1,pcols
          ut(i,k) = 0.
          vt(i,k) = 0.
       end do
    end do

!---------------------------------------------------------------------------
! Compute the stress profiles and diffusivities
!---------------------------------------------------------------------------

! Loop from bottom to top to get stress profiles      

    do k = kbot-1, ktop, -1

! Determine the absolute value of the saturation stress and the diffusivity
! for each wave.
! Define critical levels where the sign of (u-c) changes between interfaces.

       d(:ncol) = dback
       do l = -ngwv, ngwv
          do i = 1, ncol
             ubmc(i,l) = ubi(i,k) - c(l)
             tausat(i,l) = abs (effkwv * rhoi(i,k) * ubmc(i,l)**3 / (2.*ni(i,k)) )
             if (tausat(i,l) .le. taumin) tausat(i,l) = 0.0
             if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0) tausat(i,l) = 0.0
             dsat(i,l) = (ubmc(i,l) / ni(i,k))**2 * &
                  (effkwv * ubmc(i,l)**2 / (rog * ti(i,k) * ni(i,k)) - alpha(k))
             dscal = min (1.0_r8, tau(i,l,k+1) / (tausat(i,l)+taumin))
             d(i) = max( d(i), dscal * dsat(i,l))
          end do
       end do

! Compute stress for each wave. The stress at this level is the min of 
! the saturation stress and the stress at the level below reduced by damping.
! The sign of the stress must be the same as at the level below.

       do l = -ngwv, ngwv
          do i = 1, ncol
             ubmc2 = max(ubmc(i,l)**2, ubmc2mn)
             mi = ni(i,k) / (2. * kwv * ubmc2) * (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
             taudmp = tau(i,l,k+1) * exp(-2.*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k)))
             if (taudmp .le. taumin) taudmp = 0.
             tau(i,l,k) = min (taudmp, tausat(i,l))
          end do
       end do

! The orographic stress term does not change across the source region
! Note that k ge ksrcmn cannot occur without an orographic source term

       if (k .ge. ksrcmn) then
          do i = 1, ncol
             if (k .ge. ksrc(i)) then
                tau(i,0,k) = tau(i,0,pver) 
             end if
          end do
       end if

! Require that the orographic term decrease linearly (with pressure) 
! within the low level stress divergence region. This supersedes the above
! requirment of constant stress within the source region.
! Note that k ge kldvmn cannot occur without an orographic source term, since
! kldvmn=pver then and k<=pver-1

       if (k .ge. kldvmn) then
          do i = 1, ncol
             if (k .ge. kldv(i)) then
                tau(i,0,k) = min (tau(i,0,k), tau(i,0,pver)  * &
                     (1. - fracldv * (pi(i,k)-pi(i,pver)) * rdpldv(i)))
             end if
          end do
       end if

! Apply lower bounds to the stress if ngwv > 0.

       if (ngwv .ge. 1) then

! Determine the max value of tau for any l

          do i = 1, ncol
             taumax(i) = tau(i,-ngwv,k)
          end do
          do l = -ngwv+1, ngwv
             do i = 1, ncol
                taumax(i) = max(taumax(i), tau(i,l,k))
             end do
          end do
          do i = 1, ncol
             taumax(i) = mxrange * taumax(i)
          end do

! Set the min value of tau for each wave to the max of mxrange*taumax or
! mxasym*tau(-c)

          do l = 1, ngwv
             do i = 1, ncol
                tau(i, l,k) = max(tau(i, l,k), taumax(i))
                tau(i, l,k) = max(tau(i, l,k), mxasym*tau(i,-l,k))
                tau(i,-l,k) = max(tau(i,-l,k), taumax(i))
                tau(i,-l,k) = max(tau(i,-l,k), mxasym*tau(i, l,k))
             end do
          end do
          l = 0
          do i = 1, ncol
             tau(i,l,k) = max(tau(i,l,k), mxasym * 0.5 * (tau(i,l-1,k) + tau(i,l+1,k)))
          end do

       end if

    end do

! Put an upper bound on the stress at the top interface to force some stress
! divergence in the top layer. This prevents the easterlies from running
! away in the summer mesosphere, since most of the gravity wave activity
! will pass through a top interface at 75--80 km under strong easterly
! conditions. 
!++BAB fix to match ccm3.10
!!$    do l = -ngwv, ngwv
!!$       do i = 1, ncol
!!$          tau(i,l,0) = min(tau(i,l,0), 0.5*tau(i,l,1))
!!$       end do
!!$    end do
!--BAB fix to match ccm3.10

!---------------------------------------------------------------------------
! Compute the tendencies from the stress divergence.
!---------------------------------------------------------------------------

! Loop over levels from top to bottom
    do k = ktop+1, kbot

! Accumulate the mean wind tendency over wavenumber.
       do i = 1, ncol
          ubt (i,k) = 0.0
       end do
       do l = -ngwv, ngwv
          do i = 1, ncol

! Determine the wind tendency including excess stress carried down from above.
             ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)

! Require that the tendency be no larger than the analytic solution for
! a saturated region [proportional to (u-c)^3].
             ubtlsat = effkwv * abs((c(l)-ubm(i,k))**3) /(2.*rog*t(i,k)*nm(i,k))
             ubtl = min(ubtl, ubtlsat)

! Apply tendency limits to maintain numerical stability.
! 1. du/dt < |c-u|/dt  so u-c cannot change sign (u^n+1 = u^n + du/dt * dt)
! 2. du/dt < tndmax    so that ridicuously large tendency are not permitted
             ubtl = min(ubtl, umcfac * abs(c(l)-ubm(i,k)) / dt)
             ubtl = min(ubtl, tndmax)

! Accumulate the mean wind tendency over wavenumber.
             ubt (i,k) = ubt (i,k) + sign(ubtl, c(l)-ubm(i,k))

! Redetermine the effective stress on the interface below from the wind 
! tendency. If the wind tendency was limited above, then the new stress
! will be small than the old stress and will cause stress divergence in
! the next layer down. This has the effect of smoothing large stress 
! divergences downward while conserving total stress.
             tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g
          end do
       end do

! Project the mean wind tendency onto the components and scale by "efficiency".
       do i = 1, ncol
          ut(i,k) = ubt(i,k) * xv(i) * effgw
          vt(i,k) = ubt(i,k) * yv(i) * effgw
       end do

! End of level loop
    end do

!-----------------------------------------------------------------------
! Project the c=0 stress (scaled) in the direction of the source wind for recording
! on the output file.
!-----------------------------------------------------------------------
    do i = 1, ncol
       tau0x(i) = tau(i,0,kbot) * xv(i) * effgw
       tau0y(i) = tau(i,0,kbot) * yv(i) * effgw
    end do

    return
  end subroutine gw_drag_prof

end module gw_drag

