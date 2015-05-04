#include <misc.h>
#include <params.h>
subroutine grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
                  lam     ,phi     ,dphi    ,gw      ,sinlam  , &
                  coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
                  lbasiz  ,lbassi  ,detam   ,detai   ,kdpmpf  , &
                  kdpmph  ,cwava   ,phigs   )
!-----------------------------------------------------------------------
!
! Purpose:
! Initialize model and extended grid parameters
! Initialize weights for Lagrange cubic derivative estimates
! Initialize weights for Lagrange cubic interpolant
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: grdini.F90,v 1.2.42.2 2002/06/15 13:48:24 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use rgrid
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pmap                 ! dimension of artificial vert. grid
  real(r8), intent(in)   :: etamid(plev)         ! full-level model vertical grid
  real(r8), intent(in)   :: etaint(plevp)        ! half-level model vertical grid
  real(r8), intent(in)   :: gravit               ! gravitational constant
  real(r8), intent(out)  :: dlam  (platd)        ! longitudinal grid interval (radians)
  real(r8), intent(out)  :: lam   (plond,platd)  ! longitudinal coords of extended grid
  real(r8), intent(out)  :: phi   (platd)        ! latitudinal  coords of extended grid
  real(r8), intent(out)  :: dphi  (platd)        ! latitude intervals (radians)
  real(r8), intent(out)  :: gw    (plat)         ! Gaussian weights
  real(r8), intent(out)  :: sinlam(plond,platd)  ! sin(lam) model domain only
  real(r8), intent(out)  :: coslam(plond,platd)  ! cos(lam) model domain only
  real(r8), intent(out)  :: lbasdy(4,2,platd)    ! latitude derivative weights
  real(r8), intent(out)  :: lbasdz(4,2,plev)     ! vertical (full levels) deriv weights
  real(r8), intent(out)  :: lbassd(4,2,plevp)    ! vertical (half levels) deriv weights
  real(r8), intent(out)  :: lbasiy(4,2,platd)    ! Lagrange cubic interp weights (lat.)
  real(r8), intent(out)  :: lbasiz(4,2,plev)     ! Lagrange cubic interp wghts (full lev)
  real(r8), intent(out)  :: lbassi(4,2,plevp)    ! Lagrange cubic interp wghts (half lev)
  real(r8), intent(out)  :: detam (plev)         ! intervals between vertical full levs.
  real(r8), intent(out)  :: detai (plevp)        ! intervals between vertical half levs.
  integer , intent(out)  :: kdpmpf(pmap)         ! artificial full vertical grid indices
  integer , intent(out)  :: kdpmph(pmap)         ! artificial half vertical grid indices
  real(r8), intent(out)  :: cwava(plat)          ! weight applied to global integrals
  real(r8), intent(out)  :: phigs                ! Cutoff latitude for local geodesic
!
!---------------------------Local variables-----------------------------
!
  integer j                 ! index
  integer k                 ! index
  real(r8) etamln (plev)    ! log(etamid)
  real(r8) etailn (plevp)   ! log(etaint)
  real(r8) detamln(plev)    ! dlog(etamid)
  real(r8) detailn(plevp)   ! dlog(etaint)
!
!-----------------------------------------------------------------------
!
! Initialize extended horizontal grid coordinates.
!
  call grdxy(dlam    ,lam     ,phi     ,gw      ,sinlam  , &
             coslam  )
!
! Basis functions for computing Lagrangian cubic derivatives
! on unequally spaced latitude and vertical grids.
!
  call basdy(phi     ,lbasdy  )
  call basdz(plev    ,etamid  ,lbasdz  )
  call basdz(plevp   ,etaint  ,lbassd  )
!
! Basis functions for computing weights for Lagrangian cubic
! interpolation on unequally spaced latitude and height grids.
!
  call basiy(phi     ,lbasiy  )
  call basiz(plev    ,etamid  ,lbasiz  )
  call basiz(plevp   ,etaint  ,lbassi  )
!
! Compute interval lengths in latitudinal grid
!
  do j = 1,platd-1
     dphi(j) = phi(j+1) - phi(j)
  end do
!
! Compute interval lengths in vertical grids.
!
  do k = 1,plev
     etamln(k) = log(etamid(k))
  end do
  do k = 1,plevp
     etailn(k) = log(etaint(k))
  end do
  do k = 1,plev-1
     detam  (k) = etamid(k+1) - etamid(k)
     detamln(k) = etamln(k+1) - etamln(k)
  end do
  do k = 1,plev
     detai  (k) = etaint(k+1) - etaint(k)
     detailn(k) = etailn(k+1) - etailn(k)
  end do
!
! Build artificial evenly spaced vertical grid for use in determining
! vertical position of departure point.
! Build one grid for full model levels and one for half levels.
!
  call vrtmap(plev    ,pmap    ,etamln  ,detamln ,kdpmpf  )
  call vrtmap(plevp   ,pmap    ,etailn  ,detailn ,kdpmph  )
!
! Compute tracer integral constant
!
  do j=1,plat
     cwava(j) = 1./(nlon(j)*gravit)
  end do
!
! Initialize cutoff latitude (poleward of which local geodesic
! coordinates will be used); radians
!
  phigs = 1.221730
!
  return
end subroutine grdini

