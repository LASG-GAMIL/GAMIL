#include <misc.h>
#include <params.h>
module comslt
!
! Semi-Lagrangian transport
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plond, plev, plat, beglat, endlat
  use constituents, only: pcnst
  use infnan

  implicit none

  real(r8) hw1    (pcnst)          ! Pre-SLT global integral of constituent
  real(r8) hw2    (pcnst)          ! Post-SLT global integral of const.
  real(r8) hw3    (pcnst)          ! Global integral for denom. of expr. for alpha
  real(r8) alpha  (pcnst)          ! alpha(m) = ( hw1(m) - hw2(m) )/hw3(m)
  real(r8) hw1lat (pcnst,plat)     ! lat contribution to const. mass integral
  real(r8) engy1lat(plat)          ! lat contribution to total energy integral
  real(r8) levknt (plev,plev)      ! counter for departure point binning statistics
  real(r8) levkntl(plev,plev,plat) ! counter for departure point binning statistics
  real(r8) epssld                  ! "epsilon" for SLD decentering algorithm
  real(r8) phigs                   ! latitude cutoff for using local geodesic algorithm
  real(r8) gamma  (plev,plev)      ! SLD coefficient
  real(r8) hortalc(plev)           ! analytic Hortal temperature correction term
  real(r8) hdel   (0:plev-1)       ! del-"hortalc" in vertical

  real(r8), allocatable :: qfcst(:,:,:,:)

contains

  subroutine initialize_comslt

    allocate (qfcst(plond,plev,pcnst,beglat:endlat))
    qfcst (:,:,:,:) = inf

    return
  end subroutine initialize_comslt

end module comslt

