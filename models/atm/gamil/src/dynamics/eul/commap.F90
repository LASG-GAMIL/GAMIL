module commap

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plond, plev, plat, beglatexdyn, endlatexdyn, numbnd

  real(r8) w(plat)            ! gaussian weights 
  real(r8) clat  (plat)       ! model latitudes (radians)
  real(r8) latdeg(plat)       ! model latitudes (degrees)
  real(r8) clon  (plond,plat) ! model longitudes (radians)
  real(r8) londeg(plond,plat) ! model longitudes (degrees)

  real(r8) dy
  real(r8) dx
  real(r8) ythu(1-numbnd:plat+numbnd)
  real(r8) ythv(1-numbnd:plat+numbnd)
  real(r8) wtgu(1-numbnd:plat+numbnd) 
  real(r8) wtgv(1-numbnd:plat+numbnd)

  real(r8), allocatable :: sinu(:)
  real(r8), allocatable :: sinv(:)
  real(r8), allocatable :: oux (:)
  real(r8), allocatable :: ouy (:)
  real(r8), allocatable :: ovx (:)
  real(r8), allocatable :: ovy (:)
  real(r8), allocatable :: ff  (:)
  real(r8), allocatable :: cur (:)

  real(r8) latmesh_b ! Controls grid density around Poles along meridional direction.

contains

  subroutine initialize_hpar

    allocate(sinu(beglatexdyn:endlatexdyn))
    allocate(sinv(beglatexdyn:endlatexdyn))
    allocate(oux (beglatexdyn:endlatexdyn))
    allocate(ouy (beglatexdyn:endlatexdyn))
    allocate(ovx (beglatexdyn:endlatexdyn))
    allocate(ovy (beglatexdyn:endlatexdyn))
    allocate(ff  (beglatexdyn:endlatexdyn))
    allocate(cur (beglatexdyn:endlatexdyn))

    sinu(:) = inf
    sinv(:) = inf
    oux (:) = inf
    ouy (:) = inf
    ovx (:) = inf
    ovy (:) = inf
    ff  (:) = inf
    cur (:) = inf

  end subroutine initialize_hpar

end module commap
