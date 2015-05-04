module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond, plat

   real(r8) gw(plat)         ! weights used in physics, history
   real(r8) w(plat)          ! non-gaussian weights (one hemisphere)
   real(r8) w_staggered(plat-1)      ! non-gaussian weights for the staggered wind arrays
!
   real(r8) clat(plat)             ! model latitudes (radians)
   real(r8) clat_staggered (plat-1)  ! model latitudes on staggered grid (radians)
   real(r8) clon(plond,plat)             ! model longitudes (radians)
   real(r8) latdeg(plat)           ! model latitudes (degrees)
   real(r8) londeg(plond,plat)           ! model longitudes (degrees)
end module commap
