module commap

!! (wh 2003.05.14)
!! (wh 2003.10.24)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond, plev, plat, beglatexdyn, endlatexdyn  &  !!(wh 2003.10.24)
                                                   , numbnd          !!(wh 2003.11.14)

   real(r8) :: w(plat)            ! gaussian weights 
   real(r8) :: clat  (plat)       ! model latitudes (radians)
   real(r8) :: latdeg(plat)       ! model latitudes (degrees)
   real(r8) :: clon  (plond,plat) ! model longitudes (radians)
   real(r8) :: londeg(plond,plat) ! model longitudes (degrees)
   real(r8) :: dy                 ! DY in fm2003
   real(r8) :: dx                 ! DX in fm2003

   real(r8) :: ythu(1-numbnd:plat+numbnd)  ! fm2003
   real(r8) :: ythv(1-numbnd:plat+numbnd)  ! fm2003
   real(r8) :: wtgu(1-numbnd:plat+numbnd)  ! fm2003 weights 
   real(r8) :: wtgv(1-numbnd:plat+numbnd)  ! fm2003 weights

   real(r8),allocatable :: sinu(:)         !
   real(r8),allocatable :: sinv(:)         !
   real(r8),allocatable :: oux (:)         !
   real(r8),allocatable :: ouy (:)         !
   real(r8),allocatable :: ovx (:)         !  fm2003 HPAR variables
   real(r8),allocatable :: ovy (:)         !
   real(r8),allocatable :: ff  (:)         !
   real(r8),allocatable :: cur (:)         !

CONTAINS

   subroutine initialize_hpar

       allocate (sinu (beglatexdyn:endlatexdyn))
       allocate (sinv (beglatexdyn:endlatexdyn))
       allocate (oux  (beglatexdyn:endlatexdyn))
       allocate (ouy  (beglatexdyn:endlatexdyn))
       allocate (ovx  (beglatexdyn:endlatexdyn))
       allocate (ovy  (beglatexdyn:endlatexdyn))
       allocate (ff   (beglatexdyn:endlatexdyn))
       allocate (cur  (beglatexdyn:endlatexdyn))

       sinu(:) = inf
       sinv(:) = inf
       oux (:) = inf
       ouy (:) = inf
       ovx (:) = inf
       ovy (:) = inf
       ff  (:) = inf
       cur (:) = inf

       return
   end subroutine initialize_hpar

end module commap
