#include <misc.h>
#include <params.h>

subroutine dynpkg (pmap    ,t2      ,fu      ,fv      ,etamid  , &
                   etaint  ,cwava   ,detam   ,dlam    ,lam     , &
                   phi     ,dphi    ,sinlam  ,coslam  ,lbasdy  , &
                   lbasdz  ,lbasiy  ,lbassi  ,lbasiz  ,detai   , &
                   kdpmpf  ,kdpmph  ,flx_net , ztodt   )
!-----------------------------------------------------------------------
!
! Purpose:
! driving routines dynamics,and transport
!
! Note that this routine has many "#if ..." constructs.  There are 2 code 
! branches based on the token PVP to address the fact that storage order for 
! Fourier coefficients is different to enable optimal performance 
! characteristics on both shared-memory and distributed-memory computing 
! platforms.  The message-passing model needs to test SPMD since many
! arrays have their space allocated dynamically.
! COUP_CSM must be checked in order to invoke the proper calling
! sequence for running the CSM model
!
! Original version:  CCM3
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pmap                 ! max dimension of evenly spaced vert.
!                                                ! grid used to map the departure pts
!                                                ! into true model levels.
  real(r8), intent(inout):: t2(plond,plev,beglat:endlat)  ! temperature tendency
  real(r8), intent(inout):: fu(plond,plev,beglat:endlat)  ! u wind tendency
  real(r8), intent(inout):: fv(plond,plev,beglat:endlat)  ! v wind tendency
  real(r8), intent(in)   :: etamid (plev)        ! vertical coords at midpoints 
  real(r8), intent(in)   :: etaint (plevp)       ! vertical coords at interfaces
  real(r8), intent(in)   :: cwava  (plat)        ! weight applied to global integrals
  real(r8), intent(in)   :: detam  (plev)        ! intervals between vert full levs.
  real(r8), intent(in)   :: dlam   (platd)       ! longitudinal grid interval (radians)
  real(r8), intent(in)   :: lam    (plond,platd) ! longitude coords of extended grid
  real(r8), intent(in)   :: phi    (platd)       ! latitude  coords of extended grid
  real(r8), intent(in)   :: dphi   (platd)       ! latitude intervals (radians)
  real(r8), intent(in)   :: sinlam (plond,platd) ! sin(lam) model domain only
  real(r8), intent(in)   :: coslam (plond,platd) ! cos(lam) model domain only
  real(r8), intent(in)   :: lbasdy (4,2,platd)   ! latitude derivative weights
  real(r8), intent(in)   :: lbasdz (4,2,plev)    ! vert (full levels) deriv wghts 
  real(r8), intent(in)   :: lbasiy (4,2,platd)   ! Lagrange cubic interp wghts (lat.) 
  real(r8), intent(in)   :: lbassi (4,2,plevp)   ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: lbasiz (4,2,plev)    ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: detai  (plevp)       ! intervals between vert half levs.
  integer , intent(in)   :: kdpmpf (pmap)        ! artificial full vert grid indices
  integer , intent(in)   :: kdpmph (pmap)        ! artificial half vert grid indices
  real(r8), intent(in)   :: flx_net(plond,beglat:endlat)  ! net flux from physics
  real(r8), intent(in)   :: ztodt                ! time step (s)
!
!---------------------------Local workspace-----------------------------
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSAC and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
#if ( defined PVP )
  real(r8) grlps1(2*pmmax,plat/2)      ! ------------------------------
  real(r8) grlps2(2*pmmax,plat/2)      ! |
  real(r8) grt1  (2*pmmax,plev,plat/2) ! |
  real(r8) grt2  (2*pmmax,plev,plat/2) ! |
  real(r8) grq1  (2*pmmax,plev,plat/2) ! |
  real(r8) grq2  (2*pmmax,plev,plat/2) ! |
  real(r8) grz1  (2*pmmax,plev,plat/2) ! |
  real(r8) grz2  (2*pmmax,plev,plat/2) ! |
  real(r8) grd1  (2*pmmax,plev,plat/2) ! |
  real(r8) grd2  (2*pmmax,plev,plat/2) ! |
  real(r8) grfu1 (2*pmmax,plev,plat/2) ! |- see linemsac and quad for
  real(r8) grfu2 (2*pmmax,plev,plat/2) ! |  definitions
  real(r8) grfv1 (2*pmmax,plev,plat/2) ! |
  real(r8) grfv2 (2*pmmax,plev,plat/2) ! |
#else
  real(r8) grlps1(     2*pmmax,plat/2) ! ------------------------------
  real(r8) grlps2(     2*pmmax,plat/2) ! |
  real(r8) grt1  (plev,2*pmmax,plat/2) ! |
  real(r8) grt2  (plev,2*pmmax,plat/2) ! |
  real(r8) grq1  (plev,2*pmmax,plat/2) ! |
  real(r8) grq2  (plev,2*pmmax,plat/2) ! |
  real(r8) grz1  (plev,2*pmmax,plat/2) ! |
  real(r8) grz2  (plev,2*pmmax,plat/2) ! |
  real(r8) grd1  (plev,2*pmmax,plat/2) ! |
  real(r8) grd2  (plev,2*pmmax,plat/2) ! |
  real(r8) grfu1 (plev,2*pmmax,plat/2) ! |- see linemsac and quad for
  real(r8) grfu2 (plev,2*pmmax,plat/2) ! |  definitions
  real(r8) grfv1 (plev,2*pmmax,plat/2) ! |
  real(r8) grfv2 (plev,2*pmmax,plat/2) ! |
#endif
  real(r8) vcour(plev,plat)            ! maximum Courant number in slice
  real(r8) vmax2d (plev,plat)          ! max. wind at each level, latitude
  real(r8) vmax2dt(plev,plat)          ! max. truncated wind at each lvl,lat
!
!----------------------------------------------------------
! SCANDYN Dynamics scan
!----------------------------------------------------------
!
  call t_startf ('scandyn')
  call scandyn (ztodt   ,pmap    ,kdpmpf  ,kdpmph  ,lam     , &
                phi     ,dphi    ,sinlam  ,coslam  ,lbasdy  , &
                lbasdz  ,lbasiy  ,lbasiz  ,lbassi  ,detam   , &
                detai   ,dlam    ,cwava   ,etamid  ,etaint  , &
                grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                fu      ,fv      ,t2      ,flx_net , &
                vcour   ,vmax2d  ,vmax2dt )
  call t_stopf ('scandyn')
!
!----------------------------------------------------------
! Accumulate spectral coefficients
!----------------------------------------------------------
!
  call t_startf('dyndrv')
  call dyndrv(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
              grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
              grz2    ,grd2    ,grfu2   ,grfv2   ,vmax2d  , &
              vmax2dt ,vcour   )
  call t_stopf('dyndrv')
!
!----------------------------------------------------------
! Second gaussian scan (spectral -> grid)
!----------------------------------------------------------
!
  call t_startf('scan2')
  call scan2 (ztodt   ,cwava   ,etamid  )
  call t_stopf('scan2')

  return
end subroutine dynpkg

