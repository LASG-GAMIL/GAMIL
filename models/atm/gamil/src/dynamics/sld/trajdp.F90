#include <misc.h>
#include <params.h>
subroutine trajdp(lat     ,jcen    ,dt      ,locgeo  ,lvert   , &
                  lam     ,phi     ,etamid  ,upr     ,vpr     , &
                  lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
                  sigdp   ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Determine trajectory departure point for each arrival point based upon
! departure point increments
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: trajdp.F90,v 1.4.22.1 2002/06/15 13:48:33 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: lat               ! index of lat slice (model)
  integer , intent(in)   :: jcen              ! index of lat slice(extend)
  real(r8), intent(in)   :: dt                ! time step (seconds)
  logical , intent(in)   :: locgeo            ! local geodesic flag
  logical , intent(in)   :: lvert             ! flag to compute vertical trajectory
  real(r8), intent(in)   :: lam   (plond)     ! long. coord of model grid
  real(r8), intent(in)   :: phi   (platd)     ! lat.  coord of model grid
  real(r8), intent(in)   :: etamid(plev)      ! vertical full levels
  real(r8), intent(in)   :: upr   (plon,plev) ! interpolated u field (local geodesic)
  real(r8), intent(in)   :: vpr   (plon,plev) ! interpolated v field (local geodesic)
  real(r8), intent(in)   :: lampr (plon,plev) ! trajectory increment (x-direction)
  real(r8), intent(in)   :: phipr (plon,plev) ! trajectory increment (y-direction)
  real(r8), intent(in)   :: sigpr (plon,plev) ! trajectory increment (vertical)
  real(r8), intent(out)  :: lamdp (plon,plev) ! zonal      departure pt. coord.
  real(r8), intent(out)  :: phidp (plon,plev) ! meridional departure pt. coord.
  real(r8), intent(out)  :: sigdp (plon,plev) ! vertical   departure pt. coord.
  integer , intent(in)   :: nlon              ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) fac              ! 1 - eps
  real(r8) botlim           ! bottom limit of the trajectory
  real(r8) phi0             ! Current latitude (radians)
  real(r8) cphi0            ! cos latitude
  real(r8) sphi0            ! sin latitude
  real(r8) dist             ! computed distance**2
  real(r8) distmx           ! max distance**2
  real(r8) pi2              ! pi/2
  real(r8) phipi2           ! pi/2
  real(r8) sgnphi0          ! holds sign of phi0
  real(r8) pi               ! pi
  real(r8) twopi            ! 2*pi
!
  real(r8) clamgc           ! |
  real(r8) cphigc           ! |
  real(r8) dlamsc           ! | 
  real(r8) slamgc           ! | -- temporary variables
  real(r8) slam2            ! |
  real(r8) sphigc           ! |
!
  integer i                 ! |
  integer ii                ! |
  integer k                 ! | -- indices
  integer ibad              ! |
  integer kbad              ! |
!
  logical lstop             ! flag to stop run if departure point is over the pole
!
!-----------------------------------------------------------------------
!
  fac     = 1. - 10.*epsilon(fac)
  pi      = 4.*atan(1.)
  twopi   = 2.*pi
  pi2     = pi/2.
  phi0    = phi(jcen)
  cphi0   = cos( phi0 )
  sphi0   = sin( phi0 )
  sgnphi0 = sign( 1._r8, phi0 )
  distmx  = (sign(pi2,phi0) - phi0)/(1.1*dt)
  distmx  = distmx*distmx
!
! Compute coordinates of departure point.
!
  lstop = .false.
  do k = 1,plev
!
! if near pole, convert from g.c. to spherical
!
     if (locgeo) then
        do i = 1,nlon
           ii = i + i1 - 1
           sphigc = sin( phipr(i,k) )
           cphigc = cos( phipr(i,k) )
           slamgc = sin( lampr(i,k) )
           clamgc = cos( lampr(i,k) )
           phidp(i,k)  = asin((sphigc*cphi0 + cphigc*sphi0*clamgc)*fac)
           if ( abs(phidp(i,k)) .ge. phi(j1+plat)*fac ) &
                                          phidp(i,k) = sign( phi(j1+plat),phidp(i,k) )*fac
           dlamsc = asin((slamgc*cphigc/cos(phidp(i,k)))*fac)
!
! If traj is over pole, check for proper branch of arcsin
!
           dist   = upr(i,k)**2 + vpr(i,k)**2
           if( dist .gt. distmx) then
              slam2  = slamgc**2
              phipi2 = asin((sqrt((slam2 - 1.)/(slam2 - 1./cphi0**2)))*fac)
              if(sgnphi0*phipr(i,k) .gt. phipi2) dlamsc = sign(pi,lampr(i,k)) - dlamsc
           end if
           lamdp(i,k) = lam(ii) + dlamsc
        end do
     else
        do i = 1,nlon
           ii = i + i1 - 1
           lamdp(i,k) = lam(ii) + lampr(i,k)
           phidp(i,k) = phi0    + phipr(i,k)
!
! If traj is over pole, STOP
!
           if (phidp(i,k) >= phi(j1+plat)) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
           if (phidp(i,k) < phi(j1-1   )) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
        end do
     end if
#if ( defined SPMD )
!
! If traj goes off-processor (out-of-bounds), STOP
!
        do i = 1,nlon
           if (phidp(i,k) >= phi(endlatex-nxpt) ) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
           if (phidp(i,k) < phi(beglatex+nxpt) ) then
              lstop = .true.
              ibad  = i
              kbad  = k
           end if
        end do
#endif
!
! Apply appropriate limiters
!
     do i = 1,nlon
        if(lamdp(i,k) .ge. twopi) lamdp(i,k) = lamdp(i,k) - twopi + 10.*epsilon(lamdp)
        if(lamdp(i,k) .lt.   0.0) lamdp(i,k) = lamdp(i,k) + twopi - 10.*epsilon(lamdp)
     end do
!
! Compute vertical departure points and apply appropriate limiters
!
     if (lvert) then
        botlim  = etamid(plev)*fac
        do i = 1,nlon
           sigdp(i,k) = etamid(k) + sigpr(i,k)
           if(sigdp(i,k) .lt. etamid(   1)) sigdp(i,k) = etamid(1)
           if(sigdp(i,k) .ge. etamid(plev)) sigdp(i,k) = botlim
        end do
     end if
  end do
!
! Test that the latitudinal extent of trajectory is NOT over the poles
!
  if (lstop) then
     write(6,*)'ERROR IN TRAJDP: ****** MODEL IS BLOWING UP *********'
     write(6,9000) ibad,kbad,lat
     call endrun
  end if
!
  return
!
! Formats
!
9000 format(//'Departure point out of bounds.  Parcel associated with longitude ' &
              ,i5,', level ',i5, ' and latitude ',i5/' is outside the model domain. ')
end subroutine trajdp

