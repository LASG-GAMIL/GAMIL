#include <misc.h>
#include <params.h>
subroutine sltwgts(limdrh  ,limdrv  ,lhrzwgt ,lvrtwgt ,kdim    , &
                   idp     ,jdp     ,kdp     ,lam     ,phi     , &
                   z       ,dphi    ,dz      ,lamdp   ,phidp   , &
                   sigdp   ,lbasiy  ,wiz     ,kkdp    ,xl      , &
                   xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
                   hl      ,hr      ,dhl     ,dhr     ,ys      , &
                   yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
                   hs      ,hn      ,dhs     ,dhn     ,rdphi   , &
                   wgt1z   ,wgt2z   ,wgt3z   ,wgt4z   ,hb      , &
                   ht      ,dhb     ,dht     ,rdz     ,zt      , &
                   zb      ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose: 
! Compute weights for SLT interpolation
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: sltwgts.F90,v 1.4.28.1 2002/06/15 13:48:31 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use srchutil
  implicit none

!------------------------------Arguments--------------------------------
!
  logical , intent(in)   :: limdrh              ! horizontal derivative limiter flag
  logical , intent(in)   :: limdrv              ! vertical   derivative limiter flag
  logical , intent(in)   :: lhrzwgt             ! flag to compute horizontal weights
  logical , intent(in)   :: lvrtwgt             ! flag to compute vertical   weights

  integer , intent(in)   :: kdim                ! vertical coordinate
  integer , intent(in)   :: idp   (plon,plev,4) ! index of x-coordinate of dep pt
  integer , intent(in)   :: jdp   (plon,plev)   ! index of y-coordinate of dep pt
  integer , intent(in)   :: kdp   (plon,plev)   ! index of z-coordinate of dep pt

  real(r8), intent(in)   :: lam   (plond,platd) ! longitude coordinates of model grid
  real(r8), intent(in)   :: phi   (platd)       ! latitude  coordinates of model grid
  real(r8), intent(in)   :: z     (kdim)        ! vertical  coordinates of model grid
  real(r8), intent(in)   :: dphi  (platd)       ! latitudinal grid increments
  real(r8), intent(in)   :: dz    (kdim)        ! vertical grid increments
  real(r8), intent(in)   :: lamdp (plon,plev)   ! x-coordinates of dep pts.
  real(r8), intent(in)   :: phidp (plon,plev)   ! y-coordinates of dep pt
  real(r8), intent(in)   :: sigdp (plon,plev)   ! z-coordinates of dep pt
  real(r8), intent(in)   :: lbasiy(4,2,platd)   ! y-interpolation weights for Lag. cubic
  real(r8), intent(in)   :: wiz   (4,2,kdim)    ! z-interpolation weights for Lag. cubic
  integer , intent(out)  :: kkdp  (plon,plev)   ! index of z-coordinate of dep pt (alt)

  real(r8), intent(out)  :: xl    (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(out)  :: xr    (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(out)  :: wgt1x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hl    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: hr    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: dhl   (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(out)  :: dhr   (plon,plev,4) ! weight for x-interpolants (Hermite)

  real(r8), intent(out)  :: ys    (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(out)  :: yn    (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(out)  :: wgt1y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hs    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: hn    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: dhs   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: dhn   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(out)  :: rdphi (plon,plev)   ! reciprocal of y-interval

  real(r8), intent(out)  :: wgt1z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt2z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt3z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: wgt4z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(out)  :: hb    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: ht    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: dhb   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: dht   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(out)  :: rdz   (plon,plev)   ! reciprocal of z-interval
  real(r8), intent(out)  :: zt    (plon,plev)   ! linear interpolation weight
  real(r8), intent(out)  :: zb    (plon,plev)   ! linear interpolation weight
  integer , intent(in)   :: nlon                ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i                 ! |
  integer ii,jj(plon,plev,4)! |
  integer k,n,j             ! |
  integer icount            ! |
  integer jdpval            ! |
  integer kdpval            ! | -- indices
  integer jmin              ! |
  integer jmax              ! |
  integer kdimm2            ! |
  integer nval              ! |
  integer indx (plond)      ! |
!
  real(r8) dx  (platd)      ! |
  real(r8) rdx (platd)      ! |
  real(r8) dyj              ! |
  real(r8) tmp1             ! |
  real(r8) tmp2             ! |
  real(r8) tmp3             ! |
  real(r8) tmp4             ! | -- tmp variables
  real(r8) dzk              ! |
  real(r8) denom1           ! |
  real(r8) denom2           ! |
  real(r8) denom3           ! |
  real(r8) denom4           ! |
  real(r8) coef12           ! |
  real(r8) coef34           ! |
!
!-----------------------------------------------------------------------
!
  denom1 = -1._r8/6._r8
  denom2 =  0.5_r8
  denom3 = -0.5_r8
  denom4 =  1._r8/6._r8
!
! HORIZONTAL weights
!
  if (lhrzwgt) then
!
! Determine N/S extent of all departure points in this latitude slice
!
     jmin =  1000000
     jmax = -1000000
     do k=1,plev
        do i=1,nlon
           if(jdp(i,k) .lt. jmin) jmin = jdp(i,k)
           if(jdp(i,k) .gt. jmax) jmax = jdp(i,k)
        end do
     end do
!
! Compute weights for x-direction
!
     do j = 1,platd
        dx (j) = lam(nxpt+2,j) - lam(nxpt+1,j)
        rdx(j) = 1./dx(j)
     end do
!
     do n=1,4
        do k=1,plev
           do i=1,nlon
              jj(i,k,n) = jdp(i,k) - 2 + n
              xl(i,k,n) = (lam(idp(i,k,n)+1,jj(i,k,n)) - lamdp(i,k))*rdx(jj(i,k,n))
              xr(i,k,n) = 1. - xl(i,k,n)
           end do
        end do
     end do
     do n=2,3
!!#ifndef HADVTEST
!!        if (limdrh) then
!!#endif
           do k=1,plev
              do i=1,nlon
                 hl (i,k,n)   = ( 3.0 - 2.0*xl(i,k,n) )*xl(i,k,n)**2
                 hr (i,k,n)   = ( 3.0 - 2.0*xr(i,k,n) )*xr(i,k,n)**2
                 dhl(i,k,n)   = -dx(jj(i,k,n))*( xl(i,k,n) - 1. )*xl(i,k,n)**2
                 dhr(i,k,n)   =  dx(jj(i,k,n))*( xr(i,k,n) - 1. )*xr(i,k,n)**2
              end do
           end do
!!#ifndef HADVTEST
!!        else
!!#endif
           do k=1,plev
              do i=1,nlon
                 tmp1         =  xr(i,k,n) + 1.
                 tmp4         =  xr(i,k,n) - 2.
                 coef12       = -xl(i,k,n)*tmp4
                 coef34       =  xr(i,k,n)*tmp1
                 wgt1x(i,k,n) =  denom1*coef12*xr(i,k,n)
                 wgt2x(i,k,n) =  denom2*coef12*tmp1
                 wgt3x(i,k,n) =  denom3*coef34*tmp4
                 wgt4x(i,k,n) = -denom4*coef34*xl(i,k,n)
              end do
           end do
!!#ifndef HADVTEST
!!        endif
!!#endif
     end do
!
! Compute weights for y-direction
!
     icount = 0
     do jdpval=jmin,jmax
        do k=1,plev
           call wheneq(nlon    ,jdp(1,k),1       ,jdpval  ,indx    ,nval    )
           icount = icount + nval
           dyj    = dphi(jdpval)
           do ii = 1,nval
              i = indx(ii)
              ys(i,k) = ( phi(jdpval+1) - phidp(i,k) )/dyj
              yn(i,k) = 1. - ys(i,k)
           end do
!!#ifndef HADVTEST
!!           if (limdrh) then
!!#endif
              do ii = 1,nval
                 i = indx(ii)
                 rdphi(i,k) = 1./dyj
                 hs   (i,k) = ( 3.0 - 2.0*ys(i,k) )*ys(i,k)**2
                 hn   (i,k) = ( 3.0 - 2.0*yn(i,k) )*yn(i,k)**2
                 dhs  (i,k) = -dyj*( ys(i,k) - 1. )*ys(i,k)**2
                 dhn  (i,k) =  dyj*( yn(i,k) - 1. )*yn(i,k)**2
              end do
!!#ifndef HADVTEST
!!           else
!!#endif
              do ii = 1,nval
                 i          = indx(ii)
                 tmp1       = phidp(i,k) - lbasiy(1,1,jdpval)
                 tmp2       = phidp(i,k) - lbasiy(2,1,jdpval)
                 tmp3       = phidp(i,k) - lbasiy(3,1,jdpval)
                 tmp4       = phidp(i,k) - lbasiy(4,1,jdpval)
                 coef12     = tmp3*tmp4   
                 coef34     = tmp1*tmp2   
                 wgt1y(i,k) = coef12*tmp2*lbasiy(1,2,jdpval)
                 wgt2y(i,k) = coef12*tmp1*lbasiy(2,2,jdpval)
                 wgt3y(i,k) = coef34*tmp4*lbasiy(3,2,jdpval)
                 wgt4y(i,k) = coef34*tmp3*lbasiy(4,2,jdpval)
              end do
!!#ifndef HADVTEST
!!           endif
!!#endif
        end do
     end do
     if (icount.ne.nlon*plev) then
        write(6,*)'SLTWGTS:  Did not complete computations for all departure points'
        call endrun
     end if
  end if
!
! VERTICAL weights
!
  if (lvrtwgt) then
!
! Limit kdp to between "2" and "kdim-2" when computing weights and
! derivatives.
!
     kdimm2 = kdim - 2
     do k=1,plev
        do i=1,nlon
           kkdp(i,k) = min0( kdimm2,max0( 2,kdp(i,k) ) )
           dzk       = dz(kdp(i,k))
           rdz(i,k)  = 1./dzk
           zt(i,k)   = ( z  (kdp(i,k)+1) - sigdp(i,k) )/dzk
           zb(i,k)   = 1. - zt(i,k)
           ht (i,k)  = ( 3.0 - 2.0*zt(i,k) )*zt(i,k)**2
           hb (i,k)  = ( 3.0 - 2.0*zb(i,k) )*zb(i,k)**2
           dht(i,k)  = -dzk*( zt(i,k) - 1. )*zt(i,k)**2
           dhb(i,k)  =  dzk*( zb(i,k) - 1. )*zb(i,k)**2
        end do
     end do
!
     if(.not. limdrv) then
        icount = 0
        do kdpval=2,kdimm2
           do k=1,plev
              call wheneq(nlon    ,kkdp(1,k),1       ,kdpval  ,indx    ,nval    )
              icount = icount + nval
              do ii = 1,nval
                 i          = indx(ii)
                 tmp1       = sigdp(i,k) -  wiz(1,1,kdpval)
                 tmp2       = sigdp(i,k) -  wiz(2,1,kdpval)
                 tmp3       = sigdp(i,k) -  wiz(3,1,kdpval)
                 tmp4       = sigdp(i,k) -  wiz(4,1,kdpval)
                 coef12     = tmp3*tmp4   
                 coef34     = tmp1*tmp2   
                 wgt1z(i,k) = coef12*tmp2*wiz(1,2,kdpval)
                 wgt2z(i,k) = coef12*tmp1*wiz(2,2,kdpval)
                 wgt3z(i,k) = coef34*tmp4*wiz(3,2,kdpval)
                 wgt4z(i,k) = coef34*tmp3*wiz(4,2,kdpval)
              end do
           end do
        end do
        if (icount.ne.nlon*plev) then
           write(6,*)'SLTWGTS:  Did not complete computations for all departure points'
           call endrun
        end if
     end if
  end if
!
  return
end subroutine sltwgts

