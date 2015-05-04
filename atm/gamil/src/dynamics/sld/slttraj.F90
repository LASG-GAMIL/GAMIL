#include <misc.h>
#include <params.h>
subroutine slttraj(pmap    ,jcen    ,lat     ,ztodt   ,ra      , &
                   iter    ,lam     ,phi     ,dphi    ,etamid  , &
                   etaint  ,detam   ,detai   ,lbasiy  ,lbasiz  , &
                   lbassi  ,kdpmpf  ,kdpmph  ,idp     ,jdp     , &
                   kdp     ,kkdp    ,xl      ,xr      ,wgt1x   , &
                   wgt2x   ,wgt3x   ,wgt4x   ,hl      ,hr      , &
                   dhl     ,dhr     ,ys      ,yn      ,wgt1y   , &
                   wgt2y   ,wgt3y   ,wgt4y   ,hs      ,hn      , &
                   dhs     ,dhn     ,rdphi   ,wgt1z   ,wgt2z   , &
                   wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
                   dht     ,rdz     ,lampr   ,phipr   ,upr     , &
                   vpr     ,lamdp   ,phidp   ,sigdp   ,u3      , &
                   v3      ,u3sld   ,v3sld   ,etadot  ,n3      , &
                   n3m1    ,dlam    ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Determine trajectory departure point for each arrival point.
! Assume that the first guess for the dep point coordinates are the
! arrival point coordinates
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: slttraj.F90,v 1.4.8.2 2002/06/15 13:48:31 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use comslt
  use prognostics, only: ptimelevels
#if (!defined CRAY)
  use srchutil
#endif
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pmap                ! artificial vert grid dim.
  integer , intent(in)   :: jcen                ! index of lat slice(extend)
  integer , intent(in)   :: lat                 ! index of lat slice (model)
  real(r8), intent(in)   :: ztodt               ! time step (seconds)
  real(r8), intent(in)   :: ra                  ! 1./(radius of earth)
  integer , intent(in)   :: iter                ! iteration count
  real(r8), intent(in)   :: lam   (plond,platd) ! long. coord of model grid
  real(r8), intent(in)   :: phi   (platd)       ! lat.  coord of model grid
  real(r8), intent(in)   :: dphi  (platd)       ! increment between lats.
  real(r8), intent(in)   :: etamid(plev)        ! vertical full levels
  real(r8), intent(in)   :: etaint(plevp)       ! vertical half levels
  real(r8), intent(in)   :: detam (plev)        ! increment between full levs
  real(r8), intent(in)   :: detai (plevp)       ! increment between half levs
  real(r8), intent(in)   :: lbasiy(4,2,platd)   ! lat interp wts(lagrng)
  real(r8), intent(in)   :: lbasiz(4,2,plev)    ! vert interp wghts (lagrng)
  real(r8), intent(in)   :: lbassi(4,2,plevp)   ! vert interp wghts (lagrng)

  integer , intent(in)   :: kdpmpf(pmap)        ! artificial vert grid indices
  integer , intent(in)   :: kdpmph(pmap)        ! artificial vert grid indices
  integer , intent(out)  :: idp   (plon,plev,4) ! zonal      dep point index
  integer , intent(out)  :: jdp   (plon,plev)   ! meridional dep point index
  integer , intent(out)  :: kdp   (plon,plev)   ! vertical   dep point index
  integer , intent(in)   :: kkdp  (plon,plev)   ! index of z-coordinate of dep pt (alt)

  real(r8), intent(out)  :: xl    (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(out)  :: xr    (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(in)   :: wgt1x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hl    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: hr    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhl   (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhr   (plon,plev,4) ! weight for x-interpolants (Hermite)

  real(r8), intent(out)  :: ys    (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(out)  :: yn    (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(in)   :: wgt1y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hs    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: hn    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhs   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhn   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: rdphi (plon,plev)   ! reciprocal of y-interval

  real(r8), intent(in)   :: wgt1z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hb    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: ht    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dhb   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dht   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: rdz   (plon,plev)   ! reciprocal of z-interval

  real(r8), intent(out)  :: lampr (plon,plev)   ! trajectory increment (x-direction)
  real(r8), intent(out)  :: phipr (plon,plev)   ! trajectory increment (y-direction)
  real(r8), intent(out)  :: upr   (plon,plev)   ! interpolated u field (local geodesic)
  real(r8), intent(out)  :: vpr   (plon,plev)   ! interpolated v field (local geodesic)
  real(r8), intent(out)  :: lamdp (plon,plev)   ! zonal      departure pt. coord.
  real(r8), intent(inout):: phidp (plon,plev)   ! meridional departure pt. coord.
  real(r8), intent(inout):: sigdp (plon,plev)   ! vertical   departure pt. coord.
  real(r8), intent(in)   :: u3    (plond, plev, beglatex:endlatex,ptimelevels) ! u wind component
  real(r8), intent(in)   :: v3    (plond, plev, beglatex:endlatex,ptimelevels) ! v wind component
  real(r8), intent(in)   :: u3sld (plond,plev ,beglatex:endlatex)    ! u3 inpt for SLD int
  real(r8), intent(in)   :: v3sld (plond,plev ,beglatex:endlatex)    ! v3 inpt for SLD int

  real(r8), intent(in)   :: etadot(plond,plevp,beglatex:endlatex,ptimelevels)  ! Vertical motion
  integer , intent(in)   :: n3                                       ! time indicies
  integer , intent(in)   :: n3m1                                     ! time indicies
  real(r8), intent(in)   :: dlam  (platd)                            ! d-long for each lat
  integer , intent(in)   :: nlon                                     ! # of longitudes
#if (defined CRAY)
  integer,external :: ismin
#endif
!
!---------------------------Local workspace-----------------------------
!
  real(r8) fac                  ! 1-eps
  real(r8) tmp                  ! temp variable
  real(r8) ump    (plon,plev)   ! interpolated u field
  real(r8) vmp    (plon,plev)   ! interpolated v field
  real(r8) wmp    (plon,plev)   ! interpolated w field
  real(r8) wmpa   (plon,plev)   ! interpolated w field (arrival level)
  real(r8) sigpr  (plon,plev)   ! trajectory increment (vertical)
  real(r8) etamin               ! minimum eta
  real(r8) etamax               ! maximum eta
  real(r8) delsig (plev)        ! dist between dep point and model levs
  real(r8) zt     (plon,plev)   ! top vertical interpolation weight 
  real(r8) zb     (plon,plev)   ! bot vertical interpolation weight 
  real(r8) zttmp                ! top vertical interpolation weight 
  real(r8) zbtmp                ! bot vertical interpolation weight 
  real(r8) rdx    (platd)       ! reciprocal of del-x

  logical locgeo                ! local geodesic flag
  logical lvert                 ! flag to compute vertical trajectory
  logical larrival(plon)        ! flag to indicate whether or not to
!                               ! place the alternative dep point at the
!                               ! arrival point.

  integer kount                 ! index counter
  integer kkmin                 ! index of minimum "delsig"
  integer i                     ! |
  integer ii                    ! |
  integer j                     ! |
  integer jj2                   ! | - indices
  integer jj3                   ! |
  integer k                     ! |
  integer kk                    ! |
  integer n                     ! |
!
!-----------------------------------------------------------------------
!
  fac     = 1. - 10.*epsilon (fac)
  locgeo  = .false.
  if(abs(phi(jcen)) .ge. phigs) locgeo = .true.
!
! Interpolate arrival and departure vertical velocities
!
  do k = 1,plev
     do i = 1,nlon
        zttmp = ( etaint(k+1) - sigdp(i,k) )/detai(k)
        zbtmp = 1. - zttmp
        wmpa(i,k) = etadot(i1+i-1,k  ,jcen,n3  )*zttmp &
                  + etadot(i1+i-1,k+1,jcen,n3  )*zbtmp
        wmp (i,k) = etadot(i1+i-1,k  ,jcen,n3m1)*zttmp &
                  + etadot(i1+i-1,k+1,jcen,n3m1)*zbtmp
     end do
  end do
!
! Set up computation of trajectory
!
  do k = 1,plev
     do i = 1,nlon
        ii = i + i1 - 1
!
! Place u/v on unit sphere
!
        ump(i,k) = 0.5*(u3sld(ii,k,jcen) + u3(ii,k,jcen,n3))*ra
        vmp(i,k) = 0.5*(v3sld(ii,k,jcen) + v3(ii,k,jcen,n3))*ra
        wmp(i,k) = 0.5*(wmp(i,k)         + wmpa(i,k)       )
!
! Estimate departure point of parcel trajectory.
!
        lampr(i,k) = -ztodt*ump(i,k)
        phipr(i,k) = -ztodt*vmp(i,k)
        sigpr(i,k) = -ztodt*wmp(i,k)
        if (.not. locgeo) lampr(i,k) = lampr(i,k)/cos( phidp(i,k) )
!
! Initialize winds for use in g.c. calculations.
!
        if (locgeo) then
           upr(i,k) = ump(i,k)
           vpr(i,k) = vmp(i,k)
        end if
     end do
  end do
!
! Estimate initial trajectory departure points
!
  lvert = .true.
  call trajdp(lat     ,jcen    ,ztodt   ,locgeo  ,lvert   , &
              lam(1,jcen),phi  ,etamid  ,upr     ,vpr     , &
              lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
              sigdp   ,nlon    )
!
! Loop over departure point iterates.
!
  do n = 1,iter
!
! Determine departure point indicies and interpolate u,v,w
!
     call bandij (dlam    ,phi     ,lamdp   ,phidp   ,idp     , &
                  jdp     ,nlon    )
     call kdpfnd (plev    ,pmap    ,etamid  ,sigdp   ,kdpmpf  , &
                  kdp     ,nlon    )
!
! Compute weights for x,y,z dimensions
!
     do j = 1,platd
        rdx(j) = 1./(lam(nxpt+2,j) - lam(nxpt+1,j))
     end do
     do k = 1,plev
        do i = 1,nlon
           jj2 = jdp(i,k)
           jj3 = jdp(i,k) + 1
           xl(i,k,2) = (lam(idp(i,k,2)+1,jj2) - lamdp(i,k))*rdx(jj2)
           xl(i,k,3) = (lam(idp(i,k,3)+1,jj3) - lamdp(i,k))*rdx(jj3)
           xr(i,k,2) = 1. - xl(i,k,2)
           xr(i,k,3) = 1. - xl(i,k,3)
           ys(i,k)   = ( phi(jdp(i,k)+1) - phidp(i,k) )/dphi(jdp(i,k))
           yn(i,k)   = 1. - ys(i,k)
           zt(i,k)   = (etamid(kdp(i,k)+1)-sigdp(i,k))/detam(kdp(i,k))
           zb(i,k)   = 1. - zt(i,k)
        end do
     end do
!
     call sltlinint(plev    ,1       ,u3sld(1,1,beglatex),xl      ,xr      , &
                    ys      ,yn      ,zt                 ,zb      ,idp     , &
                    jdp     ,kdp     ,ump                ,nlon    )
     call sltlinint(plev    ,1       ,v3sld(1,1,beglatex),xl      ,xr      , &
                    ys      ,yn      ,zt                 ,zb      ,idp     , &
                    jdp     ,kdp     ,vmp                ,nlon    )
     call kdpfnd (plevp   ,pmap    ,etaint  ,sigdp   ,kdpmph  , &
                  kdp     ,nlon    )
     do k = 1,plev
        do i = 1,nlon
           zt(i,k)   = (etaint(kdp(i,k)+1)-sigdp(i,k))/detai(kdp(i,k))
           zb(i,k)   = 1. - zt(i,k)
        end do
     end do
     call sltlinint(plevp   ,1       ,etadot(1,1,beglatex,n3m1),xl      ,xr      , &
                    ys      ,yn      ,zt                       ,zb      ,idp     , &
                    jdp     ,kdp     ,wmp                      ,nlon    )
!
! Add arrival point velocities (n) to interpolated velocities
!
     do k = 1,plev
        do i = 1,nlon
           ii = i + i1 - 1
           ump(i,k) = 0.5*(ump(i,k) + u3(ii,k,jcen,n3))
           vmp(i,k) = 0.5*(vmp(i,k) + v3(ii,k,jcen,n3))
           wmp(i,k) = 0.5*(wmp(i,k) + wmpa(i,k)       )
        end do
     end do
!
! Compute new trajectory departure points
!
     lvert = .true.
     call depinc (jcen    ,ztodt   ,ra      ,locgeo  ,lvert   , &
                  lam(1,jcen),phi  ,ump     ,vmp     ,wmp     , &
                  upr     ,vpr     ,lamdp   ,phidp   ,lampr   , &
                  phipr   ,sigpr   ,nlon    )
     call trajdp (lat     ,jcen    ,ztodt   ,locgeo  ,lvert   , &
                  lam(1,jcen),phi  ,etamid  ,upr     ,vpr     , &
                  lampr   ,phipr   ,sigpr   ,lamdp   ,phidp   , &
                  sigdp   ,nlon    )
  end do
!
! Compile departure point binning statistics
!
  do k = 1,plev
     if(k .eq. 1) then
        etamin = etamid(k  )
        etamax = etamid(k+1)
     elseif(k .eq. plev) then
        etamin = etamid(k-1)
        etamax = etamid(k  )
     else
        etamin = etamid(k-1)
        etamax = etamid(k+1)
     end if
     do i = 1,nlon
        larrival(i) = sigdp(i,k) .ge. etamin .and. sigdp(i,k) .le. etamax
     end do
!
! If:  Departure point is within one grid interval from the arrival
! point, bin it in the ARRIVAL point bin
!
     tmp = 0.
     do i = 1,nlon
        if (larrival(i)) then
           tmp           = tmp + 1.
        end if
     end do
     levkntl(k,k,lat) = levkntl(k,k,lat) + tmp
!
! Else:  find departure point bin
!
     kount = 0
     do i=1,nlon
        if (larrival(i)) kount = kount + 1
     end do
     if (kount .ne. nlon) then
        do i = 1,nlon
           if (.not. larrival(i)) then
              delsig(1) = 1.e+35
              do kk = 2,plev
                 delsig(kk) = abs( sigdp(i,k) - etaint(kk) )
              end do
              kkmin = ismin( plev,delsig,1 )
              if(kkmin .gt. k) kkmin = kkmin-1
              levkntl(kkmin,k,lat) = levkntl(kkmin,k,lat) + 1.
           end if
        end do
     end if
  end do
!
  return
end subroutine slttraj

