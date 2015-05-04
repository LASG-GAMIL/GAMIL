#include <misc.h>
#include <params.h>

subroutine scanslt (ztodt   ,lat     ,dtr     ,iter    ,pmap    , &
                    kdpmpf  ,kdpmph  ,lam     ,phi     ,dphi    , &
                    lbasdy  ,lbasdz  ,lbasiy  ,lbasiz  ,lbassi  , &
                    detam   ,detai   ,dlam    ,cwava   ,etamid  , &
                    etaint  ,grfu    ,grfv    ,grlps1  ,grlps2  , &
                    grt1    ,grt2    ,grq1    ,grq2    ,grfu1   , &
                    grfu2   ,grfv1   ,grfv2   ,ps      ,u3      , &
                    v3      ,t3      ,q3      ,lnpssld ,prhssld , &
                    tarrsld ,parrsld ,n3      ,n3m1    ,u3sld   , &
                    v3sld   ,etadot  ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Interpolate terms for semi-lagrangian transport and SLD dynamics.
! One latitude slice only
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: scanslt.F90,v 1.14.2.3 2002/06/21 05:37:40 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,      only: plon, plond, plev, plevp, plat, platd, beglat, endlat, beglatex, &
                         endlatex, plndlv, i1, j1
  use constituents,only: pcnst, pnats
  use comslt,      only: qfcst, gamma, hw1lat
  use rgrid,       only: nmmax
  use pspect,      only: pmmax
  use commap,      only: w, clat, t0
  use prognostics, only: ptimelevels
  use physconst,   only: cappa
  use dynconst,    only: ra

  implicit none

#include <comctl.h>
#include <comfft.h>
#include <comhyb.h>

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                 ! twice the time step unless nstep = 0
  integer , intent(in)   :: lat                   ! latitude index
  real(r8), intent(in)   :: dtr                   ! 1/dt
  integer , intent(in)   :: iter                  ! number of iterations for trajectory
  integer , intent(in)   :: pmap                  ! dimension of artificial array
!                                                 ! used to locate vertical interval
!                                                 ! in which departure point falls
  integer , intent(in)   :: kdpmpf  (pmap)        ! mapping from artificial array to
!                                                 ! model levels
  integer , intent(in)   :: kdpmph  (pmap)        ! mapping from artificial array to
!                                                 ! model interfaces
  real(r8), intent(in)   :: lam     (plond,platd) ! longitude coordinates of model grid
  real(r8), intent(in)   :: phi     (platd)       ! latitude  coordinates of model grid
  real(r8), intent(in)   :: dphi    (platd)       ! latitudinal grid increments
  real(r8), intent(in)   :: lbasdy  (4,2,platd)   ! basis functions for lat deriv est.
  real(r8), intent(in)   :: lbasdz  (4,2,plev)    ! basis functions for vert deriv est.
  real(r8), intent(in)   :: lbasiy  (4,2,platd)   ! basis functions for Lagrange interp
  real(r8), intent(in)   :: lbasiz  (4,2,plev)    ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: lbassi  (4,2,plevp)   ! Lagrange cubic interp wghts (vert)
  real(r8), intent(in)   :: detam   (plev)        ! delta eta at levels
  real(r8), intent(in)   :: detai   (plevp)       ! delta eta at interfaces
  real(r8), intent(in)   :: dlam    (platd)       ! longitudinal grid increment
  real(r8), intent(in)   :: cwava   (plat)        ! weight for global water vapor int.
  real(r8), intent(in)   :: etamid  (plev)        ! eta at levels
  real(r8), intent(in)   :: etaint  (plevp)       ! eta at interfaces
  real(r8), intent(in)   :: grfu    (plond,plev,beglat:endlat) ! nonlinear term - u momentum eqn
  real(r8), intent(in)   :: grfv    (plond,plev,beglat:endlat) ! nonlinear term - v momentum eqn
#if ( defined PVP )
  real(r8), intent(out)  :: grlps1(2*pmmax     ,plat/2) ! ------------------------------
  real(r8), intent(out)  :: grlps2(2*pmmax     ,plat/2) ! |
  real(r8), intent(out)  :: grt1  (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grt2  (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grq1  (2*pmmax,plev,plat/2) ! |- see quad for definitions
  real(r8), intent(out)  :: grq2  (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grfu1 (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grfu2 (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grfv1 (2*pmmax,plev,plat/2) ! |
  real(r8), intent(out)  :: grfv2 (2*pmmax,plev,plat/2) ! ------------------------------
#else
  real(r8), intent(out)  :: grlps1(     2*pmmax,plat/2) ! ------------------------------
  real(r8), intent(out)  :: grlps2(     2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grt1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grt2  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grq1  (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grq2  (plev,2*pmmax,plat/2) ! |- see quad for definitions
  real(r8), intent(out)  :: grfu1 (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grfu2 (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grfv1 (plev,2*pmmax,plat/2) ! |
  real(r8), intent(out)  :: grfv2 (plev,2*pmmax,plat/2) ! ------------------------------
#endif
  real(r8), intent(in)   :: ps      (plond,beglat:endlat,ptimelevels)
  real(r8), intent(in)   :: u3      (plond,plev,beglatex:endlatex,ptimelevels) ! u-wind com
  real(r8), intent(in)   :: v3      (plond,plev,beglatex:endlatex,ptimelevels) ! v-wind comp
  real(r8), intent(in)   :: t3      (plond,plev,beglatex:endlatex,ptimelevels) ! temperature
  real(r8), intent(in)   :: q3      (plond,plev,pcnst+pnats,beglatex:endlatex,ptimelevels)
!                                                                    ! q and const
  real(r8), intent(in)   :: lnpssld (plond,plev,beglatex:endlatex)   ! RHS Ps term for SLD
  real(r8), intent(in)   :: prhssld (plond,plev,beglatex:endlatex)   ! RHS Ps term for SLD
  real(r8), intent(in)   :: tarrsld (plond,plev,beglatex:endlatex)   ! T  at arr. pt.(SLD)
  real(r8), intent(inout):: parrsld (plond,plev,beglatex:endlatex)   ! Ps at arr. pt.(SLD)
  integer , intent(in)   :: n3                                       ! time index
  integer , intent(in)   :: n3m1                                     ! time index
  real(r8), intent(in)   :: u3sld   (plond,plev ,beglatex:endlatex)  ! u3 inpt for SLD int
  real(r8), intent(in)   :: v3sld   (plond,plev ,beglatex:endlatex)  ! v3 inpt for SLD int
  real(r8), intent(in)   :: etadot  (plond,plevp,beglatex:endlatex,ptimelevels)! Vertical motion

  integer , intent(in)   :: nlon                                     ! # of longitudes
!
!---------------------------Local workspace-----------------------------
!
  integer i                    ! index
  integer k                    ! index
  integer l                    ! index
  integer m                    ! constituent index
  integer inc                  ! increment for fft991
  integer ntr                  ! number of fft's to perform
  integer isign                ! flag indicates fft transform directn
  integer irow                 ! N/S latitude pair index
  integer jcen                 ! lat index (extended grid)
!                              ! of forecast
  real(r8) qtmp (plond,plev,beglatex:endlatex)  ! Temporary for q array
  real(r8) fdp  (plon,plev,2)  ! interpolant
  real(r8) pmid (plond,plev)   ! pressure at model levels
  real(r8) pint (plond,plevp)  ! pressure at interfaces
  real(r8) pdel (plond,plev)   ! pressure difference between
  real(r8) lamdp(plon,plev)    ! x-coord of dep pt
  real(r8) phidp(plon,plev)    ! y-coord of dep pt
  real(r8) sigdp(plon,plev)    ! z-coord of dep pt

  integer idp   (plon,plev,4)  ! zonal      dep point index
  integer jdp   (plon,plev)    ! meridional dep point index
  integer kdp   (plon,plev)    ! vertical   dep point index
  integer kkdp  (plon,plev)    ! index of z-coordinate of dep pt (alt)

  real(r8) xl   (plon,plev,4)  ! weight for x-interpolants (left)
  real(r8) xr   (plon,plev,4)  ! weight for x-interpolants (right)
  real(r8) wgt1x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt2x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt3x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) wgt4x(plon,plev,4)  ! weight for x-interpolants (Lag Cubic)
  real(r8) hl   (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) hr   (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) dhl  (plon,plev,4)  ! weight for x-interpolants (Hermite)
  real(r8) dhr  (plon,plev,4)  ! weight for x-interpolants (Hermite)

  real(r8) ys   (plon,plev)    ! weight for y-interpolants (south)
  real(r8) yn   (plon,plev)    ! weight for y-interpolants (north)
  real(r8) wgt1y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt2y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt3y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) wgt4y(plon,plev)    ! weight for y-interpolants (Lag Cubic)
  real(r8) hs   (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) hn   (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) dhs  (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) dhn  (plon,plev)    ! weight for y-interpolants (Hermite)
  real(r8) rdphi(plon,plev)    ! reciprocal of y-interval

  real(r8) wgt1z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt2z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt3z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) wgt4z(plon,plev)    ! weight for z-interpolants (Lag Cubic)
  real(r8) hb   (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) ht   (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) dhb  (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) dht  (plon,plev)    ! weight for z-interpolants (Hermite)
  real(r8) rdz  (plon,plev)    ! reciprocal of z-interval
  real(r8) zt   (plon,plev)    ! top vertical interpolation weight 
  real(r8) zb   (plon,plev)    ! bot vertical interpolation weight 

  real(r8) lampr(plon,plev)    ! trajectory increment (x-direction)
  real(r8) phipr(plon,plev)    ! trajectory increment (y-direction)
  real(r8) upr  (plon,plev)    ! interpolated u field (local geodesic)
  real(r8) vpr  (plon,plev)    ! interpolated v field (local geodesic)

#if ( ! defined USEFFTLIB )
  real(r8) work((plon+1)*5*plev)   ! workspace array for fft991
#else
  real(r8) work((plon+1)*pcray)    ! workspace array for fft991
#endif
  real(r8) xnlin(plndlv*4 + plond) ! non-linear terms (equivalence 
!                                  ! region for following arrays to 
!                                  ! optimize fft performance)
  real(r8) grfulat(plond,plev)     ! non-linear terms for u-momentum 
  real(r8) grfvlat(plond,plev)     ! non-linear terms for u-momentum 
  real(r8) grtlat (plond,plev)     ! RHS of T-eqn
  real(r8) grqlat (plond,plev)     ! q
  real(r8) grpslat(plond)          ! RHS of Ps-eqn
!
! The following equivalences are for optimal fft performance
!
#if ( defined SUNOS )
!
! Stupid Sun compiler hoop-jumping again
!
  save grfulat, grfvlat, grtlat, grqlat, grpslat, xnlin
#endif
  equivalence   (grfulat,xnlin(1))
  equivalence   (grfvlat,xnlin(1+1*plndlv))
  equivalence   (grtlat ,xnlin(1+2*plndlv))
  equivalence   (grqlat ,xnlin(1+3*plndlv))
  equivalence   (grpslat,xnlin(1+4*plndlv))
!
  real(r8) pd   (plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pdsum(plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pd1  (plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pdsm1(plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pa   (plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) pasum(plond)            ! RHS term for Ps and (1/ps)etadot(dp/deta) 
  real(r8) coslat                  ! cos(latitude)
  real(r8) tmp1                    ! temp space
!
  logical limdrh                   ! horizontal derivative limiter flag
  logical limdrv                   ! vertical   derivative limiter flag
  logical lhrzint                  ! horizontal interp flag
  logical lvrtint                  ! vertical   interp flag
  logical lhrzwgt                  ! flag to compute horizontal weights
  logical lvrtwgt                  ! flag to compute vertical   weights
!
!-----------------------------------------------------------------------
!
  if(lat.le.plat/2) then
     irow = lat
  else
     irow = plat + 1 - lat
  end if
  jcen = j1 - 1 + lat
  coslat = cos(clat(lat))
!
! Initial guess for trajectory midpoints in spherical coords.
! Use arrival points as initial guess for trajectory midpoints.
!
  do k=1,plev
     do i=1,nlon
        phidp(i,k) = clat(lat)
        sigdp(i,k) = etamid(k)
     end do
  end do
!        
! Offset bottom level departure point first guess by epsilon
!
  do i = 1,nlon
     sigdp(i,plev) = sigdp(i,plev)*(1. - 10.*epsilon(sigdp))
  end do
!
! Loop through latitudes producing departure point calculation
!
  call slttraj(pmap    ,jcen    ,lat     ,ztodt   ,ra      , &
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
! Calculate mass of moisture in field being advected by slt.
!
  call plevs0(nlon    ,plond   ,plev    ,ps(1,lat,n3),pint    ,pmid    ,pdel)
  call qmassa(cwava(lat)  ,w(irow) ,q3(i1,1,1,jcen,n3),pdel    ,hw1lat(1,lat), &
              nlon)
!
! Compute constituent forecast
!
  lhrzwgt = .true.
  lvrtwgt = .true.
  lhrzint = .true.
  lvrtint = .true.
  limdrh  = .true.
  limdrv  = .true.
  call bandij (dlam    ,phi     ,lamdp   ,phidp   ,idp     , &
               jdp     ,nlon    )
  call kdpfnd (plev    ,pmap    ,etamid  ,sigdp   ,kdpmpf  , &
               kdp     ,nlon    )
  call sltwgts(limdrh  ,limdrv  ,lhrzwgt ,lvrtwgt ,plev    , &
               idp     ,jdp     ,kdp     ,lam     ,phi     , &
               etamid  ,dphi    ,detam   ,lamdp   ,phidp   , &
               sigdp   ,lbasiy  ,lbasiz  ,kkdp    ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,rdphi   , &
               wgt1z   ,wgt2z   ,wgt3z   ,wgt4z   ,hb      , &
               ht      ,dhb     ,dht     ,rdz     ,zt      , &
               zb      ,nlon    )
  do m = 1,pcnst
     qtmp(:plond,:plev,beglatex:endlatex) = q3(:plond,:plev,m,beglatex:endlatex,n3)
     call sltint (plev    ,jcen    ,qtmp    ,lam     ,rdphi   , &
                  rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
                  wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
                  hr      ,dhl     ,dhr     ,ys      ,yn      , &
                  wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
                  hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
                  wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
                  dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
                  lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp(1,1,1), &
                  nlon  )
!
! Pass q-interpolants into 3-D array
!
     do k = 1,plev
        do i = 1,nlon
           qfcst(i1-1+i,k,m,lat) = fdp(i,k,1)
        end do
        if(m .eq. 1) then
           do i = 1,nlon
              grqlat(i,k) = fdp(i,k,1)
           end do
        endif
     end do
  end do
!
! Accumulate P-interpolants into T equation
!
  do i = 1,nlon
     grpslat(i) = 0.
     pasum  (i) = 0.
     pa     (i) = 0.
  end do
  do k = 1,plev
     do i = 1,nlon
        grtlat(i,k) = tarrsld (i1-1+i,k,jcen)
     end do
  end do
  do k = 1,plev
     do l = 1,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + parrsld(i1-1+i,l,jcen)*gamma(l,k)
        end do
     end do
  end do
!
! Accumulate Ps interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + parrsld(i1-1+i,k,jcen)
        pasum  (i) = pasum  (i) + parrsld(i1-1+i,k,jcen)
     end do
  end do
!
! Compute first part of (1/ps)etadot(dp/deta)
!
  do k = 1,plev-1
     do i = 1,nlon
        pa(i) = pa(i) + parrsld(i1-1+i,k,jcen)
        parrsld(i1-1+i,k,jcen) = pa(i)
     end do
!
     if(k.ge.nprlev) then
        do i = 1,nlon
           parrsld(i1-1+i,k,jcen) = parrsld(i1-1+i,k,jcen) - hybi(k+1)*( pasum(i) )
        end do
     end if
  end do
!
! Compute U, V interpolants:  Non-monotonic, 3-D interpolation
!
  limdrh  = .false.
  limdrv  = .true.
  lhrzint = .true.
  lvrtint = .true.
  call sltint (plev    ,jcen   ,u3(1,1,beglatex,n3m1)  ,lam , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp(1,1,1),nlon  )
  call sltint (plev    ,jcen    ,v3(1,1,beglatex,n3m1)  ,lam     , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp(1,1,2),nlon  )
!
! Evaluate last half of grfu and grfv (Nu,Nv)
!
  call nunv1(lam(i1,jcen) ,phi(jcen),lamdp        ,phidp              ,fdp(1,1,1), &
             fdp(1,1,2)   ,coslat   ,grfu(1,1,lat),grfv(1,1,lat)      ,grfulat   , &
             grfvlat      ,nlon)
!
! Compute T interpolants:  Non-monotonic, 3-D interpolation
!
  limdrh  = .false.
  limdrv  = .true.
  lhrzint = .true.
  lvrtint = .true.
  call sltint (plev    ,jcen    ,t3(1,1,beglatex,n3m1)  ,lam     , &
               rdphi   ,rdz     ,lbasdy  ,lbasdz  ,xl      , &
               xr      ,wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   , &
               hl      ,hr      ,dhl     ,dhr     ,ys      , &
               yn      ,wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   , &
               hs      ,hn      ,dhs     ,dhn     ,wgt1z   , &
               wgt2z   ,wgt3z   ,wgt4z   ,hb      ,ht      , &
               dhb     ,dht     ,idp     ,jdp     ,kdp     , &
               kkdp    ,lhrzint ,lvrtint ,limdrh  ,limdrv  , &
               fdp(1,1,1),nlon  )
!
! Accumulate T interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
#ifdef HADVTEST
        grtlat(i,k) = fdp(i,k,1)
#else
        grtlat(i,k) = grtlat(i,k) + fdp(i,k,1)
#endif
     end do
  end do
!
! Reset "kdp" to arrival indices everywhere so that we can do true
! horizontal interpolation rather than "vertical non-interpolation"
!
  do k = 1,plev
     do i = 1,nlon
        kdp(i,k) = k
     end do
  end do
!
! Compute Ps and remaining T interpolants:
! Non-monotonic, 2-D interpolation
!
  limdrh  = .false.
  limdrv  = .false.
  lhrzint = .true.
  lvrtint = .false.
  call sltint (plev    ,jcen    ,lnpssld ,lam     ,rdphi   , &
               rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
               wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
               hr      ,dhl     ,dhr     ,ys      ,yn      , &
               wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
               hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
               wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
               dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
               lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp(1,1,1), &
               nlon  )
  call sltint (plev    ,jcen    ,prhssld ,lam     ,rdphi   , &
               rdz     ,lbasdy  ,lbasdz  ,xl      ,xr      , &
               wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
               hr      ,dhl     ,dhr     ,ys      ,yn      , &
               wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
               hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
               wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
               dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
               lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp(1,1,2), &
               nlon )
  do i = 1,nlon
     pdsum(i) = 0.
     pd   (i) = 0.
     pdsm1(i) = 0.
     pd1  (i) = 0.
  end do
!
! Accumulate P-interpolants into T equation
!
#if ( ! defined HADVTEST )
  do k = nprlev,plev
     tmp1 = cappa*t0(k)*hypi(plevp)/hypm(k)
     do i = 1,nlon
        grtlat(i,k) = grtlat(i,k) - fdp(i,k,1)*tmp1*hybm(k) 
     end do
  end do
  do k = nprlev,plev
     do l = nprlev,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + fdp(i,l,1)*hybd(l)*gamma(l,k)
        end do
     end do
  end do
  do k = 1,plev
     do l = 1,k
        do i = 1,nlon
           grtlat(i,k) = grtlat(i,k) + fdp(i,l,2)*gamma(l,k)
        end do
     end do
  end do
#endif
!
! Accumulate Ps interpolants in 3-D array
!
  do k = 1,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + fdp(i,k,2)
        pdsum(i) = pdsum(i) + fdp(i,k,2)
     end do
  end do
  do k = nprlev,plev
     do i = 1,nlon
        grpslat(i) = grpslat(i) + fdp(i,k,1)*hybd(k)
        pdsm1(i) = pdsm1(i) + fdp(i,k,1)*hybd(k)
     end do
  end do
!
! Compute remainder of (1/ps)etadot(dp/deta)
!
  do k = 1,plev-1
     do i = 1,nlon
        pd (i) = pd (i) + fdp(i,k,2)
        parrsld(i1-1+i,k,jcen) = parrsld(i1-1+i,k,jcen) + pd (i)
     end do
!
     if(k.ge.nprlev) then
        do i=1,nlon
           pd1(i) = pd1(i) + fdp(i,k,1)*hybd(k)
           parrsld(i1-1+i,k,jcen) = parrsld(i1-1+i,k,jcen) + pd1(i) - &
                hybi(k+1)*(pdsum(i) + pdsm1(i))
        end do
     end if
     do i = 1,nlon
        parrsld(i1-1+i,k,jcen) = parrsld(i1-1+i,k,jcen)*dtr
     end do
  end do
!
! Begin FFT
!
! fu,fv,T,q,vort,Ps: note contiguity assumptions
!
  inc   = 1
  isign = -1
  ntr   = 4*plev + 1
  call fft991(xnlin   ,work    ,trig(1,irow),ifax(1,irow),inc     , &
              plond   ,nlon    ,ntr         ,isign       )
!
  if (lat.gt.plat/2) then
!
! First of latitude pair (N. hemisphere), save Fourier coefficients
!
     do i = 1,2*nmmax(irow)
        grlps2(i,irow) = grpslat(i)
     end do
     do k = 1,plev
        do i = 1,2*nmmax(irow)
#if ( defined PVP )
           grfu2(i,k,irow) = grfulat(i,k)
           grfv2(i,k,irow) = grfvlat(i,k)
           grt2 (i,k,irow) = grtlat (i,k)
           grq2 (i,k,irow) = grqlat (i,k)
#else
           grfu2(k,i,irow) = grfulat(i,k)
           grfv2(k,i,irow) = grfvlat(i,k)
           grt2 (k,i,irow) = grtlat (i,k)
           grq2 (k,i,irow) = grqlat (i,k)
#endif
        end do
     end do
  else
!
! Second of latitude pair (S. hemisphere), save Fourier coefficients
!
     do i = 1,2*nmmax(irow)
        grlps1(i,irow) = grpslat(i)
     end do
     do k = 1,plev
        do i = 1,2*nmmax(irow)
#if ( defined PVP )
           grfu1(i,k,irow) = grfulat(i,k)
           grfv1(i,k,irow) = grfvlat(i,k)
           grt1 (i,k,irow) = grtlat (i,k)
           grq1 (i,k,irow) = grqlat (i,k)
#else
           grfu1(k,i,irow) = grfulat(i,k)
           grfv1(k,i,irow) = grfvlat(i,k)
           grt1 (k,i,irow) = grtlat (i,k)
           grq1 (k,i,irow) = grqlat (i,k)
#endif
        end do
     end do
  end if
!
  return
end subroutine scanslt

