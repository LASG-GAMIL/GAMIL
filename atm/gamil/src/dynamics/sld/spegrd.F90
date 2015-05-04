#include <misc.h>
#include <params.h>

subroutine spegrd (ztodt   ,lat     ,cwava   ,qfcst   ,q3      , &
                   etamid  ,ps      ,u3      ,v3      ,t3      , &
                   div     ,hw2al   ,hw2bl   ,hw3al   ,hw3bl   , &
                   hwxal   ,hwxbl   ,grts    ,grqs    ,grths   , &
                   grds    ,grus    ,gruhs   ,grvs    ,grvhs   , &
                   grpss   ,grdps   ,grpms   ,grpls   ,grtms   , &
                   grtls   ,grqms   ,grqls   ,grta    ,grqa    , &
                   grtha   ,grda    ,grua    ,gruha   ,grva    , &
                   grvha   ,grpsa   ,grdpa   ,grpma   ,grpla   , &
                   grtma   ,grtla   ,grqma   ,grqla   ,dps     , &
                   dpsl    ,dpsm    ,tl      ,tm      ,ql      , &
                   qm      ,t3m1    ,engy2alat,engy2blat,difftalat, &
                   difftblat,phis   ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Transfrom variables from spherical harmonic coefficients 
! to grid point values during second gaussian latitude scan (scan2)
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: spegrd.F90,v 1.13.4.3 2002/06/21 05:37:41 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use constituents, only: pcnst, pnats
   use pspect
   use comspe
   use commap
   use history, only: outfld
   use physconst, only: rga

   implicit none

#include <comctl.h>
#include <comfft.h>
#include <comhd.h>
#include <comhyb.h>
#include <comlun.h>
#include <comqfl.h>
!
! Arguments
!
   real(r8), intent(in)   :: ztodt                ! timestep
   real(r8), intent(in)   :: cwava                ! normalization factor (1/g*plon)
   real(r8), intent(in)   :: qfcst(plond,plev,pcnst)       ! fcst q + consts
   real(r8), intent(in)   :: q3(plond,plev,pcnst+pnats) ! q + consts
   real(r8), intent(in)   :: etamid(plev)                     ! vertical coords at midpts 
   real(r8), intent(inout)  :: ps(plond)      ! surface pressure
   real(r8), intent(inout)  :: u3(plond,plev) ! u-wind
   real(r8), intent(inout)  :: v3(plond,plev) ! v-wind
   real(r8), intent(inout)  :: t3(plond,plev) ! temperature
   real(r8), intent(inout) :: div(plond,plev) ! divergence

   real(r8), intent(out)  :: hw2al(pcnst)               ! -
   real(r8), intent(out)  :: hw2bl(pcnst)               !  | lat contributions to
   real(r8), intent(out)  :: hw3al(pcnst)               !  | components of slt global
   real(r8), intent(out)  :: hw3bl(pcnst)               !  | mass integrals
   real(r8), intent(out)  :: hwxal(pcnst,4)             !  |
   real(r8), intent(out)  :: hwxbl(pcnst,4)             ! -

   real(r8), intent(out) :: dps(plond)
   real(r8), intent(out) :: dpsl(plond)
   real(r8), intent(out) :: dpsm(plond)

   real(r8) :: tl(plond,plev)
   real(r8) :: tm(plond,plev)
   real(r8) :: ql(plond,plev)
   real(r8) :: qm(plond,plev)
   real(r8), intent(in)  :: t3m1(plond,plev) ! temperature
   real(r8), intent(out) :: engy2alat
   real(r8), intent(out) :: engy2blat
   real(r8), intent(out) :: difftalat
   real(r8), intent(out) :: difftblat
   real(r8), intent(in)  :: phis(plond)

   integer lat                   ! latitude index
   integer nlon            ! number of longitudes
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: duh(plond,plev) ! 
   real(r8) :: dvh(plond,plev) ! 
   real(r8) :: dth(plond,plev) ! 

   real(r8) pmid (plond,plev)    ! pressure at model levels
   real(r8) pint (plond,plevp)   ! pressure at model interfaces
   real(r8) pdel (plond,plev)    ! pdel(k) = pint(k+1) - pint(k)
   real(r8) pdelb(plond,plev)    ! pressure diff between interfaces
!                                ! (press defined using the "B" part 
!                                ! of the hybrid grid only)
!                                
! Symmetric fourier coefficient arrays for all variables transformed 
! from spherical harmonics (see subroutine grcalc)
!                                
   real(r8) grdps(plond)         ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grds (plond,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruhs(plond,plev)    ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvhs(plond,plev)    ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grths(plond,plev)    ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grpss(plond)         ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grus (plond,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvs (plond,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grts (plond,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8) grqs (plond,plev)    ! sum(n) of q(n,m)*P(n,m)
   real(r8) grpls(plond)         ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grtms(plond,plev)
   real(r8) grtls(plond,plev)
   real(r8) grqms(plond,plev)
   real(r8) grqls(plond,plev)
   real(r8) grpms(plond)         ! sum(n) of lnps(n,m)*H(n,m)
!
! Antisymmetric fourier coefficient arrays for all variables
! transformed from spherical harmonics (see grcalc)
!
   real(r8) grdpa(plond)         ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grda (plond,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruha(plond,plev)    ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvha(plond,plev)    ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grtha(plond,plev)    ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grpsa(plond)         ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grua (plond,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grva (plond,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grta (plond,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8) grqa (plond,plev)    ! sum(n) of q(n,m)*P(n,m)
   real(r8) grpla(plond)         ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grtma(plond,plev)
   real(r8) grtla(plond,plev)
   real(r8) grqma(plond,plev)
   real(r8) grqla(plond,plev)
   real(r8) grpma(plond)         ! sum(n) of lnps(n,m)*H(n,m)
!
! Local workspace
!
   real(r8) qfcst1(plond,plev,pcnst) ! workspace to please lf95 compiler
   real(r8) hcwavaw              ! 0.5*cwava*w(irow)
   real(r8) sum

#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! workspace needed by fft991
#else
   real(r8) work((plon+1)*pcray) ! workspace needed by fft991
#endif
   integer isign                 ! +1 => transform spectral to grid
   integer inc                   ! distance between transform elements
   integer i,k,m                 ! longitude, level, constituent indices
   integer ihem                  ! hemisphere index
   integer klev                  ! top level where hybrid coordinates apply
   
   real(r8) rcoslat              ! 1./cosine(latitude)
   real(r8) dotproda             ! dot product
   real(r8) dotprodb             ! dot product
!
!-----------------------------------------------------------------------
!
   qfcst1(1:nlon,:,:) = qfcst(i1:nlon+i1-1,:,:)
!
! Set local integer pointers into part of buffer which is not written
! to SSD.  Contiguity is still needed for speed of FFT
!
   inc = 1
   isign = +1
!
! Assemble northern and southern hemisphere grid values from the
! symmetric and antisymmetric fourier coefficients. 
! 1. Determine the fourier coefficients for the northern or southern
!    hemisphere latitude. 
! 2. Transform to gridpoint values
! 3. Clean up
!
   if (lat > plat/2) then                       ! Northern hemisphere
      do k=1,plev
         do i=1,nlon+2
            div(i,k) = grds(i,k) + grda(i,k)
            duh(i,k) = gruhs(i,k) + gruha(i,k)
            dvh(i,k) = grvhs(i,k) + grvha(i,k)
            dth(i,k) = grths(i,k) + grtha(i,k)
            tl (i,k) = grtls(i,k) + grtla(i,k)
            tm (i,k) = grtms(i,k) + grtma(i,k)
            ql (i,k) = grqls(i,k) + grqla(i,k)
            qm (i,k) = grqms(i,k) + grqma(i,k)
            u3(i,k)  = grus(i,k) + grua(i,k)
            v3(i,k)  = grvs(i,k) + grva(i,k)
            t3(i,k)  = grts(i,k) + grta(i,k)
         end do
      end do

      do i=1,nlon+2
         dps(i) = grdps(i) + grdpa(i)
         ps(i) = grpss(i) + grpsa(i)
         dpsl(i) = grpls(i) + grpla(i)
         dpsm(i) = grpms(i) + grpma(i)
      end do

   else                                          ! Southern hemisphere

      do k=1,plev
         do i=1,nlon+2
            div(i,k) = grds(i,k) - grda(i,k)
            duh(i,k) = gruhs(i,k) - gruha(i,k)
            dvh(i,k) = grvhs(i,k) - grvha(i,k)
            dth(i,k) = grths(i,k) - grtha(i,k)
            tl (i,k) = grtls(i,k) - grtla(i,k)
            tm (i,k) = grtms(i,k) - grtma(i,k)
            ql (i,k) = grqls(i,k) - grqla(i,k)
            qm (i,k) = grqms(i,k) - grqma(i,k)
            u3(i,k) = grus(i,k) - grua(i,k)
            v3(i,k) = grvs(i,k) - grva(i,k)
            t3(i,k) = grts(i,k) - grta(i,k)
         end do
      end do

      do i=1,nlon+2
         dps(i) = grdps(i) - grdpa(i)
         ps(i) = grpss(i) - grpsa(i)
         dpsl(i) = grpls(i) - grpla(i)
         dpsm(i) = grpms(i) - grpma(i)
      end do
   end if
!
! ps,div,tl,tm,dpsl,dpsm,ql,qm,dth,duh,dvh,dps
!
   call fft991 (ps,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, 1,    isign)
   call fft991 (div,  work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (tl,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (tm,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (dpsl, work, trig(1,lat), ifax(1,lat), inc, plond, nlon, 1,    isign)
   call fft991 (dpsm, work, trig(1,lat), ifax(1,lat), inc, plond, nlon, 1,    isign)
   call fft991 (ql,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (qm,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (dth,  work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (duh,  work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (dvh,  work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (dps,  work, trig(1,lat), ifax(1,lat), inc, plond, nlon, 1,    isign)
!
! u,v,t
!
   call fft991 (u3,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (v3,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
   call fft991 (t3,   work, trig(1,lat), ifax(1,lat), inc, plond, nlon, plev, isign)
!
! Remove cosine(latitude) from momentum variables
!
   rcoslat = 1./cos(clat(lat))
   do k=1,plev
      do i=1,nlon
         u3(i,k) = u3(i,k)*rcoslat
         v3(i,k) = v3(i,k)*rcoslat
         duh(i,k) = duh(i,k)*rcoslat
         dvh(i,k) = dvh(i,k)*rcoslat
      end do
   end do
!
! Copy transformed surface pressure back to in-core array, converting
! from log(ps) to ps.
!
   do i=1,nlon
      ps(i) = exp(ps(i))
   end do
!
! Diagnose pressure arrays needed by DIFCOR
!
   call plevs0 (nlon, plond, plev, ps, pint, pmid, pdel)
   call pdelb0 (ps, pdelb, nlon)
!
! Accumulate mass integrals
!
   sum = 0.
   do i=1,nlon
      sum = sum + ps(i)
   end do
   tmass(lat) = w(lat)*rga*sum/nlon
!
! Finish horizontal diffusion: add pressure surface correction term to
! t and q diffusions; add kinetic energy dissipation to internal energy
! (temperature)
!
   klev = max(kmnhd4,nprlev)
   call difcor (klev,   ztodt,  dps,    u3,     v3, &
                q3,     pdel,   pint,   t3,     dth, &
                duh,    dvh,    nlon)
!
! Calculate SLT moisture, constituent, energy, and temperature integrals
!
   hcwavaw   = 0.5*cwava*w(lat)
   engy2alat = 0.
   engy2blat = 0.
   difftalat = 0.
   difftblat = 0.
   do m=1,pcnst
      hw2al(m) = 0.
      hw2bl(m) = 0.
      hw3al(m) = 0.
      hw3bl(m) = 0.
      hwxal(m,1) = 0.
      hwxal(m,2) = 0.
      hwxal(m,3) = 0.
      hwxal(m,4) = 0.
      hwxbl(m,1) = 0.
      hwxbl(m,2) = 0.
      hwxbl(m,3) = 0.
      hwxbl(m,4) = 0.
      do k=1,plev
         dotproda = 0.
         dotprodb = 0.
         do i=1,nlon
            dotproda = dotproda + qfcst1(i,k,m)*pdela(i,k)
            dotprodb = dotprodb + qfcst1(i,k,m)*pdelb(i,k)
         end do
         hw2al(m) = hw2al(m) + hcwavaw*dotproda
         hw2bl(m) = hw2bl(m) + hcwavaw*dotprodb
      end do
   end do

   call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdela, engy2alat ,nlon)
   call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdelb, engy2blat ,nlon)
   call engy_tdif(cwava ,w(lat) ,t3  ,t3m1             ,pdela, difftalat ,nlon)
   call engy_tdif(cwava ,w(lat) ,t3  ,t3m1             ,pdelb, difftblat ,nlon)

   call qmassd (cwava, etamid, w(lat), q3, qfcst1, &
                pdela, hw3al, nlon)

   call qmassd (cwava, etamid, w(lat), q3, qfcst1, &
                pdelb, hw3bl, nlon)

   if (pcnst.gt.1) then
      call xqmass (cwava, etamid, w(lat), q3, qfcst1, &
                   q3, qfcst1, pdela, pdelb, hwxal, &
                   hwxbl, nlon)
   end if

   call outfld ('DTH     ',dth     ,plond   ,lat     )

   return
end subroutine spegrd
