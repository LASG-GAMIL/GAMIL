#include <misc.h>
#include <params.h>

subroutine spetru (ps      ,phis    ,u3      ,v3      ,t3      , &
                   q3      ,div     ,dpsl    ,dpsm    ,tl      , &
                   tm      ,ql      ,qm      ,phi     ,phisl   , &
                   phism   ,phis_hires)
!-----------------------------------------------------------------------
!
! Purpose:
! Spectrally truncate input fields which have already been transformed into 
! fourier space.  Some arrays are dimensioned (2,...), where (1,...) is the
! real part of the complex fourier coefficient, and (2,...) is the imaginary.
! Any array dimensioned (plond,...) *cannot* be dimensioned (2,plond/2,...) 
! because plond *may* be (and in fact currently is) odd. In these cases 
! reference to real and imaginary parts is by (2*m-1,...) and (2*m  ,...) 
! respectively.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use constituents, only: pcnst, pnats
  use pspect
  use comspe
  use rgrid,        only: nlon, nmmax
  use commap,       only: w, xm, rsq, cs
  use dynconst,     only: ez, ra, rearth

  implicit none

#include <comctl.h>
#include <comfft.h>

!------------------------------Arguments--------------------------------
!
  real(r8), intent(inout):: ps   (plond,plat)      ! Fourier -> spec. coeffs. for ln(ps)
  real(r8), intent(inout):: phis (plond,plat)      ! Fourier -> spec. coeffs. for sfc geo.
  real(r8), intent(inout):: u3   (plond,plev,plat) ! Fourier -> spec. coeffs. for u-wind
  real(r8), intent(inout):: v3   (plond,plev,plat) ! Fourier -> spec. coeffs. for v-wind
  real(r8), intent(inout):: t3   (plond,plev,plat) ! Fourier -> spec. coeffs. for temp
  real(r8), intent(inout):: q3   (plond,plev*(pcnst+pnats),plat)
!                                                  ! Fourier -> spec. coeffs. for q
  real(r8), intent(out)  :: div  (plond,plev,plat) ! Spectrally truncated divergence
  real(r8), intent(out)  :: dpsl (plond,plat)      ! Spectrally trunc d(ln(ps))/d(long)
  real(r8), intent(out)  :: dpsm (plond,plat)      ! Spectrally trunc d(ln(ps))/d(lat )
  real(r8), intent(out)  :: tl   (plond,plev,plat) ! Spectrally trunc d(T)/d(longitude)
  real(r8), intent(out)  :: tm   (plond,plev,plat) ! Spectrally trunc d(T)/d(latitude)
  real(r8), intent(out)  :: ql   (plond,plev,plat) ! Spectrally trunc d(q)/d(longitude)
  real(r8), intent(out)  :: qm   (plond,plev,plat) ! Spectrally trunc d(q)/d(latitude)
  real(r8), intent(out)  :: phi  (2,psp/2)         ! used in spectral truncation of phis
  real(r8), intent(out)  :: phisl(plond,plat)      ! Spectrally trunc d(phis)/d(longitude)
  real(r8), intent(out)  :: phism(plond,plat)      ! Spectrally trunc d(phis)/d(latitude)
  logical , intent(in)   :: phis_hires             ! true => PHIS came from hi res topog
!
!---------------------------Local workspace-----------------------------
!
  real(r8) alpn (pspt)          ! alp*rsq*xm*ra
  real(r8) dalpn(pspt)          ! dalp*rsq*ra
  real(r8) tmp1                 ! vector temporary
  real(r8) tmp2                 ! vector temporary
  real(r8) tmpr                 ! vector temporary (real)
  real(r8) tmpi                 ! vector temporary (imaginary)
  real(r8) phialpr,phialpi      ! phi*alp (real and imaginary)
  real(r8) phdalpr,phdalpi      ! phi*dalp (real and imaginary)
  real(r8) zwalp                ! zw*alp
  real(r8) zwdalp               ! zw*dalp
  real(r8) psdalpr,psdalpi      ! alps (real and imaginary)*dalp
  real(r8) psalpr,psalpi        ! alps (real and imaginary)*alp
  real(r8) zrcsj                ! ra/(cos**2 latitude)
  real(r8) zw                   ! w**2
  real(r8) filtlim              ! filter function
  real(r8) ft                   ! filter multiplier for spectral coefficients

#if ( ! defined USEFFTLIB )
  real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
  real(r8) work((plon+1)*pcray) ! Workspace for fft
#endif
  real(r8) zsqcs

  integer ir,ii                 ! indices complex coeffs. of spec. arrs.
  integer i,k                   ! longitude, level indices
  integer irow                  ! latitude pair index
  integer latm,latp             ! symmetric latitude indices
  integer lat                   ! index
  integer m                     ! longitudinal wavenumber index (non-PVP)
!                               ! along-diagonal index (PVP)
  integer n                     ! latitudinal wavenumber index (non-PVP)
!                               ! diagonal index (PVP)
  integer nspec                 ! index
#if ( defined PVP )              
  integer ne                    ! index into rsq
  integer ic                    ! complex coeff. index
  integer ialp                  ! index into legendre polynomials
#else
  integer mr,mc                 ! spectral indices
#endif
!
!-----------------------------------------------------------------------
!
! Zero spectral arrays
!
  alps(:)   = 0.
  vz  (:,:) = 0.
  d   (:,:) = 0.
  t   (:,:) = 0.
  q   (:,:) = 0.
  phi (:,:) = 0.
!
! Compute the quantities which are transformed to spectral space:
!   1. u = u*sqrt(1-mu*mu),   u * cos(phi)
!   2. v = v*sqrt(1-mu*mu),   v * cos(phi)
!   3. t = t                  T
!   4. ps= ln(ps). 
!
  do lat=1,plat
     irow = lat
     if (lat.gt.plat/2) irow = plat - lat + 1
     zsqcs = sqrt(cs(irow))

     do k=1,plev
        do i=1,nlon(lat)
           u3(i,k,lat) = u3(i,k,lat)*zsqcs
           v3(i,k,lat) = v3(i,k,lat)*zsqcs
        end do
     end do

     do i=1,nlon(lat)
        ps(i,lat) = log(ps(i,lat))
     end do
!
! Transform grid -> fourier
! 1st transform: U,V,T and Q
! 2nd transform: LN(PS).  3rd transform: surface geopotential
!
     call fft991 (u3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond   ,nlon(lat),plev        ,-1         )
     call fft991 (v3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,plev        ,-1         )
     call fft991 (t3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,plev        ,-1         )
     call fft991 (q3(1,1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,plev        ,-1         )
     call fft991 (ps(1,lat),work       ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,1           ,-1         )
     call fft991 (phis(1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond, nlon(lat)  ,1           ,-1         )
  end do                    ! lat=1,plat
!
! Loop over latitude pairs
!
  do irow=1,plat/2
     latp = irow
     latm = plat - irow + 1
     zrcsj = ra/cs(irow)
     zw = w(irow)*2.
     do i=1,2*nmmax(irow)
!
! Compute symmetric and antisymmetric components: ps first, then phis
!
        tmp1 = 0.5*(ps(i,latm) - ps(i,latp))
        tmp2 = 0.5*(ps(i,latm) + ps(i,latp))
        ps(i,latm) = tmp1
        ps(i,latp) = tmp2

        tmp1 = 0.5*(phis(i,latm) - phis(i,latp))
        tmp2 = 0.5*(phis(i,latm) + phis(i,latp))
        phis(i,latm) = tmp1
        phis(i,latp) = tmp2
     end do
!
! Multi-level fields: U, V, T
!
     do k=1,plev
        do i=1,2*nmmax(irow)

           tmp1 = 0.5*(u3(i,k,latm) - u3(i,k,latp))
           tmp2 = 0.5*(u3(i,k,latm) + u3(i,k,latp))
           u3(i,k,latm) = tmp1
           u3(i,k,latp) = tmp2

           tmp1 = 0.5*(v3(i,k,latm) - v3(i,k,latp))
           tmp2 = 0.5*(v3(i,k,latm) + v3(i,k,latp))
           v3(i,k,latm) = tmp1
           v3(i,k,latp) = tmp2

           tmp1 = 0.5*(t3(i,k,latm) - t3(i,k,latp))
           tmp2 = 0.5*(t3(i,k,latm) + t3(i,k,latp))
           t3(i,k,latm) = tmp1
           t3(i,k,latp) = tmp2

           tmp1 = 0.5*(q3(i,k,latm) - q3(i,k,latp))
           tmp2 = 0.5*(q3(i,k,latm) + q3(i,k,latp))
           q3(i,k,latm) = tmp1
           q3(i,k,latp) = tmp2
        end do
     end do
!     
! Compute vzmn,dmn and ln(p*)mn and also phi*mn,tmn and qmn
!
#if ( defined PVP )
     do n=1,pmax,2
        ic = ncoefi(n) - 1
        ialp = nalp(n)
        do m=1,nmreduced(n,irow)
           zwalp = zw*alp(ialp+m,irow)
           phi(1,ic+m) = phi(1,ic+m) + zwalp*phis(2*m-1,latp)
           phi(2,ic+m) = phi(2,ic+m) + zwalp*phis(2*m  ,latp)
           ir = 2*(ic+m) - 1
           ii = ir + 1
           alps(ir) = alps(ir) + zwalp*ps(2*m-1,latp)
           alps(ii) = alps(ii) + zwalp*ps(2*m  ,latp)
        end do
     end do
!
     do n=2,pmax,2
        ic = ncoefi(n) - 1
        ialp = nalp(n)
        do m=1,nmreduced(n,irow)
           zwalp = zw*alp(ialp+m,irow)
           phi(1,ic+m) = phi(1,ic+m) + zwalp*phis(2*m-1,latm)
           phi(2,ic+m) = phi(2,ic+m) + zwalp*phis(2*m  ,latm)
           ir = 2*(ic+m) - 1
           ii = ir + 1
           alps(ir) = alps(ir) + zwalp*ps(2*m-1,latm)
           alps(ii) = alps(ii) + zwalp*ps(2*m  ,latm)
        end do
     end do
#else
     do m=1,nmmax(irow)
        mr = nstart(m)
        mc = 2*mr
        do n=1,nlen(m),2
           zwalp = zw*alp(mr+n,irow)
           phi(1,mr+n) = phi(1,mr+n) + zwalp*phis(2*m-1,latp)
           phi(2,mr+n) = phi(2,mr+n) + zwalp*phis(2*m  ,latp)
           ir = mc + 2*n - 1
           ii = ir + 1
           alps(ir) = alps(ir) + zwalp*ps(2*m-1,latp)
           alps(ii) = alps(ii) + zwalp*ps(2*m  ,latp)
        end do

        do n=2,nlen(m),2
           zwalp = zw*alp(mr+n,irow)
           phi(1,mr+n) = phi(1,mr+n) + zwalp*phis(2*m-1,latm)
           phi(2,mr+n) = phi(2,mr+n) + zwalp*phis(2*m  ,latm)
           ir = mc + 2*n - 1
           ii = ir + 1
           alps(ir) = alps(ir) + zwalp*ps(2*m-1,latm)
           alps(ii) = alps(ii) + zwalp*ps(2*m  ,latm)
        end do
     end do
#endif
     do k=1,plev
#if ( defined PVP )
        do n=1,pmax,2
           ic = ncoefi(n) - 1
           ialp = nalp(n)
           do m=1,nmreduced(n,irow)
              zwdalp   = zw*dalp(ialp+m,irow)
              zwalp    = zw*alp (ialp+m,irow)
              ir       = 2*(ic+m) - 1
              ii       = ir + 1
              d(ir,k)  = d(ir,k) - (zwdalp*v3(2*m-1,k,latm) + &
                         xm(m)*zwalp*u3(2*m  ,k,latp))*zrcsj
              d(ii,k)  = d(ii,k) - (zwdalp*v3(2*m  ,k,latm) - &
                         xm(m)*zwalp*u3(2*m-1,k,latp))*zrcsj
              t(ir,k)  = t(ir,k) + zwalp*t3(2*m-1,k,latp)
              t(ii,k)  = t(ii,k) + zwalp*t3(2*m  ,k,latp)
              q(ir,k)  = q(ir,k) + zwalp*q3(2*m-1,k,latp)
              q(ii,k)  = q(ii,k) + zwalp*q3(2*m  ,k,latp)
              vz(ir,k) = vz(ir,k) + (zwdalp*u3(2*m-1,k,latm) - &
                         xm(m)*zwalp*v3(2*m  ,k,latp))*zrcsj
              vz(ii,k) = vz(ii,k) + (zwdalp*u3(2*m  ,k,latm) + &
                         xm(m)*zwalp*v3(2*m-1,k,latp))*zrcsj
           end do
        end do
!
        do n=2,pmax,2
           ic = ncoefi(n) - 1
           ialp = nalp(n)
           do m=1,nmreduced(n,irow)
              zwdalp   = zw*dalp(ialp+m,irow)
              zwalp    = zw*alp (ialp+m,irow)
              ir       = 2*(ic+m) - 1
              ii       = ir + 1
              d(ir,k)  = d(ir,k) - (zwdalp*v3(2*m-1,k,latp) + &
                         xm(m)*zwalp*u3(2*m  ,k,latm))*zrcsj
              d(ii,k)  = d(ii,k) - (zwdalp*v3(2*m  ,k,latp) - &
                         xm(m)*zwalp*u3(2*m-1,k,latm))*zrcsj
              t(ir,k)  = t(ir,k) + zwalp*t3(2*m-1,k,latm)
              t(ii,k)  = t(ii,k) + zwalp*t3(2*m  ,k,latm)
              q(ir,k)  = q(ir,k) + zwalp*q3(2*m-1,k,latm)
              q(ii,k)  = q(ii,k) + zwalp*q3(2*m  ,k,latm)
              vz(ir,k) = vz(ir,k) + (zwdalp*u3(2*m-1,k,latp) - &
                         xm(m)*zwalp*v3(2*m  ,k,latm))*zrcsj
              vz(ii,k) = vz(ii,k) + (zwdalp*u3(2*m  ,k,latp) + &
                         xm(m)*zwalp*v3(2*m-1,k,latm))*zrcsj
           end do
        end do
#else
        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=1,nlen(m),2
              zwdalp   = zw*dalp(mr+n,irow)
              zwalp    = zw*alp (mr+n,irow)
              ir       = mc + 2*n - 1
              ii       = ir + 1
              d(ir,k)  = d(ir,k) - (zwdalp*v3(2*m-1,k,latm) + &
                         xm(m)*zwalp*u3(2*m  ,k,latp))*zrcsj
              d(ii,k)  = d(ii,k) - (zwdalp*v3(2*m  ,k,latm) - &
                         xm(m)*zwalp*u3(2*m-1,k,latp))*zrcsj
              t(ir,k)  = t(ir,k) + zwalp*t3(2*m-1,k,latp)
              t(ii,k)  = t(ii,k) + zwalp*t3(2*m  ,k,latp)
              q(ir,k)  = q(ir,k) + zwalp*q3(2*m-1,k,latp)
              q(ii,k)  = q(ii,k) + zwalp*q3(2*m  ,k,latp)
              vz(ir,k) = vz(ir,k) + (zwdalp*u3(2*m-1,k,latm) - &
                         xm(m)*zwalp*v3(2*m  ,k,latp))*zrcsj
              vz(ii,k) = vz(ii,k) + (zwdalp*u3(2*m  ,k,latm) + &
                         xm(m)*zwalp*v3(2*m-1,k,latp))*zrcsj
           end do
        end do

        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=2,nlen(m),2
              zwdalp   = zw*dalp(mr+n,irow)
              zwalp    = zw*alp (mr+n,irow)
              ir       = mc + 2*n - 1
              ii       = ir + 1
              d(ir,k)  = d(ir,k) - (zwdalp*v3(2*m-1,k,latp) + &
                         xm(m)*zwalp*u3(2*m  ,k,latm))*zrcsj
              d(ii,k)  = d(ii,k) - (zwdalp*v3(2*m  ,k,latp) - &
                         xm(m)*zwalp*u3(2*m-1,k,latm))*zrcsj
              t(ir,k)  = t(ir,k) + zwalp*t3(2*m-1,k,latm)
              t(ii,k)  = t(ii,k) + zwalp*t3(2*m  ,k,latm)
              q(ir,k)  = q(ir,k) + zwalp*q3(2*m-1,k,latm)
              q(ii,k)  = q(ii,k) + zwalp*q3(2*m  ,k,latm)
              vz(ir,k) = vz(ir,k) + (zwdalp*u3(2*m-1,k,latp) - &
                         xm(m)*zwalp*v3(2*m  ,k,latm))*zrcsj
              vz(ii,k) = vz(ii,k) + (zwdalp*u3(2*m  ,k,latp) + &
                         xm(m)*zwalp*v3(2*m-1,k,latm))*zrcsj
           end do
        end do
#endif
     end do
  end do
!
  if (phis_hires) then
!
! Apply spectral filter to phis
!     filter is a function of n 
!        if n < filter limit then
!           spectral_coeff = spectral_coeff * (1. - (float(n)/filtlim)**2)
!        else         
!           spectral_coeff = 0.
!        endif
!     where filter limit = 1.4*PTRN
!     
     filtlim = float(int(1.4*float(ptrn)))
#if ( defined PVP )
     do n=1,pmax
        ic = ncoefi(n) - 1
        do m=1,nm(n)
           nspec=m-1+n
           ft = 1. - (float(nspec)/filtlim)**2
           if (float(nspec) .ge. filtlim) ft = 0.
           phi(1,ic+m) = phi(1,ic+m)*ft
           phi(2,ic+m) = phi(2,ic+m)*ft
        end do
     end do
#else    
     do m=1,pmmax
        mr = nstart(m)
        do n=1,nlen(m)
           nspec=m-1+n
           ft = 1. - (float(nspec)/filtlim)**2
           if (float(nspec) .ge. filtlim) ft = 0.
           phi(1,mr+n) = phi(1,mr+n)*ft 
           phi(2,mr+n) = phi(2,mr+n)*ft 
        end do
     end do
#endif   
     call hordif1(rearth  ,phi     )
  end if
!
! Compute grid point values of:phi*,u,v,ln(p*),t,q,vz,d and grad(ln(p*)).
!
  do irow=1,plat/2
     latp = irow
     latm = plat - irow + 1
!
! Zero fourier fields
!
     phis(:,latm) = 0.
     phis(:,latp) = 0.

     phisl(:,latm) = 0.
     phisl(:,latp) = 0.

     phism(:,latm) = 0.
     phism(:,latp) = 0.

     ps(:,latm) = 0.
     ps(:,latp) = 0.

     dpsl(:,latm) = 0.
     dpsl(:,latp) = 0.

     dpsm(:,latm) = 0.
     dpsm(:,latp) = 0.

     u3(:,:plev,latm) = 0.
     u3(:,:plev,latp) = 0.

     v3(:,:plev,latm) = 0.
     v3(:,:plev,latp) = 0.

     t3(:,:plev,latm) = 0.
     t3(:,:plev,latp) = 0.

     tl(:,:plev,latm) = 0.
     tl(:,:plev,latp) = 0.
!
     tm(:,:plev,latm) = 0.
     tm(:,:plev,latp) = 0.
!
     ql(:,:plev,latm) = 0.
     ql(:,:plev,latp) = 0.
!
     qm(:,:plev,latm) = 0.
     qm(:,:plev,latp) = 0.

     div(:,:plev,latm) = 0.
     div(:,:plev,latp) = 0.
!
! Compute(phi*,grad(phi*),u,d(u)/d(lamda),v,d(v)/d(lamda),
! ln(p*),t,grad(T),q,vz,d,grad(ln(p*)))m
!
#if ( defined PVP )
     do n=1,pmax,2
        ic = ncoefi(n) - 1
        ialp = nalp(n)
        do m=1,nmreduced(n,irow)
           phialpr = phi(1,ic+m)*alp(ialp+m,irow)
           phialpi = phi(2,ic+m)*alp(ialp+m,irow)
!
           phis(2*m-1,latm) = phis(2*m-1,latm) + phialpr
           phis(2*m  ,latm) = phis(2*m  ,latm) + phialpi
!
           phisl(2*m-1,latm) = phisl(2*m-1,latm) - phialpi*ra
           phisl(2*m  ,latm) = phisl(2*m  ,latm) + phialpr*ra
!
           phdalpr = phi(1,ic+m)*dalp(ialp+m,irow)
           phdalpi = phi(2,ic+m)*dalp(ialp+m,irow)
!
           phism(2*m-1,latp) = phism(2*m-1,latp) + phdalpr*ra
           phism(2*m  ,latp) = phism(2*m  ,latp) + phdalpi*ra
!
           ir = 2*(ic+m) - 1
           ii = ir + 1
           psalpr = alps(ir)*alp(ialp+m,irow)
           psalpi = alps(ii)*alp(ialp+m,irow)
!
           ps(2*m-1,latm) = ps(2*m-1,latm) + psalpr
           ps(2*m  ,latm) = ps(2*m  ,latm) + psalpi
           dpsl(2*m-1,latm) = dpsl(2*m-1,latm) - psalpi*ra
           dpsl(2*m  ,latm) = dpsl(2*m  ,latm) + psalpr*ra
!
           psdalpr = alps(ir)*dalp(ialp+m,irow)
           psdalpi = alps(ii)*dalp(ialp+m,irow)
!
           dpsm(2*m-1,latp) = dpsm(2*m-1,latp) + psdalpr*ra
           dpsm(2*m  ,latp) = dpsm(2*m  ,latp) + psdalpi*ra
        end do
     end do
!
     do n=2,pmax,2
        ic = ncoefi(n) - 1
        ialp = nalp(n)
        do m=1,nmreduced(n,irow)
           phialpr = phi(1,ic+m)*alp(ialp+m,irow)
           phialpi = phi(2,ic+m)*alp(ialp+m,irow)
!
           phis(2*m-1,latp) = phis(2*m-1,latp) + phialpr
           phis(2*m  ,latp) = phis(2*m  ,latp) + phialpi
           phisl(2*m-1,latp) = phisl(2*m-1,latp) - phialpi*ra
           phisl(2*m  ,latp) = phisl(2*m  ,latp) + phialpr*ra
!
           phdalpr = phi(1,ic+m)*dalp(ialp+m,irow)
           phdalpi = phi(2,ic+m)*dalp(ialp+m,irow)
!
           phism(2*m-1,latm) = phism(2*m-1,latm) + phdalpr*ra
           phism(2*m  ,latm) = phism(2*m  ,latm) + phdalpi*ra
!
           ir = 2*(ic+m) - 1
           ii = ir + 1
           psalpr = alps(ir)*alp(ialp+m,irow)
           psalpi = alps(ii)*alp(ialp+m,irow)
!
           ps(2*m-1,latp) = ps(2*m-1,latp) + psalpr
           ps(2*m  ,latp) = ps(2*m  ,latp) + psalpi
           dpsl(2*m-1,latp) = dpsl(2*m-1,latp) - psalpi*ra
           dpsl(2*m  ,latp) = dpsl(2*m  ,latp) + psalpr*ra
!
           psdalpr = alps(ir)*dalp(ialp+m,irow)
           psdalpi = alps(ii)*dalp(ialp+m,irow)
!     
           dpsm(2*m-1,latm) = dpsm(2*m-1,latm) + psdalpr*ra
           dpsm(2*m  ,latm) = dpsm(2*m  ,latm) + psdalpi*ra
        end do
     end do
!
     do n=1,pmax
        ne = n - 1
        ialp = nalp(n)
        do m=1,nmreduced(n,irow)
           alpn (ialp+m) =  alp(ialp+m,irow)*rsq(ne+m)*xm(m)*ra
           dalpn(ialp+m) = dalp(ialp+m,irow)*rsq(ne+m)      *ra
        end do
     end do
#else
     do m=1,nmmax(irow)
        mr = nstart(m)
        mc = 2*mr
        do n=1,nlen(m),2
           phialpr = phi(1,mr+n)*alp(mr+n,irow)
           phialpi = phi(2,mr+n)*alp(mr+n,irow)
!     
           phis(2*m-1,latm) = phis(2*m-1,latm) + phialpr
           phis(2*m  ,latm) = phis(2*m  ,latm) + phialpi
!
           phisl(2*m-1,latm) = phisl(2*m-1,latm) - phialpi*ra
           phisl(2*m  ,latm) = phisl(2*m  ,latm) + phialpr*ra
!
           phdalpr = phi(1,mr+n)*dalp(mr+n,irow)
           phdalpi = phi(2,mr+n)*dalp(mr+n,irow)
!
           phism(2*m-1,latp) = phism(2*m-1,latp) + phdalpr*ra
           phism(2*m  ,latp) = phism(2*m  ,latp) + phdalpi*ra
!
           ir = mc + 2*n - 1
           ii = ir + 1
           psalpr = alps(ir)*alp(mr+n,irow)
           psalpi = alps(ii)*alp(mr+n,irow)
!     
           ps(2*m-1,latm) = ps(2*m-1,latm) + psalpr
           ps(2*m  ,latm) = ps(2*m  ,latm) + psalpi
           dpsl(2*m-1,latm) = dpsl(2*m-1,latm) - psalpi*ra
           dpsl(2*m  ,latm) = dpsl(2*m  ,latm) + psalpr*ra
!
           psdalpr = alps(ir)*dalp(mr+n,irow)
           psdalpi = alps(ii)*dalp(mr+n,irow)
!
           dpsm(2*m-1,latp) = dpsm(2*m-1,latp) + psdalpr*ra
           dpsm(2*m  ,latp) = dpsm(2*m  ,latp) + psdalpi*ra
        end do
     end do

     do m=1,nmmax(irow)
        mr = nstart(m)
        mc = 2*mr
        do n=2,nlen(m),2
           phialpr = phi(1,mr+n)*alp(mr+n,irow)
           phialpi = phi(2,mr+n)*alp(mr+n,irow)
!     
           phis(2*m-1,latp) = phis(2*m-1,latp) + phialpr
           phis(2*m  ,latp) = phis(2*m  ,latp) + phialpi
           phisl(2*m-1,latp) = phisl(2*m-1,latp) - phialpi*ra
           phisl(2*m  ,latp) = phisl(2*m  ,latp) + phialpr*ra
!
           phdalpr = phi(1,mr+n)*dalp(mr+n,irow)
           phdalpi = phi(2,mr+n)*dalp(mr+n,irow)
!
           phism(2*m-1,latm) = phism(2*m-1,latm) + phdalpr*ra
           phism(2*m  ,latm) = phism(2*m  ,latm) + phdalpi*ra
!
           ir = mc + 2*n - 1
           ii = ir + 1
           psalpr = alps(ir)*alp(mr+n,irow)
           psalpi = alps(ii)*alp(mr+n,irow)
!     
           ps(2*m-1,latp) = ps(2*m-1,latp) + psalpr
           ps(2*m  ,latp) = ps(2*m  ,latp) + psalpi
           dpsl(2*m-1,latp) = dpsl(2*m-1,latp) - psalpi*ra
           dpsl(2*m  ,latp) = dpsl(2*m  ,latp) + psalpr*ra
!
           psdalpr = alps(ir)*dalp(mr+n,irow)
           psdalpi = alps(ii)*dalp(mr+n,irow)
!
           dpsm(2*m-1,latm) = dpsm(2*m-1,latm) + psdalpr*ra
           dpsm(2*m  ,latm) = dpsm(2*m  ,latm) + psdalpi*ra
        end do
     end do

     do m=1,nmmax(irow)
        mr = nstart(m)
        do n=1,nlen(m)
!
! These statements will likely not be bfb since xm*ra is now a scalar
!
           alpn (mr+n) =  alp(mr+n,irow)*rsq(n+m-1)*xm(m)*ra
           dalpn(mr+n) = dalp(mr+n,irow)*rsq(n+m-1)      *ra
        end do
     end do
#endif
     do m=1,nmmax(irow)
        dpsl(2*m-1,latm) = xm(m)*dpsl(2*m-1,latm)
        dpsl(2*m  ,latm) = xm(m)*dpsl(2*m  ,latm)
        dpsl(2*m-1,latp) = xm(m)*dpsl(2*m-1,latp)
        dpsl(2*m  ,latp) = xm(m)*dpsl(2*m  ,latp)
        phisl(2*m-1,latm) = xm(m)*phisl(2*m-1,latm)
        phisl(2*m  ,latm) = xm(m)*phisl(2*m  ,latm)
        phisl(2*m-1,latp) = xm(m)*phisl(2*m-1,latp)
        phisl(2*m  ,latp) = xm(m)*phisl(2*m  ,latp)
     end do

     do k=1,plev
#if ( defined PVP )
        do n=1,pmax,2
           ic = ncoefi(n) - 1
           ialp = nalp(n)
!DIR$ IVDEP
           do m=1,nmreduced(n,irow)
              ir = 2*(ic+m) - 1
              ii = ir + 1
!
              tmpr = d(ir,k)*alpn(ialp+m)
              tmpi = d(ii,k)*alpn(ialp+m)
              u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpi
              u3(2*m  ,k,latm) = u3(2*m  ,k,latm) - tmpr
!
              tmpr = d(ir,k)*dalpn(ialp+m)
              tmpi = d(ii,k)*dalpn(ialp+m)
              v3(2*m-1,k,latp) = v3(2*m-1,k,latp) - tmpr
              v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpi
!
              tmpr = vz(ir,k)*dalpn(ialp+m)
              tmpi = vz(ii,k)*dalpn(ialp+m)
              u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpr
              u3(2*m  ,k,latp) = u3(2*m  ,k,latp) + tmpi
!
              tmpr = vz(ir,k)*alpn(ialp+m)
              tmpi = vz(ii,k)*alpn(ialp+m)
              v3(2*m-1,k,latm) = v3(2*m-1,k,latm) + tmpi
              v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpr
!
              tmpr = t(ir,k)*alp(ialp+m,irow)
              tmpi = t(ii,k)*alp(ialp+m,irow)
              t3(2*m-1,k,latm) = t3(2*m-1,k,latm) + tmpr
              t3(2*m  ,k,latm) = t3(2*m  ,k,latm) + tmpi
              tl(2*m-1,k,latm) = tl(2*m-1,k,latm) - tmpi*ra
              tl(2*m  ,k,latm) = tl(2*m  ,k,latm) + tmpr*ra
!
              tmpr = t(ir,k)*dalp(ialp+m,irow)
              tmpi = t(ii,k)*dalp(ialp+m,irow)
              tm(2*m-1,k,latp) = tm(2*m-1,k,latp) + tmpr*ra
              tm(2*m  ,k,latp) = tm(2*m  ,k,latp) + tmpi*ra
!
              tmpr = q(ir,k)*alp(ialp+m,irow)
              tmpi = q(ii,k)*alp(ialp+m,irow)
              ql(2*m-1,k,latm) = ql(2*m-1,k,latm) - tmpi*ra
              ql(2*m  ,k,latm) = ql(2*m  ,k,latm) + tmpr*ra
!
              tmpr = q(ir,k)*dalp(ialp+m,irow)
              tmpi = q(ii,k)*dalp(ialp+m,irow)
              qm(2*m-1,k,latp) = qm(2*m-1,k,latp) + tmpr*ra
              qm(2*m  ,k,latp) = qm(2*m  ,k,latp) + tmpi*ra
!
              tmpr = d(ir,k)*alp(ialp+m,irow)
              tmpi = d(ii,k)*alp(ialp+m,irow)
              div(2*m-1,k,latm) = div(2*m-1,k,latm) + tmpr
              div(2*m  ,k,latm) = div(2*m  ,k,latm) + tmpi
           end do
        end do
!
        do n=2,pmax,2
           ic = ncoefi(n) - 1
           ialp = nalp(n)
!DIR$ IVDEP
           do m=1,nmreduced(n,irow)
              ir = 2*(ic+m) - 1
              ii = ir + 1
!
              tmpr = d(ir,k)*alpn(ialp+m)
              tmpi = d(ii,k)*alpn(ialp+m)
              u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpi
              u3(2*m  ,k,latp) = u3(2*m  ,k,latp) - tmpr
!     
              tmpr = d(ir,k)*dalpn(ialp+m)
              tmpi = d(ii,k)*dalpn(ialp+m)
              v3(2*m-1,k,latm) = v3(2*m-1,k,latm) - tmpr
              v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpi
!
              tmpr = vz(ir,k)*dalpn(ialp+m)
              tmpi = vz(ii,k)*dalpn(ialp+m)
              u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpr
              u3(2*m  ,k,latm) = u3(2*m  ,k,latm) + tmpi
!
              tmpr = vz(ir,k)*alpn(ialp+m)
              tmpi = vz(ii,k)*alpn(ialp+m)
              v3(2*m-1,k,latp) = v3(2*m-1,k,latp) + tmpi
              v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpr
!
              tmpr = t(ir,k)*alp(ialp+m,irow)
              tmpi = t(ii,k)*alp(ialp+m,irow)
              t3(2*m-1,k,latp) = t3(2*m-1,k,latp) + tmpr
              t3(2*m  ,k,latp) = t3(2*m  ,k,latp) + tmpi
              tl(2*m-1,k,latp) = tl(2*m-1,k,latp) - tmpi*ra
              tl(2*m  ,k,latp) = tl(2*m  ,k,latp) + tmpr*ra
!
              tmpr = t(ir,k)*dalp(ialp+m,irow)
              tmpi = t(ii,k)*dalp(ialp+m,irow)
              tm(2*m-1,k,latm) = tm(2*m-1,k,latm) + tmpr*ra
              tm(2*m  ,k,latm) = tm(2*m  ,k,latm) + tmpi*ra
!
              tmpr = q(ir,k)*alp(ialp+m,irow)
              tmpi = q(ii,k)*alp(ialp+m,irow)
              ql(2*m-1,k,latp) = ql(2*m-1,k,latp) - tmpi*ra
              ql(2*m  ,k,latp) = ql(2*m  ,k,latp) + tmpr*ra
!
              tmpr = q(ir,k)*dalp(ialp+m,irow)
              tmpi = q(ii,k)*dalp(ialp+m,irow)
              qm(2*m-1,k,latm) = qm(2*m-1,k,latm) + tmpr*ra
              qm(2*m  ,k,latm) = qm(2*m  ,k,latm) + tmpi*ra
!
              tmpr = d(ir,k)*alp(ialp+m,irow)
              tmpi = d(ii,k)*alp(ialp+m,irow)
              div(2*m-1,k,latp) = div(2*m-1,k,latp) + tmpr
              div(2*m  ,k,latp) = div(2*m  ,k,latp) + tmpi
           end do
        end do
#else
        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=1,nlen(m),2
              ir = mc + 2*n - 1
              ii = ir + 1
!
              tmpr = d(ir,k)*alpn(mr+n)
              tmpi = d(ii,k)*alpn(mr+n)
              u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpi
              u3(2*m  ,k,latm) = u3(2*m  ,k,latm) - tmpr
!
              tmpr = d(ir,k)*dalpn(mr+n)
              tmpi = d(ii,k)*dalpn(mr+n)
              v3(2*m-1,k,latp) = v3(2*m-1,k,latp) - tmpr
              v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpi
!
              tmpr = vz(ir,k)*dalpn(mr+n)
              tmpi = vz(ii,k)*dalpn(mr+n)
              u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpr
              u3(2*m  ,k,latp) = u3(2*m  ,k,latp) + tmpi
!
              tmpr = vz(ir,k)*alpn(mr+n)
              tmpi = vz(ii,k)*alpn(mr+n)
              v3(2*m-1,k,latm) = v3(2*m-1,k,latm) + tmpi
              v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpr
!
              tmpr = t(ir,k)*alp(mr+n,irow)
              tmpi = t(ii,k)*alp(mr+n,irow)
              t3(2*m-1,k,latm) = t3(2*m-1,k,latm) + tmpr
              t3(2*m  ,k,latm) = t3(2*m  ,k,latm) + tmpi
              tl(2*m-1,k,latm) = tl(2*m-1,k,latm) - tmpi*ra
              tl(2*m  ,k,latm) = tl(2*m  ,k,latm) + tmpr*ra
!
              tmpr = t(ir,k)*dalp(mr+n,irow)
              tmpi = t(ii,k)*dalp(mr+n,irow)
              tm(2*m-1,k,latp) = tm(2*m-1,k,latp) + tmpr*ra
              tm(2*m  ,k,latp) = tm(2*m  ,k,latp) + tmpi*ra
!
              tmpr = q(ir,k)*alp(mr+n,irow)
              tmpi = q(ii,k)*alp(mr+n,irow)
              ql(2*m-1,k,latm) = ql(2*m-1,k,latm) - tmpi*ra
              ql(2*m  ,k,latm) = ql(2*m  ,k,latm) + tmpr*ra
!
              tmpr = q(ir,k)*dalp(mr+n,irow)
              tmpi = q(ii,k)*dalp(mr+n,irow)
              qm(2*m-1,k,latp) = qm(2*m-1,k,latp) + tmpr*ra
              qm(2*m  ,k,latp) = qm(2*m  ,k,latp) + tmpi*ra
!
              tmpr = d(ir,k)*alp(mr+n,irow)
              tmpi = d(ii,k)*alp(mr+n,irow)
              div(2*m-1,k,latm) = div(2*m-1,k,latm) + tmpr
              div(2*m  ,k,latm) = div(2*m  ,k,latm) + tmpi
           end do
        end do

        do m=1,nmmax(irow)
           mr = nstart(m)
           mc = 2*mr
           do n=2,nlen(m),2
              ir = mc + 2*n - 1
              ii = ir + 1
!     
              tmpr = d(ir,k)*alpn(mr+n)
              tmpi = d(ii,k)*alpn(mr+n)
              u3(2*m-1,k,latp) = u3(2*m-1,k,latp) + tmpi
              u3(2*m  ,k,latp) = u3(2*m  ,k,latp) - tmpr
!
              tmpr = d(ir,k)*dalpn(mr+n)
              tmpi = d(ii,k)*dalpn(mr+n)
              v3(2*m-1,k,latm) = v3(2*m-1,k,latm) - tmpr
              v3(2*m  ,k,latm) = v3(2*m  ,k,latm) - tmpi
!
              tmpr = vz(ir,k)*dalpn(mr+n)
              tmpi = vz(ii,k)*dalpn(mr+n)
              u3(2*m-1,k,latm) = u3(2*m-1,k,latm) + tmpr
              u3(2*m  ,k,latm) = u3(2*m  ,k,latm) + tmpi
!
              tmpr = vz(ir,k)*alpn(mr+n)
              tmpi = vz(ii,k)*alpn(mr+n)
              v3(2*m-1,k,latp) = v3(2*m-1,k,latp) + tmpi
              v3(2*m  ,k,latp) = v3(2*m  ,k,latp) - tmpr
!
              tmpr = t(ir,k)*alp(mr+n,irow)
              tmpi = t(ii,k)*alp(mr+n,irow)
              t3(2*m-1,k,latp) = t3(2*m-1,k,latp) + tmpr
              t3(2*m  ,k,latp) = t3(2*m  ,k,latp) + tmpi
              tl(2*m-1,k,latp) = tl(2*m-1,k,latp) - tmpi*ra
              tl(2*m  ,k,latp) = tl(2*m  ,k,latp) + tmpr*ra
!
              tmpr = t(ir,k)*dalp(mr+n,irow)
              tmpi = t(ii,k)*dalp(mr+n,irow)
              tm(2*m-1,k,latm) = tm(2*m-1,k,latm) + tmpr*ra
              tm(2*m  ,k,latm) = tm(2*m  ,k,latm) + tmpi*ra
!
              tmpr = q(ir,k)*alp(mr+n,irow)
              tmpi = q(ii,k)*alp(mr+n,irow)
              ql(2*m-1,k,latp) = ql(2*m-1,k,latp) - tmpi*ra
              ql(2*m  ,k,latp) = ql(2*m  ,k,latp) + tmpr*ra
!
              tmpr = q(ir,k)*dalp(mr+n,irow)
              tmpi = q(ii,k)*dalp(mr+n,irow)
              qm(2*m-1,k,latm) = qm(2*m-1,k,latm) + tmpr*ra
              qm(2*m  ,k,latm) = qm(2*m  ,k,latm) + tmpi*ra
!
              tmpr = d(ir,k)*alp(mr+n,irow)
              tmpi = d(ii,k)*alp(mr+n,irow)
              div(2*m-1,k,latp) = div(2*m-1,k,latp) + tmpr
              div(2*m  ,k,latp) = div(2*m  ,k,latp) + tmpi
           end do
        end do
#endif
!
! d(T)/d(lamda)
! d(U)/d(lamda)
! d(V)/d(lamda)
!
!DIR$ IVDEP
        do m=1,nmmax(irow)
           tl(2*m-1,k,latm) = xm(m)*tl(2*m-1,k,latm)
           tl(2*m  ,k,latm) = xm(m)*tl(2*m  ,k,latm)
           tl(2*m-1,k,latp) = xm(m)*tl(2*m-1,k,latp)
           tl(2*m  ,k,latp) = xm(m)*tl(2*m  ,k,latp)
           ql(2*m-1,k,latm) = xm(m)*ql(2*m-1,k,latm)
           ql(2*m  ,k,latm) = xm(m)*ql(2*m  ,k,latm)
           ql(2*m-1,k,latp) = xm(m)*ql(2*m-1,k,latp)
           ql(2*m  ,k,latp) = xm(m)*ql(2*m  ,k,latp)
        end do
     end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
     do i=1,nlon(latm)+2
        tmp1 = phis(i,latm) + phis(i,latp)
        tmp2 = phis(i,latm) - phis(i,latp)
        phis(i,latm) = tmp1
        phis(i,latp) = tmp2
!
        tmp1 = phisl(i,latm) + phisl(i,latp)
        tmp2 = phisl(i,latm) - phisl(i,latp)
        phisl(i,latm) = tmp1
        phisl(i,latp) = tmp2
!
        tmp1 = phism(i,latm) + phism(i,latp)
        tmp2 = phism(i,latm) - phism(i,latp)
        phism(i,latm) = tmp1
        phism(i,latp) = tmp2
!
        tmp1 = ps(i,latm) + ps(i,latp)
        tmp2 = ps(i,latm) - ps(i,latp)
        ps(i,latm) = tmp1
        ps(i,latp) = tmp2
!
        tmp1 = dpsl(i,latm) + dpsl(i,latp)
        tmp2 = dpsl(i,latm) - dpsl(i,latp)
        dpsl(i,latm) = tmp1
        dpsl(i,latp) = tmp2
!
        tmp1 = dpsm(i,latm) + dpsm(i,latp)
        tmp2 = dpsm(i,latm) - dpsm(i,latp)
        dpsm(i,latm) = tmp1
        dpsm(i,latp) = tmp2
     end do
!
     do k=1,plev
        do i=1,nlon(latm)+2
           tmp1 = u3(i,k,latm) + u3(i,k,latp)
           tmp2 = u3(i,k,latm) - u3(i,k,latp)
           u3(i,k,latm) = tmp1
           u3(i,k,latp) = tmp2
!
           tmp1 = v3(i,k,latm) + v3(i,k,latp)
           tmp2 = v3(i,k,latm) - v3(i,k,latp)
           v3(i,k,latm) = tmp1
           v3(i,k,latp) = tmp2
!
           tmp1 = t3(i,k,latm) + t3(i,k,latp)
           tmp2 = t3(i,k,latm) - t3(i,k,latp)
           t3(i,k,latm) = tmp1
           t3(i,k,latp) = tmp2
!
           tmp1 = tl(i,k,latm) + tl(i,k,latp)
           tmp2 = tl(i,k,latm) - tl(i,k,latp)
           tl(i,k,latm) = tmp1
           tl(i,k,latp) = tmp2
!
           tmp1 = tm(i,k,latm) + tm(i,k,latp)
           tmp2 = tm(i,k,latm) - tm(i,k,latp)
           tm(i,k,latm) = tmp1
           tm(i,k,latp) = tmp2
!
           tmp1 = ql(i,k,latm) + ql(i,k,latp)
           tmp2 = ql(i,k,latm) - ql(i,k,latp)
           ql(i,k,latm) = tmp1
           ql(i,k,latp) = tmp2
!
           tmp1 = qm(i,k,latm) + qm(i,k,latp)
           tmp2 = qm(i,k,latm) - qm(i,k,latp)
           qm(i,k,latm) = tmp1
           qm(i,k,latp) = tmp2
!
           tmp1 = div(i,k,latm) + div(i,k,latp)
           tmp2 = div(i,k,latm) - div(i,k,latp)
           div(i,k,latm) = tmp1
           div(i,k,latp) = tmp2
        end do
     end do
  end do
!
! 2nd pass through initial data to obtain and merge all untruncated
! fields:read in and store the data which do not need to be spectrally
! truncated, skipping over header records first.  Also complete
! initialization of prognostics.F
!     
!
  do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
! 1st transform: U,V,T
! 2nd: ln(PS). 3rd: PHIS. 4th: longitudinal derivative of ln(PS)
! 5th: meridional derivative of ln(PS)
! 6th: divergence
!
     irow = lat
     if (lat.gt.plat/2) irow = plat - lat + 1

     call fft991 (u3   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (v3   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (t3   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (ps   (1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
     call fft991 (phis (1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
     call fft991 (dpsl (1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
     call fft991 (dpsm (1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
     call fft991 (div  (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
!
! Still more fft's
!
! 1st: zonal t derivative
! 2nd: meridional t derivative
! 3rd: zonal phis derivative
! 4th: meridional phis derivative
!
     call fft991 (tl   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (tm   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (ql   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (qm   (1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),plev        ,+1      )
     call fft991 (phisl(1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
     call fft991 (phism(1,lat)   ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plond       ,nlon(lat),1           ,+1      )
!
! Convert U,V to u,v
!
     zsqcs = sqrt(cs(irow))
     do k=1,plev
        do i=1,nlon(lat)
           u3(i,k,lat) = u3(i,k,lat)/zsqcs
           v3(i,k,lat) = v3(i,k,lat)/zsqcs
        end do
     end do
!
! Convert from ln(ps) to ps
!
     do i=1,nlon(lat)
        ps(i,lat) = exp(ps(i,lat))
     end do
  end do

  return
end subroutine spetru
