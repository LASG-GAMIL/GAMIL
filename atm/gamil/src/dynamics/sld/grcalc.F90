#include <misc.h>
#include <params.h>

! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  Make sure to make appropriate coding changes where
! necessary.

#if ( defined PVP )

subroutine grcalcs (irow    ,ztodt   ,grts    ,grqs    ,grths   , &
                    grds    ,grus    ,gruhs   ,grvs    ,grvhs   , &
                    grpss   ,grdpss  ,grpms   ,grpls   ,grtms   , &
                    grtls   ,grqms   ,grqls   )
!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.8.1 2002/06/15 13:48:23 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra

   implicit none

#include <comhd.h>
!
! Input arguments
!
   integer , intent(in)   :: irow              ! latitude pair index
   real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: symmetric fourier coefficients
!
   real(r8), intent(out) :: grts(plond,plev)  ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grqs(plond,plev)  ! sum(n) of q(n,m)*P(n,m)
   real(r8), intent(out) :: grths(plond,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grds(plond,plev)  ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grus(plond,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruhs(plond,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvs(plond,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvhs(plond,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpss(plond)      ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpss(plond)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpms(plond)    ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpls(plond)    ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(out) :: grtms (plond,plev)
   real(r8), intent(out) :: grtls (plond,plev)
   real(r8), intent(out) :: grqms (plond,plev)
   real(r8), intent(out) :: grqls (plond,plev)
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1s (plond,plev)   ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1s(plond,plev)   ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1s (plond,plev)   ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1s(plond,plev)   ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) zdfac (2*pnmax,plev) ! horiz. diffusion factor (vort,div) (complex)
   real(r8) tqfac (2*pnmax,plev) ! horiz. diffusion factor (t,q) (complex)
   real(r8) alp2  (2*pspt)       ! Legendre functions (complex)
   real(r8) dalp2 (2*pspt)       ! derivative of Legendre functions (complex)
   real(r8) alpn2 (2*pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn2(2*pspt)       ! (a/(n(n+1)))*derivative of Legendre functions (complex)
!                                ! betwn absolute and relative vorticity
   integer k                     ! level index
   integer m                     ! diagonal element(index) of spec. array
   integer n                     ! meridional wavenumber index
   integer ne                    ! index into spectral arrays
   integer mn                    ! index into spectral arrays
   integer mnc                   ! index into spectral arrays
   integer mnev                  ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
! Expand polynomials and derivatives to complex form to allow largest 
! possible vector length and multiply by appropriate factors
!
   do n=1,pmax
      ne = n - 1
!dir$ ivdep
      do m=1,nmreduced(n,irow)
         mnc = 2*(m+nalp(n))
         mn = m + nalp(n)
         alp2(mnc-1) = alp(mn,irow)
         alp2(mnc  ) = alp(mn,irow)
         dalp2(mnc-1) = dalp(mn,irow)*ra
         dalp2(mnc  ) = dalp(mn,irow)*ra
         alpn2(mnc-1) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         alpn2(mnc  ) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         dalpn2(mnc-1) = dalp(mn,irow)*(rsq(m+ne)*ra)
         dalpn2(mnc  ) = dalp(mn,irow)*(rsq(m+ne)*ra)
      end do
   end do
!
! Initialize sums
!
   grts(:,:) = 0.
   grqs(:,:) = 0.
   grths(:,:) = 0.
   grds(:,:)  = 0.
   grus(:,:)  = 0.
   gruhs(:,:) = 0.
   grvs(:,:)  = 0.
   grvhs(:,:) = 0.
   grpss(:)   = 0.
   grdpss(:)   = 0.
   grpms(:)   = 0.
   grpls(:)   = 0.
   grtms(:,:)   = 0.
   grtls(:,:)   = 0.
   grqms(:,:)   = 0.
   grqls(:,:)   = 0.
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
   do k=1,plev
!
! Diffusion factors: expand for longest possible vectors
!
!dir$ ivdep
      do n = 1,pnmax
         zdfac(n*2-1,k) = -hdifzd(n,k)
         zdfac(n*2  ,k) = -hdifzd(n,k)
         tqfac(n*2-1,k) = -hdiftq(n,k)
         tqfac(n*2  ,k) = -hdiftq(n,k)
      end do
!
! Initialize local sums
!
      gru1s(:,k) = 0.
      gruh1s(:,k) = 0.
      grv1s(:,k) = 0.
      grvh1s(:,k) = 0.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n = 1,ncutoff,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            grts  (m,k) = grts  (m,k) + t (mnev,k)*alp2  (mn)
            grths (m,k) = grths (m,k) + t (mnev,k)*alp2  (mn)*tqfac(m+ne,k)
            grqs  (m,k) = grqs  (m,k) + q (mnev,k)*alp2  (mn)
            grds  (m,k) = grds  (m,k) + d (mnev,k)*alp2  (mn)
            gru1s (m,k) = gru1s (m,k) + d (mnev,k)*alpn2 (mn)
            gruh1s(m,k) = gruh1s(m,k) + d (mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
            grv1s (m,k) = grv1s (m,k) + vz(mnev,k)*alpn2 (mn)
            grvh1s(m,k) = grvh1s(m,k) + vz(mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
         end do
      end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n = 2,ncutoff,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            grtms (m,k) = grtms (m,k) + t (mnev,k)*dalp2 (mn)
            grqms (m,k) = grqms (m,k) + q (mnev,k)*dalp2 (mn)
            grus  (m,k) = grus  (m,k) + vz(mnev,k)*dalpn2(mn)
            gruhs (m,k) = gruhs (m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            grvs  (m,k) = grvs  (m,k) - d (mnev,k)*dalpn2(mn)
            grvhs (m,k) = grvhs (m,k) - d (mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
         end do
      end do
   end do                    ! k=1,plev
!
! For short diagonals, repeat above loops with vectorization in
! vertical, instead of along diagonals, to keep vector lengths from
! getting too short.
!
   if (ncutoff.lt.pmax) then
      do n = ncutoff+1,pmax,2 ! ncutoff guaranteed even
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            do k = 1,plev
               grts  (m,k) = grts  (m,k) + t (mnev,k)*alp2  (mn)
               grths (m,k) = grths (m,k) + t (mnev,k)*alp2  (mn)*tqfac(m+ne,k)
               grqs  (m,k) = grqs  (m,k) + q (mnev,k)*alp2  (mn)
               grds  (m,k) = grds  (m,k) + d (mnev,k)*alp2  (mn)
               gru1s (m,k) = gru1s (m,k) + d (mnev,k)*alpn2 (mn)
               gruh1s(m,k) = gruh1s(m,k) + d (mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
               grv1s (m,k) = grv1s (m,k) + vz(mnev,k)*alpn2 (mn)
               grvh1s(m,k) = grvh1s(m,k) + vz(mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
            end do
         end do
      end do
      
      do n = ncutoff+2,pmax,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            do k = 1,plev
               grtms (m,k) = grtms (m,k) + t (mnev,k)*dalp2 (mn)
               grqms (m,k) = grqms (m,k) + q (mnev,k)*dalp2 (mn)
               grus  (m,k) = grus  (m,k) + vz(mnev,k)*dalpn2(mn)
               gruhs (m,k) = gruhs (m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
               grvs  (m,k) = grvs  (m,k) - d (mnev,k)*dalpn2(mn)
               grvhs (m,k) = grvhs (m,k) - d (mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do
   end if                    ! ncutoff.lt.pmax

   do k=1,plev
!
! Combine the two parts of u(m) and v(m) and compute derivatives
!
!dir$ ivdep
      do m=1,nmmax(irow)
         grus (2*m-1,k) = grus (2*m-1,k) + gru1s (2*m  ,k)
         gruhs(2*m-1,k) = gruhs(2*m-1,k) + gruh1s(2*m  ,k)
         grus (2*m  ,k) = grus (2*m  ,k) - gru1s (2*m-1,k)
         gruhs(2*m  ,k) = gruhs(2*m  ,k) - gruh1s(2*m-1,k)
         grvs (2*m-1,k) = grvs (2*m-1,k) + grv1s (2*m  ,k)
         grvhs(2*m-1,k) = grvhs(2*m-1,k) + grvh1s(2*m  ,k)
         grvs (2*m  ,k) = grvs (2*m  ,k) - grv1s (2*m-1,k)
         grvhs(2*m  ,k) = grvhs(2*m  ,k) - grvh1s(2*m-1,k)
!
! Derivatives
!
         grtls(2*m-1,k) = -grts(2*m  ,k)*ra*xm(m)
         grtls(2*m  ,k) =  grts(2*m-1,k)*ra*xm(m)
         grqls(2*m-1,k) = -grqs(2*m  ,k)*ra*xm(m)
         grqls(2*m  ,k) =  grqs(2*m-1,k)*ra*xm(m)
      end do
   end do
!
!-----------------------------------------------------------------------
! Computation for single level variables.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=1,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpss (m) = grpss (m) + alps(mnev)*alp2(mn)
         grdpss(m) = grdpss(m) + alps(mnev)*alp2(mn)*hdfst4(ne+(m+1)/2)*ztodt
      end do
   end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=2,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpms(m) = grpms(m) + alps(mnev)*dalp2(mn)
      end do
   end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
   do m=1,nmmax(irow)
      grpls(2*m-1) = -grpss(2*m  )*ra*xm(m)
      grpls(2*m  ) =  grpss(2*m-1)*ra*xm(m)
   end do

   return
end subroutine grcalcs


subroutine grcalca (irow    ,ztodt   ,grta    ,grqa    ,grtha   , &
                    grda    ,grua    ,gruha   ,grva    ,grvha   , &
                    grpsa   ,grdpsa  ,grpma   ,grpla   ,grtma   , &
                    grtla   ,grqma   ,grqla   )

!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.8.1 2002/06/15 13:48:23 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use dynconst, only: ra

   implicit none

#include <comhd.h>
!
! Input arguments
!
   integer , intent(in)   :: irow              ! latitude pair index
   real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: anti-symmetric fourier coefficients
!
   real(r8), intent(out) :: grta(plond,plev)  ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grqa(plond,plev)  ! sum(n) of q(n,m)*P(n,m)
   real(r8), intent(out) :: grtha(plond,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grda(plond,plev)  ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grua(plond,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruha(plond,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grva(plond,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvha(plond,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpsa(plond)      ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpsa(plond)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt
!                                              ! *lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpma(plond)      ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpla(plond)      ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(out) :: grtma (plond,plev)
   real(r8), intent(out) :: grtla (plond,plev)
   real(r8), intent(out) :: grqma (plond,plev)
   real(r8), intent(out) :: grqla (plond,plev)
!
!---------------------------Local workspace-----------------------------
!
   real(r8) gru1a (plond,plev)   ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) gruh1a(plond,plev)   ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grv1a (plond,plev)   ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) grvh1a(plond,plev)   ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
   real(r8) zdfac (2*pnmax,plev) ! horiz. diffusion factor (vort,div) (complex)
   real(r8) tqfac (2*pnmax,plev) ! horiz. diffusion factor (t,q) (complex)
   real(r8) alp2  (2*pspt)       ! Legendre functions (complex)
   real(r8) dalp2 (2*pspt)       ! derivative of Legendre functions (complex)
   real(r8) alpn2 (2*pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
   real(r8) dalpn2(2*pspt)       ! (a/(n(n+1)))*derivative of Legendre 

   integer k                     ! level index
   integer m                     ! diagonal element(index) of spec. array
   integer n                     ! meridional wavenumber index
   integer ne                    ! index into spectral arrays
   integer mn                    ! index into spectral arrays
   integer mnc                   ! index into spectral arrays
   integer mnev                  ! index into spectral arrays
!
!-----------------------------------------------------------------------
! Compute alpn and dalpn
! Expand polynomials and derivatives to complex form to allow largest 
! possible vector length and multiply by appropriate factors
!
   do n=1,pmax
      ne = n - 1
!dir$ ivdep
      do m=1,nmreduced(n,irow)
         mnc = 2*(m+nalp(n))
         mn = m + nalp(n)
         alp2(mnc-1) = alp(mn,irow)
         alp2(mnc  ) = alp(mn,irow)
         dalp2(mnc-1) = dalp(mn,irow)*ra
         dalp2(mnc  ) = dalp(mn,irow)*ra
         alpn2(mnc-1) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         alpn2(mnc  ) = alp(mn,irow)*(rsq(m+ne)*ra)*xm(m)
         dalpn2(mnc-1) = dalp(mn,irow)*(rsq(m+ne)*ra)
         dalpn2(mnc  ) = dalp(mn,irow)*(rsq(m+ne)*ra)
      end do
   end do
!
! Initialize sums
!
   grta(:,:) = 0.
   grqa(:,:) = 0.
   grtha(:,:) = 0.
   grda(:,:)  = 0.
   grua(:,:)  = 0.
   gruha(:,:) = 0.
   grva(:,:)  = 0.
   grvha(:,:) = 0.
   grpsa(:)   = 0.
   grdpsa(:)   = 0.
   grpma(:)   = 0.
   grpla(:)   = 0.
   grtma(:,:)   = 0.
   grtla(:,:)   = 0.
   grqma(:,:)   = 0.
   grqla(:,:)   = 0.
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
   do k=1,plev
!
! Diffusion factors: expand for longest possible vectors
!
!dir$ ivdep
      do n = 1,pnmax
         zdfac(n*2-1,k) = -hdifzd(n,k)
         zdfac(n*2  ,k) = -hdifzd(n,k)
         tqfac(n*2-1,k) = -hdiftq(n,k)
         tqfac(n*2  ,k) = -hdiftq(n,k)
      end do
!
! Initialize local sums
!
      gru1a(:,k) = 0.
      gruh1a(:,k) = 0.
      grv1a(:,k) = 0.
      grvh1a(:,k) = 0.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n = 1,ncutoff,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            grtma (m,k) = grtma (m,k) + t (mnev,k)*dalp2 (mn)
            grqma (m,k) = grqma (m,k) + q (mnev,k)*dalp2 (mn)
            grua  (m,k) = grua  (m,k) + vz(mnev,k)*dalpn2(mn)
            gruha (m,k) = gruha (m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            grva  (m,k) = grva  (m,k) - d (mnev,k)*dalpn2(mn)
            grvha (m,k) = grvha (m,k) - d (mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
         end do
      end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H. Loop over n for t(m), q(m), d(m),and the two parts of u(m) and v(m).
! The inner (vector) loop accumulates sums over n along the diagonals
! of the spectral truncation to obtain the maximum length vectors.
!
! "ncutoff" is used to switch to vectorization in the vertical when the 
! length of the diagonal is less than the number of levels.
!
      do n = 2,ncutoff,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            grta  (m,k) = grta  (m,k) + t (mnev,k)*alp2  (mn)
            grtha (m,k) = grtha (m,k) + t (mnev,k)*alp2  (mn)*tqfac(m+ne,k)
            grqa  (m,k) = grqa  (m,k) + q (mnev,k)*alp2  (mn)
            grda  (m,k) = grda  (m,k) + d (mnev,k)*alp2  (mn)
            gru1a (m,k) = gru1a (m,k) + d (mnev,k)*alpn2 (mn)
            gruh1a(m,k) = gruh1a(m,k) + d (mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
            grv1a (m,k) = grv1a (m,k) + vz(mnev,k)*alpn2 (mn)
            grvh1a(m,k) = grvh1a(m,k) + vz(mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
         end do
      end do
   end do                    ! k=1,plev
!
! For short diagonals, repeat above loops with vectorization in
! vertical, instead of along diagonals, to keep vector lengths from
! getting too short.
!
   if (ncutoff.lt.pmax) then
      do n = ncutoff+1,pmax,2 ! ncutoff guaranteed even
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            do k = 1,plev
               grtma (m,k) = grtma (m,k) + t (mnev,k)*dalp2 (mn)
               grqma (m,k) = grqma (m,k) + q (mnev,k)*dalp2 (mn)
               grua  (m,k) = grua  (m,k) + vz(mnev,k)*dalpn2(mn)
               gruha (m,k) = gruha (m,k) + vz(mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
               grva  (m,k) = grva  (m,k) - d (mnev,k)*dalpn2(mn)
               grvha (m,k) = grvha (m,k) - d (mnev,k)*dalpn2(mn)*zdfac(m+ne,k)
            end do
         end do
      end do
      do n = ncutoff+2,pmax,2
         ne = 2*(n-1)
         do m = 1,2*nmreduced(n,irow)
            mnev = m + nco2(n) - 2
            mn   = m + 2*nalp(n)
            do k = 1,plev
               grta  (m,k) = grta  (m,k) + t (mnev,k)*alp2  (mn)
               grtha (m,k) = grtha (m,k) + t (mnev,k)*alp2  (mn)*tqfac(m+ne,k)
               grqa  (m,k) = grqa  (m,k) + q (mnev,k)*alp2  (mn)
               grda  (m,k) = grda  (m,k) + d (mnev,k)*alp2  (mn)
               gru1a (m,k) = gru1a (m,k) + d (mnev,k)*alpn2 (mn)
               gruh1a(m,k) = gruh1a(m,k) + d (mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
               grv1a (m,k) = grv1a (m,k) + vz(mnev,k)*alpn2 (mn)
               grvh1a(m,k) = grvh1a(m,k) + vz(mnev,k)*alpn2 (mn)*zdfac(m+ne,k)
            end do
         end do
      end do
   end if                    ! ncutoff.lt.pmax

   do k=1,plev
!
! Combine the two parts of u(m) and v(m) and compute derivatives
!
!dir$ ivdep
      do m=1,nmmax(irow)
         grua (2*m-1,k) = grua (2*m-1,k) + gru1a (2*m  ,k)
         gruha(2*m-1,k) = gruha(2*m-1,k) + gruh1a(2*m  ,k)
         grua (2*m  ,k) = grua (2*m  ,k) - gru1a (2*m-1,k)
         gruha(2*m  ,k) = gruha(2*m  ,k) - gruh1a(2*m-1,k)
         grva (2*m-1,k) = grva (2*m-1,k) + grv1a (2*m  ,k)
         grvha(2*m-1,k) = grvha(2*m-1,k) + grvh1a(2*m  ,k)
         grva (2*m  ,k) = grva (2*m  ,k) - grv1a (2*m-1,k)
         grvha(2*m  ,k) = grvha(2*m  ,k) - grvh1a(2*m-1,k)
!
! Derivatives
!
         grtla(2*m-1,k) = -grta(2*m  ,k)*ra*xm(m)
         grtla(2*m  ,k) =  grta(2*m-1,k)*ra*xm(m)
         grqla(2*m-1,k) = -grqa(2*m  ,k)*ra*xm(m)
         grqla(2*m  ,k) =  grqa(2*m-1,k)*ra*xm(m)
      end do
   end do
!
!-----------------------------------------------------------------------
! Computation for single level variables.
!
! Evaluate symmetric components involving P and antisymmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=1,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpma(m) = grpma(m) + alps(mnev)*dalp2(mn)
      end do
   end do
!
! Evaluate antisymmetric components involving P and symmetric involving 
! H.  Loop over n for lnps(m) and derivatives.
! The inner loop accumulates over n along diagonal of the truncation.
!
   do n=2,pmax,2
      ne = n - 1
      do m=1,2*nmreduced(n,irow)
         mnev = m + nco2(n) - 2
         mn = m + 2*nalp(n)
         grpsa (m) = grpsa (m) + alps(mnev)*alp2(mn)
         grdpsa(m) = grdpsa(m) + alps(mnev)*alp2(mn)*hdfst4(ne+(m+1)/2)*ztodt
      end do
   end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
   do m=1,nmmax(irow)
      grpla(2*m-1) = -grpsa(2*m  )*ra*xm(m)
      grpla(2*m  ) =  grpsa(2*m-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalca

#else

subroutine grcalcs (irow    ,ztodt   ,grts    ,grqs    ,grths   , &
                    grds    ,grus    ,gruhs   ,grvs    ,grvhs   , &
                    grpss   ,grdpss  ,grpms   ,grpls   ,grtms   , &
                    grtls   ,grqms   ,grqls   )
!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.8.1 2002/06/15 13:48:23 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use pspect
    use comspe
    use rgrid
    use commap
    use dynconst, only: ra

    implicit none

#include <comhd.h>
!
! Input arguments
!
    integer , intent(in)   :: irow              ! latitude pair index
    real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: symmetric fourier coefficients
!
    real(r8), intent(out) :: grts(plond,plev)  ! sum(n) of t(n,m)*P(n,m)
    real(r8), intent(out) :: grqs(plond,plev)  ! sum(n) of q(n,m)*P(n,m)
    real(r8), intent(out) :: grths(plond,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
    real(r8), intent(out) :: grds(plond,plev)  ! sum(n) of d(n,m)*P(n,m)
    real(r8), intent(out) :: grus(plond,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: gruhs(plond,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvs(plond,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvhs(plond,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grpss(plond)      ! sum(n) of lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grdpss(plond)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grpms(plond)    ! sum(n) of lnps(n,m)*H(n,m)
    real(r8), intent(out) :: grpls(plond)    ! sum(n) of lnps(n,m)*P(n,m)*m/a
    real(r8), intent(out) :: grtms (plond,plev)
    real(r8), intent(out) :: grtls (plond,plev)
    real(r8), intent(out) :: grqms (plond,plev)
    real(r8), intent(out) :: grqls (plond,plev)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) gru1s (plond)      ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) gruh1s(plond)      ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grv1s (plond)      ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grvh1s(plond)      ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) alpn  (pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
    real(r8) dalpn (pspt)       ! (a/(n(n+1)))*derivative of Legendre functions (complex)

    integer k                   ! level index
    integer m                   ! Fourier wavenumber index of spectral array
    integer n                   ! meridional wavenumber index
    integer ir,ii               ! spectral indices
    integer mr,mc               ! spectral indices
    real(r8) tmp,tmpr,tmpi,raxm ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
    do m=1,nmmax(irow)
       mr = nstart(m)
       raxm = ra*xm(m)
       do n=1,nlen(m)
          alpn(mr+n) = alp(mr+n,irow)*rsq(m+n-1)*raxm
          dalpn(mr+n) = dalp(mr+n,irow)*rsq(m+n-1)*ra
       end do
    end do
!
! Initialize sums
!
    grts(:,:) = 0.
    grqs(:,:) = 0.
    grths(:,:) = 0.
    grds(:,:)  = 0.
    grus(:,:)  = 0.
    gruhs(:,:) = 0.
    grvs(:,:)  = 0.
    grvhs(:,:) = 0.
    grpss(:)   = 0.
    grdpss(:)   = 0.
    grpms(:)   = 0.
    grpls(:)   = 0.
    grtms(:,:)   = 0.
    grtls(:,:)   = 0.
    grqms(:,:)   = 0.
    grqls(:,:)   = 0.
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
    do k=1,plev
!
! Initialize local sums
!
       gru1s(:) = 0.
       gruh1s(:) = 0.
       grv1s(:) = 0.
       grvh1s(:) = 0.
!
! Loop over n for t,q,d,and end of u and v
!
       do m=1,nmmax(irow)
          mr = nstart(m)
          mc = 2*mr
          do n=1,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1

             grts (2*m-1,k) = grts (2*m-1,k) + t(ir,k)*alp(mr+n,irow)
             grts (2*m  ,k) = grts (2*m  ,k) + t(ii,k)*alp(mr+n,irow)
!
             grqs (2*m-1,k) = grqs (2*m-1,k) + q(ir,k)*alp(mr+n,irow)
             grqs (2*m  ,k) = grqs (2*m  ,k) + q(ii,k)*alp(mr+n,irow)
!
             tmp = alp(mr+n,irow)*hdiftq(n+m-1,k)
             grths(2*m-1,k) = grths(2*m-1,k) - t(ir,k)*tmp
             grths(2*m  ,k) = grths(2*m  ,k) - t(ii,k)*tmp
!
             grds(2*m-1,k) = grds(2*m-1,k) + d(ir,k)*alp(mr+n,irow)
             grds(2*m  ,k) = grds(2*m  ,k) + d(ii,k)*alp(mr+n,irow)
!
             gru1s (2*m-1) = gru1s (2*m-1) + d(ir,k)*alpn(mr+n)
             gru1s (2*m  ) = gru1s (2*m  ) + d(ii,k)*alpn(mr+n)
!
             tmp = alpn(mr+n)*hdifzd(n+m-1,k)
             gruh1s(2*m-1) = gruh1s(2*m-1) - d(ir,k)*tmp
             gruh1s(2*m  ) = gruh1s(2*m  ) - d(ii,k)*tmp
!
             grv1s (2*m-1) = grv1s (2*m-1) + vz(ir,k)*alpn(mr+n)
             grv1s (2*m  ) = grv1s (2*m  ) + vz(ii,k)*alpn(mr+n)
!
             grvh1s(2*m-1) = grvh1s(2*m-1) - vz(ir,k)*tmp
             grvh1s(2*m  ) = grvh1s(2*m  ) - vz(ii,k)*tmp
          end do
       end do
       do m=1,nmmax(irow)
          mr = nstart(m)
          mc = 2*mr
          do n=2,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
!
             grtms(2*m-1,k) = grtms(2*m-1,k) + t(ir,k)*dalp(mr+n,irow)*ra
             grtms(2*m  ,k) = grtms(2*m  ,k) + t(ii,k)*dalp(mr+n,irow)*ra
!
             grqms(2*m-1,k) = grqms(2*m-1,k) + q(ir,k)*dalp(mr+n,irow)*ra
             grqms(2*m  ,k) = grqms(2*m  ,k) + q(ii,k)*dalp(mr+n,irow)*ra
!
             grus (2*m-1,k) = grus (2*m-1,k) + vz(ir,k)*dalpn(mr+n)
             grus (2*m  ,k) = grus (2*m  ,k) + vz(ii,k)*dalpn(mr+n)
!
             tmp = dalpn(mr+n)*hdifzd(n+m-1,k)
             gruhs(2*m-1,k) = gruhs(2*m-1,k) - vz(ir,k)*tmp
             gruhs(2*m  ,k) = gruhs(2*m  ,k) - vz(ii,k)*tmp
!
             grvs (2*m-1,k) = grvs (2*m-1,k) - d(ir,k)*dalpn(mr+n)
             grvs (2*m  ,k) = grvs (2*m  ,k) - d(ii,k)*dalpn(mr+n)
!
             grvhs(2*m-1,k) = grvhs(2*m-1,k) + d(ir,k)*tmp
             grvhs(2*m  ,k) = grvhs(2*m  ,k) + d(ii,k)*tmp
          end do
       end do
!
! Combine the two parts of u(m) and v(m)
!
       do m=1,nmmax(irow)
          grus (2*m-1,k) = grus (2*m-1,k) + gru1s (2*m  )
          gruhs(2*m-1,k) = gruhs(2*m-1,k) + gruh1s(2*m  )
          grus (2*m  ,k) = grus (2*m  ,k) - gru1s (2*m-1)
          gruhs(2*m  ,k) = gruhs(2*m  ,k) - gruh1s(2*m-1)
          grvs (2*m-1,k) = grvs (2*m-1,k) + grv1s (2*m  )
          grvhs(2*m-1,k) = grvhs(2*m-1,k) + grvh1s(2*m  )
          grvs (2*m  ,k) = grvs (2*m  ,k) - grv1s (2*m-1)
          grvhs(2*m  ,k) = grvhs(2*m  ,k) - grvh1s(2*m-1)
!
! Derivatives
!
          grtls(2*m-1,k) = -grts(2*m  ,k)*ra*xm(m)
          grtls(2*m  ,k) =  grts(2*m-1,k)*ra*xm(m)
          grqls(2*m-1,k) = -grqs(2*m  ,k)*ra*xm(m)
          grqls(2*m  ,k) =  grqs(2*m-1,k)*ra*xm(m)
       end do
    end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
    do m=1,nmmax(irow)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          tmpr = alps(ir)*alp(mr+n,irow)
          tmpi = alps(ii)*alp(mr+n,irow)
          grpss(2*m-1) = grpss(2*m-1) + tmpr
          grpss(2*m  ) = grpss(2*m  ) + tmpi
!
          grdpss(2*m-1) = grdpss(2*m-1) + tmpr*hdfst4(m+n-1)*ztodt
          grdpss(2*m  ) = grdpss(2*m  ) + tmpi*hdfst4(m+n-1)*ztodt
       end do
    end do
    do m=1,nmmax(irow)
       mr = nstart(m)
       mc = 2*mr
       do n=2,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpms(2*m-1) = grpms(2*m-1) + alps(ir)*dalp(mr+n,irow)*ra
          grpms(2*m  ) = grpms(2*m  ) + alps(ii)*dalp(mr+n,irow)*ra
       end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
       grpls(2*m-1) = -grpss(2*m  )*ra*xm(m)
       grpls(2*m  ) =  grpss(2*m-1)*ra*xm(m)
    end do
!
    return
end subroutine grcalcs


subroutine grcalca (irow    ,ztodt   ,grta    ,grqa    ,grtha   , &
                    grda    ,grua    ,gruha   ,grva    ,grvha   , &
                    grpsa   ,grdpsa  ,grpma   ,grpla   ,grtma   , &
                    grtla   ,grqma   ,grqla   )

!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: grcalc.F90,v 1.5.8.1 2002/06/15 13:48:23 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use pspect
    use comspe
    use rgrid
    use commap
    use dynconst, only: ra

    implicit none

#include <comhd.h>
!
! Input arguments
!
    integer , intent(in)   :: irow              ! latitude pair index
    real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: anti-symmetric fourier coefficients
!
    real(r8), intent(out) :: grta(plond,plev)  ! sum(n) of t(n,m)*P(n,m)
    real(r8), intent(out) :: grqa(plond,plev)  ! sum(n) of q(n,m)*P(n,m)
    real(r8), intent(out) :: grtha(plond,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
    real(r8), intent(out) :: grda(plond,plev)  ! sum(n) of d(n,m)*P(n,m)
    real(r8), intent(out) :: grua(plond,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: gruha(plond,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grva(plond,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvha(plond,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grpsa(plond)      ! sum(n) of lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grdpsa(plond)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt
!                                               ! *lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grpma(plond)      ! sum(n) of lnps(n,m)*H(n,m)
    real(r8), intent(out) :: grpla(plond)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
    real(r8), intent(out) :: grtma (plond,plev)
    real(r8), intent(out) :: grtla (plond,plev)
    real(r8), intent(out) :: grqma (plond,plev)
    real(r8), intent(out) :: grqla (plond,plev)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) gru1a (plond)      ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) gruh1a(plond)      ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grv1a (plond)      ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grvh1a(plond)      ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) alpn  (pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
    real(r8) dalpn (pspt)       ! (a/(n(n+1)))*derivative of Legendre functions (complex)

    integer k                   ! level index
    integer m                   ! Fourier wavenumber index of spectral array
    integer n                   ! meridional wavenumber index
    integer ir,ii               ! spectral indices
    integer mr,mc               ! spectral indices
    real(r8) tmp,tmpr,tmpi,raxm ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
    do m=1,nmmax(irow)
       mr = nstart(m)
       raxm = ra*xm(m)
       do n=1,nlen(m)
          alpn(mr+n) = alp(mr+n,irow)*rsq(m+n-1)*raxm
          dalpn(mr+n) = dalp(mr+n,irow)*rsq(m+n-1)*ra
       end do
    end do
!
! Initialize sums
!
    grta(:,:) = 0.
    grqa(:,:) = 0.
    grtha(:,:) = 0.
    grda(:,:)  = 0.
    grua(:,:)  = 0.
    gruha(:,:) = 0.
    grva(:,:)  = 0.
    grvha(:,:) = 0.
    grpsa(:)   = 0.
    grdpsa(:)   = 0.
    grpma(:)   = 0.
    grpla(:)   = 0.
    grtma(:,:)   = 0.
    grtla(:,:)   = 0.
    grqma(:,:)   = 0.
    grqla(:,:)   = 0.
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
    do k=1,plev
!
! Initialize local sums
!
       gru1a(:) = 0.
       gruh1a(:) = 0.
       grv1a(:) = 0.
       grvh1a(:) = 0.
!
! Loop over n for t,q,d,and end of u and v
!
       do m=1,nmmax(irow)
          mr = nstart(m)
          mc = 2*mr
          do n=1,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
!
             grtma(2*m-1,k) = grtma(2*m-1,k) + t(ir,k)*dalp(mr+n,irow)*ra
             grtma(2*m  ,k) = grtma(2*m  ,k) + t(ii,k)*dalp(mr+n,irow)*ra
!
             grqma(2*m-1,k) = grqma(2*m-1,k) + q(ir,k)*dalp(mr+n,irow)*ra
             grqma(2*m  ,k) = grqma(2*m  ,k) + q(ii,k)*dalp(mr+n,irow)*ra
!
             grua (2*m-1,k) = grua (2*m-1,k) + vz(ir,k)*dalpn(mr+n)
             grua (2*m  ,k) = grua (2*m  ,k) + vz(ii,k)*dalpn(mr+n)
!
             tmp = dalpn(mr+n)*hdifzd(n+m-1,k)
             gruha(2*m-1,k) = gruha(2*m-1,k) - vz(ir,k)*tmp
             gruha(2*m  ,k) = gruha(2*m  ,k) - vz(ii,k)*tmp
!
             grva (2*m-1,k) = grva (2*m-1,k) - d(ir,k)*dalpn(mr+n)
             grva (2*m  ,k) = grva (2*m  ,k) - d(ii,k)*dalpn(mr+n)
!
             grvha(2*m-1,k) = grvha(2*m-1,k) + d(ir,k)*tmp
             grvha(2*m  ,k) = grvha(2*m  ,k) + d(ii,k)*tmp
          end do
       end do
       do m=1,nmmax(irow)
          mr = nstart(m)
          mc = 2*mr
          do n=2,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
             grta (2*m-1,k) = grta (2*m-1,k) + t(ir,k)*alp(mr+n,irow)
             grta (2*m  ,k) = grta (2*m  ,k) + t(ii,k)*alp(mr+n,irow)
!
             grqa (2*m-1,k) = grqa (2*m-1,k) + q(ir,k)*alp(mr+n,irow)
             grqa (2*m  ,k) = grqa (2*m  ,k) + q(ii,k)*alp(mr+n,irow)
!
             tmp = alp(mr+n,irow)*hdiftq(n+m-1,k)
             grtha(2*m-1,k) = grtha(2*m-1,k) - t(ir,k)*tmp
             grtha(2*m  ,k) = grtha(2*m  ,k) - t(ii,k)*tmp
!
             grda(2*m-1,k) = grda(2*m-1,k) + d(ir,k)*alp(mr+n,irow)
             grda(2*m  ,k) = grda(2*m  ,k) + d(ii,k)*alp(mr+n,irow)
!
             gru1a (2*m-1) = gru1a (2*m-1) + d(ir,k)*alpn(mr+n)
             gru1a (2*m  ) = gru1a (2*m  ) + d(ii,k)*alpn(mr+n)
!
             tmp = alpn(mr+n)*hdifzd(n+m-1,k)
             gruh1a(2*m-1) = gruh1a(2*m-1) - d(ir,k)*tmp
             gruh1a(2*m  ) = gruh1a(2*m  ) - d(ii,k)*tmp
!
             grv1a (2*m-1) = grv1a (2*m-1) + vz(ir,k)*alpn(mr+n)
             grv1a (2*m  ) = grv1a (2*m  ) + vz(ii,k)*alpn(mr+n)
!
             grvh1a(2*m-1) = grvh1a(2*m-1) - vz(ir,k)*tmp
             grvh1a(2*m  ) = grvh1a(2*m  ) - vz(ii,k)*tmp
          end do
       end do
!
! Combine the two parts of u(m) and v(m)
!
       do m=1,nmmax(irow)
          grua (2*m-1,k) = grua (2*m-1,k) + gru1a (2*m  )
          gruha(2*m-1,k) = gruha(2*m-1,k) + gruh1a(2*m  )
          grua (2*m  ,k) = grua (2*m  ,k) - gru1a (2*m-1)
          gruha(2*m  ,k) = gruha(2*m  ,k) - gruh1a(2*m-1)
          grva (2*m-1,k) = grva (2*m-1,k) + grv1a (2*m  )
          grvha(2*m-1,k) = grvha(2*m-1,k) + grvh1a(2*m  )
          grva (2*m  ,k) = grva (2*m  ,k) - grv1a (2*m-1)
          grvha(2*m  ,k) = grvha(2*m  ,k) - grvh1a(2*m-1)
!
! Derivatives
!
          grtla(2*m-1,k) = -grta(2*m  ,k)*ra*xm(m)
          grtla(2*m  ,k) =  grta(2*m-1,k)*ra*xm(m)
          grqla(2*m-1,k) = -grqa(2*m  ,k)*ra*xm(m)
          grqla(2*m  ,k) =  grqa(2*m-1,k)*ra*xm(m)
       end do
    end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
    do m=1,nmmax(irow)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpma(2*m-1) = grpma(2*m-1) + alps(ir)*dalp(mr+n,irow)*ra
          grpma(2*m  ) = grpma(2*m  ) + alps(ii)*dalp(mr+n,irow)*ra
       end do
    end do
    do m=1,nmmax(irow)
       mr = nstart(m)
       mc = 2*mr
       do n=2,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          tmpr = alps(ir)*alp(mr+n,irow)
          tmpi = alps(ii)*alp(mr+n,irow)
          grpsa(2*m-1) = grpsa(2*m-1) + tmpr
          grpsa(2*m  ) = grpsa(2*m  ) + tmpi
!
          grdpsa(2*m-1) = grdpsa(2*m-1) + tmpr*hdfst4(m+n-1)*ztodt
          grdpsa(2*m  ) = grdpsa(2*m  ) + tmpi*hdfst4(m+n-1)*ztodt
       end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
       grpla(2*m-1) = -grpsa(2*m  )*ra*xm(m)
       grpla(2*m  ) =  grpsa(2*m-1)*ra*xm(m)
    end do
!
    return
end subroutine grcalca

#endif
