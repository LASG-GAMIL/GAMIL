#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  Make sure to make appropriate coding changes where
! necessary.
#if ( defined PVP )
subroutine quad(n       ,grlps1  ,grlps2  ,grt1    ,grq1    , &
                grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
                grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Perform gaussian quadrature to obtain the spectral coefficients of 
! ln(ps), temperature, vorticity, divergence, and specific humidity
! (to calculate derivatives only). Add the tendency terms requiring
! meridional derivatives during the transform.
!
! Computational note: This routine is multitasked over each diagonal in
! the spectral (i.e., [m,n]) data structure since each of these vector
! elements are computationally independent of each other. Care should
! be taken in interpreting the code since the diagonals are indexed
! using the variable "n" (which is often used to denote 2-dimensional
! wavenumber).
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: quad.F90,v 1.4.8.1 2002/06/15 13:48:27 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use rgrid
  use dynconst, only: rearth
  use commap

  implicit none

!------------------------------Arguments--------------------------------
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMS and and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
! Suffixes 1 and 2 refer to symmetric and antisymmetric components
! respectively.
!
  integer , intent(in)   :: n                           ! spectral space "diagonal" index
  real(r8), intent(in)   :: grlps1(2*pmmax,plat/2)      ! ln(ps) - symmetric
  real(r8), intent(in)   :: grlps2(2*pmmax,plat/2)      ! ln(ps) - antisymmetric
!
! symmetric components
!
  real(r8), intent(in)   :: grt1  (2*pmmax,plev,plat/2) ! temperature
  real(r8), intent(in)   :: grq1  (2*pmmax,plev,plat/2) ! moisture
  real(r8), intent(in)   :: grz1  (2*pmmax,plev,plat/2) ! vorticity
  real(r8), intent(in)   :: grd1  (2*pmmax,plev,plat/2) ! divergence
  real(r8), intent(in)   :: grfu1 (2*pmmax,plev,plat/2) ! partial u momentum tendency (fu)
  real(r8), intent(in)   :: grfv1 (2*pmmax,plev,plat/2) ! partial v momentum tendency (fv)
!
! antisymmetric components
!
  real(r8), intent(in)   :: grt2  (2*pmmax,plev,plat/2) ! temperature
  real(r8), intent(in)   :: grq2  (2*pmmax,plev,plat/2) ! moisture
  real(r8), intent(in)   :: grz2  (2*pmmax,plev,plat/2) ! vorticity
  real(r8), intent(in)   :: grd2  (2*pmmax,plev,plat/2) ! divergence
  real(r8), intent(in)   :: grfu2 (2*pmmax,plev,plat/2) ! partial u momentum tend (fu)
  real(r8), intent(in)   :: grfv2 (2*pmmax,plev,plat/2) ! partial v momentum tend (fv)
!
!---------------------------Local workspace-----------------------------
!
  real(r8) alp2  (2*pmax,plat/2)  ! expand polynomials to complex
  real(r8) dalp2 (2*pmax,plat/2)  ! expand polynomials to complex
  integer j                       ! latitude pair index
  integer m                       ! diagonal element(index)of sp.arr.
  integer isp                     ! index into levls of sp. arrs.
  integer k                       ! level index
  integer ialp                    ! index into polynomials

  real(r8) zcsj                   ! cos**2(lat)*radius of earth
  real(r8) zrcsj                  ! 1./(a*cos^2(lat))
  real(r8) zw    (plat/2)         ! 2*w
  real(r8) ztdtrw(plat/2)         ! 2w*2dt/(a*cos^2(lat))
!
!-----------------------------------------------------------------------
!
! Compute constants
!
  do j = 1,plat/2
     zcsj      = cs(j)*rearth
     zrcsj     = 1./zcsj
     zw(j)     = w(j)*2.
     ztdtrw(j) = zrcsj*zw(j)
  end do
  ialp = nalp(n)
  isp  = nco2(n) - 2

  if (nm(n).gt.plat/2) then ! vectorize over diagonals
     do j = 1,plat/2
!
! Expand polynomials to complex form to allow largest possible vector
! length
!
!DIR$ IVDEP
        do m = 1,nmreduced(n,j)
           alp2 (2*m-1,j) = alp(ialp+m,j)*zw(j)
           alp2 (2*m  ,j) = alp(ialp+m,j)*zw(j)
           dalp2(2*m-1,j) = dalp(ialp+m,j)*ztdtrw(j)
           dalp2(2*m  ,j) = dalp(ialp+m,j)*ztdtrw(j)
        end do
     end do
  else                      ! vectorize over latitude
     do m = 1,nm(n)
!
! Expand polynomials to complex form to allow largest possible vector
! length
!
!DIR$ IVDEP
        do j = 1,plat/2
           alp2 (2*m-1,j) = alp(ialp+m,j)*zw(j)
           alp2 (2*m  ,j) = alp(ialp+m,j)*zw(j)
           dalp2(2*m-1,j) = dalp(ialp+m,j)*ztdtrw(j)
           dalp2(2*m  ,j) = dalp(ialp+m,j)*ztdtrw(j)
        end do
     end do
  end if
!
! Accumulate contributions to spectral coefficients of ln(p*), the only
! single level field. Use symmetric or antisymmetric fourier cofficients
! depending on whether the total wavenumber is even or odd.  Vectorize
! over diagonals or latitude pair index depending on vector length.
! Comparison is 2*nm(n) vs. plat/8 instead of plat/2 since latitude
! sum is into a scalar variable.
!
  do m = 1,2*nm(n)
     alps(isp+m) = 0.
  end do
  if (2*nm(n).gt.plat/8) then ! vectorize over diagonals
     if (mod(n,2).ne.0) then
        do j = 1,plat/2
           do m = 1,2*nmreduced(n,j)
              alps(isp+m) = alps(isp+m) + grlps1(m,j)*alp2(m,j)
           end do
        end do
     else
        do j = 1,plat/2
           do m = 1,2*nmreduced(n,j)
              alps(isp+m) = alps(isp+m) + grlps2(m,j)*alp2(m,j)
           end do
        end do
     end if

  else                      ! vectorize over latitude

     if (mod(n,2).ne.0) then
        do m = 1,2*nm(n)
           do j = beglatpair((m+1)/2),plat/2
              alps(isp+m) = alps(isp+m) + grlps1(m,j)*alp2(m,j)
           end do
        end do
     else
        do m = 1,2*nm(n)
           do j = beglatpair((m+1)/2),plat/2
              alps(isp+m) = alps(isp+m) + grlps2(m,j)*alp2(m,j)
           end do
        end do
     end if
  end if
!
! Accumulate contributions to spectral coefficients of the multilevel
! fields (temperature, divergence, vorticity).  Use symmetric or
! antisymmetric fourier coefficients depending on whether the total
! wavenumber is even or odd.
!
  if (2*nm(n).gt.plev) then ! vectorize over diagonals
     do k = 1,plev
        do m = 1,2*nm(n)
           t (isp+m,k) = 0.
           q (isp+m,k) = 0.
           d (isp+m,k) = 0.
           vz(isp+m,k) = 0.
        end do
        if (mod(n,2).ne.0) then ! n is odd
           do j = 1,plat/2
              do m = 1,2*nmreduced(n,j)
                 t (isp+m,k) = t (isp+m,k) + grt1 (m,k,j)*alp2 (m,j)
                 q (isp+m,k) = q (isp+m,k) + grq1 (m,k,j)*alp2 (m,j)
                 d (isp+m,k) = d (isp+m,k) + grd1 (m,k,j)*alp2 (m,j) - &
                      grfv2(m,k,j)*dalp2(m,j)
                 vz(isp+m,k) = vz(isp+m,k) + grz1 (m,k,j)*alp2 (m,j) + &
                      grfu2(m,k,j)*dalp2(m,j)
              end do
           end do
        else                    ! n is even
           do j = 1,plat/2
              do m = 1,2*nmreduced(n,j)
                 t (isp+m,k) = t (isp+m,k) + grt2 (m,k,j)*alp2 (m,j)
                 q (isp+m,k) = q (isp+m,k) + grq2 (m,k,j)*alp2 (m,j)
                 d (isp+m,k) = d (isp+m,k) + grd2 (m,k,j)*alp2 (m,j) - &
                      grfv1(m,k,j)*dalp2(m,j)
                 vz(isp+m,k) = vz(isp+m,k) + grz2 (m,k,j)*alp2 (m,j) + &
                      grfu1(m,k,j)*dalp2(m,j)
              end do
           end do
        end if
     end do

  else                       ! vectorize over levels

     do m = 1,2*nm(n)
        do k = 1,plev
           t (isp+m,k) = 0.
           q (isp+m,k) = 0.
           d (isp+m,k) = 0.
           vz(isp+m,k) = 0.
        end do
        if (mod(n,2).ne.0) then ! n is odd
           do j = beglatpair((m+1)/2),plat/2
              do k = 1,plev
                 t (isp+m,k) = t (isp+m,k) + grt1 (m,k,j)*alp2 (m,j)
                 q (isp+m,k) = q (isp+m,k) + grq1 (m,k,j)*alp2 (m,j)
                 d (isp+m,k) = d (isp+m,k) + grd1 (m,k,j)*alp2 (m,j) - &
                      grfv2(m,k,j)*dalp2(m,j)
                 vz(isp+m,k) = vz(isp+m,k) + grz1 (m,k,j)*alp2 (m,j) + &
                      grfu2(m,k,j)*dalp2(m,j)
              end do
           end do
        else                    ! n is even
           do j=beglatpair((m+1)/2),plat/2
              do k = 1,plev
                 t (isp+m,k) = t (isp+m,k) + grt2 (m,k,j)*alp2 (m,j)
                 q (isp+m,k) = q (isp+m,k) + grq2 (m,k,j)*alp2 (m,j)
                 d (isp+m,k) = d (isp+m,k) + grd2 (m,k,j)*alp2 (m,j) - &
                      grfv1(m,k,j)*dalp2(m,j)
                 vz(isp+m,k) = vz(isp+m,k) + grz2 (m,k,j)*alp2 (m,j) + &
                      grfu1(m,k,j)*dalp2(m,j)
              end do
           end do
        end if
     end do
  end if
!
  return
end subroutine quad
#else
subroutine quad(m       ,grlps1  ,grlps2  ,grt1    ,grq1    , &
                grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
                grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Perform gaussian quadrature for 1 Fourier wavenumber (m) to obtain the
! spectral coefficients of ln(ps), temperature, vorticity, divergence
! and specific humidity.  Add the tendency terms requiring meridional
! derivatives during the transform.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use rgrid
  use commap
  use dynconst, only: rearth

  implicit none

!------------------------------Arguments--------------------------------
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMS and and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
! Suffixes 1 and 2 refer to symmetric and antisymmetric components
! respectively.
!
  integer , intent(in)   :: m                           ! Fourier wavenumber
  real(r8), intent(in)   :: grlps1(2*pmmax,plat/2)      ! ln(ps) - symmetric
  real(r8), intent(in)   :: grlps2(2*pmmax,plat/2)      ! ln(ps) - antisymmetric
!
! symmetric components
!
  real(r8), intent(in)   :: grt1  (plev,2*pmmax,plat/2) ! temperature
  real(r8), intent(in)   :: grq1  (plev,2*pmmax,plat/2) ! moisture
  real(r8), intent(in)   :: grz1  (plev,2*pmmax,plat/2) ! vorticity
  real(r8), intent(in)   :: grd1  (plev,2*pmmax,plat/2) ! divergence
  real(r8), intent(in)   :: grfu1 (plev,2*pmmax,plat/2) ! partial u momentum tendency (fu)
  real(r8), intent(in)   :: grfv1 (plev,2*pmmax,plat/2) ! partial v momentum tendency (fv)
!
! antisymmetric components
!
  real(r8), intent(in)   :: grt2  (plev,2*pmmax,plat/2) ! temperature
  real(r8), intent(in)   :: grq2  (plev,2*pmmax,plat/2) ! moisture
  real(r8), intent(in)   :: grz2  (plev,2*pmmax,plat/2) ! vorticity
  real(r8), intent(in)   :: grd2  (plev,2*pmmax,plat/2) ! divergence
  real(r8), intent(in)   :: grfu2 (plev,2*pmmax,plat/2) ! partial u momentum tend (fu)
  real(r8), intent(in)   :: grfv2 (plev,2*pmmax,plat/2) ! partial v momentum tend (fv)
!
!---------------------------Local workspace-----------------------------
!
  integer j                 ! latitude pair index
  integer n                 ! index
  integer ir,ii             ! indices
  integer mr,mc             ! indices
  integer k                 ! level index
  real(r8) zcsj             ! cos**2(lat)*radius of earth
  real(r8) zrcsj            ! 1./(a*cos^2(lat))
  real(r8) zw    (plat/2)   ! 2*w
  real(r8) ztdtrw(plat/2)   ! 2w*2dt/(a*cos^2(lat))
  real(r8) zwalp
  real(r8) zwdalp
!
!-----------------------------------------------------------------------
!
! Compute constants
!
  do j = 1,plat/2
     zcsj      = cs(j)*rearth
     zrcsj     = 1./zcsj
     zw(j)     = w(j)*2.
     ztdtrw(j) = zrcsj*zw(j)
  end do
!
! Accumulate contributions to spectral coefficients of ln(p*), the only
! single level field. Use symmetric or antisymmetric fourier cofficients
! depending on whether the total wavenumber is even or odd.
!
  mr = nstart(m)
  mc = 2*mr
  do n = 1,2*nlen(m)
     alps(mc+n) = 0.
  end do
  do j=beglatpair(m),plat/2
     do n=1,nlen(m),2
        ir = mc + 2*n - 1
        ii = ir + 1
        zwalp    = zw(j)*alp(mr+n,j)
        alps(ir) = alps(ir) + grlps1(2*m-1,j)*zwalp
        alps(ii) = alps(ii) + grlps1(2*m  ,j)*zwalp
     end do
     do n = 2,nlen(m),2
        ir = mc + 2*n - 1
        ii = ir + 1
        zwalp    = zw(j)*alp(mr+n,j)
        alps(ir) = alps(ir) + grlps2(2*m-1,j)*zwalp
        alps(ii) = alps(ii) + grlps2(2*m  ,j)*zwalp
     end do
  end do
!
! Accumulate contributions to spectral coefficients of the multilevel
! fields.  Use symmetric or antisymmetric fourier coefficients depending
! on whether the total wavenumber is even or odd.
!
  do k = 1,plev
     do n = 1,2*nlen(m)
        t (mc+n,k) = 0.
        q (mc+n,k) = 0.
        d (mc+n,k) = 0.
        vz(mc+n,k) = 0.
     end do
     do j=beglatpair(m),plat/2
        do n=1,nlen(m),2
           zwdalp = ztdtrw(j)*dalp(mr+n,j)
           zwalp  = zw(j)    *alp (mr+n,j)
           ir = mc + 2*n - 1
           ii = ir + 1
           t (ir,k) = t (ir,k) + zwalp*grt1 (k,2*m-1,j)
           t (ii,k) = t (ii,k) + zwalp*grt1 (k,2*m  ,j)
           q (ir,k) = q (ir,k) + zwalp*grq1 (k,2*m-1,j)
           q (ii,k) = q (ii,k) + zwalp*grq1 (k,2*m  ,j)
           d (ir,k) = d (ir,k) + grd1 (k,2*m-1,j)*zwalp - grfv2(k,2*m-1,j)*zwdalp
           d (ii,k) = d (ii,k) + grd1 (k,2*m  ,j)*zwalp - grfv2(k,2*m  ,j)*zwdalp
           vz(ir,k) = vz(ir,k) + grz1 (k,2*m-1,j)*zwalp + grfu2(k,2*m-1,j)*zwdalp
           vz(ii,k) = vz(ii,k) + grz1 (k,2*m  ,j)*zwalp + grfu2(k,2*m  ,j)*zwdalp
        end do
     end do
     do j=beglatpair(m),plat/2
        do n = 2,nlen(m),2
           zwdalp = ztdtrw(j)*dalp(mr+n,j)
           zwalp  = zw(j)    *alp (mr+n,j)
           ir = mc + 2*n - 1
           ii = ir + 1
           t (ir,k) = t (ir,k) + zwalp*grt2(k,2*m-1,j)
           t (ii,k) = t (ii,k) + zwalp*grt2(k,2*m  ,j)
           q (ir,k) = q (ir,k) + zwalp*grq2(k,2*m-1,j)
           q (ii,k) = q (ii,k) + zwalp*grq2(k,2*m  ,j)
           d (ir,k) = d (ir,k) + grd2 (k,2*m-1,j)*zwalp - grfv1(k,2*m-1,j)*zwdalp
           d (ii,k) = d (ii,k) + grd2 (k,2*m  ,j)*zwalp - grfv1(k,2*m  ,j)*zwdalp
           vz(ir,k) = vz(ir,k) + grz2 (k,2*m-1,j)*zwalp + grfu1(k,2*m-1,j)*zwdalp
           vz(ii,k) = vz(ii,k) + grz2 (k,2*m  ,j)*zwalp + grfu1(k,2*m  ,j)*zwdalp
        end do
     end do
  end do
!
  return
end subroutine quad
#endif

