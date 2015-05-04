#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  Make sure to make appropriate coding changes where
! necessary.
#if ( defined PVP )
subroutine vertnm(k)
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the system of semi-implicit divergence/vorticity
! equations.  Equations have been de-coupled in the vertical by
! transforming the spectral terms into vertical normal modes.
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: vertnm.F90,v 1.3.28.1 2002/06/15 13:48:35 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  implicit none
#include <comhyb.h>
use commap

!------------------------------Arguments--------------------------------
!
  integer, intent(in) :: k  ! vertical index
!
!---------------------------Local workspace-----------------------------
!
  integer m                 ! diagonal element of cmplx array
  integer ir                ! index pointing into real(r8)      element
  integer ii                ! index pointing into imaginary element
  integer n                 ! wave number
  integer nn                ! index
  integer nr                ! n   (real)
  integer ni                ! n   (imag)
  integer np1r              ! n+1 (real)
  integer np1i              ! n+1 (imag)
  integer nm1r              ! n-1 (real)
  integer nm1i              ! n-1 (imag)
  integer mr                ! real spectral index
  integer mc                ! imaginary spectral index
  integer isp               ! index into spectral arrays
  real(r8) tmp  (2)         ! real/imaginary temp spaces
  real(r8) tmp1 (2)         ! real/imaginary temp spaces
  real(r8) tmp2 (2)         ! real/imaginary temp spaces
  real(r8) btrie(2*pmax)    ! "btri" + eigenvalue of ref atmosphere
  real(r8) dtri (2*pmax)    ! RHS in the solution of the tri-diagonal matrix
  real(r8) dnmn (psp)       ! remapped "dnm"
  real(r8) dsnmn(psp)       ! remapped "dsnm"
  real(r8) hsnmn(psp)       ! remapped "hsnm"
  real(r8) vznmn(psp)       ! remapped "vznm"
  real(r8) denom            ! denominator
!
!-----------------------------------------------------------------------
!
! Remap spectral arrays from sequential "along diagonals" to sequential
! along columns ("N")
!
  do n = 1,pmax
     isp = nco2(n) - 2
     nn  = 2*n-1
     do m = 1,nm(n)
        mr = nstart(m)
        mc = 2*mr
        ir = mc + nn
        ii = ir + 1
        dsnmn(ir) = dsnm(isp+(2*m-1),k)
        dsnmn(ii) = dsnm(isp+(2*m  ),k)
        hsnmn(ir) = hsnm(isp+(2*m-1),k)
        hsnmn(ii) = hsnm(isp+(2*m  ),k)
        vznmn(ir) = vznm(isp+(2*m-1),k)
        vznmn(ii) = vznm(isp+(2*m  ),k)
     end do
  end do
!
! Complete computation of btri and compute dtri (complex arithmetic)
!
  do m = 1,pmmax
     mr = nstart(m)
     mc = 2*mr
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        btrie(2*n-1) = btri(nr) + zcr(m+n-1,k)
        btrie(2*n  ) = btri(ni)
!
        tmp1(1) = 0.
        tmp1(2) = 0. 
        tmp2(1) = 0.
        tmp2(2) = 0. 
        if(n .ne. nlen(m)) then
           tmp(1)  =  bpnm(nr  )*vznmn(np1r)
           tmp(2)  =  bpnm(nr  )*vznmn(np1i)
           denom   =  a0nm(np1r)*a0nm (np1r) + a0nm(np1i)*a0nm (np1i)
           tmp1(1) = (a0nm(np1r)*tmp(1) + a0nm(np1i)*tmp(2))/denom
           tmp1(2) = (a0nm(np1r)*tmp(2) - a0nm(np1i)*tmp(1))/denom
        endif
        if(n .ne. 1    ) then
           tmp(1)  =  bmnm(nr  )*vznmn(nm1r)
           tmp(2)  =  bmnm(nr  )*vznmn(nm1i)
           denom   =  a0nm(nm1r)*a0nm (nm1r) + a0nm(nm1i)*a0nm (nm1i)
           tmp2(1) = (a0nm(nm1r)*tmp(1) + a0nm(nm1i)*tmp(2))/denom
           tmp2(2) = (a0nm(nm1r)*tmp(2) - a0nm(nm1i)*tmp(1))/denom
        endif
        dtri(2*n-1) = dsnmn(nr) + hsnmn(nr) + tmp1(1) + tmp2(1)
        dtri(2*n  ) = dsnmn(ni) + hsnmn(ni) + tmp1(2) + tmp2(2)
     end do
!
! Solve tridiagonal matrix:  call once for ODDs and once for EVENs
!
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m)
     call trisolve(n       ,atri(mc+1),btrie(1),ctri(mc+1),dtri(1),dnmn(mc+1)    )
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m) - 2
     if(n .gt. 0) then
        call trisolve(n       ,atri(mc+3),btrie(3),ctri(mc+3),dtri(3),dnmn(mc+3)    )
     endif
!
! Solve for vorticity
!
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        tmp1(1) = 0.
        tmp1(2) = 0. 
        tmp2(1) = 0.
        tmp2(2) = 0. 
        if(n .ne. nlen(m)) then
           tmp1(1) = bpnm(nr)*dnmn(np1r)
           tmp1(2) = bpnm(nr)*dnmn(np1i)
        endif
        if(n .ne. 1      ) then
           tmp2(1) = bmnm(nr)*dnmn(nm1r)
           tmp2(2) = bmnm(nr)*dnmn(nm1i)
        endif
        vznmn(nr) = vznmn(nr) - tmp1(1) - tmp2(1)
        vznmn(ni) = vznmn(ni) - tmp1(2) - tmp2(2)
        denom  =  a0nm(nr)*a0nm (nr) + a0nm(ni)*a0nm (ni)
        tmp(1) = (a0nm(nr)*vznmn(nr) + a0nm(ni)*vznmn(ni))/denom
        tmp(2) = (a0nm(nr)*vznmn(ni) - a0nm(ni)*vznmn(nr))/denom
        vznmn(nr) = tmp(1)
        vznmn(ni) = tmp(2)
     end do
  end do
!
! Map spectral arrays back to "along diagonals"
!
  do n = 1,pmax
     isp = nco2(n) - 2
     nn  = 2*n-1
     do m = 1,nm(n)
        mr = nstart(m)
        mc = 2*mr
        ir = mc + nn
        ii = ir + 1
        dnm (isp+(2*m-1),k) = dnmn (ir)
        dnm (isp+(2*m  ),k) = dnmn (ii)
        vznm(isp+(2*m-1),k) = vznmn(ir)
        vznm(isp+(2*m  ),k) = vznmn(ii)
     end do
  end do
!
  return
end subroutine vertnm
#else
subroutine vertnm(m)
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the system of semi-implicit divergence/vorticity
! equations.  Equations have been de-coupled in the vertical by
! transforming the spectral terms into vertical normal modes.
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: vertnm.F90,v 1.3.28.1 2002/06/15 13:48:35 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use commap
  
  implicit none

#include <comhyb.h>

!------------------------------Arguments--------------------------------
!
  integer, intent(in) :: m  ! diagonal element of cmplx array
!
!---------------------------Local workspace-----------------------------
!
  integer k                 ! vertical index
  integer n                 ! wave number
  integer nr                ! n   (real)
  integer ni                ! n   (imag)
  integer np1r              ! n+1 (real)
  integer np1i              ! n+1 (imag)
  integer nm1r              ! n-1 (real)
  integer nm1i              ! n-1 (imag)
  integer mr                ! real spectral index
  integer mc                ! imaginary spectral index
  real(r8) tmp  (2)         ! real/imaginary temp spaces
  real(r8) tmp1 (2)         ! real/imaginary temp spaces
  real(r8) tmp2 (2)         ! real/imaginary temp spaces
  real(r8) btrie(2*pmax)    ! "btri" + eigenvalue of ref atmosphere
  real(r8) dtri (2*pmax)    ! RHS in the solution of the tri-diagonal matrix
  real(r8) denom            ! denominator
!
!-----------------------------------------------------------------------
!
! Complete computation of btri and compute dtri (complex arithmetic)
!
  mr = nstart(m)
  mc = 2*mr
  do k = 1,plev
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        btrie(2*n-1) = btri(nr) + zcr(m+n-1,k)
        btrie(2*n  ) = btri(ni)
!
        tmp1(1) = 0.
        tmp1(2) = 0. 
        tmp2(1) = 0.
        tmp2(2) = 0. 
        if(n .ne. nlen(m)) then
           tmp(1)  =  bpnm(nr  )*vznm(np1r,k)
           tmp(2)  =  bpnm(nr  )*vznm(np1i,k)
           denom   =  a0nm(np1r)*a0nm(np1r) + a0nm(np1i)*a0nm(np1i)
           tmp1(1) = (a0nm(np1r)*tmp(1) + a0nm(np1i)*tmp(2))/denom
           tmp1(2) = (a0nm(np1r)*tmp(2) - a0nm(np1i)*tmp(1))/denom
        endif
        if(n .ne. 1    ) then
           tmp(1)  =  bmnm(nr  )*vznm(nm1r,k)
           tmp(2)  =  bmnm(nr  )*vznm(nm1i,k)
           denom   =  a0nm(nm1r)*a0nm(nm1r) + a0nm(nm1i)*a0nm(nm1i)
           tmp2(1) = (a0nm(nm1r)*tmp(1) + a0nm(nm1i)*tmp(2))/denom
           tmp2(2) = (a0nm(nm1r)*tmp(2) - a0nm(nm1i)*tmp(1))/denom
        endif
        dtri(2*n-1) = dsnm(nr,k) + hsnm(nr,k) + tmp1(1) + tmp2(1)
        dtri(2*n  ) = dsnm(ni,k) + hsnm(ni,k) + tmp1(2) + tmp2(2)
     end do
!
! Solve tridiagonal matrix:  call once for ODDs and once for EVENs
!
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m)
     call trisolve(n       ,atri(mc+1),btrie(1),ctri(mc+1),dtri(1),dnm(mc+1,k)    )
     if(mod(nlen(m),2) .eq. 0) n = nlen(m) - 1
     if(mod(nlen(m),2) .ne. 0) n = nlen(m) - 2
     if(n .gt. 0) then
        call trisolve(n       ,atri(mc+3),btrie(3),ctri(mc+3),dtri(3),dnm(mc+3,k)    )
     endif
!
! Solve for vorticity
!
     do n = 1,nlen(m)
        nr   = mc   + 2*n     - 1
        ni   = nr   + 1
        np1r = mc   + 2*(n+1) - 1
        np1i = np1r + 1
        nm1r = mc   + 2*(n-1) - 1
        nm1i = nm1r + 1
!
        tmp1(1) = 0.
        tmp1(2) = 0. 
        tmp2(1) = 0.
        tmp2(2) = 0. 
        if(n .ne. nlen(m)) then
           tmp1(1) = bpnm(nr)*dnm(np1r,k)
           tmp1(2) = bpnm(nr)*dnm(np1i,k)
        endif
        if(n .ne. 1      ) then
           tmp2(1) = bmnm(nr)*dnm(nm1r,k)
           tmp2(2) = bmnm(nr)*dnm(nm1i,k)
        endif
        vznm(nr,k) = vznm(nr,k) - tmp1(1) - tmp2(1)
        vznm(ni,k) = vznm(ni,k) - tmp1(2) - tmp2(2)
        denom  =  a0nm(nr)*a0nm(nr)   + a0nm(ni)*a0nm(ni)
        tmp(1) = (a0nm(nr)*vznm(nr,k) + a0nm(ni)*vznm(ni,k))/denom
        tmp(2) = (a0nm(nr)*vznm(ni,k) - a0nm(ni)*vznm(nr,k))/denom
        vznm(nr,k) = tmp(1)
        vznm(ni,k) = tmp(2)
     end do
  end do
!
  return
end subroutine vertnm
#endif
