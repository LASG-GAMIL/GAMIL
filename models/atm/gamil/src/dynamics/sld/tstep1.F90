#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  This is due to the fact that spectral coefficients are
! stored consecutively along diagonals of M-N wavenumber space when the
! target architecture is PVP (optimal for vectorization), and along
! total wavenumber N otherwise (optimal for message-passing).
#if ( defined PVP )
subroutine tstep1(n       ,zdt     )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! the n(th) diagonal. (Note, n is distinct from the two dimensional
! wavenumber which is also often denoted n.) The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: tstep1.F90,v 1.3.28.1 2002/06/15 13:48:34 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use comslt, only: epssld
  implicit none
#include <comhyb.h>
use commap

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: n   ! index of spectral diagonal being calculated
!                               ! this call (not two dimensional wavenumber)
  real(r8), intent(in)   :: zdt ! timestep, dt (seconds)
!
!---------------------------Local workspace-----------------------------
!
  real(r8) z (2*pmmax,plev) ! workspace for computation of spectral array d
  real(r8) zz(2*pmmax,plev) ! workspace for computation of spectral array vz
  real(r8) ztemp            ! temporary workspace
  real(r8) onepeps          ! decentering coefficient
  integer m                 ! diagonal element (index) of complex array
  integer k,kk              ! level indices
  integer irh               ! index into levels of spectral arrays
  integer irhr,irhi         ! index into real, imaginary coefficients
  integer isp               ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Set offsets for beginning of diagonal being calculated this call
!
  isp = nco2(n) - 2
  onepeps = 1. + epssld
!
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
  do k=1,plev
     do m=1,2*pmmax
        z (m,k) = 0.
        zz(m,k) = 0.
     end do
  end do
!
! Transform back from normal mode space
!
  do k=1,plev
     irhr = nco2(n) - 3
     irhi = irhr + 1
     do kk=1,plev
        do m=1,nm(n)
           z (2*m-1,k) = z (2*m-1,k) + bm1(kk,k)*dnm (irhr+2*m,kk)
           z (2*m  ,k) = z (2*m  ,k) + bm1(kk,k)*dnm (irhi+2*m,kk)
           zz(2*m-1,k) = zz(2*m-1,k) + bm1(kk,k)*vznm(irhr+2*m,kk)
           zz(2*m  ,k) = zz(2*m  ,k) + bm1(kk,k)*vznm(irhi+2*m,kk)
        end do
     end do                  ! inner loop over levels
  end do                    ! outer loop over levels
!
! Move solution for divergence and vorticity to d and vz.
!
  irh = nco2(n) - 2
  do k=1,plev
     do m=1,2*nm(n)
        d (irh+m,k) = z (m,k)
        vz(irh+m,k) = zz(m,k)
     end do
  end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
  do k=1,plev
     ztemp = onepeps*zdt*hypd(k)/hypi(plevp)
     do m=1,2*nm(n)
        alps(isp+m) = alps(isp+m) - ztemp*d(isp+m,k)
     end do
  end do
!
! Add ln(Ps)star back in to get full ln(Ps)
!
  do m=1,2*nm(n)
     alps(isp+m) = alps(isp+m) + lnpstar(isp+m)
  end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
  do k=1,plev
     do kk=1,plev
        ztemp = onepeps*zdt*tau(kk,k)
        do m=1,2*nm(n)
           t(isp+m,k) = t(isp+m,k) - ztemp*d(isp+m,kk)
        end do
     end do
  end do
!
  return
#else
  subroutine tstep1(m       ,zdt     )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: tstep1.F90,v 1.3.28.1 2002/06/15 13:48:34 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use pspect
    use comspe
    use comslt, only: epssld
    use commap

    implicit none

#include <comhyb.h>

!------------------------------Arguments--------------------------------
!
    integer , intent(in)   :: m   ! Fourier wavenumber               
    real(r8), intent(in)   :: zdt ! timestep, dt (seconds)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) z (2*pnmax,plev) ! workspace for computation of spectral array d
    real(r8) zz(2*pnmax,plev) ! workspace for computation of spectral array vz
    real(r8) ztemp            ! temporary workspace
    real(r8) onepeps          ! decentering coefficient
    integer n,j               ! 2-d wavenumber index
    integer k,kk              ! level indices
    integer mr,mc             ! real and imaginary spectral indices
    integer ir,ii             ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
!
    mr = nstart(m)
    mc = 2*mr
    onepeps = 1. + epssld
!
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
    do k=1,plev
       do j=1,2*pnmax
          z (j,k) = 0.
          zz(j,k) = 0.
       end do
    end do
!
! Transform back from normal mode space
!
    do k=1,plev
       do kk=1,plev
          do n=1,nlen(m)
             ir = mc + 2*n - 1
             ii = ir + 1
             z (2*n-1,k) = z (2*n-1,k) + bm1(kk,k)*dnm (ir,kk)
             z (2*n  ,k) = z (2*n  ,k) + bm1(kk,k)*dnm (ii,kk)
             zz(2*n-1,k) = zz(2*n-1,k) + bm1(kk,k)*vznm(ir,kk)
             zz(2*n  ,k) = zz(2*n  ,k) + bm1(kk,k)*vznm(ii,kk)
          end do
       end do                  ! inner loop over levels
    end do                    ! outer loop over levels
!
! Move solution for divergence and vorticity to d and vz.
!
    do k=1,plev
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          d (ir,k) = z (2*n-1,k)
          d (ii,k) = z (2*n  ,k)
          vz(ir,k) = zz(2*n-1,k)
          vz(ii,k) = zz(2*n  ,k)
       end do
    end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
    do k=1,plev
       ztemp = onepeps*zdt*hypd(k)/hypi(plevp)
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          alps(ir) = alps(ir) - ztemp*d(ir,k)
          alps(ii) = alps(ii) - ztemp*d(ii,k)
       end do
    end do
!
! Add ln(Ps)star back in to get full ln(Ps)
!
    do n=1,nlen(m)
       ir = mc + 2*n - 1
       ii = ir + 1
       alps(ir) = alps(ir) + lnpstar(ir)
       alps(ii) = alps(ii) + lnpstar(ii)
    end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
    do k=1,plev
       do kk=1,plev
          ztemp = onepeps*zdt*tau(kk,k)
          do n=1,nlen(m)
             ir = mc + 2*n - 1
             ii = ir + 1
             t(ir,k) = t(ir,k) - ztemp*d(ir,kk)
             t(ii,k) = t(ii,k) - ztemp*d(ii,kk)
          end do
       end do
    end do
!
    return
#endif
  end subroutine tstep1
