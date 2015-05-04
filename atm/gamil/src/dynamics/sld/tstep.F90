#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  This is due to the fact that spectral coefficients are
! stored consecutively along diagonals of M-N wavenumber space when the
! target architecture is PVP (optimal for vectorization), and along
! total wavenumber N otherwise (optimal for message-passing).
#if ( defined PVP )
subroutine tstep(n       ,ztdtsq  )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-implicit equations for each spectral element along
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
! $Id: tstep.F90,v 1.3.28.1 2002/06/15 13:48:34 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use pmgrid
  use pspect
  use comspe
  use comslt, only: epssld
  implicit none
#include <comhyb.h>
use commap

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: n               ! index of spectral diagonal being calculated
!                                           ! this call (not two dimensional wavenumber)
  real(r8), intent(in)   :: ztdtsq(2*pnmax) ! 2*dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
  real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
  real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
  real(r8) onepeps         ! decentering coefficient
  integer ne               ! index into ztdtsq
  integer m                ! diagonal element (index) of complex array
  integer k,kk             ! level indices
  integer irh              ! index into levels of spectral arrays
  integer irhr,irhi        ! index into real, imaginary coefficients
  integer isp              ! index into spectral arrays
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
! Set offsets for beginning of diagonal being calculated this call
!
  isp = nco2(n) - 2
  ne = 2*(n-1)
  onepeps = 1. + epssld
!
  do k=1,plev
!
! Coefficients for diagonal terms
!
     hhref = onepeps*0.5*href(k,k)
     hbps  = onepeps*0.5*bps (k)
!
! Loop along current diagonal (in spectral space)
! Add lnps and diagonal (vertical space) T terms to initialize hs
!
     do m=1,2*nm(n)
        hs(isp+m,k) = ztdtsq(ne+m)*(hhref*t(isp+m,k) + hbps*alps(isp+m))
     end do
     if (k.lt.plev) then
        do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to hs
!
           hhref = onepeps*0.5*href(kk,k)
           do m=1,2*nm(n)
              hs(isp+m,k) = hs(isp+m,k) + ztdtsq(ne+m)*hhref*t(isp+m,kk)
           end do
        end do
     end if
  end do                    ! k=1,plev (calculation level)
!
! Transform semi-implicit vectors to vertical normal mode space
!
  do k = 1,plev
     do m = 1,2*nm(n)
        dsnm(isp+m,k) = 0.
        hsnm(isp+m,k) = 0.
        vznm(isp+m,k) = 0.
        do kk = 1,plev
           dsnm(isp+m,k) = dsnm(isp+m,k) + bmi(kk,k)*d (isp+m,kk)
           hsnm(isp+m,k) = hsnm(isp+m,k) + bmi(kk,k)*hs(isp+m,kk)
           vznm(isp+m,k) = vznm(isp+m,k) + bmi(kk,k)*vz(isp+m,kk)
        end do
     end do
  end do
!
  return
#else
  subroutine tstep(m       ,ztdtsq  )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-implicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
!
! $Id: tstep.F90,v 1.3.28.1 2002/06/15 13:48:34 erik Exp $
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
! Input arguments
!
    integer , intent(in)   :: m             ! Fourier wavenumber               
    real(r8), intent(in)   :: ztdtsq(pnmax) ! 2*dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
    real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
    real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
    real(r8) onepeps         ! decentering coefficient
    integer n                ! 2-d wavenumber index
    integer k,kk             ! level indices
    integer mr,mc            ! real and imaginary spectral indices
    integer ir,ii            ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
!
    mr = nstart(m)
    mc = 2*mr
    onepeps = 1. + epssld
!
    do k=1,plev
!
! Coefficients for diagonal terms
!
       hhref = onepeps*0.5*href(k,k)
       hbps = onepeps*0.5*bps(k)
!
! Loop along total wavenumber index (in spectral space)
! Add lnps and diagonal (vertical space) T terms to initialize hs
!
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          hs(ir,k) = ztdtsq(n+m-1)*(hhref*t(ir,k) + hbps*alps(ir))
          hs(ii,k) = ztdtsq(n+m-1)*(hhref*t(ii,k) + hbps*alps(ii))
       end do
       if (k.lt.plev) then
          do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to hs
!
             hhref = onepeps*0.5*href(kk,k)
             do n=1,nlen(m)
                ir = mc + 2*n - 1
                ii = ir + 1
                hs(ir,k) = hs(ir,k) + ztdtsq(n+m-1)*hhref*t(ir,kk)
                hs(ii,k) = hs(ii,k) + ztdtsq(n+m-1)*hhref*t(ii,kk)
             end do
          end do
       end if
    end do                    ! k=1,plev (calculation level)
!
! Transform semi-implicit vectors to vertical normal mode space
!
    do k = 1,plev
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          dsnm(ir,k) = 0.
          dsnm(ii,k) = 0.
          hsnm(ir,k) = 0.
          hsnm(ii,k) = 0.
          vznm(ir,k) = 0.
          vznm(ii,k) = 0.
          do kk = 1,plev
             dsnm(ir,k) = dsnm(ir,k) + bmi(kk,k)*d (ir,kk)
             dsnm(ii,k) = dsnm(ii,k) + bmi(kk,k)*d (ii,kk)
             hsnm(ir,k) = hsnm(ir,k) + bmi(kk,k)*hs(ir,kk)
             hsnm(ii,k) = hsnm(ii,k) + bmi(kk,k)*hs(ii,kk)
             vznm(ir,k) = vznm(ir,k) + bmi(kk,k)*vz(ir,kk)
             vznm(ii,k) = vznm(ii,k) + bmi(kk,k)*vz(ii,kk)
          end do
       end do
    end do
!
    return
#endif
  end subroutine tstep

