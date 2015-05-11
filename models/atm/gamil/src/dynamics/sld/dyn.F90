#include <misc.h>
#include <params.h>
! Note that this routine has 2 complete blocks of code for PVP vs.
! non-PVP.  Make sure to make appropriate coding changes where
! necessary.
#if ( defined PVP )
subroutine dyn(irow    ,grlps1  ,grt1    ,grq1    ,grz1    , &
               grd1    ,grfu1   ,grfv1   ,grlps2  ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Combine undifferentiated and longitudinally differentiated Fourier
! coefficient terms for later use in the Gaussian quadrature
!
! Computational note: Index "2*m-1" refers to the real part of the
! complex coefficient, and "2*m" to the imaginary.
!
! The naming convention is as follows:
!  - t, q, d, z refer to temperature, specific humidity, divergence
!     and vorticity
!  - "1" suffix to an array => symmetric component of current latitude
!     pair
!  - "2" suffix to an array => antisymmetric component
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: dyn.F90,v 1.5.2.1 2002/06/15 13:48:21 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use rgrid
  use commap
  use dynconst, only: rearth

  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: irow                 ! latitude pair index
  real(r8), intent(inout):: grlps1(2*pmmax)      ! sym. undifferentiated term in Ps eqn.
  real(r8), intent(inout):: grt1  (2*pmmax,plev) ! sym. undifferentiated term in t eqn.
  real(r8), intent(inout):: grq1  (2*pmmax,plev) ! sym. undifferentiated term in q
  real(r8), intent(out)  :: grz1  (2*pmmax,plev) ! sym. undifferentiated term in z eqn.
  real(r8), intent(out)  :: grd1  (2*pmmax,plev) ! sym. undifferentiated term in d eqn.
  real(r8), intent(inout):: grfu1 (2*pmmax,plev) ! sym. nonlinear terms in u eqn.
  real(r8), intent(inout):: grfv1 (2*pmmax,plev) ! sym. nonlinear terms in v eqn.
  real(r8), intent(inout):: grlps2(2*pmmax)      ! antisym. undifferentiated term in Ps eq
  real(r8), intent(inout):: grt2  (2*pmmax,plev) ! antisym. undifferentiated term in t eq
  real(r8), intent(inout):: grq2  (2*pmmax,plev) ! antisym. undifferentiated term in q 
  real(r8), intent(out)  :: grz2  (2*pmmax,plev) ! antisym. undifferentiated term in z eq
  real(r8), intent(out)  :: grd2  (2*pmmax,plev) ! antisym. undifferentiated term in d eq
  real(r8), intent(inout):: grfu2 (2*pmmax,plev) ! antisym. nonlinear terms in u eqn.
  real(r8), intent(inout):: grfv2 (2*pmmax,plev) ! antisym. nonlinear terms in v eqn.
!
!---------------------------Local workspace-----------------------------
!
  real(r8) tmp1
  real(r8) tmp2
  real(r8) zxm(pmmax)  ! m*2dt/(a*cos(lat)**2)
  real(r8) zrcsj       ! 1./(a*cos(lat)**2)
  integer m            ! Fourier index
  integer k            ! level index
!
!-----------------------------------------------------------------------
!
  do m=1,2*nmmax(irow)
     tmp1 = 0.5*(grlps2(m) + grlps1(m))
     tmp2 = 0.5*(grlps2(m) - grlps1(m))
     grlps1(m) = tmp1
     grlps2(m) = tmp2
  end do
  do k=1,plev
     do m=1,2*nmmax(irow)
        tmp1 = 0.5*(grt2(m,k) + grt1(m,k))
        tmp2 = 0.5*(grt2(m,k) - grt1(m,k))
        grt1(m,k) = tmp1
        grt2(m,k) = tmp2
!
        tmp1 = 0.5*(grq2(m,k) + grq1(m,k))
        tmp2 = 0.5*(grq2(m,k) - grq1(m,k))
        grq1(m,k) = tmp1
        grq2(m,k) = tmp2
!
        tmp1 = 0.5*(grfu2(m,k) + grfu1(m,k))
        tmp2 = 0.5*(grfu2(m,k) - grfu1(m,k))
        grfu1(m,k) = tmp1
        grfu2(m,k) = tmp2
!
        tmp1 = 0.5*(grfv2(m,k) + grfv1(m,k))
        tmp2 = 0.5*(grfv2(m,k) - grfv1(m,k))
        grfv1(m,k) = tmp1
        grfv2(m,k) = tmp2
     end do
  end do
!
! Set constants
!
  zrcsj = 1./(cs(irow)*rearth)
!
! Combine constants with Fourier wavenumber m
!
  do m=1,nmmax(irow)
     zxm(m) = zrcsj*xm(m)
  end do
!
! Combine undifferentiated and longitudinal derivative terms for
! later use in Gaussian quadrature
!
  do k=1,plev
     do m=1,nmmax(irow)
        grd1(2*m-1,k) = - zxm(m)*grfu1(2*m,k)
        grd1(2*m,k)   =   zxm(m)*grfu1(2*m-1,k)
        grz1(2*m-1,k) = - zxm(m)*grfv1(2*m,k)
        grz1(2*m,k)   =   zxm(m)*grfv1(2*m-1,k)
!
        grd2(2*m-1,k) = - zxm(m)*grfu2(2*m,k)
        grd2(2*m,k)   =   zxm(m)*grfu2(2*m-1,k)
        grz2(2*m-1,k) = - zxm(m)*grfv2(2*m,k)
        grz2(2*m,k)   =   zxm(m)*grfv2(2*m-1,k)
     end do
  end do
  return
end subroutine dyn

#else

subroutine dyn(irow    ,grlps1  ,grt1    ,grq1    ,grz1    , &
               grd1    ,grfu1   ,grfv1   ,grlps2  ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Combine undifferentiated and longitudinally differentiated Fourier
! coefficient terms for later use in the Gaussian quadrature
!
! Computational note: Index "2*m-1" refers to the real part of the
! complex coefficient, and "2*m" to the imaginary.
!
! The naming convention is as follows:
!  - t, q, d, z refer to temperature, specific humidity, divergence
!     and vorticity
!  - "1" suffix to an array => symmetric component of current latitude
!     pair
!  - "2" suffix to an array => antisymmetric component
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: dyn.F90,v 1.5.2.1 2002/06/15 13:48:21 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use rgrid
  use commap
  use dynconst, only: rearth

  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: irow                 ! latitude pair index
  real(r8), intent(inout):: grlps1(2*pmmax)      ! sym. undifferentiated term in Ps eqn.
  real(r8), intent(inout):: grt1  (plev,2*pmmax) ! sym. undifferentiated term in t eqn.
  real(r8), intent(inout):: grq1  (plev,2*pmmax) ! sym. undifferentiated term in q
  real(r8), intent(out)  :: grz1  (plev,2*pmmax) ! sym. undifferentiated term in z eqn.
  real(r8), intent(out)  :: grd1  (plev,2*pmmax) ! sym. undifferentiated term in d eqn.
  real(r8), intent(inout):: grfu1 (plev,2*pmmax) ! sym. nonlinear terms in u eqn.
  real(r8), intent(inout):: grfv1 (plev,2*pmmax) ! sym. nonlinear terms in v eqn.
  real(r8), intent(inout):: grlps2(2*pmmax)      ! antisym. undifferentiated term in Ps eq
  real(r8), intent(inout):: grt2  (plev,2*pmmax) ! antisym. undifferentiated term in t eq
  real(r8), intent(inout):: grq2  (plev,2*pmmax) ! antisym. undifferentiated term in q
  real(r8), intent(out)  :: grz2  (plev,2*pmmax) ! antisym. undifferentiated term in z eq
  real(r8), intent(out)  :: grd2  (plev,2*pmmax) ! antisym. undifferentiated term in d eq
  real(r8), intent(inout):: grfu2 (plev,2*pmmax) ! antisym. nonlinear terms in u eqn.
  real(r8), intent(inout):: grfv2 (plev,2*pmmax) ! antisym. nonlinear terms in v eqn.
!
!---------------------------Local workspace-----------------------------
!
  real(r8) tmp1
  real(r8) tmp2
  real(r8) zxm(pmmax)       ! m*2dt/(a*cos(lat)**2)
  real(r8) zrcsj            ! 1./(a*cos(lat)**2)
  integer m                 ! Fourier index
  integer k                 ! level index
!
!-----------------------------------------------------------------------
!
  do m=1,2*nmmax(irow)
     tmp1 = 0.5*(grlps2(m) + grlps1(m))
     tmp2 = 0.5*(grlps2(m) - grlps1(m))
     grlps1(m) = tmp1
     grlps2(m) = tmp2
  end do
  do m=1,2*nmmax(irow)
     do k=1,plev
        tmp1 = 0.5*(grt2(k,m) + grt1(k,m))
        tmp2 = 0.5*(grt2(k,m) - grt1(k,m))
        grt1(k,m) = tmp1
        grt2(k,m) = tmp2
!
        tmp1 = 0.5*(grq2(k,m) + grq1(k,m))
        tmp2 = 0.5*(grq2(k,m) - grq1(k,m))
        grq1(k,m) = tmp1
        grq2(k,m) = tmp2
!
        tmp1 = 0.5*(grfu2(k,m) + grfu1(k,m))
        tmp2 = 0.5*(grfu2(k,m) - grfu1(k,m))
        grfu1(k,m) = tmp1
        grfu2(k,m) = tmp2
!
        tmp1 = 0.5*(grfv2(k,m) + grfv1(k,m))
        tmp2 = 0.5*(grfv2(k,m) - grfv1(k,m))
        grfv1(k,m) = tmp1
        grfv2(k,m) = tmp2
     end do
  end do
!
! Set constants
!
  zrcsj = 1./(cs(irow)*rearth)
!
! Combine constants with Fourier wavenumber m
!
  do m=1,nmmax(irow)
     zxm(m) = zrcsj*xm(m)
  end do
!
! Combine undifferentiated and longitudinal derivative terms for
! later use in Gaussian quadrature
!
  do k=1,plev
     do m=1,nmmax(irow)
        grd1(k,2*m-1) = - zxm(m)*grfu1(k,2*m)
        grd1(k,2*m)   =   zxm(m)*grfu1(k,2*m-1)
        grz1(k,2*m-1) = - zxm(m)*grfv1(k,2*m)
        grz1(k,2*m)   =   zxm(m)*grfv1(k,2*m-1)
!
        grd2(k,2*m-1) = - zxm(m)*grfu2(k,2*m)
        grd2(k,2*m)   =   zxm(m)*grfu2(k,2*m-1)
        grz2(k,2*m-1) = - zxm(m)*grfv2(k,2*m)
        grz2(k,2*m)   =   zxm(m)*grfv2(k,2*m-1)
     end do
  end do
  return
end subroutine dyn
#endif
