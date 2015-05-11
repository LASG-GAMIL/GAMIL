#include <misc.h>
#include <params.h>

module comspe

!----------------------------------------------------------------------- 
! 
! Purpose: Spectral space arrays
! 
! Method: 
! 
! Author: CCM Core Group
! $Author: erik $
! $Id: comspe.F90,v 1.2.6.1 2002/06/15 13:47:39 erik Exp $
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use infnan
  use pmgrid, only: plev, plat
  use pspect

  implicit none
!
! $Id: comspe.F90,v 1.2.6.1 2002/06/15 13:47:39 erik Exp $
! $Author: erik $
!
! Spectral space arrays
!
  real(r8) :: vz(psp,plev) = inf      ! Vorticity spectral coefficients
  real(r8) :: d(psp,plev)  = inf      ! Divergence spectral coefficients
  real(r8) :: t(psp,plev)  = inf      ! Temperature spectral coefficients
  real(r8) :: alps(psp)    = inf      ! Log-pressure spectral coefficients

  integer :: ncutoff       = bigint   ! Break-even point for vector lengths in GRCALC
  integer :: nalp(pmax)    = bigint   ! Pointer into polynomial arrays
#if ( defined SPMD )
  integer :: begm(0:plat-1) = bigint  ! Starting Fourier wavenumber owned by MPI task
  integer :: endm(0:plat-1) = bigint  ! Ending Fourier wavenumber owned by MPI task
#else
  integer :: begm(0:0)      = 1
  integer :: endm(0:0)      = pmmax
#endif

#if ( defined PVP )
  integer :: ncoefi(pmaxp) = bigint   ! Pointer to start of coefficient diagonals
  integer :: nm(pmax)      = bigint   ! Number of coeffs stored on a given diagonal
  integer :: nco2(pmax)    = bigint   ! Complex form of ncoefi
#else
  integer :: nstart(pmmax) = bigint   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)   = bigint   ! Length vectors for spectral arrays
#endif
  real(r8) :: alp(pspt,plat/2)  = inf ! Legendre polynomials
  real(r8) :: dalp(pspt,plat/2) = inf ! Legendre polynomial derivatives

end module comspe
