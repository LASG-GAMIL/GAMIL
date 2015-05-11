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
! $Id: comspe.F90,v 1.3.6.1 2002/06/15 13:48:20 erik Exp $
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use infnan
  use pmgrid, only: plev, plat
  use pspect

  implicit none

  real(r8) :: vz(psp,plev)   = inf     ! Vorticity spectral coefficients
  real(r8) :: d(psp,plev)    = inf     ! Divergence spectral coefficients
  real(r8) :: t(psp,plev)    = inf     ! Temperature spectral coefficients
  real(r8) :: q(psp,plev)    = inf     ! Moisture     spectral coefficients
  real(r8) :: alps(psp)      = inf     ! Log-pressure spectral coefficients
  real(r8) :: hs(psp,plev)   = inf     ! hydrostatic matrix for "real" atmosphere
  real(r8) :: hsnm(psp,plev) = inf   ! vertical normal modes of "hs"
  real(r8) :: dsnm(psp,plev) = inf   ! vertical normal modes of "ds"
  real(r8) :: dnm(psp,plev)  = inf   ! vertical normal modes of "d"
  real(r8) :: vznm(psp,plev) = inf   ! vertical normal modes of "vz"
  real(r8) :: lnpstar(psp)   = inf   ! ln (Ps*) (SLD term; Ritchie & Tanguay, 1995)
  real(r8) :: a0nm(psp)      = inf   ! wave # coefs (use in vert normal mode space)
  real(r8) :: bmnm(psp)      = inf   ! wave # coefs (use in vert normal mode space)
  real(r8) :: bpnm(psp)      = inf   ! wave # coefs (use in vert normal mode space)
  real(r8) :: atri(psp)      = inf   ! wave # coefs (use in vert normal mode space)
  real(r8) :: btri(psp)      = inf   ! wave # coefs (use in vert normal mode space)
  real(r8) :: ctri(psp)      = inf   ! wave # coefs (use in vert normal mode space)

  integer :: ncutoff         = bigint   ! Break-even point for vector lengths in GRCALC
  integer :: nalp(pmax)      = bigint   ! Pointer into polynomial arrays
#if ( defined SPMD )
  integer :: begm(0:plat-1)  = bigint
  integer :: endm(0:plat-1)  = bigint
#else
  integer :: begm(0:0)       = 1
  integer :: endm(0:0)       = pmmax
#endif

  integer :: ncoefi(pmaxp)   = bigint   ! Pointer to start of coefficient diagonals
  integer :: nm(pmax)        = bigint   ! Number of coeffs stored on a given diagonal
  integer :: nco2(pmax)      = bigint   ! Complex form of ncoefi
  integer :: nstart(pmmax)   = bigint   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)     = bigint   ! Length vectors for spectral arrays

  real(r8) :: alp(pspt,plat/2)  = inf   ! Legendre polynomials
  real(r8) :: dalp(pspt,plat/2) = inf   ! Legendre polynomial derivatives

end module comspe
