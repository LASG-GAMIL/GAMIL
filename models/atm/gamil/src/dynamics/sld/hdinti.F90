#include <misc.h>
#include <params.h>
subroutine hdinti(rearth, deltat)
!-----------------------------------------------------------------------
!
! Purpose:
! Time independent initialization for the horizontal diffusion.
!
! Author:  D. Williamson
!
!-----------------------------------------------------------------------
!
! $Id: hdinti.F90,v 1.2.42.1 2002/06/15 13:48:25 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  implicit none
#include <comhd.h>
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: rearth ! radius of the earth
  real(r8), intent(in)   :: deltat ! time step
!
!---------------------------Local workspace-----------------------------
!
  integer k                        ! level index
  integer n                        ! n-wavenumber index
!
!-----------------------------------------------------------------------
!
! Top level for del**4 diffusion, set for 18-level model
!
  kmnhd4 = 5
!
! Bottom level for increased del**2 diffusion (kmxhd2 < kmnhd4)
!
  kmxhd2 = 3
!
! Initialize physical constants for courant number based spect truncation
!
  nmaxhd = ptrk
  cnlim  = 0.999          ! maximum allowable Courant number
  cnfac  = deltat*float(nmaxhd)/rearth
!
! Initialize arrays used for courant number based spectral truncation
!
  do k=1,plev
     nindex(k) = 2*nmaxhd
  end do
!
! Set the Del^2 and Del^4 diffusion coefficients for each wavenumber
!
  hdfst2(1) = 0.
  hdfsd2(1) = 0.
!
  hdfst4(1) = 0.
  hdfsd4(1) = 0.
  do n=2,pnmax
     hdfst2(n) = dif2 * (n*(n-1)  ) / rearth**2
     hdfsd2(n) = dif2 * (n*(n-1)-2) / rearth**2

     hdfst4(n) = dif4 * (n*(n-1)*n*(n-1)  ) / rearth**4
     hdfsd4(n) = dif4 * (n*(n-1)*n*(n-1)-4) / rearth**4
  end do
!
  return
end subroutine hdinti

