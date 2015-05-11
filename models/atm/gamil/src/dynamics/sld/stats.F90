#include <misc.h>
#include <params.h>
subroutine stats(lat     ,pint    ,pdel    ,pstar   , &
                 div     ,t       ,q       ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Accumulation of diagnostic statistics for 1 latitude.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id: stats.F90,v 1.5.28.1 2002/06/15 13:48:32 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use commap

  implicit none

#include <comsta.h>

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: lat                ! latitude index (S->N)

  real(r8), intent(in)   :: pint (plond,plevp) ! pressure at model interfaces
  real(r8), intent(in)   :: pdel (plond,plev)  ! pdel(k) = pint(k+1) - pint(k)
  real(r8), intent(in)   :: pstar(plond)       ! ps + psr (surface pressure)
  real(r8), intent(in)   :: div  (plond,plev)  ! divergence
  real(r8), intent(in)   :: t    (plond,plev)  ! temperature
  real(r8), intent(in)   :: q    (plond,plev)  ! moisture
  integer , intent(in)   :: nlon               ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) prat            ! pdel(i,k)/pint(i,plevp)
  integer  i,k             ! longitude, level indices
!
!-----------------------------------------------------------------------
!
! Compute statistics for current latitude line
!
  rmsz (lat) = 0.
  rmsd (lat) = 0.
  rmst (lat) = 0.
  stq  (lat) = 0.
  psurf(lat) = 0.

  do i = 1,nlon
     psurf(lat) = psurf(lat) + pstar(i)
  end do
  do k = 1,plev
     do i = 1,nlon
        prat      = pdel(i,k)/pint(i,plevp)
        rmsd(lat) = rmsd(lat) + div(i,k)*div(i,k)*prat
        rmst(lat) = rmst(lat) + (t(i,k)**2)*prat
        stq (lat) = stq(lat) + q(i,k)*pdel(i,k)
     end do
  end do
  psurf(lat) = w(lat)*psurf(lat)/nlon
  rmsd (lat) = w(lat)*rmsd(lat)/nlon
  rmst (lat) = w(lat)*rmst(lat)/nlon
  stq  (lat) = w(lat)*stq(lat)/nlon
! 
  return
end subroutine stats

