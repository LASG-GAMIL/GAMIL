#include <misc.h>
#include <params.h>

subroutine plevs0 (ncol    , ncold   ,nver    ,ps      ,pint    , &
                   pmid    ,pdel)


!! (wanhui 2003.05.17)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
! 
! Method: 
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------
!
! $Id: plevs0.F90,v 1.1.2.1 2002/06/15 13:46:59 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none
#include <comhyb.h>

!-----------------------------------------------------------------------
  integer , intent(in)  :: ncol               ! Longitude dimension
  integer , intent(in)  :: ncold              ! Declared longitude dimension
  integer , intent(in)  :: nver               ! vertical dimension
  real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
  real(r8), intent(out) :: pint(ncold,nver+1) ! Pressure at model interfaces
  real(r8), intent(out) :: pmid(ncold,nver)   ! Pressure at model levels
  real(r8), intent(out) :: pdel(ncold,nver)   ! Layer thickness (pint(k+1) - pint(k))
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set interface pressures
!
  do k=1,nver+1
     do i=1,ncol
!!      pint(i,k) = hyai(k)*ps0 + hybi(k)*ps(i)
        pint(i,k) = sig(k)*( ps(i)-pmtop*100.0 ) + pmtop*100.0
     end do
  end do
!
! Set midpoint pressures and layer thicknesses
!
  do k=1,nver
     do i=1,ncol
!!      pmid(i,k) = hyam(k)*ps0 + hybm(k)*ps(i)
        pmid(i,k) = sigl(k)*( ps(i)-pmtop*100.0 ) + pmtop*100.0
        pdel(i,k) = pint(i,k+1) - pint(i,k)
     end do
  end do

  return
end subroutine plevs0

