subroutine plevs0 (nlon, nlev, hyai, hybi, hyam, hybm, ps, pmid, pdel, ps0)
!
! $Id: plevs0.F,v 1.1 2000/09/01 22:58:44 rosinski Exp $
!
  implicit none
!-----------------------------------------------------------------------
!
! the following fields are set:
! hyai     fraction of reference pressure used for interface pressures
! hyam     fraction of reference pressure used for midpoint pressures
! hybi     fraction of surface pressure used for interface pressures
! hybm     fraction of surface pressure used for midpoint pressures
!
! Input arguments
!
  integer nlon,nlev
  real*8 ps(nlon)
  real*8 pmid(nlon,nlev)
  real*8 pdel(nlon,nlev)
  real*8 hyai(nlev+1)
  real*8 hybi(nlev+1)
  real*8 hyam(nlev)
  real*8 hybm(nlev)
  real*8 ps0
!
! Local workspace
!
  integer k
  real*8 pint(nlon,nlev+1)
!
! Set interface pressures
!
  do k=1,nlev+1
    pint(:,k) = hyai(k)*ps0 + hybi(k)*ps(:)
  end do
!
! Set midpoint pressures and layer thicknesses
!
  do k=1,nlev
    pmid(:,k) = hyam(k)*ps0 + hybm(k)*ps(:)
    pdel(:,k) = pint(:,k+1) - pint(:,k)
  end do
  
  return
end subroutine plevs0
