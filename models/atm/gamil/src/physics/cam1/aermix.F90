#include <misc.h>
#include <params.h>
subroutine aermix(lchnk, ncol, pint, sulfmix, aermmr, rh)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set global mean tropospheric aerosol
! 
! Method: 
! Specify aerosol mixing ratio and compute relative humidity for later
! adjustment of aerosol optical properties. Aerosol mass mixing ratio
! is specified so that the column visible aerosol optical depth is a
! specified global number (tauvis). This means that the actual mixing
! ratio depends on pressure thickness of the lowest three atmospheric
! layers near the surface.
!
! Optical properties and relative humidity parameterization are from:
!
! J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
! and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
! 16 April 1993
!
! Visible (vis) here means 0.5-0.7 micro-meters
! Forward scattering fraction is taken as asymmetry parameter squared
! 
! Author: J.T. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <ptrrgrid.h>
!------------------------------Commons----------------------------------
#include <crdcon.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: pint(pcols,pverrp)   ! Rad level interface press. (dynes/cm2)
   real(r8), intent(in) :: sulfmix(pcols,pver)  ! time interpolated sulfate mass mixing ratio
!
! Output arguments
!
   real(r8), intent(out) :: aermmr(pcols,pverr)  ! Rad level aerosol mass mixing ratio
   real(r8), intent(out) :: rh(pcols,pverr)      ! Rad level relative humidity (fraction)

!
!---------------------------Local variables-----------------------------
!
   integer i          ! Longitude index
   integer k          ! Level index
!
   real(r8) kaervs    ! Visible extinction coefficiant of aerosol (m2/g)
   real(r8) omgvis    ! Visible single scattering albedo
   real(r8) gvis      ! Visible scattering asymmetry parameter
   real(r8) rhcnst    ! Constant relative humidity factor
!
! Relative humidity factor
!
   real(r8) rhfac     ! Multiplication factor for kaer
   real(r8) rhpc      ! Level relative humidity in %
   real(r8) a0        ! Constant in relative humidity mult factor
   real(r8) a1        ! Constant in relative humidity mult factor
   real(r8) a2        ! Constant in relative humidity mult factor
   real(r8) a3        ! Constant in relative humidity mult factor
!
!--------------------------Data Statements------------------------------
!
   data a0 / -9.2906106183    /
   data a1 /  0.52570211505   /
   data a2 / -0.0089285760691 /
   data a3 /  5.0877212432e-05/
!
   data kaervs / 5.3012   /
   data omgvis / 0.999999 /
   data gvis   / 0.694889 /
   data rhcnst /  .80     /
!
!-----------------------------------------------------------------------
!
! Set relative humidity and factor; then aerosol amount.
!
   do i=1,ncol
      do k=1,pverr
!
         rh(i,k) = rhcnst
!
! Compute relative humidity factor for the extinction coefficiant; this
! factor accounts for the dependence of size distribution on relative
! humidity:
!
         if ( rh(i,k) > .90 ) then
            rhfac = 2.8
         else if (rh(i,k) < .60 ) then
            rhfac = 1.0
         else
            rhpc  = 100. * rh(i,k)
            rhfac = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
         endif
!
! Compute aerosol mass mixing ratio for specified levels (1.e4 factor is
! for units conversion of the extinction coefficiant from m2/g to cm2/g)
!
         if ( k >= pverrp-mxaerl ) then
            aermmr(i,k) = gravit*tauvis / &
                          (1.e4*kaervs*rhfac*(1.-omgvis*gvis*gvis) * &
                          (pint(i,pverrp)-pint(i,pverrp-mxaerl))) + &
                          sulfmix(i,k)
         else
            aermmr(i,k) = sulfmix(i,k)
         endif
!
      enddo
   enddo
!
   return
end subroutine aermix

