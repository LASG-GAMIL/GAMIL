#include <misc.h>
#include <preproc.h>

subroutine SurfaceRunoff (clm, zwice, vol_liq, s, zwt, &
                          fcov ) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                           M  available land surface process model.
!  M --COMMUNITY LAND MODEL--  C
!  C                           L
!  LMCLMCLMCLMCLMCLMCLMCLMCLMCLM
!
!-----------------------------------------------------------------------
! Purpose:
! Calculate surface runoff
!
! Method:
! The original code for total runoff was provided by Robert E. Dickinson
! based on the following properties:  exponential decrease of Ksat, a
! water table level determination level including highland and lowland 
! levels and fractional area of wetland (water table above the surface). 
! Runoff is parameterized from the lowlands in terms of precip 
! incident on wet areas and a base flow, where these are estimated 
! using ideas from TOPMODEL.
!
! The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
! o  using a new method to determine water table depth and
!    the fractional wet area (fcov)
! o  computing runoff (surface and subsurface) from this
!    fraction and the remaining fraction (i.e. 1-fcov)
! o  for the 1-fcov part, using BATS1e method to compute
!    surface and subsurface runoff.
!
! The original code on soil moisture and runoff were provided by 
! R. E. Dickinson in July 1996.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 12 November 1999:  Z.-L. Yang and G.-Y. Niu
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: SurfaceRunoff.F90,v 1.2.10.3 2002/06/15 13:50:20 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : denice, denh2o
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

  real(r8), intent(out) :: zwice              ! the sum of ice mass of soil (kg/m2)
  real(r8), intent(out) :: vol_liq(1:nlevsoi) ! partial volume of liquid water in layer [-]
  real(r8), intent(out) :: s(1:nlevsoi)       ! wetness of soil (including ice) [-]
  real(r8), intent(out) :: zwt                ! water table depth [m]
  real(r8), intent(out) :: fcov               ! fractional area with water table at surface [-]

!----Local Variables----------------------------------------------------

  integer i                   ! do loop index 

  real(r8) vol_ice(1:nlevsoi) ! partial volume of ice lens in layer [-]
  real(r8) zmean              ! The surface soil layers contributing to runoff [m]
  real(r8) wmean              ! The averaged soil wetness depth in surface soil layers [m]
  real(r8) fz                 ! coefficient for water table depth

!----End Variable List--------------------------------------------------

!
! Porosity of soil, partial volume of ice and liquid
!

  zwice = 0.
  do i = 1,nlevsoi 
     zwice = zwice + clm%h2osoi_ice(i)
     vol_ice(i) = min(clm%watsat(i), clm%h2osoi_ice(i)/(clm%dz(i)*denice))
     clm%eff_porosity(i) = clm%watsat(i)-vol_ice(i)
     vol_liq(i) = min(clm%eff_porosity(i), clm%h2osoi_liq(i)/(clm%dz(i)*denh2o))
  enddo

!
! Calculate wetness of soil
!

  do i = 1,nlevsoi
     s(i) = min(1._r8,(vol_ice(i)+vol_liq(i))/clm%watsat(i))
  end do

!
! Determine water table 
!

  wmean = 0.                                                  
  fz    = 1.0                                                
  do i  = 1, nlevsoi                                        
     wmean = wmean + s(i)*clm%dz(i)                          
  enddo
  zwt = fz * (clm%zi(nlevsoi) - wmean)                   
  
!
! Saturation fraction
!

  fcov = clm%wtfact*min(1._r8,exp(-zwt))
  
  wmean = 0.                                               
  zmean = 0.                                              
  do i = 1, 3                                          
     zmean = zmean + clm%dz(i)                          
     wmean = wmean + s(i) * clm%dz(i)                 
  enddo
  wmean = wmean / zmean                           

!
! If top soil layer is impermeable then all qflx_top_soil goes to surface runoff
!

  if (clm%eff_porosity(1) < clm%wimp) then
   clm%qflx_surf =  max(0._r8, fcov*clm%qflx_top_soil) + &
                    max(0._r8, (1.-fcov)*clm%qflx_top_soil)        
  else
   clm%qflx_surf =  max(0._r8, fcov*clm%qflx_top_soil) + &
                    max(0._r8, (1.-fcov)*min(1._r8,wmean**4)*clm%qflx_top_soil)        
  endif

end subroutine SurfaceRunoff
