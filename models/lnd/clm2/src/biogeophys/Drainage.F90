#include <misc.h>
#include <preproc.h>

subroutine Drainage (clm,  zwice, vol_liq, s,   zwt, &
                     fcov, hk,    dhkdw,   dwat ) 

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
! Calculate subsurface drainage
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
! $Id: Drainage.F90,v 1.2.10.3 2002/06/15 13:50:14 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

  real(r8), intent(in) :: vol_liq(1:nlevsoi) ! partial volume of liquid water in layer [-]
  real(r8), intent(in) :: zwice              ! the sum of ice mass of soil [kg/m2]
  real(r8), intent(in) :: s(1:nlevsoi)       ! wetness of soil (including ice) [-]
  real(r8), intent(in) :: zwt                ! water table depth [m]
  real(r8), intent(in) :: fcov               ! fractional area with water table at surface [-]
  real(r8), intent(in) :: hk(1:nlevsoi)      ! hydraulic conductivity [mm/s]
  real(r8), intent(in) :: dhkdw(1:nlevsoi)   ! d(hk)/d(vol_liq)
  real(r8), intent(in) :: dwat(1:nlevsoi)    ! change in soil water [-]

!----Local Variables----------------------------------------------------

  integer i                ! do loop index
  real(r8) xs              ! excess soil water above saturation [mm]
  real(r8) dzmm(1:nlevsoi) ! layer thickness [mm]
  real(r8) watmin          ! minimum soil moisture [mm]
  real(r8) hksum           ! summation of hydraulic cond for layers 6->9 [mm/s]

! For Z.-L. Yang & G.-Y. Niu's modification

  real(r8) zsat            ! hydraulic conductivity weighted soil thickness [m mm/s]
  real(r8) wsat            ! hydraulic conductivity weighted soil wetness [m mm/s]
  real(r8) qflx_drain_wet  ! subsurface runoff from "wet" part [mm/s]
  real(r8) qflx_drain_dry  ! subsurface runoff from "dry" part [mm/s]
  real(r8) dzksum          ! hydraulic conductivity weighted soil thickness [m mm/s]

!----End Variable List--------------------------------------------------

!
! Streamflow and total runoff
!

!
! Convert layer thicknesses from m to mm
!

  do i = 1,nlevsoi 
     dzmm(i) = clm%dz(i)*1.e3
  enddo

!
! The amount of streamflow is assumed maintained by flow from the 
! lowland water table with different levels contributing according to 
! their thickness and saturated hydraulic conductivity, i.e. a given 
! layer below the water table interface loses water at a rate per unit 
! depth given by qflx_drain*hk/(sum over all layers below this water table 
! of hk*dz). Because this is a slow smooth process, and not strongly 
! coupled to water in any one layer, it should remain stable for 
! explicit time differencing. Hence, for simplicity it is removed
! explicitly prior to the main soil water calculation.
! Another assumption: no subsurface runoff for ice mixed soil 
! Zong-Liang Yang & G.-Y. Niu                                         
!

  clm%qflx_drain = 0.
  qflx_drain_wet = 0.
  qflx_drain_dry = 0.
  
  hksum = 0.
  do i = 6,nlevsoi-1                  
     hksum = hksum + hk(i)
  enddo
  if (zwice <= 0. .AND. hksum > 0.) then
     zsat = 0.                                                    
     wsat = 0.                                                    
     dzksum = 0.                                                  
     do i = 6,nlevsoi-1                  
        zsat = zsat + clm%dz(i)*hk(i)                               
        wsat = wsat + s(i)*clm%dz(i)*hk(i)                         
        dzksum  = dzksum   + hk(i)*clm%dz(i)                       
     enddo
     wsat = wsat / zsat                                         

     qflx_drain_dry = (1.-fcov)*4.e-2* wsat ** (2.*clm%bsw(1)+3.)  ! mm/s
     qflx_drain_wet = fcov * 1.e-5 * exp(-zwt)                     ! mm/s
     clm%qflx_drain = qflx_drain_dry + qflx_drain_wet
     
     do i = 6, nlevsoi-1                 
        clm%h2osoi_liq(i) = clm%h2osoi_liq(i) &
             - clm%dtime*clm%qflx_drain*clm%dz(i)*hk(i)/dzksum                                  
     enddo
  endif

!
! Limit h2osoi_liq to be greater than or equal to watmin. 
! Get water needed to bring h2osoi_liq equal watmin from lower layer. 
!

  watmin = 0.0
  do i = 1, nlevsoi-1
     if (clm%h2osoi_liq(i) < 0.) then
        xs = watmin-clm%h2osoi_liq(i)
     else
        xs = 0.
     end if
     clm%h2osoi_liq(i  ) = clm%h2osoi_liq(i  ) + xs
     clm%h2osoi_liq(i+1) = clm%h2osoi_liq(i+1) - xs
  end do
  i = nlevsoi
  if (clm%h2osoi_liq(i) < watmin) then
     xs = watmin-clm%h2osoi_liq(i)
  else
     xs = 0.
  end if
  clm%h2osoi_liq(i) = clm%h2osoi_liq(i) + xs
  clm%qflx_drain = clm%qflx_drain - xs/clm%dtime

!
! Determine water in excess of saturation
!

  xs = max(0._r8, clm%h2osoi_liq(1)-(clm%pondmx+clm%eff_porosity(1)*dzmm(1)))
  if (xs > 0.) clm%h2osoi_liq(1) = clm%pondmx+clm%eff_porosity(1)*dzmm(1)
  
  do i = 2,nlevsoi 
     xs = xs + max(clm%h2osoi_liq(i)-clm%eff_porosity(i)*dzmm(i), 0._r8)  ! [mm]
     clm%h2osoi_liq(i) = min(clm%eff_porosity(i)*dzmm(i), clm%h2osoi_liq(i))
  enddo

!
! Sub-surface runoff and drainage 
!

  clm%qflx_drain = clm%qflx_drain + xs/clm%dtime  &
       + hk(nlevsoi) + dhkdw(nlevsoi)*dwat(nlevsoi) ! [mm/s]
  
!
! Set imbalance for snow capping
!

  clm%qflx_qrgwl = clm%qflx_snowcap
  
!
! Implicit evaporation term is now zero
!

  clm%eflx_impsoil = 0.

!
! Renew the ice and liquid mass due to condensation
!
  
  if (clm%snl+1 >= 1) then
     clm%h2osoi_liq(1) = clm%h2osoi_liq(1) + clm%qflx_dew_grnd*clm%dtime
     clm%h2osoi_ice(1) = clm%h2osoi_ice(1) + (clm%qflx_dew_snow-clm%qflx_sub_snow)* &
                         clm%dtime
  endif

end subroutine Drainage
