#include <misc.h>
#include <preproc.h>

subroutine Biogeophysics2 (clm)

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
! This is the main subroutine to execute the calculation of soil/snow and
! ground temperatures and update surface fluxes based on the new ground
! temperature 
! 
! Method:
! Calling sequence is:
! Biogeophysics2:                    surface biogeophysics driver
!    -> SoilTemperature:             soil/snow and ground temperatures      
!          -> SoilTermProp           thermal conductivities and heat 
!                                     capacities        
!          -> Tridiagonal            tridiagonal matrix solution            
!          -> PhaseChange            phase change of liquid/ice contents        
!
! (1) Snow and soil temperatures
!     o The volumetric heat capacity is calculated as a linear combination 
!       in terms of the volumetric fraction of the constituent phases. 
!     o The thermal conductivity of soil is computed from 
!       the algorithm of Johansen (as reported by Farouki 1981), and the 
!       conductivity of snow is from the formulation used in
!       SNTHERM (Jordan 1991).
!     o Boundary conditions:  
!       F = Rnet - Hg - LEg (top),  F = 0 (base of the soil column).
!     o Soil / snow temperature is predicted from heat conduction 
!       in 10 soil layers and up to 5 snow layers. 
!       The thermal conductivities at the interfaces between two 
!       neighboring layers (j, j+1) are derived from an assumption that 
!       the flux across the interface is equal to that from the node j 
!       to the interface and the flux from the interface to the node j+1. 
!       The equation is solved using the Crank-Nicholson method and 
!       results in a tridiagonal system of equations.
!
! (2) Phase change (see PhaseChange.F90)
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics2.F90,v 1.2.10.3 2002/06/15 13:50:13 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : hvap, cpair, grav, vkc, tfrz, sb
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer j                          ! do loop index
  real(r8) fact(clm%snl+1 : nlevsoi) ! used in computing tridiagonal matrix
  real(r8) egsmax ! max. evaporation which soil can provide at one time step [mm/s]
  real(r8) egidif ! the excess of evaporation over "egsmax" [mm/s]
  real(r8) xmf    ! total latent heat of phase change of ground water [W/m2]
  real(r8) tinc   ! temperature difference between two time steps [K]

!----End Variable List--------------------------------------------------

!
! Determine soil temperatures including surface soil temperature
!

  call SoilTemperature(clm      , clm%tssbef, clm%htvp, clm%emg, clm%cgrnd, &
                       clm%dlrad, clm%tg    , xmf     , fact )

!
! Correct fluxes to present soil temperature
!

  tinc = clm%t_soisno(clm%snl+1) - clm%tssbef(clm%snl+1)
  clm%eflx_sh_grnd =  clm%eflx_sh_grnd + tinc*clm%cgrnds 
  clm%qflx_evap_soi =  clm%qflx_evap_soi + tinc*clm%cgrndl

!
! egidif holds the excess energy if all water is evaporated from
! the top soil layer during the timestep.  This energy is added to
! the sensible heat flux.
!

  egsmax = (clm%h2osoi_ice(clm%snl+1)+clm%h2osoi_liq(clm%snl+1)) / clm%dtime
  egidif = max( 0._r8, clm%qflx_evap_soi - egsmax )
  clm%qflx_evap_soi = min ( clm%qflx_evap_soi, egsmax )
  clm%eflx_sh_grnd = clm%eflx_sh_grnd + clm%htvp*egidif

!
! Ground heat flux
!

  clm%eflx_soil_grnd = clm%sabg + clm%dlrad + (1-clm%frac_veg_nosno)*clm%emg*clm%forc_lwrad &
       - clm%emg*sb*clm%tssbef(clm%snl+1)**3*(clm%tssbef(clm%snl+1) + 4.*tinc) &
       - (clm%eflx_sh_grnd+clm%qflx_evap_soi*clm%htvp)

!
! Total fluxes (vegetation + ground)
!

  clm%eflx_sh_tot = clm%eflx_sh_veg + clm%eflx_sh_grnd
  clm%qflx_evap_tot = clm%qflx_evap_veg + clm%qflx_evap_soi
  clm%eflx_lh_tot= hvap*clm%qflx_evap_veg + clm%htvp*clm%qflx_evap_soi   ! (account for sublimation)


!
! Assign ground evaporation to sublimation from soil ice or to dew
! on snow or ground 
!

  clm%qflx_evap_grnd = 0.
  clm%qflx_sub_snow = 0.
  clm%qflx_dew_snow = 0.
  clm%qflx_dew_grnd = 0.

  if (clm%qflx_evap_soi >= 0.) then
     ! Do not allow for sublimation in melting (melting ==> evap. ==> sublimation)
     clm%qflx_evap_grnd = min(clm%h2osoi_liq(clm%snl+1)/clm%dtime, clm%qflx_evap_soi)
     clm%qflx_sub_snow = clm%qflx_evap_soi - clm%qflx_evap_grnd
  else
     if (clm%tg < tfrz) then
        clm%qflx_dew_snow = abs(clm%qflx_evap_soi)
     else
        clm%qflx_dew_grnd = abs(clm%qflx_evap_soi)
     endif
  endif

!
! Outgoing long-wave radiation from vegetation + ground
!

! For conservation we put the increase of ground longwave to outgoing
  clm%eflx_lwrad_out = clm%ulrad &
       + (1-clm%frac_veg_nosno)*(1.-clm%emg)*clm%forc_lwrad &
       + (1-clm%frac_veg_nosno)*clm%emg*sb * clm%tssbef(clm%snl+1)**4 &
  + 4.*clm%emg*sb*clm%tssbef(clm%snl+1)**3*tinc

!
! Radiative temperature
!

  clm%t_rad = (clm%eflx_lwrad_out/sb)**0.25

!
! Soil energy balance check
!

  clm%errsoi = 0. 
  do j = clm%snl+1, nlevsoi
     clm%errsoi = clm%errsoi - (clm%t_soisno(j)-clm%tssbef(j))/fact(j) 
  enddo
  clm%errsoi = clm%errsoi + clm%eflx_soil_grnd - xmf

!
! Variables needed by history tape
!

 clm%dt_grnd        = tinc
 clm%eflx_lh_vege   = (clm%qflx_evap_veg - clm%qflx_tran_veg) * hvap
 clm%eflx_lh_vegt   = clm%qflx_tran_veg * hvap       
 clm%eflx_lh_grnd   = clm%qflx_evap_soi * clm%htvp
 clm%eflx_lwrad_net = clm%eflx_lwrad_out -  clm%forc_lwrad  

end subroutine Biogeophysics2
