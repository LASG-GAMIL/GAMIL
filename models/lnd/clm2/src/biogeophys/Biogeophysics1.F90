#include <misc.h>
#include <preproc.h>

subroutine Biogeophysics1 (clm)

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
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
!  Biogeophysics1:                   surface biogeophysics driver
!    -> QSat:                        saturated vapor pressure, specific humidity,
!                                     and derivatives at ground surface
!    -> SurfaceRadiation:            surface solar radiation
!    -> BareGroundFluxes:            surface fluxes for bare soil or
!                                     snow-covered vegetation patches
!          -> MoninObukIni:          first-guess Monin-Obukhov length and
!                                     wind speed including stability effect
!          -> FrictionVelocity:      friction velocity and potential 
!                                     temperature and humidity profiles
!    -> CanopyFluxes:                leaf temperature and surface fluxes
!                                     for vegetated patches 
!          -> QSat                   saturated vapor pressure, specific humidity,
!                                     and derivatives at leaf surface
!          -> MoninObukIni           first-guess Monin-Obukhov length and 
!                                     wind speed including stability effect
!          -> FrictionVelocity       friction velocity and potential
!                                     temperature and humidity profiles
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for sunlit leaves
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for shaded leaves
!          -> SensibleHCond          sensible heat conductance for air, leaf,
!                                     and ground
!          -> LatentHCond            latent heat conductance for ground and
!                                     leaf
!          -> QSat                   recalculation of saturated vapor pressure,
!                                     specific humidity, and derivatives at
!                                     leaf sfc using updated leaf temperature
!  Leaf temperature
!   Foliage energy conservation is given by the foliage energy budget 
!   equation:
!                  Rnet - Hf - LEf = 0 
!   The equation is solved by Newton-Raphson iteration. The
!   iteration includes the calculation of photosynthesis and 
!   stomatal resistance and the integration of turbulent flux profiles. 
!   The sensible and latent heat transfer between foliage (f) and atmosphere (a)
!   and ground (g) is linked by the equations:  
!                  Ha = Hf + Hg and Ea = Ef + Eg
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics1.F90,v 1.4.6.3 2002/06/15 13:50:12 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : denh2o, denice, roverg, hvap, hsub, istice, istwet 
  use clm_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i,j      ! loop indices
  real(r8) qred    ! soil surface relative humidity factor [-]
  real(r8) avmuir  ! IR inverse optical depth per unit leaf area [-]
  real(r8) eg      ! water vapor pressure at temperature T [pa]
  real(r8) qsatg   ! saturated humidity [kg/kg]
  real(r8) degdT   ! d(eg)/dT
  real(r8) qsatgdT ! d(qsatg)/dT
  real(r8) fac     ! soil wetness of surface layer [-]
  real(r8) psit    ! negative potential of soil [mm]
  real(r8) hr      ! relative humidity [-]
  real(r8) wx      ! partial volume of ice and water of surface layer [-]

!----End Variable List--------------------------------------------------

!
! Initial set 
!

  clm%eflx_sh_tot    = 0.
  clm%qflx_evap_tot  = 0.
  clm%eflx_lh_tot    = 0.
  clm%eflx_sh_veg    = 0.  
  clm%qflx_evap_veg  = 0.  
  clm%qflx_tran_veg  = 0.  
  clm%cgrnd          = 0._r8
  clm%cgrnds         = 0._r8
  clm%cgrndl         = 0._r8
  clm%t_ref2m        = 0.

!
! Ground and soil temperatures from previous time step
!

  clm%tg = clm%t_soisno(clm%snl+1)
  do i = clm%snl+1, nlevsoi
     clm%tssbef(i) = clm%t_soisno(i)
  enddo

!
! Saturated vapor pressure, specific humidity and their derivatives
! at ground surface
!

  qred = 1.
  if (clm%itypwat/=istwet .AND. clm%itypwat/=istice) then
     wx   = (clm%h2osoi_liq(1)/denh2o+clm%h2osoi_ice(1)/denice)/clm%dz(1)
     fac  = min(1._r8, wx/clm%watsat(1))
     fac  = max( fac, 0.01_r8 )
     psit = -clm%sucsat(1) * fac ** (- clm%bsw(1))
     psit = max(clm%smpmin, psit)
     hr   = exp(psit/roverg/clm%tg)
     qred = (1.-clm%frac_sno)*hr + clm%frac_sno
  endif

  call QSat(clm%tg, clm%forc_pbot, eg, degdT, qsatg, &
            qsatgdT)

  clm%qg = qred*qsatg  
  clm%dqgdT = qred*qsatgdT

  if (qsatg > clm%forc_q .AND. clm%forc_q > qred*qsatg) then
     clm%qg = clm%forc_q
     clm%dqgdT = 0.
  endif

!
! Emissivity
!

  if (clm%h2osno>0. .OR.clm%itypwat==istice) then
     clm%emg = 0.97
  else
     clm%emg = 0.96
  endif
  avmuir=1.
  clm%emv=1.-exp(-(clm%elai+clm%esai)/avmuir)

!
! Latent heat of vaporization or sublimation. We arbitrarily assume that
! sublimation occurs only if h2osoi_liq = 0
!

  clm%htvp = hvap
  if (clm%h2osoi_liq(clm%snl+1) <= 0. .AND. clm%h2osoi_ice(clm%snl+1) > 0.) clm%htvp = hsub

!
! Switch between vaporization and sublimation causes rapid solution
! separation in perturbation growth test
!

#if (defined PERGRO)
  clm%htvp = hvap
#endif

!
! Roughness lengths
!

  if (clm%frac_sno > 0.) then
     clm%z0mg = clm%zsno
     clm%z0hg = clm%z0mg            ! initial set only
     clm%z0qg = clm%z0mg            ! initial set only
  else
     clm%z0mg = clm%zlnd
     clm%z0hg = clm%z0mg            ! initial set only
     clm%z0qg = clm%z0mg            ! initial set only
  endif

  clm%z0m = clm%z0mr*clm%htop
  clm%displa = clm%displar*clm%htop
  clm%z0mv = clm%z0m
  clm%z0hv = clm%z0mv
  clm%z0qv = clm%z0mv

!
! Potential, virtual potential temperature, and wind speed at the 
! reference height
!

  clm%beta=1.
  clm%zii = 1000.
  clm%thm = clm%forc_t + 0.0098*clm%forc_hgt_t              
  clm%thv = clm%forc_th*(1.+0.61*clm%forc_q)
  clm%ur = max(1.0_r8,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))

!
! Surface Radiation
!

  call SurfaceRadiation (clm)

!
! Surface Temperature and Fluxes
!

!
! BARE SOIL OR SNOW-COVERED VEGETATION
! Ground fluxes
! NOTE: in the current scheme clm%frac_veg_nosno is EITHER 1 or 0
!

  if (clm%frac_veg_nosno == 0) then

     call BareGroundFluxes (clm%tg,     clm%thm,   clm%qg,    clm%thv,   clm%z0mg,   &
                            clm%z0hg,   clm%z0qg,  clm%dqgdT, clm%htvp,  clm%beta,   &
                            clm%zii,    clm%ur,    clm%dlrad, clm%ulrad, clm%cgrnds, &
                            clm%cgrndl, clm%cgrnd, clm    )
     clm%psnsun = 0.
     clm%psnsha = 0. !put these lines here to avoid psn = NaN

!
! VEGETATION
! Calculate canopy temperature, latent and sensible fluxes from the canopy,
! and leaf water change by evapotranspiration
!

  else

     call CanopyFluxes (clm%z0mv,   clm%z0hv,  clm%z0qv,  clm%thm,   clm%forc_th, &
                        clm%thv,    clm%tg,    clm%qg,    clm%dqgdT, clm%htvp,        &
                        clm%emv,    clm%emg,   clm%dlrad, clm%ulrad, clm%cgrnds,      &
                        clm%cgrndl, clm%cgrnd, clm    )

  endif

!
! Calculate total photosynthesis (only non-zero over vegetation)
! Put here so that history output can have global extent
!

  clm%fpsn = clm%psnsun*clm%laisun + clm%psnsha*clm%laisha

end subroutine Biogeophysics1
