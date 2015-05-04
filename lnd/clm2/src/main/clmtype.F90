#include <misc.h>
#include <preproc.h>

module clmtype

!=========================================================================
! Module for defining 1-D (vertical) CLM variable specification.
!=========================================================================

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : nlevsoi, nlevsno, nlevlak, numrad, nvoc, ndst
  implicit none
  public clm1d

  type clm1d

!*************************************************************************
! time invariant patch types and grid position
!*************************************************************************

     integer  :: nstep          !time step 
     integer  :: kpatch         !patch index  
     integer  :: itypveg        !vegetation type
     integer  :: itypwat        !water type
     integer  :: itypprc        !precipitation type
     integer  :: isoicol        !color classes for soil albedos
     logical  :: lakpoi         !true = lake point
     real(r8) :: lat            !latitude  for subgrid patches (radians)
     real(r8) :: lon            !longitude for subgrid patches (radians)
     real(r8) :: latdeg         !latitude  for subgrid patches (degrees)
     real(r8) :: londeg         !longitude for subgrid patches (degrees)
     real(r8) :: dtime          !time step interval

!*************************************************************************
! time variant level information
!*************************************************************************

     integer :: snl              ! number of snow layers

!*************************************************************************
! time-invariant boundary data for each of the subgrid patches
!*************************************************************************

! level values

     real(r8) :: zi(-nlevsno+0:nlevsoi) !interface level below a "z" level (m)
     real(r8) :: dz(-nlevsno+1:nlevsoi) !layer depth (m)
     real(r8) :: z (-nlevsno+1:nlevsoi) !layer thickness (m)

! soil physical properties

     real(r8) :: bsw   (nlevsoi)       !Clapp and Hornberger "b"
     real(r8) :: watsat(nlevsoi)       !volumetric soil water at saturation (porosity)
     real(r8) :: hksat (nlevsoi)       !hydraulic conductivity at saturation (mm H2O /s)
     real(r8) :: sucsat(nlevsoi)       !minimum soil suction (mm)
     real(r8) :: csol  (nlevsoi)       !heat capacity, soil solids (J/m**3/Kelvin)
     real(r8) :: tkmg  (nlevsoi)       !thermal conductivity, soil minerals  [W/m-K]  (new)
     real(r8) :: tkdry (nlevsoi)       !thermal conductivity, dry soil       (W/m/Kelvin)
     real(r8) :: tksatu(nlevsoi)       !thermal conductivity, saturated soil [W/m-K]  (new)
     real(r8) :: rootfr(nlevsoi)       !fraction of roots in each soil layer
     real(r8) :: rootr(nlevsoi)        !effective fraction of roots in each soil layer

! leaf constants

     real(r8) :: dewmx                 !Maximum allowed dew [mm]

! hydraulic constants of soil 

     real(r8) :: wtfact                !Fraction of model area with high water table
     real(r8) :: trsmx0                !Max transpiration for moist soil+100% veg. [mm/s]

! roughness lengths

     real(r8) :: zlnd                  !Roughness length for soil [m] (new)             
     real(r8) :: zsno                  !Roughness length for snow [m] (new)             
     real(r8) :: csoilc                !Drag coefficient for soil under canopy [-] (new)

! numerical finite-difference

     real(r8) :: cnfac                 !Crank Nicholson factor (between 0 and 1) (new)
     real(r8) :: capr                  !Tuning factor to turn first layer T into surface T (new)  
     real(r8) :: ssi                   !Irreducible water saturation of snow (new)
     real(r8) :: wimp                  !Water impremeable if porosity less than wimp (new)
     real(r8) :: pondmx                !Ponding depth (mm) (new)
     real(r8) :: smpmax                !wilting point potential in mm (new)
     real(r8) :: smpmin                !restriction for min of soil potential (mm) (new)

! water and energy balance check

     real(r8) :: begwb                 !water mass begining of the time step
     real(r8) :: endwb                 !water mass end of the time step
     real(r8) :: errh2o                !water conservation error (mm H2O)
     real(r8) :: errsoi                !soil/lake energy conservation error (W/m**2)
     real(r8) :: errseb                !surface energy conservation error (W/m**2)
     real(r8) :: errsol                !solar radiation conservation error (W/m**2)
     real(r8) :: errlon                !longwave radiation conservation error (W/m**2)
     real(r8) :: acc_errh2o            !accumulation of water balance error
     real(r8) :: acc_errseb            !accumulation of surface energy balance error
     real(r8) :: acc_errsoi            !accumulation of energy balance error
     real(r8) :: acc_errsol            !accumulation of energy balance error

!*************************************************************************
! subgrid patch version of atm model input
!*************************************************************************

     real(r8) :: forc_t              !atmospheric temperature (Kelvin)
     real(r8) :: forc_u              !atmospheric wind speed in east direction (m/s)
     real(r8) :: forc_v              !atmospheric wind speed in north direction (m/s)
     real(r8) :: forc_q              !atmospheric specific humidity (kg/kg)
     real(r8) :: forc_hgt            !atmospheric reference height (m) 
     real(r8) :: forc_hgt_u          !observational height of wind [m] (new)
     real(r8) :: forc_hgt_t          !observational height of temperature [m] (new)
     real(r8) :: forc_hgt_q          !observational height of humidity [m] (new)
     real(r8) :: forc_pbot           !atmospheric pressure (Pa)
     real(r8) :: forc_th             !atmospheric potential temperature (Kelvin)
     real(r8) :: forc_vp             !atmospheric vapor pressure (Pa)
     real(r8) :: forc_rho            !density (kg/m**3)
     real(r8) :: forc_co2            !atmospheric CO2 concentration (Pa)
     real(r8) :: forc_o2             !atmospheric O2 concentration (Pa)
     real(r8) :: forc_lwrad          !downward infrared (longwave) radiation (W/m**2)
     real(r8) :: forc_psrf           !surface pressure (Pa)
     real(r8) :: forc_solad(numrad)  !direct beam radiation (vis=forc_sols , nir=forc_soll )
     real(r8) :: forc_solai(numrad)  !diffuse radiation     (vis=forc_solsd, nir=forc_solld)
     real(r8) :: forc_rain           !rain rate [mm/s]
     real(r8) :: forc_snow           !snow rate [mm/s]

!*************************************************************************
! biogeophys
!*************************************************************************

! Surface solar radiation 

     real(r8) :: rssun          !sunlit stomatal resistance (s/m)
     real(r8) :: rssha          !shaded stomatal resistance (s/m)
     real(r8) :: psnsun         !sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
     real(r8) :: psnsha         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
     real(r8) :: laisun         !sunlit leaf area
     real(r8) :: laisha         !shaded leaf area
     real(r8) :: ndvi           !Normalized Difference Vegetation Index
     real(r8) :: sabg           !solar radiation absorbed by ground (W/m**2)
     real(r8) :: sabv           !solar radiation absorbed by vegetation (W/m**2)
     real(r8) :: fsa            !solar radiation absorbed (total) (W/m**2)
     real(r8) :: fsr            !solar radiation reflected (W/m**2)

! Surface energy fluxes

     real(r8) :: taux           !wind stress: e-w (kg/m/s**2)
     real(r8) :: tauy           !wind stress: n-s (kg/m/s**2)
     real(r8) :: eflx_lwrad_out !emitted infrared (longwave) radiation (W/m**2) 
     real(r8) :: eflx_lwrad_net !net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8) :: eflx_sh_tot    !total sensible heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_sh_veg    !sensible heat flux from leaves (W/m**2) [+ to atm]
     real(r8) :: eflx_sh_grnd   !sensible heat flux from ground (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_tot    !total latent heat flux (W/m8*2)  [+ to atm] 
     real(r8) :: eflx_lh_vege   !veg evaporation heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_vegt   !veg transpiration heat flux (W/m**2) [+ to atm]
     real(r8) :: eflx_lh_grnd   !ground evaporation heat flux (W/m**2) [+ to atm]   
     real(r8) :: eflx_soil_grnd !soil heat flux (W/m**2) [+ = into soil]
     real(r8) :: eflx_snomelt   !snow melt heat flux (W/m**2)
     real(r8) :: eflx_impsoil   !implicit evaporation for soil temperature equation (W/m**2)

! velocities

     real(r8) :: u10                       !10-m wind (m/s)
     real(r8) :: fv                        !friction velocity (m/s)
     real(r8) :: fm                        !used in u10 calculation

! Temperatures

     real(r8) :: t_veg                        !vegetation temperature (Kelvin)
     real(r8) :: t_grnd                       !ground temperature (Kelvin)
     real(r8) :: t_rad                        !radiative temperature (Kelvin)
     real(r8) :: t_ref2m                      !2 m height surface air temperature (Kelvin)
     real(r8) :: t_soisno(-nlevsno+1:nlevsoi) !soil temperature (Kelvin)
     real(r8) :: t_lake(1:nlevlak)            !lake temperature (Kelvin)
     real(r8) :: t_snow                       !vertically averaged snow temperature
     real(r8) :: dt_veg                       !change in t_veg, last iteration (Kelvin)
     real(r8) :: dt_grnd                      !change in t_grnd, last iteration (Kelvin)

! Soil properties

     real(r8) :: btran          !transpiration wetness factor (0 to 1) 

! Photosynthesis

     real(r8) :: fpsn           !photosynthesis (umol CO2 /m**2 /s)

!*************************************************************************
! hydrology
!*************************************************************************

     logical  :: do_capsnow                      !true => do snow capping  

     real(r8) :: qflx_infl                       !infiltration (mm H2O /s) 
     real(r8) :: qflx_surf                       !surface runoff (mm H2O /s) 
     real(r8) :: qflx_drain                      !sub-surface runoff (mm H2O /s) 
     real(r8) :: qflx_top_soil                   !net water input into soil from top (mm/s)
     real(r8) :: qflx_evap_soi                   !soil evaporation (mm H2O/s) (+ = to atm)
     real(r8) :: qflx_evap_veg                   !vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8) :: qflx_tran_veg                   !vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8) :: qflx_snomelt                    !snow melt (mm H2O /s)
     real(r8) :: qflx_evap_tot                   !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
     real(r8) :: qflx_prec_intr                  !interception of precipitation [mm/s]
     real(r8) :: qflx_prec_grnd                  !water onto ground including canopy runoff [kg/(m2 s)]
     real(r8) :: qflx_rain_grnd                  !rain on ground after interception (mm H2O/s) [+]
     real(r8) :: qflx_snow_grnd                  !snow on ground after interception (mm H2O/s) [+]
     real(r8) :: qflx_evap_grnd                  !ground surface evaporation rate (mm H2O/s) [+]
     real(r8) :: qflx_dew_grnd                   !ground surface dew formation (mm H2O /s) [+]
     real(r8) :: qflx_sub_snow                   !sublimation rate from snow pack (mm H2O /s) [+]
     real(r8) :: qflx_dew_snow                   !surface dew added to snow pack (mm H2O /s) [+]
     real(r8) :: qflx_snowcap                    !excess precipitation due to snow capping (mm H2O /s) [+]
     real(r8) :: qflx_qrgwl                      !qflx_surf at glaciers, wetlands, lakes
     real(r8) :: h2osno                          !snow water (mm H2O)
     real(r8) :: h2ocan                          !canopy water (mm H2O)
     real(r8) :: h2osoi_liq(-nlevsno+1:nlevsoi)  !liquid water (kg/m2) (new)
     real(r8) :: h2osoi_ice(-nlevsno+1:nlevsoi)  !ice lens (kg/m2) (new)
     real(r8) :: h2osoi_vol(nlevsoi)             !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
     real(r8) :: snowdp                          !snow height (m) 
     real(r8) :: snowage                         !non dimensional snow age [-] (new)
     real(r8) :: snowice                         !average snow ice lens
     real(r8) :: snowliq                         !average snow liquid water
     real(r8) :: h2osno_old                      !snow mass for previous time step (kg/m2) (new)
     integer  :: frac_veg_nosno                  !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
     integer  :: frac_veg_nosno_alb              !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
     real(r8) :: frac_sno                        !fraction of ground covered by snow (0 to 1) 
     real(r8) :: frac_iceold(-nlevsno+1:nlevsoi) !fraction of ice relative to the total water (new)
     real(r8) :: rsw                             !soil water content for root zone
     real(r8) :: eff_porosity(nlevsoi)           !effective porosity = porosity - vol_ice
     real(r8) :: sfact                           !term for implicit correction to evaporation
     real(r8) :: sfactmax                        !maximim of "sfact"

     integer  :: imelt(-nlevsno+1:nlevsoi)       !flag for melting (=1), freezing (=2), Not=0 (new)        

!*************************************************************************
! surfacealbedo (for next time step)
!*************************************************************************

     real(r8) :: parsun         !average absorbed PAR for sunlit leaves (W/m**2)
     real(r8) :: parsha         !average absorbed PAR for shaded leaves (W/m**2)
     real(r8) :: albd(numrad)   !surface albedo (direct)                     
     real(r8) :: albi(numrad)   !surface albedo (diffuse)                    
     real(r8) :: albgrd(numrad) !ground albedo (direct)                      
     real(r8) :: albgri(numrad) !ground albedo (diffuse)                     
     real(r8) :: fabd(numrad)   !flux absorbed by veg per unit direct flux   
     real(r8) :: fabi(numrad)   !flux absorbed by veg per unit diffuse flux  
     real(r8) :: ftdd(numrad)   !down direct flux below veg per unit dir flx 
     real(r8) :: ftid(numrad)   !down diffuse flux below veg per unit dir flx
     real(r8) :: ftii(numrad)   !down diffuse flux below veg per unit dif flx
     real(r8) :: fsun           !sunlit fraction of canopy                   

!*************************************************************************
! ecosysdynamics
!*************************************************************************

     real(r8) :: hbot           !canopy bottom (m)
     real(r8) :: htop           !canopy top (m)
     real(r8) :: tlai           !one-sided leaf area index, no burying by snow
     real(r8) :: tsai           !one-sided stem area index, no burying by snow
     real(r8) :: elai           !one-sided leaf area index with burying by snow
     real(r8) :: esai           !one-sided stem area index with burying by snow
     real(r8) :: fwet           !fraction of canopy that is wet (0 to 1)
     real(r8) :: fdry           !fraction of foliage that is green and dry [-] (new)

!*************************************************************************
! terms from pft_varcon - to avoid indirect indexing
!*************************************************************************

     real(r8) :: z0mr           ! ratio of momentum roughness length to canopy top height [-]
     real(r8) :: z0m            ! momentum roughness length [m]
     real(r8) :: displar        ! ratio of displacement height to canopy top height [-]
     real(r8) :: displa         ! displacement height [m]
     real(r8) :: dleaf          ! leaf dimension [m]
     real(r8) :: xl             ! pft_varcon leaf/stem orientation index
     real(r8) :: rhol(numrad)   ! pft_varcon leaf reflectance  : 1=vis, 2=nir 
     real(r8) :: rhos(numrad)   ! pft_varcon stem reflectance  : 1=vis, 2=nir 
     real(r8) :: taul(numrad)   ! pft_varcon leaf transmittance: 1=vis, 2=nir 
     real(r8) :: taus(numrad)   ! pft_varcon stem transmittance: 1=vis, 2=nir 
     real(r8) :: qe25           ! quantum efficiency at 25c (umol co2 / umol photon)
     real(r8) :: vcmx25         ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     real(r8) :: mp             ! slope for conductance-to-photosynthesis relationship
     real(r8) :: c3psn          ! photosynthetic pathway: 0. = c4, 1. = c3

!*************************************************************************
! terms due to splitting the code into Biogeophys1 and Biogeophys2
!*************************************************************************

     real(r8) cgrnd  ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
     real(r8) cgrndl ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
     real(r8) cgrnds ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
     real(r8) tg     ! ground surface temperature [K]
     real(r8) tssbef(-nlevsno:nlevsoi)  ! soil/snow temperature before update
     real(r8) qg     ! ground specific humidity [kg/kg]
     real(r8) dqgdT  ! d(qg)/dT
     real(r8) emg    ! ground emissivity
     real(r8) emv    ! vegetation emissivity
     real(r8) htvp   ! latent heat of vapor of water (or sublimation) [j/kg]
     real(r8) z0mg   ! roughness length over ground, momentum [m]
     real(r8) z0hg   ! roughness length over ground, sensible heat [m]
     real(r8) z0qg   ! roughness length over ground, latent heat [m]
     real(r8) z0mv   ! roughness length over vegetation, momentum [m]
     real(r8) z0hv   ! roughness length over vegetation, sensible heat [m]
     real(r8) z0qv   ! roughness length over vegetation, latent heat [m]
     real(r8) beta   ! coefficient of convective velocity [-]
     real(r8) zii    ! convective boundary height [m]
     real(r8) thm    ! intermediate variable (forc_t+0.0098*forc_hgt_t)
     real(r8) thv    ! virtual potential temperature (kelvin)
     real(r8) ur     ! wind speed at reference height [m/s] ***DO WE NEED THIS???
     real(r8) dlrad  ! downward longwave radiation below the canopy [W/m2]
     real(r8) ulrad  ! upward longwave radiation above the canopy [W/m2]
     real(r8) qmelt  ! snow melt [mm/s]

! -----------------------------------------------------------------




  end type clm1d

  SAVE

end module clmtype



