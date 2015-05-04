#include <misc.h>
#include <preproc.h>

module clm_varder

! Declare clm variable

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  implicit none
  
  type (clm1d), allocatable :: clm(:)     
  SAVE
  
!=======================================================================
CONTAINS
!=======================================================================

  subroutine clm_varder_ini

    use shr_kind_mod, only: r8 => shr_kind_r8
    use infnan
    use clm_varmap, only : begpatch, endpatch
    use clm_varcon, only : spval
    implicit none

    integer :: k

! allocate memory for clm derived type

    allocate (clm(begpatch:endpatch))

! set all elements to infinity

    do k = begpatch, endpatch

!*************************************************************************
! Time-invariant boundary data for each of the subgrid patches
!*************************************************************************

       clm(k)%kpatch  = bigint        
       clm(k)%itypveg = bigint        
       clm(k)%itypwat = bigint        
       clm(k)%isoicol = bigint        

! level values

       clm(k)%dz(-nlevsno+1:0) = inf       !snow layer thickness (m)
       clm(k)%z (-nlevsno+1:0) = inf       !snow layer depth (m)
       clm(k)%zi(-nlevsno+0:0) = inf       !snow layer interfaces (m)  
       clm(k)%dz(1:nlevsoi)    = inf       !soil layer thickness (m)
       clm(k)%z (1:nlevsoi)    = inf       !soil layer depth (m)
       clm(k)%zi(1:nlevsoi)    = inf       !soil layer interfaces (m)  

! soil physical properties

       clm(k)%bsw   (1:nlevsoi) = inf      !Clapp and Hornberger "b"
       clm(k)%watsat(1:nlevsoi) = inf      !volumetric soil water at saturation (porosity)
       clm(k)%hksat (1:nlevsoi) = inf      !hydraulic conductivity at saturation (mm H2O /s)
       clm(k)%sucsat(1:nlevsoi) = inf      !minimum soil suction (mm)
       clm(k)%csol  (1:nlevsoi) = inf      !heat capacity, soil solids (J/m**3/Kelvin)
       clm(k)%tkmg  (1:nlevsoi) = inf      !thermal conductivity, soil minerals  [W/m-K]  (new)
       clm(k)%tkdry (1:nlevsoi) = inf      !thermal conductivity, dry soil       (W/m/Kelvin)
       clm(k)%tksatu(1:nlevsoi) = inf      !thermal conductivity, saturated soil [W/m-K]  (new)
       clm(k)%rootfr(1:nlevsoi) = inf      !fraction of roots in each soil layer
       clm(k)%rootr (1:nlevsoi) = inf      !effective fraction of roots in each soil layer

! leaf constants

       clm(k)%dewmx = inf            !Maximum allowed dew [mm]

! hydraulic constants of soil 	     

       clm(k)%wtfact = inf           !Fraction of model area with high water table
       clm(k)%trsmx0 = inf           !Max transpiration for moist soil+100% veg. [mm/s]

! roughness lengths		     

       clm(k)%zlnd   = inf           !Roughness length for soil [m] (new)             
       clm(k)%zsno   = inf           !Roughness length for snow [m] (new)             
       clm(k)%csoilc = inf           !Drag coefficient for soil under canopy [-] (new)

! numerical finite-difference

       clm(k)%cnfac   = inf          !Crank Nicholson factor (between 0 and 1) (new)
       clm(k)%capr    = inf          !Tuning factor to turn first layer T into surface T (new)  
       clm(k)%ssi     = inf          !Irreducible water saturation of snow (new)
       clm(k)%wimp    = inf          !Water impremeable if porosity less than wimp (new)
       clm(k)%pondmx  = inf          !Ponding depth (mm) (new)
       clm(k)%smpmax  = inf          !wilting point potential in mm (new)
       clm(k)%smpmin  = inf          !restriction for min of soil potential (mm) (new)

! water and energy balance check

       clm(k)%begwb      = inf       !water mass begining of the time step
       clm(k)%endwb      = inf       !water mass end of the time step
       clm(k)%errh2o     = inf       !water conservation error (mm H2O)
       clm(k)%errsoi     = inf       !soil/lake energy conservation error (W/m**2)
       clm(k)%errseb     = inf       !surface energy conservation error (W/m**2)
       clm(k)%errsol     = inf       !solar radiation conservation error (W/m**2)
       clm(k)%errlon     = inf       !longwave radiation conservation error (W/m**2)
       clm(k)%acc_errseb = 0.        !accumulation of surface energy balance error
       clm(k)%acc_errh2o = 0.        !accumulation of water balance error

!*************************************************************************
! subgrid patch version of atm model input
!*************************************************************************

       clm(k)%forc_t     = inf         !atmospheric temperature (Kelvin)
       clm(k)%forc_u     = inf         !atmospheric wind speed in east direction (m/s)
       clm(k)%forc_v     = inf         !atmospheric wind speed in north direction (m/s)
       clm(k)%forc_q     = inf         !atmospheric specific humidity (kg/kg)
       clm(k)%forc_hgt   = inf         !atmospheric reference height (m) 
       clm(k)%forc_hgt_u = inf         !observational height of wind [m] (new)
       clm(k)%forc_hgt_t = inf         !observational height of temperature [m] (new)
       clm(k)%forc_hgt_q = inf         !observational height of humidity [m] (new)
       clm(k)%forc_pbot  = inf         !atmospheric pressure (Pa)
       clm(k)%forc_th    = inf         !atmospheric potential temperature (Kelvin)
       clm(k)%forc_vp    = inf         !atmospheric vapor pressure (Pa)
       clm(k)%forc_rho   = inf         !density (kg/m**3)
       clm(k)%forc_co2   = inf         !atmospheric CO2 concentration (Pa)
       clm(k)%forc_o2    = inf         !atmospheric O2 concentration (Pa)
       clm(k)%forc_lwrad = inf         !downward infrared (longwave) radiation (W/m**2)
       clm(k)%forc_psrf  = inf         !surface pressure (Pa)
       clm(k)%forc_solad(1:numrad) = inf !direct beam radiation (vis=forc_sols , nir=forc_soll )
       clm(k)%forc_solai(1:numrad) = inf !diffuse radiation     (vis=forc_solsd, nir=forc_solld)
       clm(k)%forc_rain  = inf         !rain rate [mm/s]
       clm(k)%forc_snow  = inf         !snow rate [mm/s]

!*************************************************************************
! biogeophys
!*************************************************************************

! Surface solar radiation 

       clm(k)%rssun  = inf        !sunlit stomatal resistance (s/m)
       clm(k)%rssha  = inf        !shaded stomatal resistance (s/m)
       clm(k)%psnsun = inf        !sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
       clm(k)%psnsha = inf        !shaded leaf photosynthesis (umol CO2 /m**2/ s)
       clm(k)%laisun = inf        !sunlit leaf area
       clm(k)%laisha = inf        !shaded leaf area
       clm(k)%ndvi   = inf        !Normalized Difference Vegetation Index
       clm(k)%sabg   = inf        !solar radiation absorbed by ground (W/m**2)
       clm(k)%sabv   = inf        !solar radiation absorbed by vegetation (W/m**2)
       clm(k)%fsa    = inf        !solar radiation absorbed (total) (W/m**2)
       clm(k)%fsr    = inf        !solar radiation reflected (W/m**2)

! Surface energy fluxes

       clm(k)%taux           = inf !wind stress: e-w (kg/m/s**2)
       clm(k)%tauy           = inf !wind stress: n-s (kg/m/s**2)
       clm(k)%eflx_lwrad_out = inf !emitted infrared (longwave) radiation (W/m**2) 
       clm(k)%eflx_lwrad_net = inf !net infrared (longwave) rad (W/m**2) [+  = to atm]
       clm(k)%eflx_sh_tot    = inf !total sensible heat flux (W/m**2) [+ to atm]
       clm(k)%eflx_sh_veg    = inf !sensible heat flux from leaves (W/m**2) [+ to atm]
       clm(k)%eflx_sh_grnd   = inf !sensible heat flux from ground (W/m**2) [+ to atm]
       clm(k)%eflx_lh_tot    = inf !total latent heat flux (W/m8*2)  [+ to atm] 
       clm(k)%eflx_lh_vege   = inf !veg evaporation heat flux (W/m**2) [+ to atm]
       clm(k)%eflx_lh_vegt   = inf !veg transpiration heat flux (W/m**2) [+ to atm]
       clm(k)%eflx_lh_grnd   = inf !ground evaporation heat flux (W/m**2) [+ to atm]   
       clm(k)%eflx_soil_grnd = inf !soil heat flux (W/m**2) [+  = into soil]
       clm(k)%eflx_snomelt   = inf !snow melt heat flux (W/m**2)

! Velocities

       clm(k)%u10 = inf            !10-m wind (m/s)
       clm(k)%fv  = inf            !friction velocity (m/s)
       clm(k)%fm  = inf            !used in u10 calculation

! Temperatures

       clm(k)%t_veg   = inf                     !vegetation temperature (Kelvin)
       clm(k)%t_grnd  = inf                     !ground temperature (Kelvin)
       clm(k)%t_rad   = inf                     !radiative temperature (Kelvin)
       clm(k)%t_ref2m = inf                     !2 m height surface air temperature (Kelvin)
       clm(k)%t_soisno(-nlevsno+1:0) = inf      !snow temperature (Kelvin)
       clm(k)%t_soisno(1:nlevsoi)    = inf      !soil temperature (Kelvin)
       clm(k)%t_lake(1:nlevlak)      = inf      !lak temperature (Kelvin)
       clm(k)%dt_veg  = spval                   !change in t_veg, last iteration (Kelvin)
       clm(k)%dt_grnd = spval                   !change in t_grnd, last iteration (Kelvin)

! Soil properties   

       clm(k)%btran = inf                       !transpiration wetness factor (0 to 1) 

! Photosynthesis

       clm(k)%fpsn = inf                        !photosynthesis (umol CO2 /m**2 /s)


!*************************************************************************
! Hydrology
!*************************************************************************

       clm(k)%qflx_infl       = inf                !Infiltration (mm H2O /s) 
       clm(k)%qflx_surf       = inf                !surface runoff (mm H2O /s) 
       clm(k)%qflx_drain      = inf                !sub-surface runoff (mm H2O /s) 
       clm(k)%qflx_top_soil   = inf                !net water input into soil from top (mm/s)
       clm(k)%qflx_evap_soi   = inf                !soil evaporation (mm H2O/s) (+ = to atm)
       clm(k)%qflx_evap_veg   = inf                !vegetation evaporation (mm H2O/s) (+ = to atm)
       clm(k)%qflx_tran_veg   = inf                !vegetation transpiration (mm H2O/s) (+ = to atm)
       clm(k)%qflx_snomelt    = inf                !snow melt (mm H2O /s)
       clm(k)%qflx_evap_tot   = inf                !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
       clm(k)%qflx_prec_intr  = inf                !interception of precipitation [mm/s]
       clm(k)%qflx_prec_grnd  = inf                !water onto ground including canopy runoff [kg/(m2 s)]
       clm(k)%qflx_rain_grnd  = inf                !rain on ground after interception (mm H2O/s) [+]
       clm(k)%qflx_snow_grnd  = inf                !snow on ground after interception (mm H2O/s) [+]
       clm(k)%qflx_evap_grnd  = inf                !ground surface evaporation rate (mm H2O/s) [+]
       clm(k)%qflx_dew_grnd   = inf                !ground surface dew formation (mm H2O /s) [+]
       clm(k)%qflx_sub_snow   = inf                !sublimation rate from snow pack (mm H2O /s) [+]
       clm(k)%qflx_dew_snow   = inf                !surface dew added to snow pack (mm H2O /s) [+]
       clm(k)%qflx_snowcap    = inf                !excess precipitation due to snow capping (mm H2O /s) [+]
       clm(k)%qflx_qrgwl      = 0                  !qflx_surf at glaciers, wetlands, lakes
       clm(k)%h2osno          = inf                !snow water (mm H2O / m**2)
       clm(k)%h2ocan          = inf                !canopy water (mm H2O / m**2)
       clm(k)%h2osoi_liq(-nlevsno+1:0) = inf       !snow liquid water (kg/m2) (new)
       clm(k)%h2osoi_ice(-nlevsno+1:0) = inf       !snow ice lens (kg/m2) (new)
       clm(k)%h2osoi_liq(1:nlevsoi) = inf          !soil liquid water (kg/m2) (new)
       clm(k)%h2osoi_ice(1:nlevsoi) = inf          !soil ice lens (kg/m2) (new)
       clm(k)%h2osoi_vol(1:nlevsoi) = inf          !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
       clm(k)%snowdp          = inf                !snow height (m) 
       clm(k)%snowage         = inf                !non dimensional snow age [-] (new)
       clm(k)%t_snow          = inf                !average snow temperature
       clm(k)%snowice         = inf                !average snow ice lens
       clm(k)%snowliq         = inf                !average snow liquid water
       clm(k)%h2osno_old      = inf                !snow mass for previous time step (kg/m2) (new)
       clm(k)%frac_veg_nosno  = bigint             !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
       clm(k)%frac_veg_nosno_alb = bigint          !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
       clm(k)%frac_sno        = inf                !fraction of ground covered by snow (0 to 1) 
       clm(k)%frac_iceold     = inf                !fraction of ice relative to the total water (new)
       clm(k)%rsw             = inf                !soil water content for root zone
       clm(k)%eff_porosity    = inf                !effective porosity
       clm(k)%sfact           = inf                !term for implicit correction to evaporation
       clm(k)%sfactmax        = inf                !maximim of "sfact"

       clm(k)%imelt(-nlevsno+1:nlevsoi) = bigint   !flag for melting (=1), freezing (=2), Not=0 (new)        

!*************************************************************************
! Surfacealbedo (for next time step)
!*************************************************************************

       clm(k)%parsun          = inf !average absorbed PAR for sunlit leaves (W/m**2)
       clm(k)%parsha          = inf !average absorbed PAR for shaded leaves (W/m**2)
       clm(k)%albd(1:numrad)  = inf !surface albedo (direct)                     
       clm(k)%albi(1:numrad)  = inf !surface albedo (diffuse)                    
       clm(k)%albgrd(1:numrad)= inf !ground albedo (direct)                      
       clm(k)%albgri(1:numrad)= inf !ground albedo (diffuse)                     
       clm(k)%fabd(1:numrad)  = inf !flux absorbed by veg per unit direct flux   
       clm(k)%fabi(1:numrad)  = inf !flux absorbed by veg per unit diffuse flux  
       clm(k)%ftdd(1:numrad)  = inf !down direct flux below veg per unit dir flx 
       clm(k)%ftid(1:numrad)  = inf !down diffuse flux below veg per unit dir flx
       clm(k)%ftii(1:numrad)  = inf !down diffuse flux below veg per unit dif flx
       clm(k)%fsun            = inf !sunlit fraction of canopy                   

!*************************************************************************
! Ecosysdynamics
!*************************************************************************

       clm(k)%displa        = inf !displacement height [m]
       clm(k)%z0m           = inf !roughness length, momentum [m]
       clm(k)%tlai          = inf !one-sided leaf area index, no burying by snow
       clm(k)%tsai          = inf !one-sided stem area index, no burying by snow
       clm(k)%elai          = inf !one-sided leaf area index with burying by snow
       clm(k)%esai          = inf !one-sided stem area index with burying by snow
       clm(k)%fwet          = inf !fraction of canopy that is wet (0 to 1)
       clm(k)%fdry          = inf !fraction of foliage that is green and dry [-] (new)
       clm(k)%hbot          = inf !canopy bottom height [m]
       clm(k)%htop          = inf !canopy top height [m]

!*************************************************************************
! Terms due to splitting the code into Biogeophys1 and Biogeophys2
!*************************************************************************

       clm(k)%cgrnd  = inf ! deriv of soil energy flux wrt to soil temp [w/m2/k]
       clm(k)%cgrndl = inf ! deriv of soil sensible heat flux wrt soil temp [w/m2/k]
       clm(k)%cgrnds = inf ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       clm(k)%tg     = inf ! ground surface temperature [K]
       clm(k)%tssbef(-nlevsno:nlevsoi) = inf ! soil/snow temperature before update
       clm(k)%qg     = inf ! ground specific humidity [kg/kg]
       clm(k)%dqgdT  = inf ! d(qg)/dT
       clm(k)%emg    = inf ! ground emissivity
       clm(k)%emv    = inf ! vegetation emissivity
       clm(k)%htvp   = inf ! latent heat of vapor of water (or sublimation) [j/kg]
       clm(k)%z0mg   = inf ! roughness length over ground, momentum [m]
       clm(k)%z0hg   = inf ! roughness length over ground, sensible heat [m]
       clm(k)%z0qg   = inf ! roughness length over ground, latent heat [m]
       clm(k)%z0mv   = inf ! roughness length over vegetation, momentum [m]
       clm(k)%z0hv   = inf ! roughness length over vegetation, sensible heat [m]
       clm(k)%z0qv   = inf ! roughness length over vegetation, latent heat [m]
       clm(k)%beta   = inf ! coefficient of convective velocity [-]
       clm(k)%zii    = inf ! convective boundary height [m]
       clm(k)%thm    = inf ! intermediate variable (forc_t+0.0098*forc_hgt_t)
       clm(k)%thv    = inf ! virtual potential temperature (kelvin)
       clm(k)%ur     = inf ! wind speed at reference height [m/s]
       clm(k)%dlrad  = inf ! downward longwave radiation below the canopy [W/m2]
       clm(k)%ulrad  = inf ! upward longwave radiation above the canopy [W/m2]
       clm(k)%qmelt  = inf ! snow melt [mm/s]

    end do  ! end of patch loop

  end subroutine clm_varder_ini

end module clm_varder


