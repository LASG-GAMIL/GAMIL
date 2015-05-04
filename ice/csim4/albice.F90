#include <misc.h>
#include <params.h>

subroutine albice(lchnk   ,ncol    ,Tair,  snowh   ,coszrs  ,&
                  asdir   ,aldir   ,asdif   ,aldif   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute surface albedos
!
! Method: 
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: albice.F90,v 1.1.4.3 2002/06/15 13:50:09 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use comsrf, only :icefrac,sicthk
  use ice_constants
  implicit none
#include <albedo.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: lchnk            ! chunk identifier
  integer , intent(in) :: ncol             ! number of atmospheric columns

  real(r8), intent(in) :: Tair(pcols)      ! bottom level air temp
  real(r8), intent(in) :: snowh(pcols)     ! Snow depth (liquid water equivalent)
  real(r8), intent(in) :: coszrs(pcols)    ! Cosine solar zenith angle
  real(r8), intent(out):: asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r8), intent(out):: aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r8), intent(out):: asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r8), intent(out):: aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
      ! albedos for ice in each category
  real(r8) :: alvdrn (pcols) ! visible, direct   (fraction)
  real(r8) :: alidrn (pcols) ! near-ir, direct   (fraction)
  real(r8) :: alvdfn (pcols) ! visible, diffuse  (fraction)
  real(r8) :: alidfn (pcols) ! near-ir, diffuse  (fraction)
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! Longitude index
!-----------------------------------------------------------------------
      real (r8), parameter :: albocn = 0.06_dbl_kind  ! ocean albedo
      real (r8), parameter :: &
           ahmax    = 0.5_dbl_kind,   &! thickns above which ice alb is const,
           albicev  = 0.78_dbl_kind,  &! visible ice albedo for h > ahmax
           albicei  = 0.36_dbl_kind,  &! near-ir ice albedo for h > ahmax
           albsnowv = 0.98_dbl_kind,  &! cold snow albedo, visible
           albsnowi = 0.70_dbl_kind,  &! cold snow albedo, near IR
           dT_mlt   = 1._dbl_kind,    &! change in temp to give dalb_mlt change
           dalb_mlt = -0.075_dbl_kind,&! albedo change per dT_mlt change
           dalb_mltv= -0.100_dbl_kind,&! albedo vis change per dT_mlt change in temp for snow
           dalb_mlti= -0.150_dbl_kind,&! albedo nir change per dT_mlt change in temp for snow
           Tf       = -1.8_dbl_kind    ! hardwired

      ! parameter for fractional snow area 
      real(r8)  fhtan ! factor used in dependence of albedo on ice thickness
      real(r8)  vicen(pcols),vsnon(pcols),aicen(pcols),tsfcn(pcols)
      real (r8) hi ! ice thickness  (m)
      real (r8) hs ! snow thickness (m)
      real (r8) snw !
      real (r8) albo  ! effective ocean albedo, function of ice thickness
      real (r8) asnow ! snow-covered area fraction
      real (r8) asnwv ! snow albedo, visible 
      real (r8) asnwi ! snow albedo, near IR
      real (r8) fh ! piecewise linear function of thickness 
      real (r8) fT ! piecewise linear function of surface temperature
      real (r8) dTs ! difference of Tsfc and Timelt
!-----------------------------------------------------------------------
!
! Initialize all sea ice surface albedos to zero
!
      asdir(:) = 0.
      aldir(:) = 0.
      asdif(:) = 0.
      aldif(:) = 0.
      alvdrn(:) = 0.
      alidrn(:) = 0.
      alvdfn(:) = 0.
      alidfn(:) = 0.

  fhtan = atan(ahmax*5._dbl_kind) 

  do i=1,ncol
     if (icefrac(i,lchnk) > 0._r8 .and. coszrs(i)>0.0) then
        hi  = sicthk(i,lchnk)
        snw = snowh(i)*rhofresh/rhos
        aicen(i) = icefrac(i,lchnk)
        vicen(i) = hi*aicen(i) 
        !---------------------------------------------------------
        ! keep snow/ice boundary above sea level by reducing snow
        !---------------------------------------------------------
        vsnon(i) = min(snw*aicen(i),p33*vicen(i))
        Tsfcn(i) = min(Tair(i)-Tffresh,-p2)   ! deg C       
        !---------------------------------------------------------
        ! make linear temp profile and compute enthalpy
        !---------------------------------------------------------

        hi = vicen(i) / aicen(i)
        hs = vsnon(i) / aicen(i)
        
        ! bare ice, thickness dependence
        fh = min(atan(hi*5.)/fhtan,c1)
        albo = albocn*(c1-fh)
        alvdfn(i) = albicev*fh + albo
        alidfn(i) = albicei*fh + albo

        ! bare ice, temperature dependence
        dTs = Timelt - Tsfcn(i)
        fT = min(dTs/dT_mlt-c1,c0)
        alvdfn(i) = alvdfn(i) - dalb_mlt*fT
        alidfn(i) = alidfn(i) - dalb_mlt*fT

        if( hs .gt. 0._r8 ) then
           ! fractional area of snow on ice (thickness dependent)
           asnow = hs / ( hs + snowpatch ) 
           asnwv = albsnowv
           asnwi = albsnowi
           ! snow on ice, temperature dependence
           asnwv = asnwv - dalb_mltv*fT
           asnwi = asnwi - dalb_mlti*fT
           
           ! combine ice and snow albedos
           alvdfn(i) = alvdfn(i)*(c1-asnow) + asnwv*asnow
           alidfn(i) = alidfn(i)*(c1-asnow) + asnwi*asnow
        endif
        alvdrn(i) = alvdfn(i)
        alidrn(i) = alidfn(i)
     endif  ! aicen > 0._r8

     if (icefrac(i,lchnk) > 0._r8 .and. coszrs(i)>0.0) then
        asdir(i)  = alvdrn(i)
        aldir(i)  = alidrn(i)
        asdif(i) = alvdfn(i)
        aldif(i) = alidfn(i)
     end if
  enddo
!
  return
end subroutine albice

