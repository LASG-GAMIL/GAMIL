#include <misc.h>
#include <preproc.h>

subroutine Hydrology2 (clm)

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
! This is the main subroutine to execute the calculation of soil/snow 
! hydrology
!
! Method:
! Calling sequence is: 
!  Hydrology2:                       surface hydrology driver
!    -> SnowWater:                   change of snow mass and snow water
!                                    onto soil
!    -> SurfaceRunoff:               surface runoff 
!    -> Infiltration:                infiltration into surface soil layer (minus the
!                                    evaporation)
!    -> SoilWater:                   soil water movement between layers
!          -> Tridiagonal            tridiagonal matrix solution
!    -> Drainage:                    subsurface runoff  
!    -> SnowCompaction:              compaction of snow layers
!    -> CombineSnowLayers:           combine snow layers that are thinner
!                                     than minimum
!    -> DivideSnowLayers:            subdivide snow layers that are thicker
!                                     than maximum
!    -> WetIceHydrology:             calculate hydrology for wetland and land
!                                     ice 
! Author:
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: Hydrology2.F90,v 1.4.6.4.6.1 2002/10/03 20:07:28 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : denh2o, denice, istice, istwet, istsoil, spval
  use clm_varpar, only : nlevsoi, nlevsno
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer j                   ! do loop index
  real(r8) zwice              ! the sum of ice mass of soil [kg/m2]
  real(r8) vol_liq(1:nlevsoi) ! partial volume of liquid water in layer [-]
  real(r8) s(1:nlevsoi)       ! wetness of soil (including ice) [-]
  real(r8) zwt                ! water table depth [m]
  real(r8) fcov               ! fractional area with water table at surface [-]
  real(r8) dwat(1:nlevsoi)    ! change in soil water [-]
  real(r8) hk(1:nlevsoi)      ! hydraulic conductivity [mm/s]
  real(r8) dhkdw(1:nlevsoi)   ! d(hk)/d(vol_liq)

!----End Variable List--------------------------------------------------

!
! Determine the change of snow mass and the snow water onto soil
!

     call SnowWater (clm)

!
! Determine soil hydrology
!

     if (clm%itypwat == istsoil) then
      call SurfaceRunoff  (clm, zwice, vol_liq, s, zwt, &
                           fcov)
      call Infiltration   (clm)
      call SoilWater      (clm, vol_liq, dwat, hk, dhkdw)
      call Drainage       (clm,  zwice, vol_liq, s,   zwt, &
                           fcov, hk,    dhkdw,   dwat )
     endif

     if (clm%snl < 0) then

!
! Natural compaction and metamorphosis.
!

        call SnowCompaction (clm)

!
! Combine thin snow elements
!

        call CombineSnowLayers (clm)

!
! Divide thick snow elements
!

        call DivideSnowLayers (clm)

!
! Set empty snow layers to zero
!

        if (clm%snl > -nlevsno) then
           clm%snowage = 0.
           do j = -nlevsno+1, clm%snl
              clm%h2osoi_ice(j) = 0.
              clm%h2osoi_liq(j) = 0.
              clm%t_soisno(j) = 0.
              clm%dz(j) = 0.
              clm%z(j) = 0.
              clm%zi(j-1) = 0.
           enddo
        endif

     endif

!
! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice 
! over all snow layers for the given patch - these will be written out
! to the history tape
!
        
     if (clm%snl < 0) then
        clm%t_snow  = 0.
        clm%snowice = 0.
        clm%snowliq = 0.
        do j = clm%snl+1, 0
           clm%t_snow  = clm%t_snow  + clm%t_soisno(j)
           clm%snowice = clm%snowice + clm%h2osoi_ice(j)
           clm%snowliq = clm%snowliq + clm%h2osoi_liq(j)
        end do
        clm%t_snow = clm%t_snow/abs(clm%snl)
     else
        clm%t_snow  = spval
        clm%snowice = spval
        clm%snowliq = spval
     endif

!
! Update ground temperature
!

     clm%t_grnd = clm%t_soisno(clm%snl+1)

!
! Determine volumetric soil water
!

     do j = 1,nlevsoi
        clm%h2osoi_vol(j) = clm%h2osoi_liq(j)/(clm%dz(j)*denh2o) &
                          + clm%h2osoi_ice(j)/(clm%dz(j)*denice)
     end do

!
! Determine ending water balance
!

     clm%endwb=clm%h2ocan+clm%h2osno
     do j = 1, nlevsoi
        clm%endwb = clm%endwb + clm%h2osoi_ice(j) + clm%h2osoi_liq(j)
     enddo

!
! Determine wetland and land ice hydrology
!

     if (clm%itypwat==istwet .or. clm%itypwat==istice) call WetIceHydrology (clm)

end subroutine Hydrology2




