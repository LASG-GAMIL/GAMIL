#include <misc.h>
#include <preproc.h>

subroutine EcosystemDyn (clm, doalb, endofyr)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Ecosystem dynamics: phenology, vegetation
! 
! Method: 
! Calling sequence:
!    EcosystemDynamics:             !ecosystem dynamics driver
! This subroutine calculates:
!    o leaf areas (tlai, elai)
!    o stem areas (tsai, esai)
!    o height (htop)
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: EcosystemDyn.F90,v 1.6.6.6.6.1 2002/10/03 20:07:29 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use mvegFileMod , only : mlai1, msai1, mlai2, msai2, &
                           mhvt1, mhvt2, mhvb1, mhvb2, timwt 
  use time_manager, only : get_curr_date
  implicit none

! ------------------------ arguments ------------------------------
  type (clm1d), intent(inout) :: clm   !CLM 1-D Module
  logical     , intent(in)    :: doalb !true = surface albedo calculation time step
  logical     , intent(in)    :: endofyr
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
  real(r8) :: ol      !thickness of canopy layer covered by snow (m)
  real(r8) :: fb      !fraction of canopy layer covered by snow
! -----------------------------------------------------------------

  if (doalb) then

! need to update elai and esai only every albedo time step so do not 
! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e., 
! if albedos are not done every time step).

! leaf phenology
! Set leaf and stem areas based on day of year
! Interpolate leaf area index, stem area index, and vegetation heights
! between two monthly 
! The weights below (timwt(1) and timwt(2)) were obtained by a call to 
! routine InterpMonthlyVeg in subroutine NCARlsm. 
!                 Field   Monthly Values
!                -------------------------
! leaf area index LAI  <- mlai1 and mlai2
! leaf area index SAI  <- msai1 and msai2
! top height      HTOP <- mhvt1 and mhvt2
! bottom height   HBOT <- mhvb1 and mhvb2

     clm%tlai = timwt(1)*mlai1(clm%kpatch) + timwt(2)*mlai2(clm%kpatch)
     clm%tsai = timwt(1)*msai1(clm%kpatch) + timwt(2)*msai2(clm%kpatch)
     clm%htop = timwt(1)*mhvt1(clm%kpatch) + timwt(2)*mhvt2(clm%kpatch)
     clm%hbot = timwt(1)*mhvb1(clm%kpatch) + timwt(2)*mhvb2(clm%kpatch)

! adjust lai and sai for burying by snow. if exposed lai and sai
! are less than 0.05, set equal to zero to prevent numerical 
! problems associated with very small lai and sai.

     ol = min( max(clm%snowdp-clm%hbot,0._r8), clm%htop-clm%hbot)
     fb = 1. - ol / max(1.e-06_r8, clm%htop-clm%hbot)
     clm%elai = clm%tlai*fb
     clm%esai = clm%tsai*fb
     if (clm%elai < 0.05) clm%elai = 0._r8
     if (clm%esai < 0.05) clm%esai = 0._r8

! Fraction of vegetation free of snow

     if ((clm%elai + clm%esai) >= 0.05) then
        clm%frac_veg_nosno_alb = 1
     else
        clm%frac_veg_nosno_alb = 0
     endif
  
  endif

  return
end subroutine EcosystemDyn
