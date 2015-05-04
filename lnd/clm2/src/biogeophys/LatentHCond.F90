#include <misc.h>
#include <preproc.h>

subroutine LatentHCond (raw,   rbw,   rdw,   rpp,   wtaq,  &
                        wtlq,  wtgq,  wtaq0, wtlq0, wtgq0, &
                        wtalq, wtgaq, wtglq, clm    )

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
! Provides dimensional and non-dimensional latent heat 
! conductances for canopy and soil flux calculations.  Latent flux
! differs from the sensible heat flux due to stomatal resistance.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: LatentHCond.F90,v 1.1.10.3 2002/06/15 13:50:16 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

  real(r8), intent(in) :: raw    ! aerodynamic resistance [s/m]
  real(r8), intent(in) :: rbw    ! leaf boundary layer resistance [s/m]
  real(r8), intent(in) :: rdw    ! latent heat resistance between ground and bottom of canopy [s/m]
  real(r8), intent(in) :: rpp    ! fraction of potential evaporation from leaf [-]

  real(r8), intent(out) :: wtaq  ! latent heat conductance for air [m/s]
  real(r8), intent(out) :: wtlq  ! latent heat conductance for leaf [m/s]
  real(r8), intent(out) :: wtgq  ! latent heat conductance for ground [m/s]
  real(r8), intent(out) :: wtaq0 ! normalized latent heat conductance for air [-]
  real(r8), intent(out) :: wtlq0 ! normalized latent heat conductance for leaf [-]
  real(r8), intent(out) :: wtgq0 ! normalized heat conductance for ground [-]
  real(r8), intent(out) :: wtalq ! normalized latent heat cond. for air and leaf [-]
  real(r8), intent(out) :: wtglq ! normalized latent heat cond. for leaf and ground [-]
  real(r8), intent(out) :: wtgaq ! normalized latent heat cond. for air and ground [-]

!----Local Variables----------------------------------------------------

  real(r8) wtsqi                 ! latent heat resistance for air, grd and leaf [s/m]

!----End Variable List--------------------------------------------------

  wtaq  = clm%frac_veg_nosno/raw                                ! air
  wtlq  = clm%frac_veg_nosno*(clm%elai+clm%esai)/rbw * rpp      ! leaf
  wtgq  = clm%frac_veg_nosno/rdw                                ! ground
  wtsqi = 1./(wtaq+wtlq+wtgq)

  wtgq0 = wtgq*wtsqi                    ! ground
  wtlq0 = wtlq*wtsqi                    ! leaf
  wtaq0 = wtaq*wtsqi                    ! air

  wtglq = wtgq0+wtlq0                   ! ground + leaf
  wtgaq = wtaq0+wtgq0                   ! air + ground
  wtalq = wtaq0+wtlq0                   ! air + leaf

end subroutine LatentHCond
