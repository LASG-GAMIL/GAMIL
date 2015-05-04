#include <misc.h>
#include <preproc.h>

subroutine WetIceHydrology (clm)

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
! Calculate hydrology for ice and wetland
!
! Method:
!
! Author:
! 7 November 2000: Mariana Vertenstein; Initial code
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: WetIceHydrology.F90,v 1.2.6.3 2002/06/15 13:50:21 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : istwet, istice
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------
! None
!----End Variable List--------------------------------------------------

!
! Wetland and land ice runoff
!

  clm%qflx_drain  = 0. 
  clm%qflx_surf   = 0.
  clm%qflx_infl   = 0.
  clm%qflx_qrgwl  = clm%forc_rain + clm%forc_snow - clm%qflx_evap_tot - &
                    (clm%endwb - clm%begwb)/clm%dtime

end subroutine WetIceHydrology
