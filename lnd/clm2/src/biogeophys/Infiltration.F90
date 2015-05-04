#include <misc.h>
#include <preproc.h>

subroutine Infiltration (clm) 

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
! Calculate infiltration into surface soil layer (minus the evaporation)
!
! Method:
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
! $Id: Infiltration.F90,v 1.1.10.3 2002/06/15 13:50:16 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm	 !CLM 1-D Module

!----Local Variables----------------------------------------------------
!None
!----End Variable List--------------------------------------------------

!
! Infiltration into surface soil layer (minus the evaporation)
!

  if (clm%snl+1 >= 1) then
     clm%qflx_infl = clm%qflx_top_soil - clm%qflx_surf - clm%qflx_evap_grnd
  else
     clm%qflx_infl = clm%qflx_top_soil - clm%qflx_surf
  endif
  
end subroutine Infiltration
