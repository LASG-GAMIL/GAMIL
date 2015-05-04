#include <misc.h>
#include <preproc.h>

function StabilityFunc(k, zeta)

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
! Stability function for rib < 0.
!
! Method:
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: StabilityFunc.F90,v 1.1.10.4 2002/06/15 13:50:19 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_PI
  implicit none

!----Local Variables----------------------------------------------------

  integer k               ! stability function form (wind or temperature/humidity) [-]
  real(r8) zeta           ! dimensionless height used in Monin-Obukhov theory [-]
  real(r8) StabilityFunc  ! stability function for unstable case
  real(r8) chik           ! temporary variable [-]

!=== End Variable List ===================================================

  chik = (1.-16.*zeta)**0.25
  if (k == 1) then
     StabilityFunc = 2.*log((1.+chik)*0.5) &
          + log((1.+chik*chik)*0.5)-2.*atan(chik)+SHR_CONST_PI*0.5
  else
     StabilityFunc = 2.*log((1.+chik*chik)*0.5)
  endif

end function StabilityFunc
