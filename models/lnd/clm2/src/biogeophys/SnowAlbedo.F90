#include <misc.h>
#include <preproc.h>

subroutine SnowAlbedo (clm, coszen, nband, ind, alb)

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
! Determine snow albedos
! 
! Method: 
! 
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! April 2002: Vertenstein/Oleson/Levis; Final form
! 
!-----------------------------------------------------------------------
! $Id: SnowAlbedo.F90,v 1.2.10.3 2002/06/15 13:50:17 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm       !CLM 1-D Module

  real(r8), intent(in) :: coszen ! cosine solar zenith angle for next time step
  integer, intent(in) :: nband   ! number of solar radiation waveband classes
  integer, intent(in) :: ind     ! 0=direct beam, 1=diffuse radiation

  real(r8), intent(out):: alb(numrad) !snow albedo by waveband

!----Local Variables----------------------------------------------------

  integer  :: ib          !waveband class

! Variables and constants for snow albedo calculation

  real(r8) :: snal0 = 0.95 !vis albedo of new snow for sza<60 [-]
  real(r8) :: snal1 = 0.65 !nir albedo of new snow for sza<60 [-]
  real(r8) :: conn  = 0.5  !constant for visible snow alb calculation [-]
  real(r8) :: cons  = 0.2  !constant (=0.2) for nir snow albedo calculation [-]
  real(r8) :: sl    = 2.0  !factor that helps control alb zenith dependence [-]
  real(r8) :: age          !factor to reduce vis snow alb due to snow age [-]
  real(r8) :: albs         !temporary vis snow albedo [-]
  real(r8) :: albl         !temporary nir snow albedo [-]
  real(r8) :: cff          !snow alb correction factor for zenith angle > 60 [-]
  real(r8) :: czf          !solar zenith correction for new snow albedo [-]

!----End Variable List--------------------------------------------------

!
! Zero albedos
!

  do ib = 1, nband
     alb(ib) = 0._r8
  end do

!
! Albedo for snow cover.
! Snow albedo depends on snow-age, zenith angle, and thickness of snow,
! age gives reduction of visible radiation
!

!
! Correction for snow age
!

  age = 1.-1./(1.+clm%snowage)
  albs = snal0*(1.-cons*age)
  albl = snal1*(1.-conn*age)

  if (ind == 0) then

!
! czf corrects albedo of new snow for solar zenith
!

    cff    = ((1.+1./sl)/(1.+max(0.001_r8,coszen)*2.*sl )- 1./sl)
    cff    = max(cff,0._r8)
    czf    = 0.4*cff*(1.-albs)
    albs = albs+czf
    czf    = 0.4*cff*(1.-albl)
    albl = albl+czf

  endif

  alb(1) = albs
  alb(2) = albl

  return
end subroutine SnowAlbedo
