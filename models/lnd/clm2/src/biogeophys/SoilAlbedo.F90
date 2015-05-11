#include <misc.h>
#include <preproc.h>

subroutine SoilAlbedo (clm, coszen, nband, albsnd, albsni)

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
! Determine ground surface albedo, accounting for snow
! 
! Method: 
! 
! Author: 
! Gordon Bonan
! April 2002: Vertenstein/Oleson/Levis; Final form
! 
!-----------------------------------------------------------------------
! $Id: SoilAlbedo.F90,v 1.2.10.3 2002/06/15 13:50:18 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varpar, only : numrad
  use clm_varcon, only : albsat, albdry, alblak, albice, tfrz, istice, istsoil
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm     !CLM 1-D Module

  integer, intent(in) :: nband           ! number of solar radiation waveband classes
  real(r8), intent(in) :: coszen         ! cosine solar zenith angle for next time step
  real(r8), intent(in) :: albsnd(numrad) ! snow albedo (direct)
  real(r8), intent(in) :: albsni(numrad) ! snow albedo (diffuse)

!----Local Variables----------------------------------------------------

  integer  ib      !waveband number (1=vis, 2=nir)
  real(r8) inc     !soil water correction factor for soil albedo
  real(r8) albsod  !soil albedo (direct)
  real(r8) albsoi  !soil albedo (diffuse)

!----End Variable List--------------------------------------------------

  do ib = 1, nband
     if (clm%itypwat == istsoil)  then               !soil
        inc    = max(0.11-0.40*clm%h2osoi_vol(1), 0._r8)
        albsod = min(albsat(clm%isoicol,ib)+inc, albdry(clm%isoicol,ib))
        albsoi = albsod
     else if (clm%itypwat == istice)  then           !land ice
        albsod = albice(ib)
        albsoi = albsod
     else if (clm%t_grnd > tfrz) then                !unfrozen lake, wetland
        albsod = 0.05/(max(0.001_r8,coszen) + 0.15)
        albsoi = albsod
     else                                            !frozen lake, wetland
        albsod = alblak(ib)
        albsoi = albsod
     end if

     clm%albgrd(ib) = albsod*(1.-clm%frac_sno) + albsnd(ib)*clm%frac_sno
     clm%albgri(ib) = albsoi*(1.-clm%frac_sno) + albsni(ib)*clm%frac_sno

  end do

  return
end subroutine SoilAlbedo
