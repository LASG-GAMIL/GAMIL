#include <misc.h>
#include <preproc.h>

module pft_varcon

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module of vegetation constants
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: pft_varcon.F90,v 1.4.6.2.6.1 2002/10/03 20:07:41 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar
  implicit none

! Vegetation type constants

  character(len=40) pftname(0:numpft) !PFT description

  integer ncorn                  !value for corn
  integer nwheat                 !value for wheat
  integer noveg                  !value for not vegetated 
  integer ntree                  !value for last type of tree

  real(r8) dleaf(0:numpft)       !characteristic leaf dimension (m) 
  real(r8) c3psn(0:numpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(r8) vcmx25(0:numpft)      !max rate of carboxylation at 25C (umol CO2/m**2/s)
  real(r8) mp(0:numpft)          !slope of conductance-to-photosynthesis relationship
  real(r8) qe25(0:numpft)        !quantum efficiency at 25C (umol CO2 / umol photon)
  real(r8) xl(0:numpft)          !leaf/stem orientation index
  real(r8) rhol(0:numpft,numrad) !leaf reflectance: 1=vis, 2=nir 
  real(r8) rhos(0:numpft,numrad) !stem reflectance: 1=vis, 2=nir 
  real(r8) taul(0:numpft,numrad) !leaf transmittance: 1=vis, 2=nir 
  real(r8) taus(0:numpft,numrad) !stem transmittance: 1=vis, 2=nir 
  real(r8) z0mr(0:numpft)        !ratio of momentum roughness length to canopy top height (-)
  real(r8) displar(0:numpft)     !ratio of displacement height to canopy top height (-)
  real(r8) roota_par(0:numpft)   !CLM rooting distribution parameter [1/m]
  real(r8) rootb_par(0:numpft)   !CLM rooting distribution parameter [1/m]

end module pft_varcon

