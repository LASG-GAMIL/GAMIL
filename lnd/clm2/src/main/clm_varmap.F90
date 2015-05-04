#include <misc.h>
#include <preproc.h>

module clm_varmap

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module of mapping arrays
! 
! Method: 
! The land surface model works by gathering all the land points on a
! [lon]x[lat] grid into a vector of [numland] land points. This
! is then expanded into a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch] subgrid patches per land point. 
! [ixy], [jxy], [patch], and [land] are indices for the mapping: 
! [lsmlon]x[lsmlat] grid <-> [numland] vector of land points <-> [numpatch] vector 
! of subgrid points. [landvec%wtxy] are the weights to obtain the 
! grid average from the subgrid patches.
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm_varmap.F90,v 1.1.22.4 2002/06/15 13:50:29 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  integer :: numland                 !number of land points
  integer :: begland                 !beginning land index (minimum value is 1)
  integer :: endland                 !ending land index (maximum value is  numpatch)
  integer :: numpatch                !total number of patches allowing for subgrid patches
  integer :: begpatch                !beginning patch index (minimum value is 1)
  integer :: endpatch                !ending patch index (maximum value is  numpatch)

  type land1d
     integer ,pointer :: patch(:,:)  !patch vector index: 1 to numpatch
     integer, pointer :: ixy(:)      !longitude index for each land point: 1 to lsmlon
     integer, pointer :: jxy(:)      !latitude index for each land point: 1 to lsmlat
     real(r8),pointer :: wtxy(:,:)   !subgrid weights 
  end type land1d

  type subgrid1d
     integer , pointer :: ixy(:)     !longitude index for each patch point (1 to lsmlon)
     integer , pointer :: jxy(:)     !latitude index for each patch point (1 to lsmlat)
     integer , pointer :: mxy(:)     !subgrid patch type index (1 to maxpatch)
     real(r8), pointer :: wtxy(:)    !subgrid weight for each patch point 
     integer , pointer :: land(:)    !land index for each patch point
  end type subgrid1d

  type (land1d)    :: landvec
  type (subgrid1d) :: patchvec

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine mapvar_ini

    use shr_kind_mod, only: r8 => shr_kind_r8
    use infnan
    use clm_varpar, only : maxpatch
    implicit none

! Initializes mapping vectors

    allocate (landvec%ixy(numland))
    allocate (landvec%jxy(numland))
    allocate (landvec%wtxy(numland,maxpatch))
    allocate (landvec%patch(numland,maxpatch))

    allocate (patchvec%ixy(numpatch))
    allocate (patchvec%jxy(numpatch))
    allocate (patchvec%mxy(numpatch))
    allocate (patchvec%wtxy(numpatch))
    allocate (patchvec%land(numpatch))

! Initialize appropriate components as infinity

    landvec%wtxy(:,:) = inf
    patchvec%wtxy(:)  = inf

  end subroutine mapvar_ini

!=======================================================================

end module clm_varmap
