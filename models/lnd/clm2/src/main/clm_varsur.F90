#include <misc.h>
#include <preproc.h>

module clm_varsur

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 2-d surface boundary data 
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm_varsur.F90,v 1.6.6.3.6.1 2002/10/03 20:07:35 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar
  implicit none

! land model grid

  integer  numlon(lsmlat)                !longitude points for each latitude strip
  real(r8) latixy(lsmlon,lsmlat)         !latitude of grid cell (degrees)
  real(r8) longxy(lsmlon,lsmlat)         !longitude of grid cell (degrees)
  real(r8) area(lsmlon,lsmlat)           !grid cell area (km**2)
  real(r8) lats(lsmlat+1)                !grid cell latitude, southern edge (degrees)
  real(r8) lonw(lsmlon+1,lsmlat)         !grid cell longitude, western edge (degrees)
  real(r8) lsmedge(4)                    !North,East,South,West edges of grid (deg)
  logical :: pole_points                 !true => grid has pole points
  logical :: fullgrid  = .true.          !true => no grid reduction towards poles
  logical :: offline_rdgrid              !true => read offline grid rather than creating it

! fractional land and mask

  integer  landmask(lsmlon,lsmlat)       !land mask: 1 = land. 0 = ocean
  real(r8) landfrac(lsmlon,lsmlat)       !fractional land

! surface boundary data 

  integer , allocatable :: soic2d(:,:)   !soil color
  real(r8), allocatable :: sand3d(:,:,:) !soil texture: percent sand
  real(r8), allocatable :: clay3d(:,:,:) !soil texture: percent clay
  real(r8), allocatable :: pctgla(:,:)   !percent of grid cell that is glacier
  real(r8), allocatable :: pctlak(:,:)   !percent of grid cell that is lake
  real(r8), allocatable :: pctwet(:,:)   !percent of grid cell that is wetland
  real(r8), allocatable :: pcturb(:,:)   !percent of grid cell that is urbanized

! lake and soil levels

  real(r8) :: zlak(1:nlevlak)            !lake z  (layers) 
  real(r8) :: dzlak(1:nlevlak)           !lake dz (thickness)
  real(r8) :: zsoi(1:nlevsoi)            !soil z  (layers)
  real(r8) :: dzsoi(1:nlevsoi)           !soil dz (thickness)
  real(r8) :: zisoi(0:nlevsoi)           !soil zi (interfaces)  


!=======================================================================
CONTAINS
!=======================================================================

  subroutine varsur_alloc
    allocate (soic2d(lsmlon,lsmlat))   
    allocate (sand3d(lsmlon,lsmlat,nlevsoi))
    allocate (clay3d(lsmlon,lsmlat,nlevsoi))
    allocate (pctgla(lsmlon,lsmlat))   
    allocate (pctlak(lsmlon,lsmlat))   
    allocate (pctwet(lsmlon,lsmlat))   
    allocate (pcturb(lsmlon,lsmlat))   
  end subroutine varsur_alloc

  subroutine varsur_dealloc
    deallocate (soic2d)   
    deallocate (sand3d)
    deallocate (clay3d)
    deallocate (pctgla)   
    deallocate (pctlak)   
    deallocate (pctwet)   
    deallocate (pcturb)   
  end subroutine varsur_dealloc

end module clm_varsur
