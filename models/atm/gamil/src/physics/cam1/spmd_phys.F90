#include <misc.h>
#include <params.h>

module spmd_phys

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM for physics.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: beglat, endlat
   use ppgrid, only: begchunk, endchunk
   implicit none

   private
   public spmdinit_phys

CONTAINS

!========================================================================

   subroutine spmdinit_phys ()
      begchunk = beglat
      endchunk = endlat

      return
   end subroutine spmdinit_phys

end module spmd_phys
