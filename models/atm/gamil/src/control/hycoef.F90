#include <misc.h>
#include <params.h>

subroutine hycoef 

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Defines the locations of model interfaces from input data in the
! hybrid coordinate scheme.  Actual pressure values of model level
! interfaces are determined elsewhere from the fields set here.
! 
! Method: 
! the following fields are set:
! hyai     fraction of reference pressure used for interface pressures
! hyam     fraction of reference pressure used for midpoint pressures
! hybi     fraction of surface pressure used for interface pressures
! hybm     fraction of surface pressure used for midpoint pressures
! hybd     difference of hybi's
! hypi     reference state interface pressures
! hypm     reference state midpoint pressures
! hypd     reference state layer thicknesses
! hypdln   reference state layer thicknesses (log p)
! hyalph   distance from interface to level (used in integrals)
! prsfac   log pressure extrapolation factor (used to compute psl)
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid

  implicit none

#include <comhyb.h>

  integer k

  ps0 = 1.0e5            ! Base state surface pressure (pascals)
!
! Set layer locations
!
! Reference state interface pressures
!
  do k = 1, plevp
     hypi(k) = pmtop + sig(k) * (ps0 - pmtop)
  end do

  do k = 1, plev
     hypm(k) = 0.5 * (hypi(k) + hypi(k+1)) 
  end do

end subroutine hycoef

