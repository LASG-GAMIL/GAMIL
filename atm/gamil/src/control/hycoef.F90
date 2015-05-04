#include <misc.h>
#include <params.h>

subroutine hycoef 

!! (wanhui 2003.05.16)
!! (wanhui 2003.10.23)
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
!
! $Id: hycoef.F90,v 1.2.4.1 2002/06/15 13:47:27 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use stdatm, only: p00   !!(wh 2003.10.23)

  implicit none

#include <comhyb.h>

!---------------------------Local workspace-----------------------------
    integer k                 ! Level index
!!  real(r8) amean,bmean,atest,btest,eps
!-----------------------------------------------------------------------
!
!!  eps    = 1.e-05
!!  nprlev = 0
    ps0    = 1.0e5            ! Base state surface pressure (pascals)
!!  psr    = 1.0e5            ! Reference surface pressure (pascals)
!
! Set layer locations
!
! Reference state interface pressures
!
  do k=1,plevp
     hypi(k) = hyai(k)*ps0 + hybi(k)*p00
  end do
!
  do k=1,plev
     hypm(k) = 0.5*( hypi(k) + hypi(k+1) ) 
  end do

  return
end subroutine hycoef

