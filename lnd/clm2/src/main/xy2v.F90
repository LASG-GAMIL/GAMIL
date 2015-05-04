#include <misc.h>
#include <preproc.h>

subroutine xy2v (nx, ny, fldxy, ki, kf, fldv)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! convert a grid-average field to subgrid patch vector
! 
! Method: 
! This code converts a grid-average field [fldxy] dimensioned
! [lsmlon] x [lsmlat] to a subgrid patch vector [fldv] for 
! [numpatch] subgrid patches. 
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: xy2v.F90,v 1.1.22.2 2002/06/15 13:50:38 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varmap, only : patchvec
  implicit none

! ------------------------ arguments----------------------------------
  integer , intent(in)  :: nx, ny          !x-y dimension	
  integer , intent(in)  :: ki, kf          !beginning and end patch indices
  real(r8), intent(in)  :: fldxy(nx,ny)    !gridded input
  real(r8), intent(out) :: fldv(ki:kf)     !subgrid vector output
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
  integer i,j,k             !indices
! ---------------------------------------------------------------

  do k = ki,kf
     i = patchvec%ixy(k)
     j = patchvec%jxy(k)
     fldv(k) = fldxy(i,j)
  end do

  return
end subroutine xy2v





