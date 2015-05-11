#include <misc.h>
#include <params.h>
subroutine basiz(pkdim   ,sig     ,lbasiz  )
!-----------------------------------------------------------------------
!
! Purpose:
! Compute weights used in Lagrange cubic polynomial interpolation in 
! the central interval of a four point stencil.  Done for each interval
! in the unequally spaced height grid.
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: basiz.F90,v 1.2.42.2 2002/06/15 13:48:20 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: pkdim             ! vertical dimension
  real(r8), intent(in)   :: sig(pkdim)        ! sigma levels
  real(r8), intent(out)  :: lbasiz(4,2,pkdim) ! vertical interpolation weights
!
!---------------------------Local workspace-----------------------------
!
  integer  kk                                 ! index
!
!-----------------------------------------------------------------------
!
  do kk = 2,pkdim-2
     call lcbas(sig(kk-1),lbasiz(1,1,kk),lbasiz(1,2,kk) )
  end do
!
  return
end subroutine basiz




