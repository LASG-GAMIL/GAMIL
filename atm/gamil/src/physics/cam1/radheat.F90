#include <misc.h>
#include <params.h>

module radheat
!-----------------------------------------------------------------------
! Purpose: to provide an interface to convert shortwave and longwave 
! radiative heating terms from radctl into net heating. This module exists
! primarily to provide a hook for waccm, to allow incorporating additional
! radiative terms (eUV heating and nonLTE longwave cooling).
! Public functions/subroutines:
!     radheat_net
! 
! Author: B.A. Boville
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver
  use physics_types, only: physics_state, physics_ptend

  implicit none
  private          ! Make default type private to the module
!
! Public interfaces
!
  public radheat_net   ! return net radiative heating

!===============================================================================
contains
!===============================================================================
  subroutine radheat_net(state, ptend, qrl, qrs)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and return in ptend.
! This routine provides the waccm hook for computing nonLTE cooling and 
! eUV heating. 
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: qrl(pcols,pver)        ! LTE longwave heating
    real(r8), intent(in) :: qrs(pcols,pver)        ! shortwave heating (>200 nm)

    type(physics_state), intent(in)    :: state   ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend   ! indivdual parameterization tendencie

!---------------------------Local storage-------------------------------
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
!-----------------------------------------------------------------------
    lchnk = state%lchnk
    ncol = state%ncol

    ptend%name       = 'radheat'
    ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))
    ptend%ls         = .TRUE.

    return
  end subroutine radheat_net
end module radheat
