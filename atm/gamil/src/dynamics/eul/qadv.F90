!! (wanhui 2003.10.24)
!!---------------------

module qadv
!-----------------------------------------------------------------------
! Purpose:  variables for the H2O advection scheme
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plond, beglatexdyn, endlatexdyn, plev
  use infnan

  implicit none

    integer,parameter :: nonos = 0
    integer,parameter :: isor  = 3
    real(r8),   parameter :: ep    = 1.0d-10

    integer  :: iord
    integer  :: ipq (plond)
    real(r8) :: dtdsg(plev)
    real(r8) :: dsnp, dssp

    real(r8),allocatable ::   gc  (:)
    real(r8),allocatable :: dtdlt (:)
    real(r8),allocatable :: dtdln (:)

CONTAINS

  subroutine initialize_qadv
!
! Purpose:  Allocate and initialize the arrays of the standard atmosphere.
!
    allocate (gc    (beglatexdyn:endlatexdyn))
    allocate (dtdlt (beglatexdyn:endlatexdyn))
    allocate (dtdln (beglatexdyn:endlatexdyn))

    iord     = 3
    dsnp     = inf
    dssp     = inf

    ipq  (:) = bigint
    dtdsg(:) = inf

    gc   (:) = inf
    dtdlt(:) = inf
    dtdln(:) = inf

    return
  end subroutine initialize_qadv


end module qadv
