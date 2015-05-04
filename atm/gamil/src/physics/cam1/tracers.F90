#include <misc.h>
#include <params.h>

module tracers

!-----------------------------------------------------------------------
!
! Purpose: Contains data and functions for manipulating advected and non-advected tracers
!
! Public functions/subroutines:
!   initindx, add_cnst
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------
    use constituents, only: pcnst, pnats, ppcnst

    implicit none

    public

    save

    integer nusr_adv                       ! Number of user defined advected tracers
    integer nusr_nad                       ! Number of user defined non-advected tracers
    !
    ! Offsets for both advected and non-advected tracer species
    !
    integer ixuadv   ! user advected tracer beginning index
    integer ixunad   ! user non-advected tracer beginning index
    integer ixmoist  ! moisture beginning index
    integer ixcldw   ! cloud water beginning index
    integer ixtrct   ! test tracers beginning index
    !
    ! Lists of tracer names and diagnostics
    !
    character(len=8) :: hadvnam(pcnst)        ! names of horizontal advection tendencies
    character(len=8) :: vadvnam(pcnst)        ! names of vertical advection tendencies
    character(len=8) :: dcconnam(ppcnst)      ! names of convection tendencies
    character(len=8) :: fixcnam(pcnst)        ! names of species slt fixer tendencies
    character(len=8) :: tendnam(pcnst)        ! names of total tendencies of species
    character(len=8) :: sflxnam(ppcnst)       ! names of surface fluxes of species
    character(len=8) :: tottnam(pcnst)        ! names for horz + vert + fixer tendencies

end module tracers
