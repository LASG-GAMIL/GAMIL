#include <misc.h>
#include <preproc.h>

module ioUnitMod

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  logical :: lsmiou(99)  !I/O file unit numbers (1 to 99): true if active

!=======================================================================
contains
!=======================================================================

  integer function getavu()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! get next available Fortran unit number
!
! Method: 
! Get next available Fortran unit number itst. Set lsmiou(itst), in 
! lsmio common block, true. If coupled to CAM, use CAM function navu
! to get available unit number, in which case lsmiou is not needed.
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

#if (defined COUP_CAM)
    use units     !CAM units module
#endif

! ------------------------ local variables ------------------------
    integer itst  !Fortran unit number
! -----------------------------------------------------------------

#if (defined COUP_CAM)
    getavu = getunit()
    RETURN
#else
    do itst = 1, 99
       if (.not.lsmiou(itst)) then
          getavu = itst
          lsmiou(itst) = .true.
          RETURN
       end if
    end do
    write (6,*) 'GETAVU error: ran out of Fortran unit numbers'
    call endrun
#endif
  end function getavu

!=======================================================================

  subroutine relavu (iunit)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! close and release Fortran unit no longer in use
!
! Method: 
! Close and release Fortran unit number iunit. Set lsmiou(iunit) to 
! false. If coupled to cam, use cam function relunit to close/release 
! unit number.
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

#if (defined COUP_CAM)
    use units     !CAM units module
#endif

! ------------------------ arguments ------------------------------
    integer, intent(in) :: iunit    !Fortran unit number
! -----------------------------------------------------------------

#if (defined COUP_CAM)
    close(iunit)
    call freeunit(iunit)
#else
    if (.not.lsmiou(iunit)) then
       write (6,*) 'RELAVU eror: unit ',iunit,' is not flagged as in use'
       call endrun
    end if
    if (iunit<1 .or. iunit>99) then
       write (6,*) 'RELAVU error: attempt to return out of range unit'
       call endrun
    end if
    close(iunit)
    lsmiou(iunit) = .false.
#endif
    return
  end subroutine relavu

!=======================================================================

end module ioUnitMod







