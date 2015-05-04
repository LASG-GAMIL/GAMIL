#include <misc.h>
#include <params.h>

module units

implicit none

PRIVATE

   logical :: lutag(0:99) = .false.      ! list of flags marking logical units in use

   public :: getunit, freeunit

CONTAINS

   integer function getunit (iu)
!
! Arguments
!
   integer, intent(in), optional :: iu   ! desired unit number
!
! Local workspace
!
   integer :: n                          ! loop index

   if (present (iu)) then
      if (iu < 0 .or. iu > 99) then
         write(6,*)'GETUNIT: invalid unit number request:', iu
         call endrun
      else if (lutag(iu) .or. iu == 0 .or. iu == 5 .or. iu == 6) then
         write(6,*)'GETUNIT: unit number ', iu, ' is already in use'
         call endrun
      else
         getunit = iu
         lutag (iu) = .true.
         return
      end if

   else
!
! Choose first available unit other than 0, 5, or 6
!
      do n=1,99
         if (n == 5 .or. n == 6) then
            cycle
         end if
         if (.not.lutag(n)) then
            getunit = n
            lutag(n) = .true.
            return
         end if
      end do
   end if

   write(6,*)'GETUNIT: an available unit could not be found'
   call endrun

   end function getunit

!#######################################################################

   subroutine freeunit (iu)
!
! Arguments
!
   integer, intent(in) :: iu       ! unit number to be freed

   if (iu < 0 .or. iu > 99) then
      write(6,*)'FREEUNIT: invalid unit number request:', iu
   else if (iu == 0 .or. iu == 5 .or. iu == 6) then
      write(6,*)'FREEUNIT: units 0, 5, and 6 must not be freed'
      call endrun
   else if (lutag(iu)) then
      lutag (iu) = .false.
   else
      write(6,*)'FREEUNIT: unit ', iu, ' was not in use'
   end if

   return
   end subroutine freeunit

end module units
