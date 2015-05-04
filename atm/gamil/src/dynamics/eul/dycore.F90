module dycore
!
! Data and utility routines related to the dycore
!
   implicit none

PRIVATE

   public :: dycore_is

CONTAINS

   logical function dycore_is (name)
!
! Input arguments
!
      character(len=*) :: name
      
      if (name == 'eul' .or. name == 'EUL') then
         dycore_is = .true.
      else
         dycore_is = .false.
      end if
      
      return
   end function dycore_is

end module dycore


