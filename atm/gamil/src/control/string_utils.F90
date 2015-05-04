module string_utils


   implicit none
   private

! Public interface methods

   public ::&
      to_upper   ! Convert character string to upper case

contains

function to_upper(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to upper case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id: string_utils.F90,v 1.1 2001/10/19 17:50:28 eaton Exp $
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to upper case
   character(len=len(str))      :: to_upper

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= 97  .and.  aseq <= 122 ) ctmp = achar(aseq - 32)
      to_upper(i:i) = ctmp
   end do

end function to_upper

end module string_utils
