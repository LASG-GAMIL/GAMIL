module dycore
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: dycore
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC     dycore_is

!
! !DESCRIPTION:
!
! Utility routines related to the dycore
!
!      \begin{tabular}{|l|l|} \hline \hline
!        dycore         & dycore being used \\ \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.01   Rosinski    Creation
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: dycore_is --- determine dynamical core in use
!
! !INTERFACE: 
   logical function dycore_is (name)

      implicit none

! !INPUT PARAMETERS:
      character(len=*) :: name

! !DESCRIPTION:
!   Determine the dynamical core in use.  
!   True if the dycore name is "lr" or "LR"
!
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC


      if (name == 'lr' .or. name == 'LR') then
         dycore_is = .true.
      else
         dycore_is = .false.
      end if

      return
!EOC
   end function dycore_is
!-----------------------------------------------------------------------
end module dycore
