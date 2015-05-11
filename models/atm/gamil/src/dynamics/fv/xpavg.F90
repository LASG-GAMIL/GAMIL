!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  xpaxg --- Average a scalar latitude field
!
! !INTERFACE:
      subroutine xpavg(p, im)

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none

! !INPUT PARAMETERS:
      integer im

! !INPUT/OUTPUT PARAMETERS:
      real(r8) p(im)

!
! !DESCRIPTION:
!   This routine determines the average of the scalar latitude field p
!   and then sets all the values in p to the average.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i
      real(r8) sum1

      sum1 = 0.
      do i=1,im
         sum1 = sum1 + p(i)
      enddo
      sum1 = sum1 / im

      do i=1,im
         p(i) = sum1
      enddo

      return
!EOC
      end subroutine xpavg
!-----------------------------------------------------------------------

