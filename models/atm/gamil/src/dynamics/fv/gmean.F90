#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean --- Calculate the mean of a 2D field
!
! !INTERFACE:

subroutine gmean(im,  jm,  jfirst,  jlast,  q, qmean)

! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8
  use dynamics_vars, only : gw
  use pmgrid, only: npr_y

#if defined( SPMD )
  use parutilitiesmodule, only : parcollective, sumop
  use spmd_dyn, only: comm_y
#endif

  implicit none

! !INPUT PARAMETERS:

  integer im, jm                       ! Horizontal dimensions
  integer jfirst, jlast                ! Latitude strip
  real(r8), intent(in) :: q(im,jfirst:jlast)              ! 2D field 

  real(r8) qmean

! !DESCRIPTION:
!     Calculate the mean of a 2D field
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.10   Lin     Revised
!   01.06.27   Mirin   Use y communicator
!
!EOP
!-----------------------------------------------------------------------
!BOC

  real(r8)  xsum(jm)
  integer i, j

  do j=1,jm
     xsum(j) = 0.
  enddo
  do j=jfirst,jlast
     do i=1,im
        xsum(j) = xsum(j) + q(i,j)
     enddo
     xsum(j) = xsum(j)*gw(j)
  enddo

#if defined( SPMD )
  if (npr_y .ne. 1) then
     call parcollective( comm_y, sumop, jm, xsum )
  endif
#endif

  qmean = 0.0
  do j=1,jm
     qmean = qmean + xsum(j)
  enddo
  qmean = qmean / (2*im)

  return
!EOC
end subroutine gmean
!-----------------------------------------------------------------------
