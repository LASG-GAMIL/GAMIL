#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum --- Calculate x-sum bit-wise consistently
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum(a, ifirst, ilast, im, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      use pmgrid
#if defined ( SPMD )
      use spmd_dyn, only : commxy_x
      use parutilitiesmodule, only : parexchangevector
#endif
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

! !INPUT PARAMETERS:
      integer im                     ! global longitudes
      integer ifirst                 ! first longitude on this PE
      integer ilast                  ! last longitude on this PE
      integer ltot                   ! number of quantities to be summed
      real (r8) a(ifirst:ilast,ltot) ! input vector to be summed

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)           ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !CALLED FROM:
!     te_map
!
! !REVISION HISTORY:
!
!     AAM 00.11.01 : Created
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real (r8) quan_all(im*ltot)
      integer i,l,icount,ipe

#if defined ( SPMD )
      real (r8) quan_send(nprxy_x*ltot*(ilast-ifirst+1))
      integer  sendcount(nprxy_x)
      integer  recvcount(nprxy_x)
#endif

#if defined ( SPMD ) 
      icount=0
      do ipe=1,nprxy_x
        sendcount(ipe) = ltot*(ilast-ifirst+1)
        do i=ifirst,ilast
        do l=1,ltot
          icount=icount+1    
          quan_send(icount)=a(i,l)
        enddo
        enddo
      enddo
      call parexchangevector( commxy_x, sendcount, quan_send,     &
                                          recvcount, quan_all )
#else
      do l=1,ltot
        do i=i,im
          quan_all((i-1)*ltot+l)=a(i,l)
        enddo
      enddo
#endif

      do l=1,ltot
        sum(l)=0.
        do i=1,im
          sum(l) = sum(l) + quan_all((i-1)*ltot+l)
        enddo
      enddo

      return
!EOC
      end
!-----------------------------------------------------------------------

