#include <misc.h>
!---------------------------------------------------------------------------
!
! Purpose:
!
! 	Wrapper routines for the MPI (Message Passing) library for the
!	distributed memory (SPMD) version of the code. Also data with
!	"shorthand" names for the MPI data types.
!
! Author: Jim Rosinski
!
!---------------------------------------------------------------------------
!
! Compile these routines only when SPMD is defined
!
#if (defined SPMD)

!****************************************************************

   subroutine mpibarrier (comm)
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
!
! MPI barrier, have threads wait until all threads have reached this point
!
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_barrier (comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_barrier failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpibarrier
 
!****************************************************************
 
   subroutine mpifinalize
!
! End of all MPI communication
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   integer ier   !MP error code
 
   call mpi_finalize (ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_finalize failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpifinalize
 
!****************************************************************
 
   subroutine mpipack_size (incount, datatype, comm, size)
!
! Returns the size of the packed data
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(in):: comm
   integer, intent(out):: size
 
   integer ier   !MP error code
 
   call mpi_pack_size (incount, datatype, comm, size, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_pack_size failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack_size
 
!****************************************************************
 
   subroutine mpipack (inbuf, incount, datatype, outbuf, outsize,    &
                       position, comm)
!
! Pack the data and send it.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(out):: outsize
   integer, intent(inout):: position
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_pack (inbuf, incount, datatype, outbuf, outsize,         &
     &            position, comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_pack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack
 
!****************************************************************
 
   subroutine mpiunpack (inbuf, insize, position, outbuf, outcount,  &
                         datatype, comm)
!
! Un-packs the data from the packed receive buffer
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: insize
   integer, intent(inout):: position
   integer, intent(in):: outcount
   integer, intent(in):: datatype
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_unpack (inbuf, insize, position, outbuf, outcount,       &
     &              datatype, comm, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_unpack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpiunpack
 
!****************************************************************
 
   subroutine mpisendrecv (sendbuf, sendcount, sendtype, dest, sendtag,  &
                           recvbuf, recvcount, recvtype, source,recvtag, &
                           comm)
!
! Blocking send and receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real(r8), intent(in):: sendbuf(*)
   real(r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcount
   integer, intent(in):: sendtype
   integer, intent(in):: dest
   integer, intent(in):: sendtag
   integer, intent(in):: recvcount
   integer, intent(in):: recvtype
   integer, intent(in):: source
   integer, intent(in):: recvtag
   integer, intent(in):: comm
 
   integer :: status(MPI_STATUS_SIZE)
   integer ier   !MP error code
 
   call t_startf ('mpi_sendrecv')
   call mpi_sendrecv (sendbuf, sendcount, sendtype, dest, sendtag,   &
     &                recvbuf, recvcount, recvtype, source, recvtag, &
     &                comm, status, ier)
   if (ier.ne.mpi_success) then
      write(6,*)'mpi_sendrecv failed ier=',ier
      call endrun
   end if
!
! ASSUME nrecv = nsend for stats gathering purposes.  This is not actually
! correct, but its the best we can do since recvcount is a Max number
!
   nsend = nsend + 1
   nrecv = nrecv + 1
   nwsend = nwsend + sendcount
   nwrecv = nwrecv + sendcount

   call t_stopf ('mpi_sendrecv')
 
   return
   end subroutine mpisendrecv
 
!****************************************************************
 
   subroutine mpiisend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_isend')
   call mpi_isend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_isend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_isend')
 
   return
   end subroutine mpiisend
 
!****************************************************************
 
   subroutine mpiirecv (buf, count, datatype, source, tag, comm, request)
!
! Does a non-blocking receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
   call t_startf ('mpi_irecv')
   call mpi_irecv (buf, count, datatype, source, tag, comm, request, ier )
   if (ier/=mpi_success) then
      write(6,*)'mpi_irecv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
   call t_stopf ('mpi_irecv')
 
   return
   end subroutine mpiirecv
 
!****************************************************************
 
   subroutine mpiwaitall (count, array_of_requests, array_of_statuses)
!
! Waits for a collection of nonblocking operations to complete.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   integer, intent(in):: count
   integer, intent(inout):: array_of_requests(*)
   integer, intent(out):: array_of_statuses(*)
 
   integer ier   !MP error code
 
   call t_startf ('mpi_waitall')
   call mpi_waitall (count, array_of_requests, array_of_statuses, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_waitall failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_waitall')
 
   return
   end subroutine mpiwaitall
 
!****************************************************************
 
   subroutine mpisend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_send')
   call mpi_send (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_send failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
   call t_stopf ('mpi_send')
 
   return
   end subroutine mpisend
 
!****************************************************************
 
   subroutine mpirecv (buf, count, datatype, source, tag, comm)
!
! Does a blocking receive
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ier   !MP error code
 
   call t_startf ('mpi_recv')
   call mpi_recv (buf, count, datatype, source, tag, comm, status, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_recv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
   call t_stopf ('mpi_recv')
 
   return
   end subroutine mpirecv
 
!****************************************************************
 
   subroutine mpigather (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                         recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcnt
   integer, intent(in):: sendtype
   integer, intent(in):: recvcnt
   integer, intent(in):: recvtype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_gather')
   call mpi_gather (sendbuf, sendcnt, sendtype,                      &
     &              recvbuf, recvcnt, recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_gather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_gather')
 
   return
   end subroutine mpigather
 
!****************************************************************
 
   subroutine mpigatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
                          displs, recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   ! MPI error code
 
   call t_startf ('mpi_gather')
   call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
                     root, comm, ier)
   if (ier /= mpi_success) then
      write(6,*)'mpi_gather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_gather')
 
   return
   end subroutine mpigatherv
 
!****************************************************************
 
   subroutine mpisum (sendbuf, recvbuf, cnt, datatype, root, comm)
!
! Sums sendbuf across all processors on communicator, returning 
! result to root.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: cnt
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_reduce')
   call mpi_reduce (sendbuf, recvbuf, cnt, datatype, mpi_sum, &
                    root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_reduce failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_reduce')
 
   return
   end subroutine mpisum
 
!****************************************************************
 
   subroutine mpiscatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                          recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8),intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer,intent(in):: sendcnt
   integer,intent(in):: sendtype
   integer,intent(in):: recvcnt
   integer,intent(in):: recvtype
   integer,intent(in):: root
   integer,intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_scatter')
   call mpi_scatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
   &                 recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_scatter')
 
   return
   end subroutine mpiscatter
 
!****************************************************************
 
   subroutine mpiscatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, &
                           recvcnt, recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnt
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_scatter')
   call mpi_scatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, &
                      recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_scatter')
 
   return
   end subroutine mpiscatterv
 
!****************************************************************
 
   subroutine mpibcast (buffer, count, datatype, root, comm )
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(inout):: buffer(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_bcast')
   call mpi_bcast (buffer, count, datatype, root, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_bcast failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_bcast')
 
   return
   end subroutine mpibcast
!****************************************************************
 
   subroutine mpialltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                            recvbuf, recvcnts, rdispls, recvtype, &
                            comm)
!
! All-to-all scatter/gather
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: sdispls(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: rdispls(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
   call t_startf ('mpi_alltoallv')
   call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                       recvbuf, recvcnts, rdispls, recvtype, &
                       comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_alltoallv failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_alltoallv')
 
   return
   end subroutine mpialltoallv
!****************************************************************
 
   subroutine mpiallgatherint (sendbuf, scount, recvbuf, rcount, comm)
!
! Collects integer data from each task and broadcasts resulting
! vector to all tasks
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   implicit none
 
   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in) :: scount
   integer, intent(in) :: rcount
   integer, intent(in) :: comm
 
   integer ier   !MP error code

   call t_startf ('mpi_allgather')
   call mpi_allgather (sendbuf, scount, mpiint, recvbuf, rcount, &
                       mpiint, comm, ier)
   if (ier/=mpi_success) then
      write(6,*)'mpi_allgather failed ier=',ier
      call endrun
   end if
   call t_stopf ('mpi_allgather')
 
   return
   end subroutine mpiallgatherint
!
! If SPMD is not turned on
!
#else
   subroutine wrap_mpi
   implicit none
!
! A unused stub routine to make the compiler happy when SPMD is
! turned off (which means you don't need anything in this file).
!
   write(6,*) '(WRAP_MPI): This should not be called at all'
   call endrun
   end subroutine wrap_mpi
#endif


