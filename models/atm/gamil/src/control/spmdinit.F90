#include <misc.h>
#include <params.h>

subroutine spmdinit
!-----------------------------------------------------------------------
!
! Purpose: MPI initialization routine:
!
! Method: get number of cpus, processes, tids, etc
!
! Author: CCM Core Group
!
!-----------------------------------------------------------------------

#if ( defined SPMD )
   use mpishorthand, only: mpiint, mpichar, mpilog, mpipk, mpir8, mpir4, &
                           mpicom, mpi_max_processor_name, mpi_integer, &
                           mpi_character, mpi_logical, mpi_real8, mpi_real4, &
                           mpi_packed, mpi_comm_world
   use spmd_dyn, only:  npes, spmdinit_dyn
   use spmd_phys, only: spmdinit_phys
#endif
   use pmgrid, only: plat, masterproc, iam

#ifdef COUP_CSM
   !use shr_msg_mod, only: shr_msg_comm_atm ! For FGOALS2.0
#endif

   implicit none

#if ( defined SPMD )
#include <commpi.h>          !!(wh)
!
! Local workspace
!
   integer i                 ! indices
   integer ier               ! return error status
   integer, allocatable :: length(:)  ! length of name

   character*(mpi_max_processor_name), allocatable :: proc_name(:) ! returned processor name

   logical mpi_running       ! returned value indicates if MPI_INIT has been called
!---------------------------------------------------------------------------
!
! Initialize mpi
!
   call mpi_initialized (mpi_running, ier)
   if (.not.mpi_running) call mpi_init (ier)
!
! Set mpishorthand variables.  Need to set as variables rather than parameters since
! some MPI implementations set values for MPI tags at run time
!
   mpiint  = mpi_integer
   mpichar = mpi_character
   mpilog  = mpi_logical
   mpir4   = mpi_real4
   mpir8   = mpi_real8
   mpipk   = mpi_packed
#ifdef COUP_CSM
   !mpicom  = shr_msg_comm_atm ! For FGOALS2.0
#else
   call mpi_comm_dup(mpi_comm_world, mpicom, ier)
#endif
!
! Get my id
!
   call mpi_comm_rank (mpicom, iam, ier)
   if (iam == 0) then
      masterproc = .true.
   else
      masterproc = .false.
   end if
   myrank = iam    !!(wh)
!
! Get number of processors
!
   call mpi_comm_size (mpicom, npes, ier)
   nprocs = npes       !!(wh)
   allocate ( length(0:npes-1) )
   allocate ( proc_name(0:npes-1) )
   proc_name(:) = ' '

!
! Get processor names and send to root. "1" is the msg tag
!
   call mpi_get_processor_name (proc_name(iam), length(iam), ier)

   if (masterproc) then
      do i=1,npes-1
         call mpirecv (proc_name(i), mpi_max_processor_name, mpichar, i, 1, mpicom)
      end do
      write(6,*) npes, 'pes participating in computation'
      write(6,*) '-----------------------------------'
      write(6,*) 'NODE#  NAME'
      do i=0,npes-1
         write(6,'(i3,2x,a)') i,trim(proc_name(i))
      end do
   else
      call mpisend (proc_name(iam), mpi_max_processor_name, mpichar, 0, 1, mpicom)
   end if

   call set_mpi_params()
!
! Currently spmdinit_dyn must be called before spmdinit_phys because the latter just copies
! in data computed in the former
!
   call spmdinit_dyn ()
   call spmdinit_phys ()
   deallocate(length)
   deallocate(proc_name)
#endif
   return
end subroutine spmdinit

