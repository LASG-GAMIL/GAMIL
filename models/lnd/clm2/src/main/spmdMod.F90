#include <misc.h>
#include <preproc.h>

module spmdMod

!-----------------------------------------------------------------------
!
! Purpose:
! MPI routines for initialization and computing arguments for
! gatherv and scatterv operations
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: spmdMod.F90,v 1.7.2.2 2001/11/26 15:28:11 mvertens Exp $
!-----------------------------------------------------------------------

#ifdef COUP_CSM
  use shr_msg_mod
#endif

#if (!defined SPMD)
  logical :: masterproc = .true. ! proc 0 logical for printing msgs
  integer :: iam = 0
#endif

#if (defined SPMD)

#if (defined COUP_CAM)
  use mpishorthand
  use spmd_dyn, only: npes
  use pmgrid  , only: masterproc, iam
#endif

#if (defined OFFLINE)
  use mpishorthand
#endif

#ifdef COUP_CSM
  use mpishorthand, only : mpiint, mpichar, mpilog, mpir8, mpicom, mpipk
#endif

#if (defined OFFLINE) || (defined COUP_CSM)
  integer :: npes        !number of processors
  integer :: iam         !proc number
  logical :: masterproc  !proc 0 logical for printing msgs
#endif

  integer, public, allocatable :: proc_landi(:)
  integer, public, allocatable :: proc_landf(:)
  integer, public, allocatable :: proc_patchi(:)
  integer, public, allocatable :: proc_patchf(:)
  integer, public, allocatable :: proc_patchpts(:)
  integer, public, allocatable :: proc_landpts(:)

  SAVE

!===============================================================================
CONTAINS
!===============================================================================

#if (defined OFFLINE) || (defined COUP_CSM)

  subroutine spmd_init

!-----------------------------------------------------------------------
!
! Purpose: MPI initialization (number of cpus, processes, tids, etc)
!
!-----------------------------------------------------------------------

  implicit none

! ------------------------ local variables -----------------------------
    integer i,j        ! indices
    integer ier        ! return error status
    integer, allocatable :: length(:)
    integer, allocatable :: displ(:)
    character*(MPI_MAX_PROCESSOR_NAME), allocatable :: proc_name(:)
#if (defined OFFLINE)
    logical mpi_running
#endif
!-----------------------------------------------------------------------

#if (defined OFFLINE)

! Initialize mpi

    call mpi_initialized (mpi_running, ier)
    if (.not. mpi_running) call mpi_init (ier)

#endif

! Set mpishorthand variables.  Need to set as variables rather
! than parameters since some MPI implementations set values for
! MPI tags at run time.

    mpiint  = mpi_integer
    mpichar = mpi_character
    mpilog  = mpi_logical
    mpir8   = mpi_real8
    mpipk   = mpi_packed
#if (defined OFFLINE)
    mpicom  = MPI_COMM_WORLD
#elif (defined COUP_CSM)
    mpicom  = SHR_MSG_COMM_LND
#endif

! Get my processor id

    call mpi_comm_rank(mpicom, iam, ier)
    if (iam==0) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

! Get number of processors

    call mpi_comm_size(mpicom, npes, ier)

! Get my processor names

    allocate (length(0:npes-1))
    allocate (displ(0:npes-1))
    allocate (proc_name(0:npes-1))
    call mpi_get_processor_name (proc_name(iam),length(iam),ier)
    call mpi_allgather (length(iam),1,mpiint,length,1,mpiint,mpicom,ier)
    do i =0,npes-1
       displ(i)=i*MPI_MAX_PROCESSOR_NAME
    end do
    call mpi_gatherv (proc_name(iam),length(iam),mpichar, &
                      proc_name,length,displ,mpichar,0,mpicom,ier)
    if (masterproc) then
       write(6,100)npes
       write(6,200)
       write(6,220)
       do i=0,npes-1
          write(6,250)i,(proc_name((i))(j:j),j=1,length(i))
       end do
    endif
    deallocate (length)
    deallocate (displ)
    deallocate (proc_name)

100 format(i3," pes participating in computation")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i3,")",2x,100a1)

    return
  end subroutine spmd_init

#endif

!===============================================================================

  subroutine spmd_init_patch

!-----------------------------------------------------------------------
!
! Purpose: Initialize arrays for number of land/patch points per proc
!
!-----------------------------------------------------------------------

    allocate (proc_landi(0:npes-1))
    allocate (proc_landf(0:npes-1))
    allocate (proc_landpts(0:npes-1))
    allocate (proc_patchi(0:npes-1))
    allocate (proc_patchf(0:npes-1))
    allocate (proc_patchpts(0:npes-1))

    return
  end subroutine spmd_init_patch

!===============================================================================

  subroutine compute_mpigs_patch (nfact, numtot, numperproc, displs)

!------------------------------------------------------------------
!
! Purpose: Compute arguments for gatherv, scatterv for patche vectors
!
!------------------------------------------------------------------

    implicit none

! ------------------- arguments -----------------------------------
    integer, intent(in ) :: nfact                ! multiplicative factor for patches
    integer, intent(out) :: numtot               ! total number of elements (to send or recv)
    integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
    integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!------------------------------------------------------------------

! ---------------------- local variables --------------------------
    integer :: p                                 ! index
!------------------------------------------------------------------

    numtot = (proc_patchpts(iam))*nfact

    do p=0,npes-1
       numperproc(p) = proc_patchpts(p)*nfact
    end do

    displs(0) = 0
    do p=1,npes-1
       displs(p) = displs(p-1) + numperproc(p-1)
    end do

  end subroutine compute_mpigs_patch

!===============================================================================

  subroutine compute_mpigs_land (nfact, numtot, numperproc, displs)

!------------------------------------------------------------------
!
! Purpose: Compute arguments for gatherv, scatterv for land vectors
!
!------------------------------------------------------------------

    implicit none

! ------------------- arguments -----------------------------------
    integer, intent(in ) :: nfact                ! multiplicative factor for patches
    integer, intent(out) :: numtot               ! total number of elements (to send or recv)
    integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
    integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!------------------------------------------------------------------

! ---------------------- local variables --------------------------
    integer :: p                                 ! index
!------------------------------------------------------------------

    numtot = (proc_landpts(iam))*nfact

    do p=0,npes-1
       numperproc(p) = proc_landpts(p)*nfact
    end do

    displs(0) = 0
    do p=1,npes-1
       displs(p) = displs(p-1) + numperproc(p-1)
    end do

  end subroutine compute_mpigs_land

!===============================================================================

#endif

end module spmdMod


