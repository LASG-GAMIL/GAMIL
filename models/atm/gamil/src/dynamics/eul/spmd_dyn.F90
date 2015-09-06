#include <misc.h>
#include <params.h>

module spmd_dyn

!! (wanhui 2003.10.18)
!! (zx,wh  2003.10.21)
!! (wanhui 2003.10.21)
!! (wanhui 2003.10.24)
!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM.  Currently used for both dynamics
!          and physics, but ultimately the physics part should be broken off.
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, iam, beglatex, endlatex, numbnd, numlats, numlatsex, &
                           beglat, endlat, beglatexdyn, endlatexdyn, plev, plond !!(2003.10.24)
   use constituents, only: pcnst
   use mpishorthand, only: mpir8, mpicom
   use infnan, only: inf
   implicit none

   private
   public spmdinit_dyn, compute_gsfactors, compute_gsfactors_dyn
   save
   integer, public :: npes                 ! Total number of MPI tasks
   integer, private :: cut(2,0:plat-1)      ! partition for MPI tasks
   integer, private :: cutex(2,0:plat-1)    ! extended partition 
   integer, private :: neighs               ! number of south neighbors to comm guardcells
   integer, private :: neighn               ! number of north neighbors to comm guardcells
   integer, private :: npessp               ! number of MPI tasks in spectral space
   integer, private :: maxm                 ! max number of Fourier wavenumbers per MPI task
   integer, private :: numm(0:plat-1)       ! number of Fourier wavenumbers owned per task
   integer, private :: bsiz                 ! buffer size
   integer, private :: maxlats              ! max number of lats on any MPI task
   integer, private, allocatable :: nlat_p(:)    ! number of latitudes per MPI task

   real(r8), private, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs

CONTAINS

!========================================================================

   subroutine spmdinit_dyn (nproc)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processors
! 
! Method: Distribution is S->N for processors 0->npes
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!     use pspect, only: maxlats
!-----------------------------------------------------------------------
!
! Local workspace

#include <commpi.h>

      integer procid    ! processor id
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer isum      ! running total of work parcelled out
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer jendtmp,jendtmp1,jr,nproc
!
!-----------------------------------------------------------------------
!
! Allocate memory for number of lats per proc
!
     npes = nproc
     allocate (nlat_p (0:npes-1))

      jendtmp = PLAT/ npes
      jr = mod( PLAT,npes)
      if( iam .lt. jr ) then
        jendtmp = jendtmp + 1
      else
      endif
      numlatsex = endlatex - beglatex + 1
      numlats=jendtmp
      
      jendtmp1 = PLAT/ npes
      do procid=0,npes-1
        if (procid.lt.jr) then
          nlat_p(procid)=jendtmp1+1
        else
          nlat_p(procid)=jendtmp1
        endif
      end do
     
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

!!    call decomp_wavenumbers ()

    do procid=0,npes-1
    if (iam.eq.procid) then
    write(6,11) '*** p',iam,'beglat:endlat=',beglat,':',endlat,'numlats=',numlats
    write(6,11) '*** p',iam,'beglatex:endlatex=',beglatex,':',endlatex,'numlatsex=',numlatsex
    write(6,11) '*** p',iam,'beglatexdyn:endlatexdyn=',beglatexdyn,':',endlatexdyn
    write(6,11) '*** p',iam,'numbnd=',numbnd
11  format(1x,a6,i2,a25,i3,a2,i3,a12,i3)
    end if
    end do
    call mpibarrier(mpicom,ierr)
!!    stop

      return
   end subroutine spmdinit_dyn


!========================================================================
  

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = displs(p-1) + numperproc(p-1)
     end do
     
  end subroutine compute_gsfactors


!############################(wanhui)##################################

!!subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
  subroutine compute_gsfactors_dyn (numhor, numver, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
!!   integer, intent(in) :: numperlat    ! number of elements per latitude
     integer, intent(in) :: numhor
     integer, intent(in) :: numver
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
!!   numtot = numperlat*numlats
     numtot = numhor*numver
   
     do p=0,npes-1
!!      numperproc(p) = numperlat*nlat_p(p)
        numperproc(p) = plond*( nlat_p(p)+numlatsex-numlats )*numver
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = displs(p-1) + numperproc(p-1)
     end do
     
!!end subroutine compute_gsfactors
  end subroutine compute_gsfactors_dyn

!############################(wanhui)##################################

#endif

end module spmd_dyn
