#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use prognostics
   use ppgrid, only: pcols, pver
   use constituents, only: ppcnst

   implicit none

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: num     ! number of values
      integer :: i,j,k,m ! temporary indices
      real(r8), allocatable :: yz_tmp(:,:,:), q3_local(:,:,:,:)

      num = plond*plat
      call wrtout(nrg, strip2d, phis, num, 2)
!
! Write finite-volume dynamical core specific fields:
! [ (u3s,v3s), delp, pt, q3, ps ]
!
      num = plndlv*plat

      allocate( yz_tmp(plon,beglat:endlat,beglev:endlev) )

!
! TEMPORARY:  copy U3S to YZ_TMP
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            yz_tmp(i,j,k) = u3s(i,j,k)
         enddo
      enddo
   enddo
      call wrtout(nrg, strip3dxyz, yz_tmp, num, 3)

!
! TEMPORARY:  copy V3S to YZ_TMP
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            yz_tmp(i,j,k) = v3s(i,j,k)
         enddo
      enddo
   enddo
      call wrtout(nrg, strip3dxyz, yz_tmp, num, 3)

      call wrtout(nrg, strip3dxyz, delp, num, 3)


!
! TEMPORARY:  copy PT to YZ_TMP
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            yz_tmp(i,j,k) = pt(i,j,k)
         enddo
      enddo
   enddo

      call wrtout(nrg, strip3dxyz, yz_tmp, num, 3)
      deallocate( yz_tmp )

      num = plndlv*ppcnst*plat
      allocate( q3_local(plon,beglev:endlev,ppcnst,beglat:endlat) )

!$omp parallel do private(i,j,k,m)
      do m=1,ppcnst
         do k=beglev,endlev
            do j = beglat,endlat
               do i=1,plon
                  q3_local(i,k,m,j) = q3(i,j,k,m)
               enddo
             enddo
         enddo
      enddo
      call wrtout(nrg, strip3dq3old, q3_local, num, 3)
      deallocate( q3_local )

      num = plond*plat
      call wrtout(nrg, strip2d, ps, num, 2)
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if

      return
   end subroutine write_restart_dynamics

   subroutine wrtout(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only: decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y
      use parutilitiesmodule, only: commglobal, pargather 
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global length of array
#if defined( SPMD )
      real(r8) arr(*)            ! Array to be gathered
#else
      real(r8) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      real(r8), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
      if ( masterproc ) then
         allocate( bufres(lenarr) ) 
      else
         allocate( bufres(1) )
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call pargather( comm_y, 0, arr, decomp, bufres )
      else
         call pargather( commglobal, 0, arr, decomp, bufres )
      endif

! WS 01.01.03: This code is OK
      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      endif
      deallocate( bufres )
#else
      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'wrt ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine wrtout

   subroutine wrtouti(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to write restart binary file 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only: decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y
      use parutilitiesmodule, only: commglobal, pargather 
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global length of array
#if defined( SPMD )
      integer arr(*)            ! Array to be gathered
#else
      integer arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#if ( defined SPMD )
      integer, allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
      if ( masterproc ) then
          allocate( bufres(lenarr) ) 
      else
          allocate( bufres(1) )
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call pargather( comm_y, 0, arr, decomp, bufres )
      else
         call pargather( commglobal, 0, arr, decomp, bufres )
      endif

! WS 01.01.03: This code is OK
      if (masterproc) then
         write (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      endif
      deallocate( bufres )
#else
      write (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'wrt ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine wrtouti

   subroutine read_restart_dynamics (nrg)

      use dynamics_vars, only: dynamics_init
      use time_manager, only: get_step_size
#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
#include <comctl.h>
!
! Input arguments
!
      integer :: nrg     ! Unit number
!
! Local workspace
!
      integer :: ioerr   ! error status
      integer :: num     ! number of values
      real(r8) :: dtime  ! timestep size
      integer :: i,j,k,m ! temporary indices
      real(r8), allocatable :: yz_tmp(:,:,:), q3_local(:,:,:,:)

      dtime = get_step_size()
      call dynamics_init( dtime, iord, jord, nsplit )
!
! Initialize memory
!
      call initialize_prognostics

      num = plond*plat
      call lrreadin(nrg, strip2d, phis, num, 2)
!
! Read finite-volume dynamical core specific fields:
! [ (u3s,v3s), delp, pt, q3, ps ]
!
      num = plndlv*plat

      allocate( yz_tmp(plon,beglat:endlat,beglev:endlev) )
      call lrreadin(nrg, strip3dxyz, yz_tmp, num, 3) ! read in U3S
!
! TEMPORARY:  copy YZ_TMP to U3S
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            u3s(i,j,k) = yz_tmp(i,j,k)
         enddo
      enddo
   enddo
      call lrreadin(nrg, strip3dxyz, yz_tmp, num, 3) ! read in V3S
!
! TEMPORARY:  copy YZ_TMP to V3S
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            v3s(i,j,k) = yz_tmp(i,j,k)
         enddo
      enddo
   enddo

      call lrreadin(nrg, strip3dxyz, delp,num, 3)
      call lrreadin(nrg, strip3dxyz, yz_tmp,  num, 3)  ! read in PT
!
! TEMPORARY:  copy YZ_TMP to PT
!
!$omp parallel do private(i,j,k)
   do k=beglev,endlev
      do j = beglat,endlat
         do i=1,plon
            pt(i,j,k) = yz_tmp(i,j,k)
         enddo
      enddo
   enddo
      deallocate( yz_tmp )

      num = plndlv*ppcnst*plat

      allocate( q3_local(plon,beglev:endlev,ppcnst,beglat:endlat) )
      call lrreadin(nrg, strip3dq3old, q3_local, num, 3)

!$omp parallel do private(i,j,k,m)
      do m=1,ppcnst
         do k=beglev,endlev
            do j = beglat,endlat
               do i=1,plon
                  q3(i,j,k,m) = q3_local(i,k,m,j)
               enddo
            enddo
         enddo
      enddo
      deallocate( q3_local )

      num = plond*plat
      call lrreadin(nrg, strip2d, ps, num, 2)
!
! Read global integrals
!
      if (masterproc) then
         read(nrg, iostat=ioerr) tmass0
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun
         end if
      end if

#if ( defined SPMD )
      call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
#endif

      return
   end subroutine read_restart_dynamics

   subroutine lrreadin(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read real variable from restart binary file 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatter, parcollective, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r8) arr(*)            ! Array to be gathered
#else
      real(r8) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      real(r8), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatter( comm_y, 0, bufres, decomp, arr )
         call parcollective( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatter( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadin

   subroutine lrreadin4(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read real*4 variable from restart binary file 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatterreal4, parcollective1dreal4, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      real(r4) arr(*)            ! Array to be gathered
#else
      real(r4) arr(lenarr)       ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      real(r4), allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatterreal4( comm_y, 0, bufres, decomp, arr )
         call parcollective1dreal4( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatterreal4( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadin4

   subroutine lrreadini(iu, decomp, arr, lenarr, ndim)
!-----------------------------------------------------------------------
! Wrapper routine to read integer variable from restart binary file 
!-----------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use decompmodule, only : decomptype
#if ( defined SPMD )
      use spmd_dyn, only: comm_y, comm_z
      use parutilitiesmodule, only : commglobal, parscatter, parcollective, BCSTOP
#endif
!------------------------------Arguments--------------------------------
      integer iu                 ! Logical unit
      type (decomptype):: decomp ! Decomposition descriptor
      integer lenarr             ! Global size of array
#if defined( SPMD )
      integer arr(*)             ! Array to be gathered
#else
      integer arr(lenarr)        ! Array (SMP-only)
#endif
      integer ndim               ! dimensionality (2 or 3) of array
!---------------------------Local variables-----------------------------
      integer ioerr              ! errorcode
#ifdef SPMD
      integer, allocatable :: bufres(:) 
#endif
!-----------------------------------------------------------------------
#ifdef SPMD
      if (masterproc) then
         allocate (bufres(lenarr))
         read (iu,iostat=ioerr) bufres
         if (ioerr /= 0 ) then
            write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
         end if
      else
         allocate (bufres(1))
      endif
      if (ndim .eq. 2 .and. twod_decomp .eq. 1) then
         if (myid_z .eq. 0) call parscatter( comm_y, 0, bufres, decomp, arr )
         call parcollective( comm_z, BCSTOP, plon*(endlat-beglat+1), arr )
      else
         call parscatter( commglobal, 0, bufres, decomp, arr )
      endif
      deallocate (bufres)
#else 
      read (iu,iostat=ioerr) arr
      if (ioerr /= 0 ) then
         write (6,*) 'READIN ioerror ',ioerr,' on i/o unit = ',iu
         call endrun
      end if
#endif
      return
   end subroutine lrreadini

end module restart_dynamics
