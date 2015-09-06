#include <misc.h>
#include <params.h>

module binary_io

  use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use pmgrid
#if ( defined SPMD )
  use spmd_dyn,  only: npes, compute_gsfactors, compute_gsfactors_dyn
  use mpi_gamil, only: nlat_each_proc
  use mpishorthand
#endif

contains

  subroutine wrtout_2D_array_dyn(iu, ibeg, iend, jbeg, jend, array)

    use mpi_gamil, only: gamil_gather_2D_array, is_rootproc

    implicit none

    real(r8), intent(in):: array(ibeg:iend,jbeg:jend)
    integer, intent(in) :: iu
    integer, intent(in) :: ibeg, iend, jbeg, jend

    integer ioerr
    real(r8) global_array(plond,plat)

    call gamil_gather_2D_array(array, ibeg, iend, jbeg, jend, global_array)

    if (is_rootproc) then
      write(iu, iostat=ioerr) global_array
      if (ioerr /= 0) then
        write(6, "('Error: binary_io::wrtout_2D_array_dyn: ioerror ', i, ' on i/o unit ', i)") ioerr, iu
        call endrun
      end if
    end if

  end subroutine wrtout_2D_array_dyn

  subroutine wrtout_3D_array_dyn(iu, ibeg, iend, jbeg, jend, kbeg, kend, array)

    use mpi_gamil, only: gamil_gather_3D_array, is_rootproc, myproc_id

    implicit none

    real(r8), intent(in):: array(ibeg:iend,jbeg:jend,kbeg:kend)
    integer, intent(in) :: iu
    integer, intent(in) :: ibeg, iend, jbeg, jend, kbeg, kend

    integer ioerr
    real(r8) global_array(plond,plat,kbeg:kend)

    call gamil_gather_3D_array(array, ibeg, iend, jbeg, jend, kbeg, kend, global_array)

    if (is_rootproc) then
      write(iu, iostat=ioerr) global_array
      if (ioerr /= 0) then
        write(6, "('Error: binary_io::wrtout_3D_array_dyn: ioerror ', i, ' on i/o unit ', i)") ioerr, iu
        call endrun
      end if
    end if

  end subroutine wrtout_3D_array_dyn

  subroutine readin_3D_array_dyn(iu, ibeg, iend, jbeg, jend, kbeg, kend, array, exchange)

    use mpi_gamil

    implicit none

    real(r8), intent(out):: array(ibeg:iend,jbeg:jend,kbeg:kend)
    integer, intent(in) :: iu
    integer, intent(in) :: ibeg, iend, jbeg, jend, kbeg, kend
    logical, intent(in) :: exchange

    integer ioerr
    real(r8) global_array(plond,plat,kbeg:kend)

    if (is_rootproc) then
      read(iu,iostat=ioerr) global_array
      if (ioerr /= 0) then
        write(6, "('Error: binary_io::readin_3D_array_dyn: ioerror ', i, ' on i/o unit ', i)") ioerr, iu
        call endrun
      end if
    end if

    call gamil_scatter_3D_array(global_array, ibeg, iend, jbeg, jend, kbeg, kend, array)

    call register_comm_array(ibeg,iend,jbeg,jend,kbeg,kend,1,1,array(:,jbeg,kbeg))
    if (exchange) then
      call gamil_arrays_comm(COMM_TO_LEFT,  1, array(:,jbeg,kbeg))
      call gamil_arrays_comm(COMM_TO_RIGHT, 1, array(:,jbeg,kbeg))
      call gamil_arrays_comm(COMM_TO_BOT,   1, array(:,jbeg,kbeg))
      call gamil_arrays_comm(COMM_TO_TOP,   1, array(:,jbeg,kbeg))
    else
      call gamil_arrays_comm(COMM_ROTATE_LEFT, 2, array(:,jbeg,kbeg))
    endif
    call remove_comm_array(array(:,jbeg,kbeg))

  end subroutine readin_3D_array_dyn

  subroutine wrtout_r8(iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: numperlat          ! Length of arr

        real(r8) :: arr(numperlat*numlats)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr

#if ( defined SPMD )
        real(r8), allocatable :: bufres(:)

        integer numrecv(0:npes-1)! number of items to be received
        integer displs(0:npes-1) ! displacement array
        integer numsend          ! number of bytes to send
        integer isiz             ! size of gather array to allocate
        integer p                ! process index
        integer ier              ! allocation status

        call compute_gsfactors(numperlat, numsend, numrecv, displs)

        if (masterproc) then
            isiz = 0
            do p = 0, npes-1
                isiz = isiz+numperlat*nlat_each_proc(p+1)
            end do
        else
            isiz = 1
        end if
        allocate (bufres(isiz), stat=ier)
        if (ier /= 0) call endrun

        call mpigatherv(arr, numsend, mpir8, bufres, numrecv, &
            displs, mpir8, 0, mpicom)

        if (masterproc) then
            write (iu, iostat=ioerr) bufres
            if (ioerr /= 0) then
                write (6, "('Error: binary_io::wrtout_r8: ioerror ', i, ' on i/o unit ', i)") ioerror, iu
                call endrun
            end if
            deallocate (bufres)
        end if

#else

        write (iu, iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6, "('Error: binary_io::wrtout_r8: ioerror ', i, ' on i/o unit ', i)") ioerror, iu
            call endrun
        end if

#endif

        return
    end subroutine wrtout_r8

    !###############################(wanhui)###############################
    !!(wanhui 2003.10.18)

    subroutine wrtout_r8_dyn (iu, arr, mx, my, mz)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        ! for the variables of LASG dynamical framework (wh)
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu         ! logical unit
        integer, intent(in) :: mx, my, mz ! shape of the array
        real(r8) :: arr(mx,my,mz)         ! array to be gathered (with the boundaries)

        integer ioerr

#if ( defined SPMD )
        real(r8), allocatable :: bufres(:)

        real(r8) :: arrtmp(mx,mz,my-2)                 ! Array to be gathered

        integer :: numrecv(0:npes-1)! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend          ! number of bytes to send
        integer :: isiz             ! size of gather array to allocate
        integer :: p                ! process index
        integer :: i,j,k            ! index
        integer :: ier              ! allocation status

        do k = 1, mz
            do j = 1, my-2
                do i = 1, mx
                    arrtmp(i,k,j) = arr(i,my-j,k)
                    !arrtmp(i,j,k) = arr(i,j+1,k)
                end do
            end do
        end do

        !call compute_gsfactors(numperlat, numsend, numrecv, displs)
        call compute_gsfactors(mx*mz, numsend, numrecv, displs)

        if (masterproc) then
            isiz = 0
            do p = 0, npes-1
                isiz = isiz+mx*mz*nlat_each_proc(p+1)
            end do
        else
            isiz = 1
        end if
        allocate (bufres(isiz), stat=ier)
        if (ier/=0) call endrun

        call mpigatherv (arrtmp, numsend, mpir8, bufres, numrecv, &
            displs, mpir8, 0, mpicom)

        if (masterproc) then
            write(iu, iostat=ioerr) bufres
            if (ioerr /= 0) then
                write(6, *) 'Error: binary_io::wrtout_r8_dyn: ioerror ', ioerr, ' on i/o unit = ', iu
                call endrun
            end if
            deallocate (bufres)
        end if
#else
        write(iu, iostat=ioerr) arr
        if (ioerr /= 0) then
            write(6, *) 'Error: binary_io::wrtout_r8_dyn: ioerror ', ioerr, ' on i/o unit = ', iu
            call endrun
        end if
#endif

        return
    end subroutine wrtout_r8_dyn

    !###############################(wanhui)###############################

    subroutine wrtout_r4 (iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: numperlat          ! Length of arr

        real(r4) :: arr(numperlat*numlats)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        real(r4), allocatable :: bufres(:)

        integer :: numrecv(0:npes-1)! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend
        integer :: isiz
        integer :: p

        call compute_gsfactors (numperlat, numsend, numrecv, displs)

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                isiz = isiz + numperlat*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
        end if

        call mpigatherv (arr, numsend, mpir4, bufres, numrecv, &
            displs, mpir4, 0, mpicom)

        if (masterproc) then
            write (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'WRTOUT_R4 ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
            deallocate (bufres)
        endif

#else

        write (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return
    end subroutine wrtout_r4

    !#######################################################################

    subroutine wrtout_int (iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: numperlat          ! number of values per latitude band
        integer arr(numperlat*numlats)            ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        integer, allocatable :: bufres(:)
        integer :: isiz
        integer :: p
        integer :: numsend          ! number of items to be sent
        integer :: displs(0:npes-1) ! displacement array
        integer :: numrecv(0:npes-1)! number of items to be received

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                isiz = isiz + numperlat*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
        end if

        call compute_gsfactors (numperlat, numsend, numrecv, displs)
        call mpigatherv (arr, numsend, mpiint, bufres, numrecv, &
            displs, mpiint, 0, mpicom)

        if (masterproc) then
            write (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'WRTOUT_INT ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
            deallocate (bufres)
        endif

#else

        write (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_INT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if
#endif
        return
    end subroutine wrtout_int

    !###############################(wanhui)###############################
    !!(wanhui 2004.02.29)

    subroutine wrtout_int_dyn (iu, arr, mx, my, mz)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        ! for the variables of LASG dynamical framework
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit

        integer, intent(in) :: mx,my,mz           ! shape of the array

        integer :: arr(mx,my,mz)                 ! Array to be gathered (with the boundaries)
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        integer, allocatable :: bufres(:)

        integer :: arrtmp(mx,mz,my-2)                 ! Array to be gathered

        integer :: numrecv(0:npes-1)! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend          ! number of bytes to send
        integer :: isiz             ! size of gather array to allocate
        integer :: p                ! process index
        integer :: i,j,k            ! index
        integer :: ier              ! allocation status

        do k=1,mz
            do j=1,my-2
                do i=1,mx
                    arrtmp(i,k,j) = arr(i,my-j,k)
                    !           arrtmp(i,j,k) = arr(i,j+1,k)
                enddo
            enddo
        enddo

        !!    call compute_gsfactors (numperlat, numsend, numrecv, displs)
        call compute_gsfactors (mx*mz,     numsend, numrecv, displs)

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                !!          isiz = isiz + numperlat*nlat_each_proc(p+1)
                isiz = isiz + mx*mz*nlat_each_proc(p+1)
            end do
        else
            isiz = 1
        end if
        allocate (bufres(isiz), stat=ier)
        if (ier/=0) call endrun

        call mpigatherv (arrtmp, numsend, mpiint, bufres, numrecv, &
            displs, mpiint, 0, mpicom)

        if (masterproc) then
            write (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'WRTOUT_INT_dyn ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
            deallocate (bufres)
        endif
#else

        write (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'WRTOUT_INT_dyn ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return

    end subroutine wrtout_int_dyn

    !###############################(wanhui)###############################


    !#######################################################################

    subroutine readin_r8 (iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to read binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: numperlat          ! Length of arr

        real(r8) :: arr(numperlat*numlats)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        real(r8), allocatable :: bufres(:)

        integer :: numrecv          ! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend(0:npes-1)
        integer :: isiz
        integer :: p

        call compute_gsfactors (numperlat, numrecv, numsend, displs)

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                isiz = isiz + numperlat*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
            read (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'READIN_R8 ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
        end if

        call mpiscatterv (bufres, numsend, displs, mpir8, arr, &
            numrecv, mpir8, 0, mpicom)

        if (masterproc) then
            deallocate (bufres)
        endif

#else

        read (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'READIN_R8 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return
    end subroutine readin_r8

    !#######################################################################
    !!#############################(wanhui)#################################

    !! subroutine readin_r8 (iu, arr, numperlat)
    !! subroutine readin_r8_dyn (iu, arr, numhor, numver)
    subroutine readin_r8_dyn (iu, arr, mx,my,mz)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to read binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: mx,my,mz


        !!    real(r8) :: arr(numperlat*numlats)        ! Array to be gathered

        real(r8) :: arr(mx,my,mz)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        real(r8), allocatable :: bufres(:)

        real(r8) :: arrtmp(mx,mz,my-2)        ! Array to be gathered

        integer :: numrecv          ! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend(0:npes-1)
        integer :: isiz
        integer :: p
        integer :: i,j,k

        !call compute_gsfactors(numperlat, numrecv, numsend, displs)
        call compute_gsfactors(mx*mz, numrecv, numsend, displs)

        if (masterproc) then
            isiz = 0
            do p = 0, npes-1
                !isiz = isiz+numperlat*nlat_each_proc(p+1)
                isiz = isiz+mx*mz*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
            read (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'READIN_R8_dyn ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
        end if

        call mpiscatterv(bufres, numsend, displs, mpir8, arrtmp, &
            numrecv, mpir8, 0, mpicom)

        do k=1,mz
            do j=1,my-2
                do i=1,mx
                    arr(i,my-j,k) = arrtmp(i,k,j)
                    !!           arr(i,j+1,k) = arrtmp(i,j,k)
                enddo
            enddo
        enddo


        if (masterproc) then
            deallocate (bufres)
        endif

#else

        read (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'READIN_R8_dyn ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return
    end subroutine readin_r8_dyn
    !!#############################(wanhui)#################################

    subroutine readin_r4 (iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to read binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: numperlat          ! Length of arr

        real(r4) :: arr(numperlat*numlats)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        real(r4), allocatable :: bufres(:)

        integer :: numrecv          ! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend(0:npes-1)
        integer :: isiz
        integer :: p

        call compute_gsfactors (numperlat, numrecv, numsend, displs)

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                isiz = isiz + numperlat*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
            read (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'READIN_R4 ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
        end if

        call mpiscatterv (bufres, numsend, displs, mpir4, arr, &
            numrecv, mpir4, 0, mpicom)

        if (masterproc) then
            deallocate (bufres)
        endif

#else

        read (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'READIN_R4 ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return
    end subroutine readin_r4

    !#######################################################################

    subroutine readin_int (iu, arr, numperlat)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to write restart binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer iu                 ! Logical unit
        integer numperlat          ! number of values per latitude band
        integer arr(numperlat*numlats)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        integer, allocatable :: bufres(:)
        integer isiz
        integer :: p
        integer :: numrecv          ! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend(0:npes-1)! number of items to be sent

        call compute_gsfactors (numperlat, numrecv, numsend, displs)

        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                isiz = isiz + numperlat*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
            read (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'READIN_INT ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
        end if

        call mpiscatterv (bufres, numsend, displs, mpiint, arr, &
            numrecv, mpiint, 0, mpicom)

        if (masterproc) then
            deallocate (bufres)
        endif

#else

        read (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'READIN_INT ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if
#endif
        return
    end subroutine readin_int

    !!#############################(wanhui 2004.02.29)#################################

    subroutine readin_int_dyn (iu, arr, mx,my,mz)
        !-----------------------------------------------------------------------
        !
        ! Purpose:
        ! Wrapper routine to read binary file
        !
        ! Author:
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments
        !
        integer, intent(in) :: iu                 ! Logical unit
        integer, intent(in) :: mx,my,mz


        !!    real(r8) :: arr(numperlat*numlats)        ! Array to be gathered

        integer :: arr(mx,my,mz)        ! Array to be gathered
        !
        ! Local workspace
        !
        integer ioerr              ! errorcode

#if ( defined SPMD )
        integer, allocatable :: bufres(:)

        integer :: arrtmp(mx,mz,my-2)        ! Array to be gathered

        integer :: numrecv          ! number of items to be received
        integer :: displs(0:npes-1) ! displacement array
        integer :: numsend(0:npes-1)
        integer :: isiz
        integer :: p
        integer :: i,j,k

        !!    call compute_gsfactors (numperlat, numrecv, numsend, displs)
        call compute_gsfactors (mx*mz, numrecv, numsend, displs)


        if (masterproc) then
            isiz = 0
            do p=0,npes-1
                !!          isiz = isiz + numperlat*nlat_each_proc(p+1)
                isiz = isiz + mx*mz*nlat_each_proc(p+1)
            end do
            allocate (bufres(isiz))
            read (iu,iostat=ioerr) bufres
            if (ioerr /= 0 ) then
                write (6,*) 'READIN_INT_dyn ioerror ',ioerr,' on i/o unit = ',iu
                call endrun
            end if
        end if

        call mpiscatterv (bufres, numsend, displs, mpiint, arrtmp, &
            numrecv, mpiint, 0, mpicom)

        do k=1,mz
            do j=1,my-2
                do i=1,mx
                    arr(i,my-j,k) = arrtmp(i,k,j)
                    !!           arr(i,j+1,k) = arrtmp(i,j,k)
                enddo
            enddo
        enddo


        if (masterproc) then
            deallocate (bufres)
        endif

#else

        read (iu,iostat=ioerr) arr
        if (ioerr /= 0 ) then
            write (6,*) 'READIN_INT_dyn ioerror ',ioerr,' on i/o unit = ',iu
            call endrun
        end if

#endif

        return
    end subroutine readin_int_dyn
    !!#############################(wanhui)#################################


end module binary_io
