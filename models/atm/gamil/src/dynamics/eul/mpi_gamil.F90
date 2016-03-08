module mpi_gamil

#include <mpif.h>
#include <params.h>
  
#define PACK_LEFT           0
#define PACK_RIGHT          1
#define PACK_TOP            2
#define PACK_BOT            3
#define UNPACK_NONE         -1
#define UNPACK_LEFT         4
#define UNPACK_RIGHT        5
#define UNPACK_TOP          6
#define UNPACK_BOT          7

  integer,parameter :: COMM_TO_LEFT = 0
  integer,parameter :: COMM_TO_RIGHT = 1
  integer,parameter :: COMM_TO_TOP = 2
  integer,parameter :: COMM_TO_BOT = 3
  integer,parameter :: COMM_ROTATE_LEFT = 4
  integer,parameter :: COMM_ROTATE_RIGHT = 5
  integer,parameter :: PACK_FOR_PHYS = 1
  integer,parameter :: PACK_FOR_1D = 2
  integer,parameter :: XLON = PLON+2
  integer,parameter :: YLAT = PLAT
  integer,parameter :: ZALT = PLEV+1

  integer, public :: ierr                      ! error code
  integer, public :: gamil_comm                ! the comm group of gamil
  integer, public :: lat_comm                  ! the comm group of latitude processes
  integer, public :: myproc_id                 ! the id of the current process
  integer, public :: my_latproc_id             ! the id of the current process in lat_comm
  integer, public :: rootproc_id               ! the id of the root process
  logical, public :: is_rootproc               ! is the current process the root process
  integer, public :: nproc                     ! the number of processes of gamil
  integer, public :: top_procid                ! the id of the process at top of the current process
  integer, public :: bot_procid                ! the id of the process at bottom of the current process
  integer, public :: left_procid               ! the id of the process at left of the current process
  integer, public :: right_procid              ! the id of the process at right of the current process
  logical, public :: has_leftBorder
  logical, public :: has_rightBorder
  logical, public :: has_topBorder
  logical, public :: has_botBorder
  integer, public :: beglonex, ibeg0, ibeg1, ibeg2
  integer, public :: endlonex, iend0, iend1, iend2
  integer, public :: jbeg0, jbeg1, jbeg2
  integer, public :: jend0, jend1, jend2
  integer, public :: num_x_proc
  integer, public :: num_y_proc
  integer, public, allocatable :: nlat_each_proc(:)
  integer, public, allocatable :: lat_local_proc_id(:,:)
  integer, private,allocatable :: disp_2D_gamil_comm(:), disp_phys_2D(:), disps_1(:), disp_2D_lat_comm(:), disp_2D_lat_comm_tmp(:)
  integer, private :: max_count_2D_lat_comm
  integer, allocatable :: counts_2D_gamil_comm(:), counts_1(:), counts_2D_lat_comm(:), counts_2D_lat_comm_tmp(:)
  integer, allocatable :: ibegs(:), jbegs(:), iends(:), jends(:)
  integer, allocatable :: lat_procs_ibegs(:), lat_procs_iends(:)
  integer, allocatable :: jbegs_phys(:), jends_phys(:)
  integer :: isend, irecv, status(MPI_STATUS_SIZE)
  real, target, allocatable :: gather_scatter_buf1_2D(:)
  real, target, allocatable :: gather_scatter_buf2_2D(:)
  real, target, allocatable :: send_buf_left(:)
  real, target, allocatable :: send_buf_right(:)
  real, target, allocatable :: send_buf_top(:)
  real, target, allocatable :: send_buf_bot(:)
  real, target, allocatable :: recv_buf_left(:)
  real, target, allocatable :: recv_buf_right(:)
  real, target, allocatable :: recv_buf_top(:)
  real, target, allocatable :: recv_buf_bot(:)
  real*16, allocatable :: reduce_buf(:)

  type, private :: comm_array_type
      integer num_dims
      integer ibeg, iend
      integer jbeg, jend
      integer kbeg, kend
      integer lbeg, lend
      real, pointer :: array(:)
      logical is_phys_array
  end type comm_array_type

  integer, private, parameter :: requests_size = 32
  type, private :: mpi_icomm_request
      logical used_mark
      integer num_requests
      integer nlevels
      integer requests(16)
      integer num_arrays
      integer comm_array_ids(64)
      real, allocatable :: send_buffer(:)
      real, allocatable :: recv_buffer(:)
      integer comm_direction                   ! 0->left, 1->right, 2->top, 3->bottom
  end type mpi_icomm_request

  type(mpi_icomm_request), private :: icomm_requests(requests_size)

  integer, private, parameter :: max_num_comm_arrays = 128
  integer, private :: num_registered_comm_arrays
  type(comm_array_type), private :: registered_comm_arrays(max_num_comm_arrays)

  REAL*16, allocatable, private :: reduce_buf_real16(:)
  integer,              private :: reduce_buf_real16_size

  save

   INTERFACE

   SUBROUTINE set_aff(k)
   INTEGER k
   END SUBROUTINE

   subroutine is_the_same_array(data1,data2,match_result)
     real data1, data2
     integer match_result
   end subroutine is_the_same_array
   END INTERFACE
  
contains

subroutine gamil_comm_init
    use omp_lib
    use mpishorthand, only: mpicom
    implicit none
    integer i,j

    rootproc_id = 0
    gamil_comm=mpicom
    call MPI_COMM_RANK(gamil_comm,myproc_id,ierr)
    call MPI_COMM_SIZE(gamil_comm,nproc,ierr)
    is_rootproc = .false.
    if (rootproc_id .eq. myproc_id) then
        is_rootproc = .true.
    endif

    allocate(disp_2D_gamil_comm(nproc), disp_phys_2D(nproc), counts_2D_gamil_comm(nproc), lat_local_proc_id(XLON,YLAT))
    allocate(disp_2D_lat_comm(nproc), disp_2D_lat_comm_tmp(nproc), counts_2D_lat_comm(nproc), counts_2D_lat_comm_tmp(nproc))
    allocate(ibegs(nproc), jbegs(nproc), iends(nproc), jends(nproc), nlat_each_proc(nproc))
    allocate(lat_procs_ibegs(nproc), lat_procs_iends(nproc))
    allocate(jbegs_phys(nproc), jends_phys(nproc))
    allocate(counts_1(nproc), disps_1(nproc))


end subroutine gamil_comm_init


subroutine gamil_2D_decomp
    use spmd_dyn
    use pmgrid, only: plat, beglat, endlat, beglatex, endlatex, beglatexdyn, endlatexdyn
    implicit none
    integer :: num_Xcell, num_Xrem
    integer :: num_Ycell, num_Yrem
    integer :: X_indx
    integer :: Y_indx
    integer :: i, j, k
    integer :: count, tmp
    integer :: num_bounds
    integer :: color
    integer :: num_comp_cost_lats
    REAL    :: pole_lat_cost_ratio
  
    if (nproc .ne. num_x_proc*num_y_proc) then
        write(6,*) 'The number of processes is not right', nproc, num_x_proc, num_y_proc
        call endrun
    endif

    pole_lat_cost_ratio = 0.0
    num_comp_cost_lats = YLAT+INT(FLOAT(YLAT)*pole_lat_cost_ratio)*2

    do i = 1, nproc
       disps_1(i) = i-1
       counts_1(i) = 1
    enddo

    num_bounds = 1
    X_indx = mod(myproc_id, num_x_proc)
    Y_indx = myproc_id / num_x_proc
    num_Xcell = XLON / num_x_proc
    num_Xrem = mod(XLON, num_x_proc)
    num_Ycell = num_comp_cost_lats / num_y_proc
    num_Yrem = mod(YLAT, num_y_proc)

    if (NUM_Xcell .lt. 1) then
        if (is_rootproc) then
           write(6,*) 'The number of processes on X direction is too much!'
        endif
        call endrun
    endif

    if (NUM_Ycell .lt. 1) then
        if (is_rootproc) then
           write(6,*) 'The number of processes on Y direction is too much!'
        endif
    endif
    
    ibeg0 = X_indx*num_Xcell+1
    if (X_indx .lt. num_Xrem) then
       ibeg0 = ibeg0+X_indx
       iend0 = ibeg0+num_Xcell
    else
       ibeg0 = ibeg0+num_Xrem
       iend0 = ibeg0+num_Xcell-1
    endif
    
    jbeg0 = Y_indx*num_Ycell+1
    if (Y_indx .lt. num_Yrem) then
       jbeg0 = jbeg0+Y_indx
       jend0 = jbeg0+num_Ycell
    else
       jbeg0 = jbeg0+num_Yrem
       jend0 = jbeg0+num_Ycell-1
    endif

    jbeg0 = jbeg0 - INT(FLOAT(YLAT)*pole_lat_cost_ratio)
    jend0 = jend0 - INT(FLOAT(YLAT)*pole_lat_cost_ratio)
    if (jend0 < 1 .or. jbeg0 > YLAT .or. pole_lat_cost_ratio < 0) then
        write(6,*) 'error for pole_lat_cost_ratio'
        call endrun
    end if

    if (jbeg0 < 1) jbeg0 = 1
    if (jend0 > YLAT) jend0 = YLAT

    if (X_indx .eq. 0) then
       left_procid = myproc_id+num_x_proc-1
       ibeg1 = ibeg0+1
       ibeg2 = ibeg0+2
    else
       left_procid = myproc_id-1
       ibeg1 = ibeg0
       ibeg2 = ibeg0
    endif
    if (X_indx .eq. num_x_proc - 1) then
       right_procid = myproc_id-num_x_proc+1
       iend1 = iend0-1
       iend2 = iend0-2
    else
       right_procid = myproc_id+1
       iend1 = iend0
       iend2 = iend0
    endif
    if (Y_indx .eq. 0) then
       top_procid = -1
       jbeg1 = jbeg0+1
       jbeg2 = jbeg0+2
    else
       top_procid = myproc_id-num_x_proc
       jbeg1 = jbeg0
       jbeg2 = jbeg0
    endif
    if (Y_indx .eq. num_y_proc-1) then
       bot_procid = -1
       jend1 = jend0-1
       jend2 = jend0-2
    else
       bot_procid = myproc_id+num_x_proc
       jend1 = jend0
       jend2 = jend0
    endif

    beglonex = ibeg0 - 1
    endlonex = iend0 + 1
    if (beglonex .lt. 1) beglonex = 1
    if (endlonex .gt. XLON) endlonex = XLON
    beglatex = jbeg0 - 1
    endlatex = jend0 + 1

    beglat = jbeg0
    endlat = jend0

    has_leftBorder = .false.
    has_rightBorder = .false.
    has_topBorder = .false.
    has_botBorder = .false.


    count = (iend0-ibeg0+1)*(jend0-jbeg0+1)

    call mpi_allgatherv(ibeg0,1,MPI_INTEGER,ibegs,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
    call mpi_allgatherv(jbeg0,1,MPI_INTEGER,jbegs_phys,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
    call mpi_allgatherv(iend0,1,MPI_INTEGER,iends,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
    call mpi_allgatherv(jend0,1,MPI_INTEGER,jends_phys,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
    call mpi_allgatherv(count,1,MPI_INTEGER,counts_2D_gamil_comm,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
!    call mpi_gatherv(count,1,MPI_INTEGER,counts_2D_gamil_comm,counts_2D_gamil_comm,disps_1,MPI_INTEGER,rootproc_id,gamil_comm,ierr)

    disp_2D_gamil_comm(nproc) = 0
    do i = nproc - 1, 1, -1
       disp_2D_gamil_comm(i) = disp_2D_gamil_comm(i+1)+counts_2D_gamil_comm(i+1)
    enddo

    disp_phys_2D(1) = 0
    do i = 2, nproc
       disp_phys_2D(i) = disp_phys_2D(i-1)+counts_2D_gamil_comm(i-1)
    enddo


    if (disp_2D_gamil_comm(1)+counts_2D_gamil_comm(1) .ne. XLON*YLAT) then
       write(6,*) 'error in decomp', disp_2D_gamil_comm(1), counts_2D_gamil_comm(1), XLON*YLAT
       call endrun
    endif

    do i = 1, nproc
       if (myproc_id .eq. i-1) then
          write(6,*) 'proc', i-1, 'ibeg', ibegs(i), 'iend', iends(i), &
                 'beglonex', beglonex, 'endlonex', endlonex, 'left_proc', left_procid, &
                 'right_proc', right_procid, 'jbeg', jbegs_phys(i), 'counts', counts_2D_gamil_comm(i), 'disp', disp_2D_gamil_comm(i)
       endif
    enddo    
    
    allocate(gather_scatter_buf1_2D(XLON*YLAT*32), gather_scatter_buf2_2D(XLON*YLAT*32))
    allocate(send_buf_left((jend0-jbeg0+1)*2*ZALT*16))
    allocate(send_buf_right((jend0-jbeg0+1)*2*ZALT*16))
    allocate(recv_buf_left((jend0-jbeg0+1)*2*ZALT*16))
    allocate(recv_buf_right((jend0-jbeg0+1)*2*ZALT*16))
    allocate(send_buf_top((iend0-ibeg0+1)*2*ZALT*16))
    allocate(send_buf_bot((iend0-ibeg0+1)*2*ZALT*16))
    allocate(recv_buf_top((iend0-ibeg0+1)*2*ZALT*16))
    allocate(recv_buf_bot((iend0-ibeg0+1)*2*ZALT*16))
    allocate(reduce_buf(nproc))

    do k = 1, nproc
!        write(6, *) 'lat intervel', jbegs(k), jends(k)
        nlat_each_proc(k) = jends_phys(k) - jbegs_phys(k) + 1
        do j = jbegs_phys(k), jends_phys(k)
        do i = ibegs(k), iends(k)
            lat_local_proc_id(i,j) = k - 1
        enddo
        enddo
    enddo

!    do i = 1, YLAT
!           write(6,*) 'lat proc map', i, lat_local_proc_id(i)
!    enddo

    
!    write(6,*) 'proc physics', myproc_id, ibeg0, iend0, jbeg0, jend0, jbeg1, jend1, jbeg2, jend2
!    write(6,*) 'proc physics', myproc_id, left_procid, right_procid, top_procid, bot_procid

    tmp = jbeg0
    jbeg0 = plat + 1 - jend0
    jend0 = plat + 1 - tmp
    tmp = jbeg1
    jbeg1 = plat + 1 - jend1
    jend1 = plat + 1 - tmp
    tmp = jbeg2
    jbeg2 = plat + 1 - jend2
    jend2 = plat + 1 - tmp    
    tmp = top_procid
    top_procid = bot_procid
    bot_procid = tmp
    beglatexdyn = jbeg0 - 1
    endlatexdyn = jend0 + 1

    has_topBorder = .false.
    has_botBorder = .false.
    if (top_procid .eq. -1) has_topBorder = .true.
    if (bot_procid .eq. -1) has_botBorder = .true.    

    call mpi_allgatherv(jend0,1,MPI_INTEGER,jends,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)
    call mpi_allgatherv(jbeg0,1,MPI_INTEGER,jbegs,counts_1,disps_1,MPI_INTEGER,gamil_comm,ierr)

    color = myproc_id / num_x_proc
    call MPI_COMM_SPLIT (gamil_comm, color, myproc_id, lat_comm, ierr)
    call MPI_COMM_RANK(lat_comm,my_latproc_id,ierr)
    write(6,*) 'split ', color, myproc_id, lat_comm, bot_procid, top_procid, my_latproc_id

    count = (iend0-ibeg0+1)*(jend0-jbeg0+1)
    call mpi_allgatherv(ibeg0,1,MPI_INTEGER,lat_procs_ibegs,counts_1,disps_1,MPI_INTEGER,lat_comm,ierr)
    call mpi_allgatherv(iend0,1,MPI_INTEGER,lat_procs_iends,counts_1,disps_1,MPI_INTEGER,lat_comm,ierr)
    call mpi_allgatherv(count,1,MPI_INTEGER,counts_2D_lat_comm,counts_1,disps_1,MPI_INTEGER,lat_comm,ierr)
    disp_2D_lat_comm(1) = 0
    do i = 2, num_x_proc
       disp_2D_lat_comm(i) = disp_2D_lat_comm(i-1)+counts_2D_lat_comm(i-1)
    enddo
    max_count_2D_lat_comm = -1
    do i = 1, num_x_proc
       if (max_count_2D_lat_comm < counts_2D_lat_comm(i)) then
           max_count_2D_lat_comm = counts_2D_lat_comm(i)
       endif
    enddo

    
    write(6,*) 'proc dynamics', myproc_id, ibeg0, iend0, jbeg0, jend0, jbeg1, jend1, jbeg2, jend2
    write(6,*) 'proc dynamics', myproc_id, left_procid, right_procid, top_procid, bot_procid
    do i = 1, nproc
       if (is_rootproc) then
          write(6,*) 'ibeg', ibegs(i), 'jbeg', jbegs(i), 'counts', counts_2D_gamil_comm(i), 'disp', disp_2D_gamil_comm(i)
       endif
    enddo   

    do i=1, requests_size
      icomm_requests(i)%used_mark = .false.
      icomm_requests(i)%num_requests = 0
      allocate(icomm_requests(i)%send_buffer((jend0-jbeg0+1+iend0-ibeg0+1)*2*ZALT*16))
      allocate(icomm_requests(i)%recv_buffer((jend0-jbeg0+1+iend0-ibeg0+1)*2*ZALT*16))
    end do

    num_registered_comm_arrays = 0

   reduce_buf_real16_size = 1
   allocate(reduce_buf_real16(reduce_buf_real16_size))

end subroutine gamil_2D_decomp


 SUBROUTINE allreduce_real16(input_data, output_data, num_data, comm, num_proc)
   implicit none
   real*16           :: input_data(:), output_data(:)
   integer           :: num_data, comm
   integer           :: ierr
   integer,optional  :: num_proc
   integer           :: local_num_proc, i, k


#if ( defined  NO_MPI_REAL16 )
   if (present(num_proc)) then
       local_num_proc = num_proc
   else
       call MPI_COMM_SIZE(comm, local_num_proc, ierr)
   end if
   if (reduce_buf_real16_size < num_data*local_num_proc) then
      reduce_buf_real16_size = num_data*local_num_proc*2
      deallocate(reduce_buf_real16)
      allocate(reduce_buf_real16(reduce_buf_real16_size))
   end if
   call mpi_allgather(input_data,num_data*2,MPI_REAL8,reduce_buf_real16,num_data*2,MPI_REAL8,comm,ierr)
   do k=1,num_data
       output_data(k) = 0.0
       do i=1,local_num_proc
          output_data(k) = output_data(k) + reduce_buf_real16((i-1)*num_data+k)
       end do
    enddo
#else
   call mpi_allreduce(input_data,output_data,num_data,MPI_REAL16,MPI_SUM,comm,ierr)
#endif
 END SUBROUTINE allreduce_real16
 


integer function get_free_icomm_request()
    implicit none
    integer :: i

    do i=1, requests_size
      if (.not.(icomm_requests(i)%used_mark)) then
          icomm_requests(i)%used_mark = .true.
          icomm_requests(i)%num_requests = 0
          icomm_requests(i)%num_arrays = 0
          get_free_icomm_request = i
          return
      end if 
    end do
    call endrun()
end function get_free_icomm_request


subroutine register_comm_array(ibeg, iend, jbeg, jend, kbeg, kend, lbeg, lend, array, is_phys_array)
    implicit none
    integer :: ibeg, iend, jbeg, jend, kbeg, kend, lbeg, lend
    logical,optional :: is_phys_array
    real,target,intent(in):: array(:)
    integer i, match_result

    num_registered_comm_arrays = num_registered_comm_arrays + 1
    if (num_registered_comm_arrays .gt. max_num_comm_arrays) then
      write(6,*) 'the number of registered comm arrays exceeds maximum number'
      call endrun
    end if

    registered_comm_arrays(num_registered_comm_arrays)%ibeg = ibeg
    registered_comm_arrays(num_registered_comm_arrays)%iend = iend
    registered_comm_arrays(num_registered_comm_arrays)%jbeg = jbeg
    registered_comm_arrays(num_registered_comm_arrays)%jend = jend
    registered_comm_arrays(num_registered_comm_arrays)%kbeg = kbeg
    registered_comm_arrays(num_registered_comm_arrays)%kend = kend
    registered_comm_arrays(num_registered_comm_arrays)%lbeg = lbeg
    registered_comm_arrays(num_registered_comm_arrays)%lend = lend
    registered_comm_arrays(num_registered_comm_arrays)%is_phys_array = .false.
    if (present(is_phys_array) .and. is_phys_array) then
      registered_comm_arrays(num_registered_comm_arrays)%is_phys_array = .true.
    end if
    registered_comm_arrays(num_registered_comm_arrays)%array => array

end subroutine register_comm_array


subroutine remove_comm_array(array)
    implicit none
    real,intent(in):: array(:)
    integer match_result

    call is_the_same_array(array(1),registered_comm_arrays(num_registered_comm_arrays)%array(1),match_result)
    if (num_registered_comm_arrays .eq. 0 .or. match_result .eq. 0) then
      write(6,*) 'can not remove comm array'
      call endrun
    end if
    num_registered_comm_arrays = num_registered_comm_arrays - 1
end subroutine remove_comm_array


subroutine add_comm_array(rid, array)
    implicit none
    integer,intent(in) :: rid
    real,target,intent(in):: array(:)
    integer i, comm_array_id, match_result

    call t_startf("dyn search array")
    comm_array_id = 0
    do i = num_registered_comm_arrays, 1, -1
      call is_the_same_array(array(1),registered_comm_arrays(i)%array(1),match_result)
      if (match_result .eq. 1) then
        comm_array_id = i
        goto 20
      end if
    end do

20 if (comm_array_id .eq. 0) then
     if (is_rootproc) then
       write(6,*) 'encounter unregistered comm array'
       call endrun
     end if
    end if
    call t_stopf("dyn search array")

    icomm_requests(rid)%num_arrays = icomm_requests(rid)%num_arrays + 1
    icomm_requests(rid)%comm_array_ids(icomm_requests(rid)%num_arrays) = comm_array_id
end subroutine add_comm_array


subroutine wait_icomm_request(id)
    implicit none
    integer :: id
    integer :: i
    
    if (.not. (icomm_requests(id)%used_mark)) call endrun

    call t_startf("dyn sendrecv")
    do i=1, icomm_requests(id)%num_requests
       call MPI_WAIT(icomm_requests(id)%requests(i),status,ierr)                         
    end do
    call t_stopf("dyn sendrecv")

    call t_startf("dyn pack/unpack")
    call unpack_request_comm_arrays(id)
    call t_stopf("dyn pack/unpack")

    icomm_requests(id)%used_mark = .false.
    icomm_requests(id)%num_requests = 0
end subroutine wait_icomm_request


subroutine gamil_scatter_2D_array_phys(array_g, ibeg, iend, jbeg,jend, array)
    implicit none
    real :: array(ibeg:iend,jbeg:jend)
    real :: array_g(XLON,YLAT)
    integer :: ibeg,iend,jbeg,jend
    integer :: i, j, k
    integer :: count

    if (is_rootproc) then
       do k=1, nproc
          do j=jbegs_phys(k), jends_phys(k)
             do i=ibegs(k), iends(k)
                gather_scatter_buf1_2D((j - jbegs_phys(k))*(iends(k)-ibegs(k)+1)+i-ibegs(k)+1+disp_phys_2D(k)) = array_g(i,j)
             enddo
          enddo
       enddo
    endif

    count = (iend0-ibeg0+1)*(jends_phys(myproc_id+1)-jbegs_phys(myproc_id+1)+1)
    call mpi_scatterv(gather_scatter_buf1_2D,counts_2D_gamil_comm,disp_phys_2D,MPI_REAL8,gather_scatter_buf2_2D, &
                      count,MPI_REAL8,rootproc_id,gamil_comm,ierr)    

    
    do j = jbegs_phys(myproc_id+1), jends_phys(myproc_id+1)
       do i = ibeg0, iend0
           array(i,j) = gather_scatter_buf2_2D((j - jbegs_phys(myproc_id+1))*(iend0-ibeg0+1)+i-ibeg0+1)
       end do
    enddo

end subroutine gamil_scatter_2D_array_phys


subroutine gamil_scatter_3D_array_phys(array_g, ibeg, iend, jbeg,jend, nlevels, array)
    implicit none
    real :: array(ibeg:iend,1:nlevels,jbeg:jend)
    real :: array_g(XLON,1:nlevels,YLAT)
    real,allocatable :: array_2d_tmp_input(:,:)
    real,allocatable :: array_2d_tmp_output(:,:)
    integer :: ibeg,iend,jbeg,jend, nlevels
    integer :: i, j, k

    allocate(array_2d_tmp_input(XLON,YLAT),array_2d_tmp_output(XLON,YLAT))

    do k = 1, nlevels
      do j = 1, YLAT
      do i = 1, XLON
         array_2d_tmp_input(i,j) = array_g(i,k,j)
      enddo
      enddo

      call gamil_scatter_2D_array_phys(array_2d_tmp_input, 1, XLON, 1, YLAT, array_2d_tmp_output)
      
      do j = jbegs_phys(myproc_id+1), jends_phys(myproc_id+1)
      do i = ibeg0, iend0
         array(i,k,j) = array_2d_tmp_output(i,j)
      enddo
      enddo
    enddo

    deallocate(array_2d_tmp_input)
    deallocate(array_2d_tmp_output)

end subroutine gamil_scatter_3D_array_phys


subroutine gamil_gather_2D_array_phys(array_g, ibeg, iend, jbeg,jend, array)
    implicit none
    real :: array(ibeg:iend,jbeg:jend)
    real :: array_g(XLON-2,YLAT)
    real :: array_g_tmp(XLON,YLAT)
    integer :: ibeg,iend,jbeg,jend
    integer :: i, j, k
    integer :: count

    do j = jbegs_phys(myproc_id+1), jends_phys(myproc_id+1)
       do i = ibeg0, iend0
           gather_scatter_buf1_2D((j - jbegs_phys(myproc_id+1))*(iend0-ibeg0+1)+i-ibeg0+1) = array(i,j)
       end do
    end do

    count = (iend0-ibeg0+1)*(jends_phys(myproc_id+1)-jbegs_phys(myproc_id+1)+1)
    call mpi_gatherv(gather_scatter_buf1_2D,count,MPI_REAL8,gather_scatter_buf2_2D, &
                     counts_2D_gamil_comm,disp_phys_2D,MPI_REAL8,rootproc_id,gamil_comm,ierr)    
    if (is_rootproc) then
       do k=1, nproc
          do j=jbegs_phys(k), jends_phys(k)
             do i=ibegs(k), iends(k)
                array_g_tmp(i,j) = gather_scatter_buf2_2D((j - jbegs_phys(k))*(iends(k)-ibegs(k)+1)+i-ibegs(k)+1+disp_phys_2D(k))
             enddo
          enddo
       enddo
       do j = 1, YLAT
       do i = 2, XLON-1
          array_g(i-1,j)=array_g_tmp(i,j)
       enddo
       enddo
    endif

end subroutine gamil_gather_2D_array_phys


subroutine gamil_gather_3D_array_phys(array_g, ibeg, iend, jbeg,jend, nlevels, array)
    implicit none
    real :: array(ibeg:iend,1:nlevels,jbeg:jend)
    real :: array_g(XLON-2,1:nlevels,YLAT)
    real,allocatable :: array_2d_tmp_input(:,:)
    real,allocatable :: array_2d_tmp_output(:,:)
    integer :: ibeg,iend,jbeg,jend, nlevels
    integer :: i, j, k

    allocate(array_2d_tmp_input(XLON-2,YLAT),array_2d_tmp_output(XLON,YLAT))

    do k = 1, nlevels

      do j = jbegs_phys(myproc_id+1), jends_phys(myproc_id+1)
      do i = ibeg0, iend0
         array_2d_tmp_output(i,j) = array(i,k,j)
      enddo
      enddo

      call gamil_gather_2D_array_phys(array_2d_tmp_input, 1, XLON, 1, YLAT, array_2d_tmp_output)
      
      if (is_rootproc) then
         do j = 1, YLAT
         do i = 1, XLON-2
            array_g(i,k,j) = array_2d_tmp_input(i,j)
         enddo
         enddo
      endif

    enddo

    deallocate(array_2d_tmp_input)
    deallocate(array_2d_tmp_output)

end subroutine gamil_gather_3D_array_phys


subroutine gamil_scatter_2D_array(array_g, ibeg, iend, jbeg,jend, array)
    implicit none
    real :: array(ibeg:iend,jbeg:jend)
    real :: array_g(XLON,YLAT)
    integer :: ibeg,iend,jbeg,jend
    integer :: i, j, k
    integer :: count

    if (is_rootproc) then
       do k=1, nproc
          do j=jbegs(k), jends(k)
             do i=ibegs(k), iends(k)
                gather_scatter_buf1_2D((j - jbegs(k))*(iends(k)-ibegs(k)+1)+i-ibegs(k)+1+disp_2D_gamil_comm(k)) = array_g(i,j)
             enddo
          enddo
       enddo
    endif

    count = (iend0-ibeg0+1)*(jends(myproc_id+1)-jbegs(myproc_id+1)+1)
    call mpi_scatterv(gather_scatter_buf1_2D,counts_2D_gamil_comm,disp_2D_gamil_comm,MPI_REAL8,gather_scatter_buf2_2D, &
                      count,MPI_REAL8,rootproc_id,gamil_comm,ierr)    

    
    do j = jbegs(myproc_id+1), jends(myproc_id+1)
       do i = ibeg0, iend0
           array(i,j) = gather_scatter_buf2_2D((j - jbegs(myproc_id+1))*(iend0-ibeg0+1)+i-ibeg0+1)
       end do
    enddo

end subroutine gamil_scatter_2D_array


subroutine gamil_scatter_3D_array(array_g, ibeg, iend, jbeg, jend, kbeg, kend, array)
    implicit none
    real :: array(ibeg:iend,jbeg:jend,kbeg:kend)
    real :: array_g(XLON,YLAT,kbeg:kend)
    integer :: ibeg,iend,jbeg,jend,kbeg,kend
    integer :: k

    do k=kbeg,kend
       call gamil_scatter_2D_array(array_g(1,1,k),ibeg,iend,jbeg,jend,array(ibeg,jbeg,k))
    enddo
end subroutine gamil_scatter_3D_array


subroutine gamil_gather_3D_array(array, ibeg, iend, jbeg, jend, kbeg, kend, array_g)
    implicit none
    real :: array(ibeg:iend,jbeg:jend,kbeg:kend)
    real :: array_g(XLON,YLAT,kbeg:kend)
    integer :: ibeg,iend,jbeg,jend,kbeg,kend
    integer :: k

    do k=kbeg,kend
       call gamil_gather_2D_array(array(ibeg,jbeg,k), ibeg, iend, jbeg, jend, array_g(1,1,k))
    enddo
    
end subroutine gamil_gather_3D_array


subroutine gamil_gather_2D_array(array, ibeg, iend, jbeg, jend, array_g)
    implicit none
    real :: array(ibeg:iend,jbeg:jend)
    real :: array_g(XLON,YLAT)
    integer :: ibeg,iend,jbeg,jend
    integer :: i, j, k
    integer :: count

    do j = jbeg0, jend0
       do i = ibeg0, iend0
           gather_scatter_buf1_2D((j - jbeg0)*(iend0-ibeg0+1)+i-ibeg0+1) = array(i,j)
       end do
    enddo

    count = (iend0-ibeg0+1)*(jend0-jbeg0+1)
    call mpi_gatherv(gather_scatter_buf1_2D,count,MPI_REAL8,gather_scatter_buf2_2D, &
                     counts_2D_gamil_comm,disp_2D_gamil_comm,MPI_REAL8,rootproc_id,gamil_comm,ierr)    

    if (is_rootproc) then
       do k=1, nproc
          do j=jbegs(k), jends(k)
             do i=ibegs(k), iends(k)
                array_g(i,j) = gather_scatter_buf2_2D((j - jbegs(k))*(iends(k)-ibegs(k)+1)+i-ibegs(k)+1+disp_2D_gamil_comm(k))
             enddo
          enddo
       enddo
    endif
end subroutine gamil_gather_2D_array


subroutine gamil_3D_array_pack(pack_direction,ibeg,iend,jbeg,jend,kbeg,kend,iter,nlevels,arr,pack_mode,buffer)
    implicit none
    integer :: pack_direction
    integer :: ibeg, iend, jbeg, jend,kbeg,kend
    integer :: nlevels
    integer :: i, j, k
    real :: arr(ibeg:iend,jbeg:jend,kbeg:kend)
    integer :: iter, local_iter
    integer :: local_jbeg, local_jend
    integer :: local_ibeg, local_iend
    integer :: pack_mode
    integer :: unit_size
    real :: buffer(:)

    if (pack_direction .eq. PACK_LEFT) then
      local_jbeg = jbeg0
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)-1
        local_jend = jends_phys(myproc_id+1)+1
      endif
      if (jbeg .eq. jend) then
        local_jbeg = jbeg
        local_jend = jend
      end if
      local_ibeg = beglonex+2-nlevels
      local_iend = beglonex+1
      if (ibeg .eq. iend) then
        write(6,*) 'can not use left right communication'
        call endrun
      end if
    else if (pack_direction .eq. PACK_RIGHT) then
      local_jbeg = jbeg0
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)-1
        local_jend = jends_phys(myproc_id+1)+1
      endif
      if (jbeg .eq. jend) then
        local_jbeg = jbeg
        local_jend = jend
      end if
      local_ibeg = endlonex-nlevels
      local_iend = endlonex-1
      if (ibeg .eq. iend) then
        write(6,*) 'can not use left right communication'
        call endrun
      end if
    else if (pack_direction .eq. PACK_TOP) then
      local_jbeg = jbeg0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)
      endif
      local_jend = local_jbeg+nlevels-1
      local_ibeg = ibeg0
      local_iend = iend0
      if (ibeg0 .gt. 1) local_ibeg = local_ibeg-1
      if (iend0 .lt. XLON) local_iend = local_iend+1
      if (jbeg .eq. jend) then
        write(6,*) 'can not use top bottom communication'
        call endrun
      end if
      if (ibeg .eq. iend) then
        local_ibeg = ibeg
        local_iend = iend
      end if
    else if (pack_direction .eq. PACK_BOT) then
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jend = jends_phys(myproc_id+1)
      endif
      local_jbeg = local_jend-nlevels+1
      local_ibeg = ibeg0
      local_iend = iend0
      if (ibeg0 .gt. 1) local_ibeg = local_ibeg-1
      if (iend0 .lt. XLON) local_iend = local_iend+1
      if (jbeg .eq. jend) then
        write(6,*) 'can not use top bottom communication'
        call endrun
      end if
      if (ibeg .eq. iend) then
        local_ibeg = ibeg
        local_iend = iend
      end if
    else
      write(6,*) 'error pack direction'
      call endrun
    end if
    
    unit_size = (local_jend-local_jbeg+1)*(local_iend-local_ibeg+1)
    if (unit_size .ge. 64) then
!$OMP PARALLEL DO PRIVATE (I,J,K,local_iter)
      do k=kbeg, kend
        local_iter = iter+unit_size*(k-kbeg)
        do j=local_jbeg, local_jend
          do i=local_ibeg, local_iend
            buffer(local_iter)=arr(i,j,k)
            local_iter = local_iter+1
          enddo
        enddo    
      enddo
    else
      do k=kbeg, kend
        local_iter = iter+unit_size*(k-kbeg)
        do j=local_jbeg, local_jend
          do i=local_ibeg, local_iend
            buffer(local_iter)=arr(i,j,k)
            local_iter = local_iter+1
          enddo
        enddo    
      enddo
    end if
    iter = iter+(local_jend-local_jbeg+1)*(local_iend-local_ibeg+1)*(kend+1-kbeg)

end subroutine gamil_3D_array_pack


subroutine gamil_3D_array_unpack(unpack_direction,ibeg,iend,jbeg,jend,kbeg,kend,iter,nlevels,arr,pack_mode,buffer)
    implicit none
    integer :: ibeg, iend, jbeg, jend, kbeg, kend
    integer :: nlevels
    integer :: i, j, k
    real :: arr(ibeg:iend,jbeg:jend,kbeg:kend)
    integer :: iter, local_iter
    integer :: local_jbeg, local_jend
    integer :: local_ibeg, local_iend
    integer :: pack_mode
    real :: buffer(:)
    integer :: unpack_direction
    integer :: unit_size

    if (unpack_direction .eq. UNPACK_LEFT) then
      local_jbeg = jbeg0
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)-1
        local_jend = jends_phys(myproc_id+1)+1
      endif
      if (jbeg .eq. jend) then
        local_jbeg = jbeg
        local_jend = jend
      end if
      local_ibeg = beglonex
      local_iend = beglonex+nlevels-1
    else if (unpack_direction .eq. UNPACK_RIGHT) then
      local_jbeg = jbeg0
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)-1
        local_jend = jends_phys(myproc_id+1)+1
      end if
      if (jbeg .eq. jend) then
        local_jbeg = jbeg
        local_jend = jend
      end if
      local_ibeg = endlonex-nlevels+1
      local_iend = endlonex
    else if (unpack_direction .eq. UNPACK_TOP) then
      local_jbeg = jbeg0
      if (pack_mode .eq. 1) then
        local_jbeg = jbegs_phys(myproc_id+1)
      endif
      local_jend = local_jbeg-1
      local_jbeg = local_jbeg-nlevels
      local_ibeg = ibeg0
      local_iend = iend0
      if (ibeg0 .gt. 1) local_ibeg = local_ibeg-1
      if (iend0 .lt. XLON) local_iend = local_iend+1
      if (ibeg .eq. iend) then
        local_ibeg = ibeg
        local_iend = iend
      end if
    else if (unpack_direction .eq. UNPACK_BOT) then
      local_jend = jend0
      if (pack_mode .eq. 1) then
        local_jend = jends_phys(myproc_id+1)
      endif
      local_jbeg = local_jend+1
      local_jend = local_jend+nlevels
      local_ibeg = ibeg0
      local_iend = iend0
      if (ibeg0 .gt. 1) local_ibeg = local_ibeg-1
      if (iend0 .lt. XLON) local_iend = local_iend+1
      if (ibeg .eq. iend) then
        local_ibeg = ibeg
        local_iend = iend
      end if
    else
      write(6,*) 'error unpack direction'
      call endrun
    end if
 
    unit_size = (local_jend-local_jbeg+1)*(local_iend-local_ibeg+1)
    if (unit_size .ge. 64) then
!$OMP PARALLEL DO PRIVATE (I,J,K,local_iter)
      do k=kbeg, kend
        local_iter = iter+unit_size*(k-kbeg)
        do j=local_jbeg, local_jend
          do i=local_ibeg, local_iend
            arr(i,j,k)=buffer(local_iter)
            local_iter = local_iter+1
          enddo
        enddo    
      enddo
    else
      do k=kbeg, kend
        local_iter = iter+unit_size*(k-kbeg)
        do j=local_jbeg, local_jend
          do i=local_ibeg, local_iend
            arr(i,j,k)=buffer(local_iter)
            local_iter = local_iter+1
          enddo
        enddo    
      enddo
    endif
    iter = iter+(local_jend-local_jbeg+1)*(local_iend-local_ibeg+1)*(kend+1-kbeg)

end subroutine gamil_3D_array_unpack


subroutine pack_request_comm_arrays(rid,iter)
    implicit none
    integer :: rid
    integer :: i, iter
    integer :: pack_direction
    integer :: comm_array_id, pack_mode

    call t_startf("dyn pack/unpack")

    pack_direction = -1
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_LEFT .or. &
      icomm_requests(rid)%comm_direction .eq. COMM_ROTATE_LEFT) then
      pack_direction = PACK_LEFT
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_RIGHT) then
      pack_direction = PACK_RIGHT
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_TOP) then
      pack_direction = PACK_TOP
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_BOT) then
      pack_direction = PACK_BOT
    end if

    iter = 1
    do i = 1, icomm_requests(rid)%num_arrays
      comm_array_id = icomm_requests(rid)%comm_array_ids(i)
      pack_mode = 0
      if (registered_comm_arrays(comm_array_id)%is_phys_array) pack_mode = 1
      call gamil_3D_array_pack(pack_direction,                                   &
                               registered_comm_arrays(comm_array_id)%ibeg,       &
                               registered_comm_arrays(comm_array_id)%iend,       &
                               registered_comm_arrays(comm_array_id)%jbeg,       &
                               registered_comm_arrays(comm_array_id)%jend,       &
                               registered_comm_arrays(comm_array_id)%kbeg,       &
                               registered_comm_arrays(comm_array_id)%kend,       &
                               iter,                                             &
                               icomm_requests(rid)%nlevels,                      &
                               registered_comm_arrays(comm_array_id)%array(:),   &
                               pack_mode,                                        &
                               icomm_requests(rid)%send_buffer) 
    end do

    call t_stopf("dyn pack/unpack")
end subroutine pack_request_comm_arrays


subroutine unpack_request_comm_arrays(rid)
    implicit none
    integer :: rid
    integer :: i, iter
    integer :: unpack_direction
    integer :: comm_array_id, pack_mode

    iter = 1
    unpack_direction = UNPACK_NONE
    if ((endlonex .eq. XLON .and. icomm_requests(rid)%comm_direction .eq. COMM_ROTATE_LEFT) .or. &
       (right_procid .ne. -1 .and. icomm_requests(rid)%comm_direction .eq. COMM_TO_LEFT)) then
      unpack_direction = UNPACK_RIGHT
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_RIGHT .and. left_procid .ne. -1) then
      unpack_direction = UNPACK_LEFT
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_TOP .and. bot_procid .ne. -1) then
      unpack_direction = UNPACK_BOT
    end if
    if (icomm_requests(rid)%comm_direction .eq. COMM_TO_BOT .and. top_procid .ne. -1) then
      unpack_direction = UNPACK_TOP
    end if

    if (unpack_direction .eq. UNPACK_NONE) return

    do i = 1, icomm_requests(rid)%num_arrays
      comm_array_id = icomm_requests(rid)%comm_array_ids(i)
      pack_mode = 0
      if (registered_comm_arrays(comm_array_id)%is_phys_array) pack_mode = 1
      call gamil_3D_array_unpack(unpack_direction,                                 &
                                 registered_comm_arrays(comm_array_id)%ibeg,       &
                                 registered_comm_arrays(comm_array_id)%iend,       &
                                 registered_comm_arrays(comm_array_id)%jbeg,       &
                                 registered_comm_arrays(comm_array_id)%jend,       &
                                 registered_comm_arrays(comm_array_id)%kbeg,       &
                                 registered_comm_arrays(comm_array_id)%kend,       &
                                 iter,                                             &
                                 icomm_requests(rid)%nlevels,                      &
                                 registered_comm_arrays(comm_array_id)%array(:),   &
                                 pack_mode,                                        &
                                 icomm_requests(rid)%recv_buffer) 
    end do

end subroutine unpack_request_comm_arrays


subroutine gamil_array_send_left(send_array, recv_array,arr_size)
    implicit none
    integer :: arr_size
    real :: send_array(arr_size)
    real :: recv_array(arr_size)

    call t_startf("dyn 1D sendrecv")
    if (num_x_proc .gt. 1) then
      call MPI_ISEND(send_array,arr_size,MPI_REAL8,left_procid,1,gamil_comm,isend,ierr)
      call MPI_WAIT(isend,status,ierr)                         
    end if
    call t_stopf("dyn 1D sendrecv")
end subroutine gamil_array_send_left


subroutine gamil_array_send_right(send_array, recv_array, arr_size)
    implicit none
    integer :: arr_size
    real :: send_array(arr_size)
    real :: recv_array(arr_size)

    call t_startf("dyn 1D sendrecv")
    if (num_x_proc .gt. 1) then
      call MPI_ISEND(send_array,arr_size,MPI_REAL8,right_procid,1,gamil_comm,isend,ierr)
      call MPI_WAIT(isend,status,ierr)                         
    end if
    call t_stopf("dyn 1D sendrecv")
end subroutine gamil_array_send_right


subroutine gamil_array_recv_left(send_array,recv_array,arr_size)
    implicit none
    integer :: arr_size
    real :: send_array(arr_size)
    real :: recv_array(arr_size)

    call t_startf("dyn 1D sendrecv")
    if (num_x_proc .gt. 1) then
      call MPI_IRECV(recv_array,arr_size,MPI_REAL8,left_procid,1,gamil_comm,irecv,ierr)
      call MPI_WAIT(irecv,status,ierr)                         
    else
      recv_array = send_array
    end if
    call t_stopf("dyn 1D sendrecv")
end subroutine gamil_array_recv_left


subroutine gamil_array_recv_right(send_array,recv_array,arr_size)
    implicit none
    integer :: arr_size
    real :: send_array(arr_size)
    real :: recv_array(arr_size)

    call t_startf("dyn 1D sendrecv")
    if (num_x_proc .gt. 1) then
      call MPI_IRECV(recv_array,arr_size,MPI_REAL8,right_procid,1,gamil_comm,irecv,ierr)
      call MPI_WAIT(irecv,status,ierr)                         
    else
      recv_array = send_array
    end if
    call t_stopf("dyn 1D sendrecv")
end subroutine gamil_array_recv_right


subroutine gamil_arrays_comm(comm_direction,nlevels, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, &
                            arr9, arr10, arr11, arr12, arr13, arr14, arr15, arr16, request_id)
    implicit none

    integer :: nlevels
    integer :: comm_direction
    real,optional :: arr1(:)
    real,optional :: arr2(:)
    real,optional :: arr3(:)
    real,optional :: arr4(:)
    real,optional :: arr5(:)
    real,optional :: arr6(:)
    real,optional :: arr7(:)
    real,optional :: arr8(:)
    real,optional :: arr9(:)
    real,optional :: arr10(:)
    real,optional :: arr11(:)
    real,optional :: arr12(:)
    real,optional :: arr13(:)
    real,optional :: arr14(:)
    real,optional :: arr15(:)
    real,optional :: arr16(:)
    integer,optional :: request_id
    integer iter, local_request_id

    local_request_id = get_free_icomm_request()
    icomm_requests(local_request_id)%nlevels = nlevels
    icomm_requests(local_request_id)%comm_direction = comm_direction
      
    if (present(arr1)) call add_comm_array(local_request_id, arr1) 
    if (present(arr2)) call add_comm_array(local_request_id, arr2) 
    if (present(arr3)) call add_comm_array(local_request_id, arr3) 
    if (present(arr4)) call add_comm_array(local_request_id, arr4) 
    if (present(arr5)) call add_comm_array(local_request_id, arr5) 
    if (present(arr6)) call add_comm_array(local_request_id, arr6) 
    if (present(arr7)) call add_comm_array(local_request_id, arr7) 
    if (present(arr8)) call add_comm_array(local_request_id, arr8) 
    if (present(arr9)) call add_comm_array(local_request_id, arr9) 
    if (present(arr10)) call add_comm_array(local_request_id, arr10) 
    if (present(arr11)) call add_comm_array(local_request_id, arr11) 
    if (present(arr12)) call add_comm_array(local_request_id, arr12) 
    if (present(arr13)) call add_comm_array(local_request_id, arr13) 
    if (present(arr14)) call add_comm_array(local_request_id, arr14) 
    if (present(arr15)) call add_comm_array(local_request_id, arr15) 
    if (present(arr16)) call add_comm_array(local_request_id, arr16) 
    call pack_request_comm_arrays(local_request_id,iter)
    if (present(request_id)) request_id = local_request_id

    call t_startf("dyn sendrecv")
    
    if (comm_direction .eq. COMM_ROTATE_LEFT) then
      if (beglonex .eq. 1) then
        call MPI_ISEND(icomm_requests(local_request_id)%send_buffer,iter,MPI_REAL8,left_procid, &
                       local_request_id,gamil_comm,isend,ierr)
        icomm_requests(local_request_id)%num_requests = icomm_requests(local_request_id)%num_requests+1
        icomm_requests(local_request_id)%requests(icomm_requests(local_request_id)%num_requests) = isend
      end if
      if (endlonex .eq. XLON) then
        call MPI_IRECV(icomm_requests(local_request_id)%recv_buffer,iter,MPI_REAL8,right_procid, &
                       local_request_id,gamil_comm,irecv,ierr)
        icomm_requests(local_request_id)%num_requests = icomm_requests(local_request_id)%num_requests+1
        icomm_requests(local_request_id)%requests(icomm_requests(local_request_id)%num_requests) = irecv
      end if
    else
      if (comm_direction .eq. COMM_TO_LEFT) then
        call MPI_ISEND(icomm_requests(local_request_id)%send_buffer,iter,MPI_REAL8,left_procid, &
                       local_request_id,gamil_comm,isend,ierr)
        call MPI_IRECV(icomm_requests(local_request_id)%recv_buffer,iter,MPI_REAL8,right_procid, & 
                       local_request_id,gamil_comm,irecv,ierr)
      else if (comm_direction .eq. COMM_TO_RIGHT) then
        call MPI_ISEND(icomm_requests(local_request_id)%send_buffer,iter,MPI_REAL8,right_procid, &
                       local_request_id,gamil_comm,isend,ierr)
        call MPI_IRECV(icomm_requests(local_request_id)%recv_buffer,iter,MPI_REAL8,left_procid, &
                       local_request_id,gamil_comm,irecv,ierr)
      else if (comm_direction .eq. COMM_TO_TOP) then
        call MPI_ISEND(icomm_requests(local_request_id)%send_buffer,iter,MPI_REAL8,top_procid, &
                       local_request_id,gamil_comm,isend,ierr)
        call MPI_IRECV(icomm_requests(local_request_id)%recv_buffer,iter,MPI_REAL8,bot_procid, &
                       local_request_id,gamil_comm,irecv,ierr)
      else if (comm_direction .eq. COMM_TO_BOT) then
        call MPI_ISEND(icomm_requests(local_request_id)%send_buffer,iter,MPI_REAL8,bot_procid, &
                       local_request_id,gamil_comm,isend,ierr)
        call MPI_IRECV(icomm_requests(local_request_id)%recv_buffer,iter,MPI_REAL8,top_procid, &
                       local_request_id,gamil_comm,irecv,ierr)
      else 
        write(6,*) 'error comm direction'
        call endrun
      end if
      icomm_requests(local_request_id)%num_requests = 2
      icomm_requests(local_request_id)%requests(1) = isend
      icomm_requests(local_request_id)%requests(2) = irecv
    end if
    call t_stopf("dyn sendrecv")

    if (.not.present(request_id)) then
      call wait_icomm_request(local_request_id)
    end if
end subroutine gamil_arrays_comm


subroutine broadcast_lon_data(lon_indx, lat_indx, lon_data1, num_data1, & 
                                                  lon_data2, num_data2, &
                                                  lon_data3, num_data3, &
                                                  lon_data4, num_data4, &
                                                  lon_data5, num_data5, &
                                                  lon_data6, num_data6, &
                                                  lon_data7, num_data7, &
                                                  lon_data8, num_data8, &
                                                  lon_data9, num_data9 )
    implicit none
    real lon_data1(:)
    integer num_data1
    integer lon_indx, lat_indx
    integer bcast_rootproc_id
    integer i, nproc_pole
    real, optional :: lon_data2(:), lon_data3(:), lon_data4(:), lon_data5(:), lon_data6(:), lon_data7(:), lon_data8(:), lon_data9(:)
    integer, optional :: num_data2, num_data3, num_data4, num_data5, num_data6, num_data7, num_data8, num_data9
    real :: tmp_array(ZALT*16)
    integer :: num_data_all, iter

    if (lat_indx .ne. jbeg0 .and. lat_indx .ne. jend0) return

    bcast_rootproc_id = -1
    call MPI_COMM_SIZE(lat_comm,nproc_pole,ierr)
    do i = 1, nproc_pole
       if (lat_procs_ibegs(i) .le. lon_indx .and. lat_procs_iends(i) .ge. lon_indx) bcast_rootproc_id = i-1
    enddo

    iter = 1
    do i = 1, num_data1
       tmp_array(iter) = lon_data1(i)
       iter = iter+1
    enddo

    if (present(lon_data2)) then
    if (.not.present(num_data2)) call endrun
    do i = 1, num_data2
       tmp_array(iter) = lon_data2(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data3)) then
    if (.not.present(num_data3)) call endrun
    do i = 1, num_data3
       tmp_array(iter) = lon_data3(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data4)) then
    if (.not.present(num_data4)) call endrun
    do i = 1, num_data4
       tmp_array(iter) = lon_data4(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data5)) then
    if (.not.present(num_data5)) call endrun
    do i = 1, num_data5
       tmp_array(iter) = lon_data5(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data6)) then
    if (.not.present(num_data6)) call endrun
    do i = 1, num_data6
       tmp_array(iter) = lon_data6(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data7)) then
    if (.not.present(num_data7)) call endrun
    do i = 1, num_data7
       tmp_array(iter) = lon_data7(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data8)) then
    if (.not.present(num_data8)) call endrun
    do i = 1, num_data8
       tmp_array(iter) = lon_data8(i)
       iter = iter+1
    enddo
    endif

    if (present(lon_data9)) then
    if (.not.present(num_data9)) call endrun
    do i = 1, num_data9
       tmp_array(iter) = lon_data9(i)
       iter = iter+1
    enddo
    endif

    call t_startf("dyn bcast lon")
    call mpi_bcast(tmp_array,iter,MPI_REAL8,bcast_rootproc_id, lat_comm, ierr) 
    call t_stopf("dyn bcast lon")

    iter = 1
    do i = 1, num_data1
       lon_data1(i) = tmp_array(iter)
       iter = iter+1
    enddo

    if (present(lon_data2)) then
    if (.not.present(num_data2)) call endrun
    do i = 1, num_data2
       lon_data2(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data3)) then
    if (.not.present(num_data3)) call endrun
    do i = 1, num_data3
       lon_data3(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data4)) then
    if (.not.present(num_data4)) call endrun
    do i = 1, num_data4
       lon_data4(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data5)) then
    if (.not.present(num_data5)) call endrun
    do i = 1, num_data5
       lon_data5(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data6)) then
    if (.not.present(num_data6)) call endrun
    do i = 1, num_data6
       lon_data6(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data7)) then
    if (.not.present(num_data7)) call endrun
    do i = 1, num_data7
       lon_data7(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data8)) then
    if (.not.present(num_data8)) call endrun
    do i = 1, num_data8
       lon_data8(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

    if (present(lon_data9)) then
    if (.not.present(num_data9)) call endrun
    do i = 1, num_data9
       lon_data9(i) = tmp_array(iter)
       iter = iter+1
    enddo
    endif

end subroutine broadcast_lon_data


subroutine gamil_all_reduce(value1,value2,value3,value4)
    implicit none
    real*16 value1
    real*16,optional:: value2
    real*16,optional:: value3
    real*16,optional:: value4
    real*16 values(10), results(10)
    real*16 tmp_result
    integer i,num_values

    values(1) = value1
    num_values = 1
    if (present(value2)) then
        values(2) = value2
        num_values = num_values+1
    endif
    if (present(value3)) then
        values(3) = value3
        num_values = num_values+1
    endif
    if (present(value4)) then
        values(4) = value4
        num_values = num_values+1
    endif

    call t_startf("dyn allreduce")
    call allreduce_real16(values,results,num_values,gamil_comm)
    call t_stopf("dyn allreduce")

    value1 = results(1)
    if (present(value2)) value2 = results(2)
    if (present(value3)) value3 = results(3)
    if (present(value4)) value4 = results(4)
    
end subroutine gamil_all_reduce


subroutine gamil_average_pole_data_phys(lat_indx,ibeg,iend,lev,array1,array2)
    implicit none
    integer :: ibeg, iend, lev
    real :: array1(ibeg:iend,lev)
    real,optional :: array2(ibeg:iend,lev)
    real*16 local_sum(lev*2), local_sum_tmp(8,lev), global_sum(lev*2), tmp_sum
    integer :: local_ibeg, local_iend
    integer :: i, k, lat_indx


    if (lat_indx .ne. jbeg0 .and. lat_indx .ne. jend0) return

    local_ibeg = ibeg0
    local_iend = iend0
    if (iend0 .eq. XLON) local_iend = iend0 - 2

!$OMP PARALLEL do private(i,k,tmp_sum)
    do k=1,lev
       tmp_sum = 0.0
       do i=local_ibeg, local_iend
          tmp_sum = tmp_sum + array1(i,k)
       enddo 
       local_sum_tmp(1,k) = tmp_sum
    enddo
    do k=1,lev
       local_sum(k) = local_sum_tmp(1,k)
    enddo

    if (present(array2)) then
!$OMP PARALLEL do private(i,k,tmp_sum)
    do k=1,lev
       tmp_sum = 0.0
       do i=local_ibeg, local_iend
          tmp_sum = tmp_sum + array2(i,k)
       enddo 
       local_sum_tmp(1,k) = tmp_sum
    enddo
    do k=1,lev
       local_sum(k+lev) = local_sum_tmp(1,k)
    enddo
    endif
    
    call t_startf("dyn polereduce")
    if (present(array2)) then
       call allreduce_real16(local_sum,global_sum,lev*2,lat_comm)
    else 
       call allreduce_real16(local_sum,global_sum,lev,lat_comm)
    endif
    call t_stopf("dyn polereduce")

!$OMP PARALLEL do private(k,i,tmp_sum)
    do k=1,lev
       tmp_sum = global_sum(k) / dble(XLON-2)
       do i=ibeg, iend
          array1(i,k) = tmp_sum 
       enddo
    enddo

    if (present(array2)) then
!$OMP PARALLEL do private(k,i,tmp_sum)
    do k=1,lev
       tmp_sum = global_sum(k+lev) / dble(XLON-2)
       do i=ibeg, iend
          array2(i,k) = tmp_sum 
       enddo
    enddo
    endif

end subroutine gamil_average_pole_data_phys

subroutine gamil_sum_pole_data_phys(lat_indx,array1,num_data1,array2,num_data2)
    implicit none
    integer :: num_data1
    integer :: lat_indx
    real*16 :: array1(:)
    real*16 :: tmp_array1(ZALT*16), tmp_array2(ZALT*16)
    real*16, optional :: array2(:)
    integer, optional :: num_data2
    integer :: iter


    if (lat_indx .ne. jbeg0 .and. lat_indx .ne. jend0) return

    tmp_array1(:num_data1) = array1(:num_data1)
    iter = num_data1
    
    if (present(array2)) then
       if (.not.(present(num_data2))) call endrun
       tmp_array1(iter+1:iter+num_data2) = array2(:num_data2)
       iter = iter+num_data2
    endif
    
    call t_startf("dyn polereduce")
    call allreduce_real16(tmp_array1,tmp_array2,iter,lat_comm)
    call t_stopf("dyn polereduce")

    array1(:num_data1) = tmp_array2(:num_data1)
    iter = num_data1

    if (present(array2)) then
       array2(:num_data2) = tmp_array2(iter+1:iter+num_data2)
       iter = iter+num_data2
    endif

end subroutine gamil_sum_pole_data_phys


subroutine gamil_min_lat_row_data(array1,num_data1,array2,num_data2)
    implicit none
    integer :: num_data1
    real*8 :: array1(:)
    real*8 :: tmp_array1((ZALT+YLAT)*16), tmp_array2((ZALT+YLAT)*16)
    real*8, optional :: array2(:)
    integer, optional :: num_data2
    integer :: iter


    tmp_array1(:num_data1) = array1(:num_data1)
    iter = num_data1
    
    if (present(array2)) then
       if (.not.(present(num_data2))) call endrun
       tmp_array1(iter+1:iter+num_data2) = array2(:num_data2)
       iter = iter+num_data2
    endif
    
    call t_startf("dyn latreduce")
    call mpi_allreduce(tmp_array1,tmp_array2,iter,MPI_REAL8,MPI_MIN,lat_comm,ierr)
    call t_stopf("dyn latreduce")

    array1(:num_data1) = tmp_array2(:num_data1)
    iter = num_data1

    if (present(array2)) then
       array2(:num_data2) = tmp_array2(iter+1:iter+num_data2)
       iter = iter+num_data2
    endif

end subroutine gamil_min_lat_row_data


subroutine endrun()
   write(6,*) 'gamil call endrun'
   call mpi_abort (gamil_comm, 1)
end subroutine endrun

end module mpi_gamil
