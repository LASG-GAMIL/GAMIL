#include "misc.h" ! added by SHI Xiangjun 2008/11/10

module phys_buffer

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Implement a physics buffer to hold arrays that must persist
    !   across timesteps or between calls to different physics packages.
    !
    !   Current implementation only supports one buffer.
    !
    ! Author: B. Eaton
    !
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use infnan,       only: nan
    use pmgrid,       only: masterproc
    use ppgrid,       only: pcols, begchunk, endchunk
    use phys_grid,    only: physgrid_set, write_field_from_chunk, &
                            read_chunk_from_field, get_ncols_p
    use dyn_grid,     only: ptimelevels
    use string_utils, only: to_upper

#ifdef SPMD
    use mpishorthand, only: mpicom
#endif

    implicit none

    private

    save

    public pbuf_defaultopts    ! set namelist defaults
    public pbuf_setopts        ! set user specified namelist values
    public pbuf_init           ! initialize physics buffer
    public pbuf_add            ! add field to physics buffer
    public pbuf_get_fld_idx    ! get index of specified field in the physics buffer
    public pbuf_old_tim_idx    ! return the index for the oldest time
    public pbuf_next_tim_idx   ! return the index for the next time
    public pbuf_update_tim_idx ! update the index for the oldest time
    public pbuf_allocate       ! allocate memory for physics buffer fields
    public pbuf_deallocate     ! deallocate memory for physics buffer fields
    public pbuf_setval         ! set value for a field in the buffer
    public pbuf_write_restart  ! write physics buffer to restart file
    public pbuf_read_restart   ! read physics buffer from restart file

    type, public :: pbuf_fld
        character(16) name
        character(16) scope
        integer fdim, mdim, ldim
        real(r8), pointer :: fld_ptr(:,:,:,:,:)
    end type pbuf_fld

    integer, public, parameter :: pbuf_size_max = 100
    type(pbuf_fld), public :: pbuf(pbuf_size_max)

    integer, public :: pbuf_times  ! number of time levels in physics buffer (dycore dependent)

    integer :: pbuf_size = 0
    integer :: old_time_idx = 1
    logical :: global_allocate_all = .true.  ! allocate all buffers as global  !!!sxj

contains

    ! return default runtime options

    subroutine pbuf_defaultopts(pbuf_global_allocate_out)
        logical, intent(out), optional :: pbuf_global_allocate_out

        if (present(pbuf_global_allocate_out)) then
            pbuf_global_allocate_out = global_allocate_all
        end if

    end subroutine pbuf_defaultopts

    ! set runtime options

    subroutine pbuf_setopts(pbuf_global_allocate_in)
        logical, intent(in), optional :: pbuf_global_allocate_in

        if (present(pbuf_global_allocate_in)) then
            global_allocate_all = pbuf_global_allocate_in
        end if

    end subroutine pbuf_setopts

    ! initialize physics buffer

    subroutine pbuf_init
        integer i

        pbuf_times = ptimelevels-1

        do i = 1, pbuf_size_max
            nullify(pbuf(i)%fld_ptr)
        end do

    end subroutine pbuf_init

    ! add a field to the physics buffer

    subroutine pbuf_add(name, scope, fdim, mdim, ldim, index)
        character(*), intent(in) :: name   ! field name
        character(*), intent(in) :: scope  ! field scope, either 'global' or 'physpkg'
        integer, intent(in)  :: fdim       ! first generic field dimension
        integer, intent(in)  :: mdim       ! middle generic field dimension
        integer, intent(in)  :: ldim       ! last generic field dimension
        integer, intent(out) :: index      ! index in the physics buffer

        character(*), parameter :: sub = 'pbuf_add'
        integer i

        if (pbuf_size >= pbuf_size_max) then
            call endrun(sub//': max number physics buffer fields exceeded. Increase pbuf_size_max in phys_buffer.F90')
        end if

        do i = 1, pbuf_size
            if (pbuf(i)%name == name) then
                call endrun (sub//': ERROR: field name '//name//' is already in use.')
            end if
        end do

        if (scope /= 'global' .and. scope /= 'physpkg') then
            call endrun (sub//': scope must be set to global or physpkg.  Current value is: '//scope)
        end if

        pbuf_size = pbuf_size+1
        index = pbuf_size
        pbuf(index)%name = name
        pbuf(index)%scope = scope
        pbuf(index)%fdim = fdim
        pbuf(index)%mdim = mdim
        pbuf(index)%ldim = ldim

        write(6, "('Notice: phys_buffer::pbuf_add: add buffer for ', A10, ' at index', I3)") name, index

    end subroutine pbuf_add

    function pbuf_get_fld_idx(name, failcode)

        ! Get index of specified field in the physics buffer.  String matching is case insensitive.
        ! Call endrun if name not found

        implicit none

        character(len=*), intent(in)  :: name   ! field name
        integer, intent(in), optional :: failcode

        ! Return value
        integer :: pbuf_get_fld_idx

        ! Local variables
        integer :: i
        character(len=len(name)) :: Uname       ! =to_upper(name)
        !-----------------------------------------------------------------------------------------

        !
        !  Search for specified field in physics buffer, assuming that case of
        !  argument "name" matches definition in pbuf structure.
        !

        do i = 1, pbuf_size
            if ( pbuf(i)%name == name ) then
                pbuf_get_fld_idx = i
                return
            end if
        end do

        !
        !  Search for specified field in physics buffer, assuming that case of
        !  argument "name" does not matches definition in pbuf structure. In this
        !  instance, convert all names to upper case.
        !

        Uname = to_upper(name)
        do i = 1, pbuf_size
            if ( to_upper(pbuf(i)%name) == Uname ) then
                pbuf_get_fld_idx = i
                return
            end if
        end do

        if ( present(failcode)) then
            pbuf_get_fld_idx = failcode
            return
        else
            write (6,*) "PBUF_GET_FLD_IDX: index not found for ",name  !!--sxj
            call endrun ('PBUF_GET_FLD_IDX: index not found for '//name)
        endif

    end function pbuf_get_fld_idx

    !=========================================================================================

    function pbuf_old_tim_idx()

        ! Return index of oldest time sample in the physics buffer.

        implicit none

        ! Return value
        integer :: pbuf_old_tim_idx
        !-----------------------------------------------------------------------------------------

        pbuf_old_tim_idx = old_time_idx

    end function pbuf_old_tim_idx

    !=========================================================================================

    function pbuf_next_tim_idx(idx)

        ! Return index of next time sample in the physics buffer.

        implicit none

        integer, intent(in) :: idx

        ! Return value
        integer :: pbuf_next_tim_idx
        !-----------------------------------------------------------------------------------------

        pbuf_next_tim_idx = mod(idx, pbuf_times) + 1

    end function pbuf_next_tim_idx

    ! update index of old time sample in the physics buffer.

    subroutine pbuf_update_tim_idx()

        old_time_idx = mod(old_time_idx, pbuf_times)+1

    end subroutine pbuf_update_tim_idx

    ! Allocate storage for fields in the physics buffer with the specified scope.
    ! If global_allocate_all=.true. then storage for both global and physpkg scope
    ! is allocated just once, when scope='global'.
    !
    ! N.B. This routine must be called after phys_grid_init because that's
    !      where begchunk and endchunk are set

    subroutine pbuf_allocate(scope)
        character(*), intent(in)  :: scope

        character(*), parameter :: sub = 'pbuf_allocate'
        integer i, fdim, mdim, ldim, istat
        logical allocate_all

        if (.not. physgrid_set) then
            call endrun(sub//': ERROR: called before physics grid initialized.')
        end if

        ! allocate_all is used to force allocation of all fields at same time as allocation
        ! for global scope.
        allocate_all = .false.
        if (masterproc) then
            write(6, "('Notice: phys_buffer::pbuf_allocate: global allocate = ', l, ' buffer size = ', i)") &
                global_allocate_all, pbuf_size
        end if
        ! TAG-2
        ! check the consistency between global_allocate_all and scope
        if (global_allocate_all) then
            if (scope == 'global') then
                allocate_all = .true.
            else
                return
            end if
        end if

        do i = 1, pbuf_size
            if (pbuf(i)%scope == scope  .or.  allocate_all) then
                fdim = pbuf(i)%fdim
                mdim = pbuf(i)%mdim
                ldim = pbuf(i)%ldim
                allocate(pbuf(i)%fld_ptr(fdim,pcols,mdim,begchunk:endchunk,ldim), stat=istat)
                if (masterproc) then
                    write(6, "('Notice: phys_buffer::pbuf_allocate: allocate buffer for ', a)") pbuf(i)%name
                end if
                if (istat /= 0) then
                    call endrun(sub//': ERROR: allocate failed for '//pbuf(i)%name)
                end if
                pbuf(i)%fld_ptr = nan
            end if
        end do

    end subroutine pbuf_allocate

    subroutine pbuf_deallocate(scope)

        ! Deallocate storage for fields in the physics buffer with the specified scope.
        ! If global_allocate_all=.true. then storage for both global and physpkg scope
        ! is deallocated just once, when scope='global'.

        implicit none

        character(len=*), intent(in)  :: scope

        ! Local variables
        character(len=*), parameter :: sub = 'pbuf_deallocate'
        integer :: i, fdim, mdim, ldim
        logical :: deallocate_all
        !-----------------------------------------------------------------------------------------

        ! deallocate_all is used to force allocation of all fields at same time as allocation
        ! for global scope.
        deallocate_all = .false.
        if ( global_allocate_all ) then
            if ( scope == 'global' ) then
                deallocate_all = .true.
            else
                return
            end if
        end if

        do i = 1, pbuf_size
            if ( pbuf(i)%scope == scope   .or.  deallocate_all ) then
                if (associated(pbuf(i)%fld_ptr)) then
                    deallocate(pbuf(i)%fld_ptr)
                else
                    call endrun (sub//': ERROR: '//pbuf(i)%name//' is not allocated')
                end if
            end if
        end do

    end subroutine pbuf_deallocate

    !=========================================================================================

    subroutine pbuf_setval(name, value)

        ! Set a value for a field in the physics buffer.

        implicit none

        character(len=*), intent(in)  :: name   ! field name
        real(r8),         intent(in)  :: value

        ! Local variables
        character(len=*), parameter :: sub = 'pbuf_setval'
        integer :: idx, lchnk, ncols
        integer :: failcode = -1
        !-----------------------------------------------------------------------------------------

        ! get field index
        idx = pbuf_get_fld_idx(name, failcode)
        ! check return value
        if (idx == failcode) then
            call endrun (sub//': ERROR: name not found:'//name)
        end if

        ! check that field pointer is associated -- if so then set the value
        if (associated(pbuf(idx)%fld_ptr)) then
            do lchnk = begchunk, endchunk
                ncols = get_ncols_p(lchnk)
                pbuf(idx)%fld_ptr(:,1:ncols,:,lchnk,:) = value
            end do
        else
            call endrun (sub//': ERROR: field '//name//' is not allocated')
        end if

    end subroutine pbuf_setval


    !=========================================================================================
    subroutine pbuf_write_restart(iu)
        ! write physics buffer to restart file
        implicit none
        integer, intent(in) :: iu  ! Fortran unit number
        ! Local variables
        character(len=*), parameter :: sub = 'pbuf_write_restart'
        integer :: i, ioerr
        !-----------------------------------------------------------------------------------------
        if (masterproc) then
            write (iu,iostat=ioerr) old_time_idx
            if (ioerr /= 0 ) then
                write (6,*) sub,' ioerror ', ioerr,' on i/o unit = ',iu
                call endrun
            end if
        endif
        do i = 1, pbuf_size
            if ( pbuf(i)%scope == 'global' ) then
                call write_field_from_chunk(iu, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim, pbuf(i)%fld_ptr)
            end if
        end do
    end subroutine pbuf_write_restart



    !=========================================================================================
    subroutine pbuf_read_restart(iu)

        ! write physics buffer to restart file

        implicit none

#include<pdyn.h>              ! sxj---nprocessor display mpiint(module mpishorthand)

        integer, intent(in) :: iu  ! Fortran unit number

        ! Local variables
        character(len=*), parameter :: sub = 'pbuf_read_restart'
        integer :: i, ioerr
        !-----------------------------------------------------------------------------------------

        if (masterproc) then
            read (iu,iostat=ioerr) old_time_idx
            if (ioerr /= 0 ) then
                write (6,*) sub,' ioerror ', ioerr,' on i/o unit = ',iu
                call endrun
            end if
        endif

#if ( defined SPMD )
        call mpibcast(old_time_idx, 1, nprocessor, 0, mpicom)
#endif

        do i = 1, pbuf_size
            if ( pbuf(i)%scope == 'global' ) then
                call read_chunk_from_field(iu, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim, pbuf(i)%fld_ptr)
            end if
        end do

    end subroutine pbuf_read_restart

end module phys_buffer
