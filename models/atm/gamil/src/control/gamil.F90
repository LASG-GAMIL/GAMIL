#include <misc.h>
#include <params.h>

!-----------------------------------------------------------------------
!
! Purpose: Entry point for GAMIL
!
!-----------------------------NOTICE------------------------------------
!
! Gridpoint Atmosphere Model of IAP/LASG, version 1.0
!
! Purpose:
!
!   Call initialization, time-stepping, and finalization routines.
!
!-----------------------------------------------------------------------

program gamil

    use shr_kind_mod, only: r8 => SHR_KIND_R8
    use pmgrid
    use dycore
    use history,      only: bldfld, intht
    use units
    use restart,      only: read_restart
    use time_manager, only: get_nstep, is_first_restart_step
    use phys_buffer ! added by SHI Xiangjun and LIU Li
    use ppgrid,       only: pcols, pverp, begchunk, endchunk
    use comsrf,       only: fld_kvh ! added by LIU Li
#if (defined SPMD || defined COUP_CSM)
    use mpishorthand, only: mpicom, nsend, nrecv, nwsend, nwrecv
#endif

#ifdef COUP_CSM
    use shr_msg_mod
    use cpl_interface_mod ! For FGOALS2.0
    use cpl_fields_mod    !
#endif

    implicit none

#include <comctl.h>
#include <comlun.h>
#include <gpt.inc>
! added by SHI Xiangjun
! DONG Li: try to merge it with the namelist variable
#include <RK_or_MG.h>

#ifdef SUNOS
!#include <floatingpoint.h>
#endif

#ifdef OSF1
#include <for_fpe_flags.f>
    integer(4) old_fpe_flags   ! old settings of floating point exception flags
    integer(4) new_fpe_flags   ! new settings of floating point exception flags
    integer(4) for_set_fpe     ! function to set the floating point exceptions
#endif
    character*8 cdate          ! System date
    character*8 ctime          ! System time
    character*13 filenam
    integer iu
    integer nstep           ! Current timestep number.
    integer kvh_idx ! added by LIU Li
    integer i
    !------------------------------Externals--------------------------------
#if ( defined SUNOS )
    !integer iexcept, ieee_handler
#endif

#ifdef OSF1
    !
    ! Compaq floating point exception handler
    ! Terminate if hit invalid, divide by zero, or overflow.
    !
    new_fpe_flags = FPE_M_TRAP_INV + FPE_M_TRAP_DIV0 + FPE_M_TRAP_OVF
    old_fpe_flags = for_set_fpe(new_fpe_flags)
#endif

#if ( defined SUNOS )
    !
    ! SUN: Trap ieee exceptions for debugging purposes
    !      iexcept = ieee_handler( 'set', 'common', SIGFPE_ABORT )
    !      if ( iexcept /= 0 ) write(6,*)'ieee trapping not supported here'
    !
#endif
    !
    ! Initialize timing library.  2nd arg 0 means disable, 1 means enable
    !
    call t_setoptionf(usrsys, 0)
    call t_initializef

    call t_startf('total')
    call t_startf('initialization')
    !
    ! Initialize internal/external MPI if appropriate
    !
#ifdef COUP_CSM
    call shr_msg_stdio('atm')
    !call shr_msg_init('atm')                             ! For FGOALS2.0
    !call shr_msg_groups('atm')                           !
    call cpl_interface_init(cpl_fields_atmname, mpicom)   !
#endif
    !
    ! Initialize SPMD environment if applicable
    !
#if ( defined SPMD )
    call spmdinit
#endif
    !
    ! Print Model heading and copyright message
    !
    if (masterproc) then
        write(6, *) "----------------------------------------------------------"
        write(6, *) "    Gridpoint Atmosphere Model of IAP/LASG  (GAMIL)       "
        write(6, *) "                    version 2.0                           "
        write(6, *) "----------------------------------------------------------"
    end if
    !
    ! Fetch and print current date and time
    !
    call datetime(cdate, ctime)
    if (masterproc) then
        write(6, *) "DATE ", cdate, " TIME ", ctime
        write(6, *) "----------------------------------------------------------"
        if (dycore_is('EUL')) then
            write(6, *) 'DYCORE is EUL'
        else if (dycore_is('SLD')) then
            write(6, *) 'DYCORE is SLD'
        else if (dycore_is('LR')) then
            write(6, *) 'DYCORE is LR'
        end if
    end if
    !
    ! Set defaults then override with user-specified input
    !
    call preset
    call parse_namelist
    !
    ! Define fortran unit numbers
    !
    nsds    = getunit()
    nrg     = getunit()
    nrg2    = getunit()
    luhrest = getunit()

    if (masterproc) then
        write(6, *)
        write(6, "('=========================================')")
        write(6, "('***Summary of Logical Unit assignments***')")
        write(6, *)
        write(6, "('   Restart pointer unit (nsds)     : ', I2)") nsds
        write(6, "('   Master restart unit (nrg)       : ', I2)") nrg
        write(6, "('   Abs/ems unit for restart (nrg2) : ', I2)") nrg2
        write(6, "('   History restart unit (luhrest)  : ', I2)") luhrest
        write(6, "('=========================================')")
        write(6, *)
    end if
    !
    ! Initialize index values for advected and non-advected tracers
    !
    call initindx
    !
    ! Do appropriate dynamics and history initialization depending on whether initial, restart, or
    ! branch.  On restart run intht need not be called because all the info is on restart dataset.
    !
    select case (nsrest)
    case (0)                ! initial run
        call inital          ! dynamics (mostly) init
        call inti            ! physics init
        call bldfld          ! master field list
        call intht           ! set up history tape contents for this run
    case (1)                ! restart run
        call read_restart    ! read restart file(s)
        call inti            ! physics init
        call bldfld          ! master field list
        call intht
    case (3)                ! branch run
        call read_restart    ! read restart file(s), minus history info
        call inti            ! physics init
        call bldfld          ! master field list
        call intht           ! set up history tape contents for this run
    case default
        write(6, *)' nsrest=', nsrest, ' must be 0, 1, or 3'
        call endrun
    end select
    !
    ! Initialize external models or datasets depending upon whether coupled
    !
    call initext
    call t_stopf('initialization')
    !
    ! Invoke driving routine for time integration
    !
    if (RK_or_MG == 'MG') then
        call pbuf_allocate('global') ! added by SHI Xiangjun
        if (is_first_restart_step() ) then ! added by LIU Li
            kvh_idx = pbuf_get_fld_idx('KVH')
            do i = begchunk, endchunk
                pbuf(kvh_idx)%fld_ptr(1,1:pcols,1:pverp,i,1) = fld_kvh(1:pcols,1:pverp,i)
            end do
        end if
    end if
    call t_startf('stepon')
    call stepon
    call t_stopf('stepon')
    if (RK_or_MG == 'MG') call pbuf_deallocate('global')  ! added by SHI Xiangjun
    !
    ! End the run cleanly
    !
    call t_stopf('total')
    call t_prf(iam)

#if ( defined SPMD )
    if (.false.) then
        write(0,*)'The following stats are exclusive of initialization/boundary datasets'
        write(0,*)'Number of messages sent by proc ',iam,' is ',nsend
        write(0,*)'Number of messages recv by proc ',iam,' is ',nrecv
    end if
#endif

    if (masterproc) then
        nstep = get_nstep()
        write (6,9300) nstep-1,nstep
9300    format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
            ' partially done to provide convectively adjusted and ', &
            'time filtered values for history tape.')
        write(6,*)'------------------------------------------------------------'
        write(6,*)'******* END OF MODEL RUN *******'
    end if

#ifdef COUP_CSM
    !call shr_msg_finalize                          ! For FGOALS2.0
    call cpl_interface_finalize(cpl_fields_atmname) !
#elif ( defined SPMD )
    call mpibarrier(mpicom)
    call mpifinalize
#endif

#if ( defined SPMD )
    iu = getunit ()
    write(filenam,'(a10,i3.3)') 'spmdstats.', iam
    open (unit=iu, file=filenam, form='formatted', status='replace')
    write (iu,*)'iam ',iam,' msgs  sent =',nsend
    write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
    write (iu,*)'iam ',iam,' words sent =',nwsend
    write (iu,*)'iam ',iam,' words recvd=',nwrecv
#endif

    stop

end program gamil
