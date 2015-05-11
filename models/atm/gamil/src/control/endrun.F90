#include <misc.h>

subroutine endrun(msg)
!-----------------------------------------------------------------------
! Purpose:
!
! Abort the model for abnormal termination
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
! $Id: endrun.F90,v 1.4 2001/08/24 16:14:45 mvertens Exp $
!-----------------------------------------------------------------------
#if (defined SPMD || defined COUP_CSM)
    use mpishorthand, only: MPI_COMM_WORLD
#endif
    use shr_sys_mod,  only: shr_sys_flush

    implicit none

    character(*), intent(in), optional :: msg

    if (present(msg)) then
        write(6, "('Notice: endrun: ', A)") msg
    else
        write(6, "('Notice: endrun is being called.')")
    end if

    call shr_sys_flush(6) ! flush all output to standard output

#if (defined SPMD) || (defined COUP_CSM) 
! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
! error code of 257
    call mpi_abort(MPI_COMM_WORLD, 1)  
#else
    call abort
#endif

end subroutine endrun
