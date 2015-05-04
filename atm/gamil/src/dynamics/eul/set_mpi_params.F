#include <misc.h>
#include <params.h>




!!(wh 2003.10.28)
!!---------------

      subroutine set_mpi_params()
     
#if (defined SPMD)
 
      implicit none

#include <PARADYN>

#include <mpif.h>
#include <commpi.h>

      if (myrank.eq.0) then

        jbeg  = 2
        jend  = ny-2
        jpole = ny-1
        inc_pole = .true.

        itop = myrank + 1
        ibot = MPI_PROC_NULL

      elseif (myrank.eq.nprocs-1) then

        jpole = 2
        jbeg  = 3
        jend  = ny-1
        inc_pole = .true.

        itop = MPI_PROC_NULL
        ibot = myrank -1

      else

        jbeg = 2
        jend = ny-1
        inc_pole = .false.

        itop = myrank + 1
        ibot = myrank - 1

      endif


      write(*,10) '*** p',myrank,': jbeg=',jbeg,'jend=',jend,'jpole=',jpole,'inc_pole=',inc_pole
      write(*,10) '*** p',myrank,': itop=',itop,'ibot=',ibot
10    format(1x,a5,i1,3(a10,i3),a12,l4)
!!      stop

#endif
      return
      end 

