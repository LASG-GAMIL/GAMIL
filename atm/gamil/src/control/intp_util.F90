#include <misc.h>

subroutine chktime( time, nrec )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Make sure the time coordinate looks like calander day, and is increasing.
   ! Calendar day can either start with 1 Jan 0Z = day 1.0  or
   !  1 Jan 0Z = day 0.0
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

   integer, intent(in) ::&
      nrec                 ! size of time array
   real(r8), intent(in) ::&
      time(nrec)           ! time coordinate expressed as calendar day.
   
   ! Local varibles:
   integer :: i
   !-----------------------------------------------------------------------

   if ( time(1) .lt. 0. .or. time(1) .ge. 366. ) then
      write(6,*)'chktime: illegal time coordinate ',time(1)
      stop
   end if
   do i = 2, nrec
      if ( time(i) .lt. 0. .or. time(i) .ge. 366. ) then
         write(6,*)'chktime: illegal time coordinate ', time(i)
         stop
      end if
      if ( time(i) .le. time(i-1) ) then
         write(6,*)'chktime: time not increasing ', time(i-1), time(i)
         stop
      end if
   end do

   return

end subroutine chktime

!#######################################################################

subroutine findplb( x, nx, xval, index )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

   integer, intent(in) ::&
      nx                   ! size of x
   real(r8), intent(in) ::&
      x(nx),              &! strictly increasing array
      xval                 ! value to be searched for in x
   
   integer, intent(out) ::&
      index

   ! Local variables:
   integer i
   !-----------------------------------------------------------------------

   if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
      index = nx
      return
   end if

   do i = 2, nx
      if ( xval .lt. x(i) ) then
         index = i-1
         return
      end if
   end do

   return

end subroutine findplb

!#######################################################################

subroutine linintp( fdimd, nver, ldimd, t1, t2, tint, f1, f2, fint )

   !-----------------------------------------------------------------------
   ! Purpose:
   ! Linearly interpolate between f1(t1) and f2(t2) to fint(tint),
   ! where f1, f2, and f3 are chunked data structures
   !
   ! Author: B. Eaton
   ! Chunked by P. Worley
   !-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid, only: get_ncols_p

   implicit none

   integer, intent(in) ::&
      fdimd,                   &! declared length of first dimension
      nver,                    &! number of vertical levels
      ldimd                     ! length of last dimension
   real(r8), intent(in) ::&
      t1,                                           &! time level of f1
      t2,                                           &! time level of f2
      tint,                                         &! interpolant time
      f1(fdimd,pcols,nver,begchunk:endchunk,ldimd), &! field at time t1
      f2(fdimd,pcols,nver,begchunk:endchunk,ldimd)   ! field at time t2

   real(r8), intent(out) ::&
      fint(fdimd,pcols,nver,begchunk:endchunk,ldimd) ! field at time tint

   ! Local variables:
   integer :: f, i, k, l, lchnk, ncol
   real(r8) :: factor
   !------------------------------------------------------------------------------

   factor = ( tint - t1 )/( t2 - t1)

   do lchnk = begchunk, endchunk
      ncol = get_ncols_p(lchnk)
      do l=1,ldimd
         do k=1,nver
            do i=1,ncol
               do f=1,fdimd
                  fint(f,i,k,lchnk,l) = f1(f,i,k,lchnk,l) &
                     + ( f2(f,i,k,lchnk,l) - f1(f,i,k,lchnk,l) )*factor
               enddo
            enddo
         enddo
      enddo
   enddo

   return

end subroutine linintp

!#######################################################################
