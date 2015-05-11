#include <misc.h>
#include <params.h>
subroutine mtdlss(a       ,b       ,c       ,rhs     ,x       , &
                  n       ,ld      ,num     ,indx    ,ws      , &
                  lws     ,ier     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Solves multiple tridiagonal linear systems ax=b for the vector x
! 
! Method: 
!     the input vectors are:
!        a  : subdiagonal elements   -
!        b  : diagonal elements       |- from matrix a above (see below)
!        c  : superdiagonal elements -
!        rhs: right hand side (vector b above)
!        indx: Index array for gather
!        npts: Number of points in indx
!
!     other inputs:
!        n  : order of linear system (length of diagonal element vector)
!        ld : leading dimension of subdiagonal matrices (see below)
!        num: number of tridiagonal systems to be solved
!        ws : workspace declared in calling routine
!             must be at least of length (ld*n)
!        lws: length of workspace ws as declared in calling program.
!
!     outputs:
!        x  : solution vector x
!        ier: error status flag
!             .eq. 0 => no errors
!             .lt. 0 => insufficient storage declared for array ws
!                       or num .gt. ld (both declared storage errors)
!             .gt. 0 => numerical error (singular matrix, etc.)
!                       value is equal to location of linear system
!                       that failed where 1 .le. ier .le. num
!
!     nb: the parameter ld is the leading dimension of all input and
!         output arrays, while n, the order of the linear systems, is
!         the second dimension.  this storage strategy is adopted for
!         efficient memory (i.e., stride 1) referencing when solving
!         multiple systems on a vector machine.  thus, num, the number
!         of systems to be solved, can take any value between 1 and ld.
!         Care should be taken defining a and c in the calling program.
!         Although the diagonal information a, b, and c, all are of
!         length n, the subdiagonal coefficients, a, are only defined
!         for n = 2, 3, 4, ... n, while the superdiagonal coefficients,
!         c, are only defined for n = 1, 2, 3, ... n-1 (i.e., the index
!         refers to the "row number" of the diagonal element and as such
!         the first element of a and last element of c do not exist).
!
!         since there is no pivoting, it is possible, but unlikely that
!         this routine could fail due to numerical instability, even if
!         the matrix a is nonsingular.  there is no attempt to diagnose
!         the nature of a failure if one is encountered in the procedure
! 
! Author: J. Hack
! 
!-----------------------------------------------------------------------
!
! $Id: mtdlss.F90,v 1.1.8.1 2002/06/15 13:48:55 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)    :: n  
  integer , intent(in)    :: ld  
  integer , intent(in)    :: num  
  integer , intent(in)    :: lws 
  integer , intent(in)    :: indx(ld)
  real(r8), intent(in)    :: a(ld,n)  
  real(r8), intent(in)    :: b(ld,n)  
  real(r8), intent(in)    :: c(ld,n)  
  real(r8), intent(in)    :: rhs(ld,n)  
  real(r8), intent(inout) :: x(ld,n)  
  real(r8), intent(inout) :: ws(ld,n) 
  integer , intent(out)   :: ier
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,ii, j          !indices
!-----------------------------------------------------------------------
!
  ier = 0
!
! Check for sufficient working storage in workspace ws and for
! logically consistant storage declaration for input arrays
!
  if ((lws < (ld*n)) .or. (num > ld)) then
     ier = -1
     return
  end if
!
! Decomposition and forward substitution loops
! Note: references to ws(i,1) denote special use of workspace ws
!
!CDIR$ IVDEP
  do ii=1,num
     i = indx(ii)
     if (b(i,1) == 0.0) then
        ier = i
        return
     end if
  end do

  do ii=1,num
     i = indx(ii)
     ws(i,1) = b(i,1)
     x (i,1) = rhs(i,1)/ws(i,1)
  end do

  do j=2,n
!CDIR$ IVDEP
     do ii=1,num
        i = indx(ii)
        ws(i,j) = c(i,j-1)/ws(i,1)
        ws(i,1) = b(i,j) - a(i,j)*ws(i,j)
     end do
!CDIR$ IVDEP
     do ii=1,num
        i = indx(ii)
        if (ws(i,1) == 0.0) then
           ier = i
           return
        end if
     end do
!CDIR$ IVDEP
     do ii=1,num
        i = indx(ii)
        x(i,j) = (rhs(i,j) - a(i,j)*x(i,j-1))/ws(i,1)
     end do
  end do
!
! Backsubstitution loop
!
  do j=n-1,1,-1
!CDIR$ IVDEP
     do ii=1,num
        i = indx(ii)
        x(i,j) = x(i,j) - ws(i,j+1)*x(i,j+1)
     end do
  end do

  return
end subroutine mtdlss

