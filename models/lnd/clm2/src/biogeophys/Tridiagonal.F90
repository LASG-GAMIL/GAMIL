#include <misc.h>
#include <preproc.h>

subroutine Tridiagonal (n, a, b, c, r, u )

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                           M  available land surface process model.
!  M --COMMUNITY LAND MODEL--  C
!  C                           L
!  LMCLMCLMCLMCLMCLMCLMCLMCLMCLM
!
!-----------------------------------------------------------------------
! Purpose:
! Solve tridiagonal system of equations
!
! Method:
!
! Author
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
! April 2002: Vertenstein/Oleson/Levis; Final form
!
!-----------------------------------------------------------------------
! $Id: Tridiagonal.F90,v 1.1.10.3 2002/06/15 13:50:20 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!----Arguments----------------------------------------------------------

  integer, intent(in) :: n
  real(r8), intent(in) :: a(1:n), b(1:n), c(1:n), r(1:n)

  real(r8), intent(out) :: u(1:n)

!----Local Variables----------------------------------------------------

  integer j
  real(r8) gam(1:n)
  real(r8) bet

!----End Variable List--------------------------------------------------

  bet = b(1)
  u(1) = r(1) / bet
  do j = 2, n
     gam(j) = c(j-1) / bet
     bet = b(j) - a(j) * gam(j)
     u(j) = (r(j) - a(j)*u(j-1)) / bet
  enddo

  do j = n-1, 1, -1
     u(j) = u(j) - gam(j+1) * u(j+1)
  enddo

end subroutine Tridiagonal
