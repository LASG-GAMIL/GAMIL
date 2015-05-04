#include <misc.h>
#include <params.h>

subroutine mr2mf( lchnk, ncol, q3 )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convert mixing ratio of the non-water trace species to mass fraction
! of total atmospheric mass.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: mr2mf.F90,v 1.1.8.1 2002/06/15 13:49:17 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use dycore,  only: dycore_is
   use tracers, only: pcnst, pnats

   implicit none

!------------------------------Arguments--------------------------------
   integer, intent(in) :: lchnk                          ! chunk identifier
   integer, intent(in) :: ncol                           ! number of atmospheric columns

   real(r8), intent(inout) :: q3(pcols,pver,pcnst+pnats) ! Constituent mass fraction
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
   integer :: m             ! constituent indices
!-----------------------------------------------------------------------
!
! LR: Do nothing here as suggested by Phil Rasch                                        
!
   if (.not. dycore_is ('LR')) then
      do m = 2,pcnst+pnats
         q3(:ncol,:pver,m) = q3(:ncol,:pver,m)*(1. - q3(:ncol,:pver,1))
      end do
   end if
!
   return
end subroutine mr2mf
 
