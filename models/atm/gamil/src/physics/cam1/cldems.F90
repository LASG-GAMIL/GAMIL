#include <misc.h>
#include <params.h>
subroutine cldems(lchnk   ,ncol    ,clwp    ,fice    ,rei     ,emis    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud emissivity using cloud liquid water path (g/m**2)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J.T. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none
!
!------------------------------Parameters-------------------------------
!
   real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
   parameter (kabsl = 0.090361)
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: clwp(pcols,pver)       ! cloud liquid water path (g/m**2)
   real(r8), intent(in) :: rei(pcols,pver)        ! ice effective drop size (microns)
   real(r8), intent(in) :: fice(pcols,pver)       ! fractional ice content within cloud
!
! Output arguments
!
   real(r8), intent(out) :: emis(pcols,pver)       ! cloud emissivity (fraction)
!
!---------------------------Local workspace-----------------------------
!
   integer i,k                 ! longitude, level indices
   real(r8) kabs                   ! longwave absorption coeff (m**2/g)
   real(r8) kabsi                  ! ice absorption coefficient
!
!-----------------------------------------------------------------------
!
   do k=1,pver
      do i=1,ncol
         kabsi = 0.005 + 1./rei(i,k)
         kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
         emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
      end do
   end do
!
   return
end subroutine cldems

