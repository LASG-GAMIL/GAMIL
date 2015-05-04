#include <misc.h>
#include <params.h>
subroutine cldclw(lchnk   ,ncol    ,zi      ,clwp    ,tpw     ,hl      )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Evaluate cloud liquid water path clwp (g/m**2)
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
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: zi(pcols,pverp)      ! height at layer interfaces(m)
   real(r8), intent(in) :: tpw(pcols)           ! total precipitable water (mm)
!
! Output arguments
!
   real(r8) clwp(pcols,pver)     ! cloud liquid water path (g/m**2)
   real(r8) hl(pcols)            ! liquid water scale height
   real(r8) rhl(pcols)           ! 1/hl

!
!---------------------------Local workspace-----------------------------
!
   integer i,k               ! longitude, level indices
   real(r8) clwc0                ! reference liquid water concentration (g/m**3)
   real(r8) emziohl(pcols,pverp) ! exp(-zi/hl)
!
!-----------------------------------------------------------------------
!
! Set reference liquid water concentration
!
   clwc0 = 0.21
!
! Diagnose liquid water scale height from precipitable water
!
   do i=1,ncol
      hl(i)  = 700.0*log(max(tpw(i)+1.0_r8,1.0_r8))
      rhl(i) = 1.0/hl(i)
   end do
!
! Evaluate cloud liquid water path (vertical integral of exponential fn)
!
   do k=1,pverp
      do i=1,ncol
         emziohl(i,k) = exp(-zi(i,k)*rhl(i))
      end do
   end do
   do k=1,pver
      do i=1,ncol
         clwp(i,k) = clwc0*hl(i)*(emziohl(i,k+1) - emziohl(i,k))
      end do
   end do
!
   return
end subroutine cldclw

