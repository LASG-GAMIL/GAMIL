#include <misc.h>
#include <params.h>
#define PCWDETRAIN
subroutine q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,nstep   ,          &
                    dl      ,evp     ,cu      )

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
! 
!-----------------------------------------------------------------------


   real(r8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg
   integer, intent(in) :: nstep

   real(r8), intent(in) :: q(pcols,pver)
   real(r8), intent(in) :: qs(pcols,pver)
   real(r8), intent(in) :: qu(pcols,pver)
   real(r8), intent(in) :: su(pcols,pver)
   real(r8), intent(in) :: du(pcols,pver)
   real(r8), intent(in) :: qhat(pcols,pver)
   real(r8), intent(in) :: shat(pcols,pver)
   real(r8), intent(in) :: dp(pcols,pver)
   real(r8), intent(in) :: mu(pcols,pver)
   real(r8), intent(in) :: md(pcols,pver)
   real(r8), intent(in) :: sd(pcols,pver)
   real(r8), intent(in) :: qd(pcols,pver)
   real(r8), intent(in) :: ql(pcols,pver)
   real(r8), intent(in) :: evp(pcols,pver)
   real(r8), intent(in) :: cu(pcols,pver)
   real(r8), intent(in) :: dsubcld(pcols)

   real(r8),intent(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
   real(r8),intent(out) :: dl(pcols,pver)
   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
!
! work fields:
!
   integer i
   integer k

   real(r8) fact
   real(r8) emc
   real(r8) rl
!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0.
         dqdt(i,k) = 0.
         dl(i,k) = 0.
      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do
!
   fact = 0.
!      fact = 1.

   do k = ktm,pver-1
      do i = il1g,il2g
#ifndef PCWDETRAIN
! cludge to make it look like the standard cam version of convection
! detrain all water into environment till 80% rh
! then make remaining water fall out as precip
         fact = 1.
         if (q(i,k) > 0.8*qs(i,k) .and. k < pver-3) fact = 0.
#endif
         emc = +fact*du(i,k)*ql(i,k+1) &         ! evaporating cloud detraining to env
               -cu(i,k)                &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft
!            emc = 0

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = (1-fact)*du(i,k)*ql(i,k+1)

      end do
   end do

!
   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(i,k) = (1./dsubcld(i))* &
                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        )
            dqdt(i,k) = (1./dsubcld(i))* &
                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        )
         else if (k > mx(i)) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do
!
   return
end subroutine q1q2_pjr
