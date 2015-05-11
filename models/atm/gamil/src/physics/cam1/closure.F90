#include <misc.h>
#include <params.h>
subroutine closure(lchnk   , &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,capelmt ,nstep   ,ptendlaq  ,ptendlat  , &
                   negcape)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,    only: pcols, pver

   implicit none

#include <guang.h>
!
!-----------------------------Arguments---------------------------------
!
   integer, intent(in) :: lchnk                 ! chunk identifier

   real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
   real(r8), intent(inout) :: t(pcols,pver)        ! temperature
   real(r8), intent(inout) :: p(pcols,pver)        ! pressure (mb)
   real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
   real(r8), intent(in) :: z(pcols,pver)        ! height (m)
   real(r8), intent(in) :: s(pcols,pver)        ! normalized dry static energy
   real(r8), intent(in) :: tp(pcols,pver)       ! parcel temp
   real(r8), intent(in) :: qs(pcols,pver)       ! sat spec humidity
   real(r8), intent(in) :: qu(pcols,pver)       ! updraft spec. humidity
   real(r8), intent(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   real(r8), intent(in) :: mc(pcols,pver)       ! net convective mass flux
   real(r8), intent(in) :: du(pcols,pver)       ! detrainment from updraft
   real(r8), intent(in) :: mu(pcols,pver)       ! mass flux of updraft
   real(r8), intent(in) :: md(pcols,pver)       ! mass flux of downdraft
   real(r8), intent(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
   real(r8), intent(in) :: sd(pcols,pver)       ! dry static energy of downdraft
   real(r8), intent(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
   real(r8), intent(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
   real(r8), intent(in) :: dp(pcols,pver)       ! pressure thickness of layers
   real(r8), intent(in) :: qstp(pcols,pver)     ! spec humidity of parcel
   real(r8), intent(in) :: zf(pcols,pver+1)     ! height of interface levels
   real(r8), intent(in) :: ql(pcols,pver)       ! liquid water mixing ratio
   real(r8), intent(in) :: ptendlaq(pcols,pver)   ! q tendency due to large-scale forcing
   real(r8), intent(in) :: ptendlat(pcols,pver)   ! t tendency due to large-scale forcing

   real(r8), intent(in) :: cape(pcols)          ! available pot. energy of column
   real(r8), intent(in) :: negcape(pcols)       ! negative available pot. energy of column
   real(r8), intent(in) :: tl(pcols)
   real(r8), intent(in) :: dsubcld(pcols)       ! thickness of subcloud layer

   integer, intent(in) :: lcl(pcols)        ! index of lcl
   integer, intent(in) :: lel(pcols)        ! index of launch leve
   integer, intent(in) :: jt(pcols)         ! top of updraft
   integer, intent(in) :: mx(pcols)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real(r8) dtmdt(pcols,pver)
   real(r8) dqmdt(pcols,pver)
   real(r8) dboydt(pcols,pver)
   real(r8) thetavp(pcols,pver)
   real(r8) thetavm(pcols,pver)
   real(r8) dbdtls(pcols,pver)

   real(r8) rhmx(pcols)
   real(r8) beta
   real(r8) capelmt
   real(r8) cp
   real(r8) dadt(pcols)
   real(r8) dltaa
   real(r8) grav
   real(r8) rhcrit
   real(r8) dadtls

   integer i
   integer il1g
   integer il2g
   integer k
   integer msg
   integer nstep

   real(r8) rd
   real(r8) rl
   real(r8) tau
   save tau
!
! tau=4800. were used in canadian climate center. however, when it
! is used here in echam3 t42, convection is too weak, thus
! adjusted to 2400. i.e the e-folding time is 1 hour now.
!
!  data tau/7200./
   data tau/3600./  !!(wh 2004.03)
!-----------------------------------------------------------------------
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
      mb(i) = 0.
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(i,k) = 0.
         dqmdt(i,k) = 0.
         dbdtls(i,k) = 0.
      end do
   end do
!
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(i,k) = (1./dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
            dqmdt(i,k) = (1./dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                         qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
         end if
      end do
   end do
!
   beta = 0.
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                         dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*(1.+1.608*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*(1.+0.608*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = - (dtmdt(i,k)/t(i,k)+0.608/ &
                          (1.+0.608*q(i,k))*dqmdt(i,k))*grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*(1.+0.608*q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*(1.+0.608*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (- dtmdt(i,k)/t(i,k)-0.608/ (1.+0.608*q(i,k))*dqmdt(i,k))* &
                          grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   do i = il1g,il2g
      dadt(i) = 0.
      do k = lel(i),mx(i) - 1
!       do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         endif
!       end do
      end do
   end do
!
!
!***********************************
! test a new closure
!***********************************
      do i = il1g,il2g
         mb(i) = 0.
      end do

      do k = msg + 1,pver
         do i = il1g,il2g
            if (k.ge.lel(i) .and. k.le.lcl(i)) then
               thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)* &
                               (1.+1.608*qstp(i,k)-q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)* &
                               (1.+0.608*q(i,k))
!
               dbdtls(i,k) = &
                      -( 1/t(i,k) * ptendlat(i,k) +  &
          0.608 / (1+0.608*q(i,k))* ptendlaq(i,k)) &
                    *grav* thetavp(i,k)/thetavm(i,k)

            end if
         end do
      end do
!
           do i = il1g,il2g
            rhmx(i)=q(i,mx(i))/qs(i,mx(i))
           enddo

      rhcrit=0.90
      do i = il1g,il2g
         dadtls = 0.
         do k = lel(i),mx(i) - 1
            dadtls = dadtls + dbdtls(i,k)* (zf(i,k)-zf(i,k+1))
         end do
         if (dadtls > 0. .and. dadt(i) < 0. ) &
                             mb(i) = max(-dadtls/dadt(i),0.)
         if ( negcape(i) < -400. .or. rhmx(i) < rhcrit ) mb(i) = 0.
      end do

   return
end subroutine closure

