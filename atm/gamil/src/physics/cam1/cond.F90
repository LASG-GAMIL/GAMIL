#include <misc.h>
#include <params.h>

subroutine cond (lchnk   ,ncol    , &
                 ztodt   ,pmid    ,pdel    ,ti      ,qi      , &
                 ds      ,dq      ,precl   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate large scale condensation
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1,  CMS Contact: J. Hack
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use physconst, only: epsilo, rga, rhoh2o, latvap, rh2o, cpair
   use wv_saturation, only: aqsat

   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: ztodt                  ! Physics time step (2 delta t)
   real(r8), intent(in) :: pmid(pcols,pver)       ! Pressure at layer midpoints
   real(r8), intent(in) :: pdel(pcols,pver)       ! Delta p at each model level
   real(r8), intent(in) :: ti(pcols,pver)          ! Temperature
   real(r8), intent(in) :: qi(pcols,pver)          ! Specific humidity
!
! Output arguments
!
   real(r8), intent(out) :: ds(pcols,pver)         ! heating rate (dry static energy tendency) due to rainout
   real(r8), intent(out) :: dq(pcols,pver)         ! Moisture tendency due to rainout
   real(r8), intent(out) :: precl(pcols)           ! Large-scale precipitation rate

!
!---------------------------Local variables-----------------------------
!
   real(r8) :: t(pcols,pver)          ! Temperature
   real(r8) :: q(pcols,pver)          ! Specific humidity
   real(r8) absqs                  ! Intermediate quantity
   real(r8) denom                  ! Intermediate quantity
   real(r8) dqsdt                  ! Change of qsat with respect to temp.
   real(r8) est(pcols,pver)        ! Saturation vapor pressure
   real(r8) omeps                  ! 1 - 0.622
   real(r8) qsat(pcols,pver)       ! Saturation specific humidity
   real(r8) rain(pcols)            ! Rain (units of kg/m^2 s)
   real(r8) rhm1                   ! RH - saturation RH
   real(r8) zqcd(pcols)            ! Intermed quantity (several actually)
   real(r8) zqdt                   ! Reciprocal of ztodt
   real(r8) ke                     ! `disposable parameter' in evaporation
   real(r8) evap                   ! Water evaporation rate
   real(r8) relhum                 ! Relative humidity
   real(r8) dpovrg                 ! deltap/grav
   integer i                       ! Longitude index
   integer jiter                   ! Iteration counter
   integer k                       ! Vertical index

   real(r8) cldcp            ! Latvap/cpair (L/cp)
   real(r8) clrh2o           ! Ratio of latvap to water vapor gas const
!
!-----------------------------------------------------------------------
!
   zqdt  = 1./ztodt
   omeps = 1. - epsilo
   clrh2o = latvap/rh2o   ! Ratio of latvap to water vapor gas const
   cldcp  = latvap/cpair
!
! First diagnose condensation rate due to stable processes
! Update column T and Q (evaporation process is `time-split')
! Return tendencies determined from updated profiles.
! Condensation calculation is hard-wired for two iterations
!
   t (:,:) = ti(:,:)
   q (:,:) = qi(:,:)
   ds(:,:) = 0.
   dq(:,:) = 0.
!
   do jiter=1,2
      call aqsat(t       ,pmid    ,est     ,qsat    ,pcols   , &
                 ncol    ,pver    ,1       ,pver    )
      do k=1,pver
!
! Calculate condensation-rate and new t- and q-values
!
         do i=1,ncol
!
! Use of critical saturation vapor pressure requires coefficient on the
! term omeps*est(i,k) in the next statement (e.g. omeps*est(i,k)*escrit)
! Corresponding changes must also be incorporated into estabv.for (e.g.,
! terms est(i,k) in qsat evaluation become escrit*est(i,k))
!
            denom   = (pmid(i,k) - omeps*est(i,k))*t(i,k)**2
            dqsdt   = clrh2o*qsat(i,k)*pmid(i,k)/denom
            absqs   = abs(qsat(i,k))
            rhm1    = q(i,k)/qsat(i,k) - 1.
            zqcd(i) = max(absqs*rhm1/(1. + cldcp*dqsdt),0._r8)
            if (q(i,k) < 0.0) zqcd(i) = 0.
            t     (i,k) = t     (i,k) + zqcd(i)*cldcp
            q     (i,k) = q     (i,k) - zqcd(i)
            dq    (i,k) = dq    (i,k) - zqcd(i)*zqdt
         end do
      end do
   end do
!
! Initialize rain vector (will be updated as rain falls through column)
!
   do i=1,ncol
      rain(i) = max(-dq(i,1)*pdel(i,1)*rga,0.0_r8)
   end do
   call aqsat(t       ,pmid    ,est     ,qsat    ,pcols   , &
              ncol    ,pver    ,1       ,pver    )
!
! Evaporate condensate on the way down (see Sundqvist, 1988: Physically
! Based Modelling ..., pp 433-461, Schlesinger, Ed., Kluwer Academic)
! variable evap has units of 1/s; variable rain has units of kg/m**2/s
! rain is used to accumuluate unevaporated rain water on the way down
!
   ke = 1.0e-5                     ! set in common block in final code
   do k=2,pver
      do i=1,ncol
         dpovrg  = pdel(i,k)*rga
         relhum  = q(i,k)/qsat(i,k)
         evap    = max(ke*(1.0 - relhum)*sqrt(rain(i)), 0.0_r8)
         evap    = min(evap, (qsat(i,k)-q(i,k))/ztodt)
         evap    = min(rain(i)/dpovrg,evap)
         dq(i,k) = dq(i,k) + evap
         q(i,k)  = q (i,k) + evap*ztodt
         t(i,k)  = t (i,k) - evap*ztodt*cldcp
         rain(i) = max(rain(i) - dq(i,k)*dpovrg,0.0_r8)
      end do
   end do

   ds(:,:) = -dq(:,:)*latvap

   do i=1,ncol
      precl(i) = rain(i)/rhoh2o
   end do
!
   return
end subroutine cond

