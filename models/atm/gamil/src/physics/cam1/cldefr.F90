#include <misc.h>
#include <params.h>
subroutine cldefr(lchnk   ,ncol    , &
                  landfrac,t       ,rel     ,rei     ,fice    , &
                  ps      ,pmid    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud drop size
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
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: landfrac(pcols)      ! Land fraction
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature
   real(r8), intent(in) :: ps(pcols)            ! Surface pressure
   real(r8), intent(in) :: pmid(pcols,pver)     ! Midpoint pressures
!
! Output arguments
!
   real(r8), intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
   real(r8), intent(out) :: rei(pcols,pver)      ! Ice effective drop size (microns)
   real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
   real(r8) pirnge                               ! Nrmlzd pres range for ice particle changes
   real(r8) picemn                               ! Normalized pressure below which rei=reimax
   real(r8) rirnge                               ! Range of ice radii (reimax - 10 microns)
   real(r8) reimax                               ! Maximum ice effective radius
   real(r8) pnrml                                ! Normalized pressure
   real(r8) weight                               ! Coef. for determining rei as fn of P/PS
!
!---------------------------Local workspace-----------------------------
!
   integer i,k               ! Lon, lev indices
   real(r8) rliq                 ! Temporary liquid drop size
!
!--------------------------Statement functions--------------------------
!
   logical land
   land(i) = nint(landfrac(i)) == 1
!
!-----------------------------------------------------------------------
!
   do k=1,pver
      do i=1,ncol
!
! Define liquid drop radius
!
         if (land(i)) then
            rliq = 7.0 + 3.0* min(1.0_r8,max(0.0_r8,(263.16-t(i,k))*0.05))
         else
            rliq = 10.0  !(ljli1208)
!            rliq = 20.0 !(ljli04052007)
         end if
         rel(i,k) = rliq   !!(ljli0730)
!ljli20090830          rel(i,k) = 2.0*rliq
!
! Determine rei as function of normalized pressure
!
!         reimax   = 30.0
         reimax   = 40.0
         rirnge   = 20.0
         pirnge   = 0.4
         picemn   = 0.4
!
         pnrml    = pmid(i,k)/ps(i)
         weight   = max(min((pnrml-picemn)/pirnge,1.0_r8),0._r8)
          rei(i,k) = reimax - rirnge*weight
!ljli20090830          rei(i,k) = 2*(reimax - rirnge*weight)
!
! Define fractional amount of cloud that is ice
!
! If warmer than -10 degrees C then water phase
!
         if (t(i,k) > 263.16) fice(i,k) = 0.0
!(ljli07-17)          if (t(i,k) > 273.16) fice(i,k) = 0.0  !(ljli-07-13)
!
! If colder than -10 degrees C but warmer than -30 C mixed phase
!
!(ljli2006-07-17)          if(t(i,k) < 273.16) then
!(07-17)   fice(i,k) = 1.0-(0.0059+(1-0.0059)*exp(-0.003102*(t(i,k)-273.16)**2))
!          end if
         if (t(i,k) <= 263.16 .and. t(i,k) >= 243.16) then
            fice(i,k) =(263.16-t(i,k)) / 20.0
         end if
!
!(ljli0718)   if (t(i,k) <= 263.16 .and. t(i,k) >= 233.16) then !(ljli07-17)
!             fice(i,k) = (263.16-t(i,k)) /30.0
!           end if
! If colder than -30 degrees C then ice phase
!
         if (t(i,k) < 243.16) fice(i,k) = 1.0
!
!!(ljli2006-07-18)              if (t(i,k) < 233.16) fice(i,k) = 1.0
! Turn off ice radiative properties by setting fice = 0.0
!
!+             fice(i,k) = 0.0
!
      end do
   end do
!
   return
end subroutine cldefr

