#include <misc.h>
#include <params.h>
subroutine radinp(lchnk   ,ncol    ,                            &
                  pmid    ,pint    ,o3vmr   , pmidrd  ,&
                  pintrd  ,eccf    ,o3mmr   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
! Convert model pressures to cgs, and compute ozone mixing ratio, needed for
! the solar radiation.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use shr_orb_mod
   use time_manager, only: get_curr_calday

   implicit none

#include <crdcon.h>
#include <comsol.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)
   real(r8), intent(in) :: o3vmr(pcols,pver)   ! ozone volume mixing ratio
!
! Output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor
   real(r8), intent(out) :: o3mmr(pcols,pver)   ! Ozone mass mixing ratio

!
!---------------------------Local variables-----------------------------
!
   integer i                ! Longitude loop index
   integer k                ! Vertical loop index

   real(r8) :: calday           ! current calendar day
   real(r8) amd                 ! Effective molecular weight of dry air (g/mol)
   real(r8) amo                 ! Molecular weight of ozone (g/mol)
   real(r8) vmmr                ! Ozone volume mixing ratio
   real(r8) delta               ! Solar declination angle

   save     amd   ,amo

   data amd   /  28.9644   /
   data amo   /  48.0000   /
!
!-----------------------------------------------------------------------
!
   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

!
! Convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0
         pintrd(i,k) = pint(i,k)*10.0
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0
   end do
!
! Convert ozone volume mixing ratio to mass mixing ratio:
!
   vmmr = amo/amd
   do k=1,pver
      do i=1,ncol
         o3mmr(i,k) = vmmr*o3vmr(i,k)
      end do
   end do
!
   return
end subroutine radinp
