#include <misc.h>
#include <params.h>
subroutine srfxfer(lchnk   ,ncol    ,psm1tmp ,um1     ,vm1     ,tm1     , &
                   qm1     ,exner   ,zm      ,pmidm1  ,rpdel   )
!-----------------------------------------------------------------------
!
! Purpose:
! Transfer atmospheric fields into common block /comsrf/
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: L. Bath  CMS Contact: M. Vertenstein
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use comsrf
#ifdef COUP_CSM
   use ccsm_msg, only: rho, netsw
#endif
   implicit none

#include <comtsc.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! Chunk index
   integer, intent(in) :: ncol
!
   real(r8), intent(in) :: um1(pcols)           ! Bottom level u wind
   real(r8), intent(in) :: vm1(pcols)           ! Bottom level v wind
   real(r8), intent(in) :: tm1(pcols)           ! Bottom level temperature
   real(r8), intent(in) :: qm1(pcols)           ! Bottom level specific humidity
   real(r8), intent(in) :: exner(pcols)         ! Bottom level Exner function
   real(r8), intent(in) :: zm(pcols)            ! Bottom level height above surface
!
   real(r8), intent(in) :: psm1tmp(pcols)       ! Surface pressure
   real(r8), intent(in) :: pmidm1(pcols,pver)   ! Level pressures
   real(r8), intent(in) :: rpdel(pcols)         ! 1./(pint(k+1)-pint(k))

!
!---------------------------Local variables-----------------------------
!
   integer i                 ! Longitude index
!
!-----------------------------------------------------------------------
!
! Stuff global fluxes and state variables into common
!
   do i=1,ncol
      surface_state2d(lchnk)%tbot(i) = tm1(i)
      surface_state2d(lchnk)%thbot(i) = tm1(i) * exner(i)
      surface_state2d(lchnk)%zbot(i) = zm(i)
      surface_state2d(lchnk)%ubot(i) = um1(i)
      surface_state2d(lchnk)%vbot(i) = vm1(i)
      surface_state2d(lchnk)%qbot(i) = qm1(i)
      surface_state2d(lchnk)%pbot(i) = pmidm1(i,pver)
      psm1(i,lchnk) = psm1tmp(i)
      srfrpdel(i,lchnk) = rpdel(i)
#ifdef COUP_CSM
      rho(i,lchnk)   = pmidm1(i,pver)/(rair*surface_state2d(lchnk)%tbot(i))
      netsw(i,lchnk) = surface_state2d(lchnk)%srfrad(i) - &
	surface_state2d(lchnk)%flwds(i)
#endif
   end do

   return
end subroutine srfxfer







