#include <misc.h>
#include <params.h>

subroutine tphysac (ztodt,   pblh,    qpert,   tpert,  shf,  &
                    taux,    tauy,    cflx,    sgh,    lhf,  &
                    landfrac,snowh,   tref,    precc,  precl,  &
                    tin,     state,   tend,    ocnfrac, cldn)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics after coupling to land, sea, and ice models.
! Computes the following:
!   o Radon surface flux and decay (optional)
!   o Vertical diffusion and planetary boundary layer
!   o Dry deposition for sulfur cycle (optional)
!   o Multiple gravity wave drag
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,             only: pcols, pver
   use chemistry,          only: chem_driver
   use gw_drag,            only: gw_intr
   use vertical_diffusion, only: vd_intr
   use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use constituents,       only: ppcnst, qmin
   use tracers,            only: ixtrct
   use physconst,          only: zvir, gravit, rhoh2o, latvap

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
   real(r8), intent(in) :: landfrac(pcols)        ! Land fraction
   real(r8), intent(in) :: ocnfrac(pcols)         ! Land fraction
   real(r8), intent(in) :: snowh(pcols)           ! snow depth (liquid water equivalent)

   real(r8), intent(in) :: tref(pcols)            ! 2m air temperature
   real(r8), intent(in) :: precc(pcols)           ! convective precipitation
   real(r8), intent(in) :: precl(pcols)           ! large-scale precipitation
   real(r8), intent(in) :: cldn(pcols,pver)       ! cloud fraction [fraction]
   
   real(r8), intent(out) :: pblh(pcols)           ! Planetary boundary layer height
   real(r8), intent(out) :: qpert(pcols,ppcnst)   ! Moisture/constit. perturbation (PBL)
   real(r8), intent(out) :: tpert(pcols)          ! Temperature perturbation (PBL)
   real(r8), intent(inout) :: shf(pcols)          ! Sensible heat flux (w/m^2)
   real(r8), intent(in) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(in) :: tauy(pcols)            ! Y surface stress (meridional)
   real(r8), intent(inout) :: cflx(pcols,ppcnst)  ! Surface constituent flux (kg/m^2/s)
   real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
   real(r8), intent(inout) :: lhf(pcols)          ! Latent heat flux (w/m^2)
   real(r8), intent(in) :: tin(pcols, pver) ! input T, to compute FV output T

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend) ptend                ! indivdual parameterization tendencies

   integer lchnk                            ! chunk identifier
   integer ncol                             ! number of atmospheric columns
   integer i

   logical labort                           ! abort flag

   real(r8) tvm(pcols,pver)           ! virtual temperature
   real(r8) prect(pcols)              ! total precipitation
   real(r8) surfric(pcols)              ! surface friction velocity
   real(r8) obklen(pcols)             ! Obukhov length

   real(r8) we(pcols)                ! entrainment velocity
   real(r8) ehf(pcols)               ! entrainment heat flux
   real(r8) esf(pcols)               ! entrainment dse flux
   real(r8) euf(pcols)               ! entrainment U flux
   real(r8) evf(pcols)               ! entrainment V flux
   real(r8) eqf(pcols,ppcnst)        ! entrainment constituent flux
   real(r8) khfs_tot(pcols,pver+1)           
   real(r8) kqfs_tot(pcols,pver+1)          
!
!-----------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol
!
! accumulate fluxes into net flux array
!
   do i=1,ncol
      tend%flx_net(i) = tend%flx_net(i) + shf(i) + (precc(i) + precl(i))*latvap*rhoh2o
   end do

! Convert mixing ratio of non-water tracers to mass fraction of total
! atmospheric mass (Overwrite non-water portions of q).

   if (ppcnst > 1) then
      call mr2mf (lchnk, ncol, state%q)
   end if

! Initialize parameterization tendency structure

   call physics_ptend_init(ptend)

! Check if latent heat flux exceeds the total moisture content of the
! lowest model layer, thereby creating negative moisture.

   call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,          &
              state%q(1,pver,1),state%rpdel(1,pver) ,shf ,lhf ,cflx(1,1) )

!===================================================
! Source/sink terms for advected tracers.
!===================================================

   if ( trace_test1 .or. trace_test2 .or. trace_test3 ) then
      write(*,*) '!! trace_test'
      call rnsfwcrp( lchnk, ncol, landfrac, cflx(:,ixtrct))
      call rndecay( lchnk, ncol, state%q(:,:,ixtrct), ztodt, ptend%q(:,:,ixtrct))
      ptend%lq(ixtrct) = .TRUE.
      if (trace_test3) state%q(:ncol,pver,ixtrct+2) =  0.
   end if

! Advected greenhouse trace gases:

   if (trace_gas) call chem_driver (state, ptend, cflx, ztodt)
   if (trace_gas) write(*,*) '!!  trace_gas = .true.'

! Add tendencies to cummulative model tendencies and update profiles

   call physics_update (state, tend, ptend, ztodt)

!===================================================
! Vertical diffusion/pbl calculation
! Call vertical diffusion code (pbl, free atmosphere and molecular)
!===================================================

   call vd_intr (ztodt    ,state    ,taux     ,tauy     , shf    ,&
                 cflx     ,pblh     ,tpert    ,qpert    , we     ,&
                 ehf      ,eqf      ,euf      ,evf      , esf    ,surfric ,&
                 obklen   ,ptend    ,khfs_tot ,kqfs_tot ,cldn)

   call physics_update (state, tend, ptend, ztodt)

!===================================================
! Gravity wave drag
!===================================================

   call gw_intr(state, sgh, pblh, ztodt, ptend)
   call physics_update(state, tend, ptend, ztodt)

!*** BAB's FV kludge

   state%t(:ncol,:pver) = tin(:ncol,:pver) + ztodt*tend%dtdt(:ncol,:pver)

   if (aqua_planet) then
      labort = .false.
      do i = 1, ncol
         if (ocnfrac(i) /= 1.) labort = .true.
      end do
      if (labort) then
         write(6,*) 'ERROR:  grid contains non-ocean point'
         call endrun ()
      end if
   end if
!
! Convert mass fractions of non-water tracers back to mixing ratios.
! (Overwrite non-water portions of q).
!
   if (ppcnst > 1) then
      call mf2mr (lchnk, ncol, state%q)
   end if

   return
end subroutine tphysac
