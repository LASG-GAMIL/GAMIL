#include <misc.h>
#include <params.h>

subroutine inti 
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set constants and call initialization procedures for time independent
! physics routines
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,             only: plev, plevp, plevstd,plond,plat ! Needed for hypm passed to vd_inti
   use chemistry,          only: chem_initialize
   use ppgrid,           only: begchunk, endchunk, pcols  !! sxj--
   use physics_types,    only: physics_state,physics_state_set_grid  !1 sxj--
   use gw_drag,            only: gw_inti
   use vertical_diffusion, only: vd_inti
   use moistconvection,    only: mfinti,convection_scheme 
   use cldwat,             only: inimc
   use zm_conv,            only: zm_convi
   use zm_conv_3,          only: zm_convi_3
   use shr_const_mod,      only: shr_const_zvir, shr_const_cpwv, shr_const_rwv
   use physconst,          only: rair, cpair, cpwv, gravit, stebol, epsilo, tmelt, &
                                 latvap, latice, rh2o, zvir, cpvir, rhoh2o, pstd,  &
                                 karman, rhodair
   use cloudsimulator,     only: doisccp, cloudsimulator_init  !!(wh 2005.01.28,following cam3.0)
   use MG,                 only:stratiform_init !! sxj 2008-11-10
   use aerosol_mass_interface ,only:aerosol_mass_init !! sxj-2009-03-09
   use prescribed_aerosols,    only:aerosol_initialize !! sxj-2009-03-09

    ! *** added by DONG Li *** !
    use pmgrid,       only: masterproc
#if (defined TRIAL_RUN_FORCING)
    use SolarForcing, only: SolarForcing_Init
    use GHGForcing,   only: GHGForcing_Init
#endif
    ! ************************ !

   implicit none

#include <comctl.h>
#include <comhyb.h>
#include <RK_or_MG.h> !sxj

 ! the following code copy form cam3.5.8  
  ! !! Input/output arguments
   type(physics_state), pointer :: phys_state(:)
   ! local variables
   integer :: lchnk
   if (RK_or_MG=='MG') then
      allocate(phys_state(begchunk:endchunk))
       ! Set chunk id, number of columns, and coordinates
      do lchnk = begchunk,endchunk
         call physics_state_set_grid(lchnk, phys_state(lchnk))
      end do
   endif
!
!-----------------------------------------------------------------------
!
! Initialize physconst variables
! In adiabatic case, set zvir and cpvir explicitly to zero instead of 
! computing as (rh2o/rair - 1.) and (cpwv/cpair - 1.) respectively, in order 
! to guarantee an identical zero.
!
   if (adiabatic .or. ideal_phys) then
      rh2o  = rair
      zvir  = 0.
      cpwv  = cpair
      cpvir = 0.
   else
      rh2o  = shr_const_rwv
      zvir  = shr_const_zvir
      cpwv  = shr_const_cpwv
      cpvir = cpwv/cpair - 1.
   end if
!
! Call time independent initialization routines for parameterizations.
!
   if (trace_gas) call chem_initialize
   call gw_inti (cpair   ,cpwv    ,gravit  ,rair    ,hypi    )
   call vd_inti (cpair   ,cpwv    ,gravit  ,rair    ,zvir   , &
                 hypm    ,karman    )
   call tsinti  (tmelt   ,latvap  ,rair    ,stebol  ,latice  )
   call radini  (gravit  ,cpair   ,epsilo  ,stebol  ,pstd*10.0 )
   call esinti  (epsilo  ,latvap  ,latice  ,rh2o    ,cpair  , &
                 tmelt   )
   call mfinti  (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )
   call zm_convi( tmelt, epsilo, latvap, cpair )
   if (convection_scheme == 'Zhang_Neale') call zm_convi_3(hypi)   !!sxj
   call cldinti ()
    !
    ! initialization routine for prognostic cloud water
    !
    if (RK_or_MG=='MG') then
        call aerosol_initialize(phys_state(begchunk:endchunk))  !! sxj 2009-03-09
        if (masterproc) write(*, "('Notice: aerosol_initialize: Finished.')")
        call stratiform_init                                    !! sxj 2008-11-10
        if (masterproc) write(*, "('Notice: stratiform_init: Finshed.')")
        call aerosol_mass_init()                                !! sxj 2009-03-09
        if (masterproc) write(*, "('Notice: aerosol_mass_init: Finshed.')")
    else if (RK_or_MG=='RK') then
        call inimc(tmelt, rhodair/1000.0, gravit, rh2o)   
    end if
    !
    ! initialization for ISCCP Cloud Simulator
    !
    if (doisccp) call cloudsimulator_init             !!(wh 2005.01.28,following cam3.0)

#if (defined TRIAL_RUN_FORCING)
    call SolarForcing_Init
    call GHGForcing_Init
#endif

   return
end subroutine inti
