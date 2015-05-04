#include <misc.h>
#include <params.h>

subroutine tphysidl (ztodt   ,taux    ,tauy    ,etamid  , state   , &
                     tend)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!  algorithm 1: Held/Suarez IDEALIZED physics
!  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
!  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
!  algorithm 4: Boer/Denis  IDEALIZED physics
!
! Author: J. Olson
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid            , only: plev,plat,plevp,plevstd,plond,plat
   use ppgrid
   use phys_grid         , only: get_lat_all_p, get_rlat_all_p
   use vertical_diffusion, only: vd_intr
   use physics_types     , only: physics_state, physics_tend, physics_ptend
   use geopotential      , only: geopotential_t
   use history,            only: outfld
   use physconst,          only: gravit, cappa, rair
   use tracers,            only: pcnst, pnats

   implicit none

#include <comhyb.h>
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                   ! Two times model timestep (2 delta-t)
!
! Output arguments
!
   real(r8), intent(out) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(out) :: tauy(pcols)            ! Y surface stress (meridional)
   real(r8), intent(in)  :: etamid(pver)           ! midpoint values of eta (a+b)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns

   real(r8) clat(pcols)                        ! latitudes(radians) for columns
   real(r8) pmid(pcols,pver)                   ! mid-point pressure
   integer  i,k                                ! Longitude, level indices

   real(r8) tmp                                ! temporary
   real(r8) kf                                 ! 1./efolding_time for wind dissipation
   real(r8) ka                                 ! 1./efolding_time for temperature diss.
   real(r8) kaa                                ! 1./efolding_time for temperature diss.
   real(r8) ks                                 ! 1./efolding_time for temperature diss.
   real(r8) kv                                 ! 1./efolding_time (normalized) for wind
   real(r8) kt                                 ! 1./efolding_time for temperature diss.
   real(r8) trefa                              ! "radiative equilibrium" T
   real(r8) trefc                              ! used in calc of "radiative equilibrium" T
   real(r8) cossq(pcols)                       ! coslat**2
   real(r8) cossqsq(pcols)                     ! coslat**4
   real(r8) sinsq(pcols)                       ! sinlat**2
   real(r8) onemsig                            ! 1. - sigma_reference
   real(r8) efoldf                             ! efolding time for wind dissipation
   real(r8) efolda                             ! efolding time for T dissipation
   real(r8) efoldaa                            ! efolding time for T dissipation
   real(r8) efolds                             ! efolding time for T dissipation
   real(r8) efold_strat                        ! efolding time for T dissipation in Strat
   real(r8) efold_meso                         ! efolding time for T dissipation in Meso
   real(r8) efoldv                             ! efolding time for wind dissipation
   real(r8) p_infint                           ! effective top of model
   real(r8) constw                             ! constant
   real(r8) lapsew                             ! lapse rate
   real(r8) p0strat                            ! threshold pressure
   real(r8) phi0                               ! threshold latitude
   real(r8) dphi0                              ! del-latitude
   real(r8) a0                                 ! coefficient
   real(r8) aeq                                ! 100 mb
   real(r8) apole                              ! 2   mb
   real(r8) pi                                 ! 3.14159...
   real(r8) coslat(pcols)                      ! cosine(latitude)
   real(r8) acoslat                            ! abs(acos(coslat))
   real(r8) constc                             ! constant
   real(r8) lapsec                             ! lapse rate
   real(r8) lapse                              ! lapse rate
   real(r8) h0                                 ! scale height (7 km)
   real(r8) sigmab                             ! threshold sigma level
   real(r8) pressmb                            ! model pressure in mb
   real(r8) t000                                ! minimum reference temperature
   integer  idlflag                            ! Flag to choose which idealized physics
!
!-----------------------------------------------------------------------
!
   idlflag = 1

   lchnk = state%lchnk
   ncol  = state%ncol
!
! Copy pressures into local array
!
   call get_rlat_all_p(lchnk, ncol, clat)
   do i=1,ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
   end do
   do k=1,pver
      do i=1,ncol
         pmid(i,k) = state%pmid(i,k)
      end do
   end do

   if (idlflag == 1) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf =  1.
      efolda = 40.
      efolds =  4.
      sigmab =  0.7
      t000    = 200.
!
      onemsig = 1. - sigmab
!
      ka = 1./(86400.*efolda)
      ks = 1./(86400.*efolds)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
!!             tmp   = kt/(1.+ ztodt*kt)  !(wanhui 2003.06.16)
               tmp   = kt
               trefc   = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa    = max(t000,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
!!          tmp   = ka/(1.+ ztodt*ka)     !(wanhui 2003.06.16)
            tmp   = ka
            do i=1,ncol
               trefc   = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa    = max(t000,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
!!          tmp = -kv/(1.+ ztodt*kv)               ! (wanhui 2003.06.16)
            tmp = -kv
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 2) then
!
!-----------------------------------------------------------------------
!
! Modified Held/Suarez IDEALIZED physics algorithm
! (modified with Williamson stratosphere):
!
!   Williamson, D. L., J. G. Olson and B. A. Boville, 1998: A comparison
!   of semi--Lagrangian and Eulerian tropical climate simulations.
!   Mon. Wea. Rev., vol 126, pp. 1001-1012.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf  =  1.
      efolda  = 40.
      efoldaa = 40.
      efolds  =  4.
      sigmab  =  0.7
      t000     = 200.
!
      onemsig = 1. - sigmab
!
      ka  = 1./(86400.*efolda)
      kaa = 1./(86400.*efoldaa)
      ks  = 1./(86400.*efolds)
!
      pi     = 4.*atan(1.)
      phi0   = 60.*pi/180.
      dphi0  = 15.*pi/180.
      a0     = 2.65/dphi0
      aeq    = 10000.
      apole  = 200.
      lapsew = -3.345e-03
      constw = rair*lapsew/gravit
      lapsec =  2.00e-03
      constc = rair*lapsec/gravit
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5*(1. + tanh(a0*(acoslat - phi0)))
               tmp     = kt/(1.+ ztodt*kt)
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t000,trefa)
               if (pmid(i,k) < 10000.) then
                  trefa = t000*((pmid(i,k)/10000.))**constc
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t000*( ((pmid(i,k)/p0strat))**constw - 1. )
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
            do i=1,ncol
               acoslat = abs(acos(coslat(i)))
               p0strat = aeq - (aeq - apole)*0.5*(1. + tanh(a0*(acoslat - phi0)))
               tmp     = ka/(1.+ ztodt*ka)
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t000,trefa)
               if (pmid(i,k) < 10000.) then
                  trefa = t000*((pmid(i,k)/10000.))**constc
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               if (pmid(i,k) < p0strat) then
                  trefa = trefa + t000*( ((pmid(i,k)/p0strat))**constw - 1. )
                  tmp   = kaa/(1.+ ztodt*kaa)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         endif
      end do

   elseif (idlflag == 3) then
!
!-----------------------------------------------------------------------
!
! Held/Suarez IDEALIZED physics algorithm:
! (modified with Lin/Williamson stratosphere/mesosphere):
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!
!-----------------------------------------------------------------------
!
! Add idealized radiative heating rates to temperature tendency
!
      efoldf      =  1.
      efolda      = 40.
      efolds      =  4.
      efold_strat = 40.
      efold_meso  = 10.
      efoldv      = 0.5
      sigmab      = 0.7
      lapse       = 0.00225
      h0          = 7000.
      t000         = 200.
      p_infint    = 0.01
!
      onemsig = 1. - sigmab
!
      ka = 1./(86400.*efolda)
      ks = 1./(86400.*efolds)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            do i=1,ncol
               kt    = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig
               tmp   = kt/(1.+ ztodt*kt)
               trefc = 315. - 60.*sinsq(i)
               trefa = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
               trefa = max(t000,trefa)
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         else
            do i=1,ncol
               tmp     = ka/(1.+ ztodt*ka)
               pressmb = pmid(i,k)*0.01
               trefc   = 315. - 60.*sinsq(i)
               trefa   = (trefc - 10.*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)** &
                                                                                     cappa
               trefa   = max(t000,trefa)
               if (pressmb <= 100. .and. pressmb > 1.) then
                  trefa = t000 + lapse*h0*coslat(i)*log(100./pressmb)
                  tmp   = (efold_strat-efold_meso)*log(pressmb)/log(100.)
                  tmp   = efold_meso + tmp
                  tmp   = 1./(86400.*tmp)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               if (pressmb <= 1. .and. pressmb > 0.01) then
                  trefa = t000 + lapse*h0*coslat(i)*log(100.*pressmb)
                  tmp   = 1./(86400.*efold_meso)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               if (pressmb <= 0.01) then
                  tmp   = 1./(86400.*efold_meso)
                  tmp   = tmp/(1.+ ztodt*tmp)
               endif
               tend%dtdt (i,k) = (trefa - state%t(i,k))*tmp
            end do
         endif
      end do
!
! Add diffusion near the surface for the wind fields
!
      do k=1,pver
         do i=1,pcols
            ptend%u(i,k) = 0.
            ptend%v(i,k) = 0.
         end do
      end do
      do i=1,pcols
         taux(i) = 0.
         tauy(i) = 0.
      end do
!
      kf = 1./(86400.*efoldf)
!
      do k=1,pver
         if (etamid(k) > sigmab) then
            kv  = kf*(etamid(k) - sigmab)/onemsig
            tmp = -kv/(1.+ ztodt*kv)
            do i=1,ncol
               ptend%u(i,k) = tmp*state%u(i,k)
               ptend%v(i,k) = tmp*state%v(i,k)
               tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
               tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
            end do
         else
            do i=1,ncol
               pressmb  = pmid(i,k)*0.01
               if (pressmb <= 100.) then
                  kv       = 1./(86400.*efoldv)
                  tmp      = 1. + tanh(1.5*log10(p_infint/pressmb))
                  kv       = kv*tmp
                  tmp      = -kv/(1.+ ztodt*kv)
                  ptend%u(i,k) = tmp*state%u(i,k)
                  ptend%v(i,k) = tmp*state%v(i,k)
                  tend%dudt(i,k)  = tend%dudt(i,k) + ptend%u(i,k)
                  tend%dvdt(i,k)  = tend%dvdt(i,k) + ptend%v(i,k)
               endif
            end do
         endif
      end do

   else
      write(6,*) 'TPHYSIDL: flag for choosing desired type of idealized ', &
                 'physics ("idlflag") is set incorrectly.'
      write(6,*) 'The valid options are 1, 2, or 3.'
      write(6,*) 'idlflag is currently set to: ',idlflag
      call endrun
   endif
!
! Archive idealized temperature tendency
!
   call outfld('QRS     ',tend%dtdt      ,pcols   ,lchnk      )

   return
end subroutine tphysidl

