#include <misc.h>
#include <params.h>

subroutine stepon
!-----------------------------------------------------------------------
!
! Purpose:
! Loop over time, calling driving routines for physics, dynamics, transport
!
! Original version:  CCM1
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use history, only: wshist, wrapup
   use pmgrid
   use pspect
   use comslt
   use rgrid
   use prognostics
   use buffer
   use commap
   use restart, only: write_restart
   use physconst, only: cappa, gravit
#if ( defined SPMD )
   use mpishorthand
   use spmd_dyn, only: npes, compute_gsfactors
#endif
#ifdef COUP_CSM
   use ccsm_msg, only: csmstop, ccsmfin
#endif
   use ppgrid   ,      only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   use time_manager, only: advance_timestep, get_step_size, get_nstep, &
                           is_first_step, is_first_restart_step, &
                           is_last_step, is_end_curr_day

   implicit none

#include <comctl.h>
#include <comhyb.h>
#include <comlun.h>
#include <comqfl.h>

!------------------------------Parameters-------------------------------
!
   integer, parameter :: pmap = 20000 ! max dimension of evenly spaced vert. grid used
!                                     ! by SLT code to map the departure pts into true
!                                     ! model levels.
!
!---------------------------Local workspace-----------------------------
!
   integer kdpmpf  (pmap)             ! artificial full vert grid indices
   integer kdpmph  (pmap)             ! artificial half vert grid indices

   type(physics_state), allocatable :: phys_state(:)
   type(physics_tend ), allocatable :: phys_tend(:)

   real(r8) hyad   (plev)             ! del (A)
   real(r8) lam    (plond,platd)      ! longitude coords of extended grid
   real(r8) phi    (platd)            ! latitude  coords of extended grid
   real(r8) dphi   (platd)            ! latitude intervals (radians)
   real(r8) gw     (plat)             ! Gaussian weights
   real(r8) sinlam (plond,platd)      ! sin(lam) model domain only
   real(r8) coslam (plond,platd)      ! cos(lam) model domain only
   real(r8) lbasdy (4,2,platd)        ! latitude derivative weights
   real(r8) lbasdz (4,2,plev)         ! vert (full levels) deriv wghts
   real(r8) lbassd (4,2,plevp)        ! vert (half levels) deriv wghts
   real(r8) lbasiy (4,2,platd)        ! Lagrange cubic interp wghts (lat.)
   real(r8) lbasiz (4,2,plev)         ! Lagrange cubic interp wghts (vert)
   real(r8) lbassi (4,2,plevp)        ! Lagrange cubic interp wghts (vert)
   real(r8) detam  (plev)             ! intervals between vert full levs.
   real(r8) detai  (plevp)            ! intervals between vert half levs.
   real(r8) dlam   (platd)            ! longitudinal grid interval (radians)
   real(r8) cwava  (plat)             ! weight applied to global integrals
   real(r8) etamid (plev)             ! vertical coords at midpoints
   real(r8) etaint (plevp)            ! vertical coords at interfaces
   real(r8), allocatable :: t2(:,:,:) ! temp tendency
   real(r8), allocatable :: fu(:,:,:) ! u wind tendency
   real(r8), allocatable :: fv(:,:,:) ! v wind tendency
   real(r8) flx_net(plond,beglat:endlat)       ! net flux from physics
   real(r8) coslat(plond)
   real(r8) rcoslat(plond)
   real(r8) rpmid(plond,plev)
   real(r8) pdel(plond,plev)
   real(r8) pint(plond,plevp)
   real(r8) pmid(plond,plev)
   real(r8) ztodt                     ! time step
   real(r8) :: wcstart, wcend   ! wallclock timestamp at start, end of timestep
   real(r8) :: usrstart, usrend ! user timestamp at start, end of timestep
   real(r8) :: sysstart, sysend ! sys timestamp at start, end of timestep

   integer kk                         ! |
   integer kk1                        ! |
   integer kk2                        ! | - indices
   integer kstep                      ! |
   integer l                          ! |
   integer lvsum                      ! counter
   integer lsum                       ! counter
   real(r8) vsum                      ! accumulator for SLD binning statistics
   logical vflag                      ! logical flag to indicate that binning has been done
!                                     ! correctly
   integer i,k,lat,j,begj             ! longitude,level,latitude indices
   real(r8) tmp1                      ! temp space
   integer year,month,day,tod         ! date info
   integer ntspdy                     ! Number of timesteps per day

#ifdef SPMD
   integer :: numsend          ! number of items to be sent
   integer :: numrecv(0:npes-1)! number of items to be received
   integer :: displs(0:npes-1) ! displacement array
   integer :: numperlat        ! number of items per latitude band
#endif
!
! Externals
!
   logical, external :: rstwr  ! whether or not to write restart files
!
!-----------------------------------------------------------------------
!
! Define eta coordinates: Used for calculation etadot vertical velocity
! for slt.
!
   call t_startf ('stepon_startup')

   do k=1,plev
      etamid(k) = hyam(k) + hybm(k)
   end do
   do k=1,plevp
      etaint(k) = hyai(k) + hybi(k)
   end do
!
! Initialize matrix gamma (used in sld T computation)
!
   gamma(:,:) = 0.
   do k = 1,plev
      tmp1 = cappa*t0(k)*hypi(plevp)/hypm(k)
      gamma(k,k) = 0.5*tmp1
      do l=1,k-1
         gamma(l,k) = tmp1
      end do
   end do
!
! Set slt common block variables
!
   call grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
               lam     ,phi     ,dphi    ,gw      ,sinlam  , &
               coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
               lbasiz  ,lbassi  ,detam   ,detai   ,kdpmpf  , &
               kdpmph  ,cwava   ,phigs   )

   if (is_first_step()) then
      do lat=beglat,endlat
         j = j1 - 1 + lat
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1./coslat(i)
         end do
!
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plond, plev, ps(1,lat,n3), pint, pmid, pdel)
!
         do k=1,plev
            do i=1,nlon(lat)
               rpmid(i,k) = 1./pmid(i,k)
            end do
         end do
!
! Compute appropriate (1/ps)etadot(dp/deta)
!
         call etadt0 (lat, nlon(lat), &
                      rcoslat ,div(1,1,lat,n3), u3(i1,1,j,n3), v3(i1,1,j,n3), dpsl(1,lat), &
                      dpsm(1,lat), pdel, ps(1,lat,n3), ed1(1,1,lat))
!
! Calculate vertical motion field
!
         call omcalc (rcoslat, div(1,1,lat,n3), u3(i1,1,j,n3), v3(i1,1,j,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid, pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
      end do
   end if
!
! Compute pdel from "A" portion of hybrid vertical grid
!
   do k=1,plev
      hyad(k) = hyai(k+1) - hyai(k)
   end do
   do k=1,plev
      do i=1,plon
         pdela(i,k) = hyad(k)*ps0
      end do
   end do
!
!----------------------------------------------------------
! Initialize departure point bin arrays to 0.
!----------------------------------------------------------
!
   kstep = 0
   do lat = beglat,endlat
      do k = 1,plev
         do kk = 1,plev
            levkntl(kk,k,lat) = 0.
         end do
      end do
   end do

   allocate(phys_state(begchunk:endchunk))
   allocate(phys_tend(begchunk:endchunk))
   allocate(t2(plond,plev,beglat:endlat))
   allocate(fu(plond,plev,beglat:endlat))
   allocate(fv(plond,plev,beglat:endlat))
!
! Beginning of basic time step loop
!
   call t_stopf ('stepon_startup')

   do

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcstart, usrstart, sysstart)
      end if

      ztodt = get_step_size()
!
!----------------------------------------------------------
! PHYSPKG  Call the Physics package
!----------------------------------------------------------
!
      begj = beglatex + numbnd

      call t_startf('d_p_coupling')
      call d_p_coupling (ps(1,beglat,n3), t3(i1,1,begj,n3), u3(i1,1,begj,n3), &
                         v3(i1,1,begj,n3), q3(i1,1,1,begj,n3), &
                         omga, phis, phys_state)
      call t_stopf('d_p_coupling')

      call t_startf('phys_driver')
      if (ideal_phys) then
         call phys_idealized(phys_state, phys_tend, ztodt, etamid)
      else if (adiabatic) then
         call phys_adiabatic(phys_state, phys_tend)
      else
         call physpkg (                     &
            phys_state, gw, ztodt, &
            phys_tend, cld(1,1,begchunk,n3m1), cld(1,1,begchunk,n3), &
            tcwat(1,1,begchunk,n3m1), tcwat(1,1,begchunk,n3), &
            qcwat(1,1,begchunk,n3m1), qcwat(1,1,begchunk,n3), &
            lcwat(1,1,begchunk,n3m1), lcwat(1,1,begchunk,n3) )
      end if
      call t_stopf('phys_driver')

      call t_startf('p_d_coupling')
      call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net,q3(i1,1,1,begj,n3))
      call t_stopf('p_d_coupling')

      if (is_first_restart_step()) then
         call print_memusage
      end if
!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------

      call t_startf ('dynpkg')
      call dynpkg (pmap    ,t2      ,fu      ,fv      ,etamid  , &
                   etaint  ,cwava   ,detam   ,dlam    ,lam     , &
                   phi     ,dphi    ,sinlam  ,coslam  ,lbasdy  , &
                   lbasdz  ,lbasiy  ,lbassi  ,lbasiz  ,detai   , &
                   kdpmpf  ,kdpmph  ,flx_net , ztodt   )
      call t_stopf ('dynpkg')

      if (is_first_restart_step()) then
         call print_memusage
      end if

      call t_startf('stepon_single')

! Set end of run flag.

#ifndef COUP_CSM
      if (is_last_step()) nlend = .true.
#else
      if (csmstop) then
         if ( masterproc ) write(6,*)'atm: Stopping at the end of this day'
         if (is_end_curr_day()) nlend = .true.
      end if
#endif
!
!----------------------------------------------------------
! History and restart logic
!----------------------------------------------------------
!
! Write and/or dispose history tapes if required
!
      call t_startf ('wshist')
      call wshist ()
      call t_stopf ('wshist')
!
! Write restart file
!
      if (rstwr() .and. nrefrq /= 0) then
         call t_startf ('write_restart')
         call write_restart
         call t_stopf ('write_restart')
      end if
!
! Dispose necessary files
!
      call t_startf ('wrapup')
      call wrapup
      call t_stopf ('wrapup')

      if (masterproc .and. print_step_cost) then
         call t_stampf (wcend, usrend, sysend)
         write(6,'(a,3f8.3,a)')'Prv timestep wallclock, usr, sys=', &
                               wcend-wcstart, usrend-usrstart, sysend-sysstart, ' seconds'
      end if
!
! Increment nstep before returning to top of loop
!
      call advance_timestep()
      kstep = kstep + 1
!
! Compute some statistical info
!
      if (nlend) then
#if ( defined SPMD )
         numperlat = plev*plev
         call compute_gsfactors (numperlat, numsend, numrecv, displs)
         call mpigatherv (levkntl(1,1,beglat), numsend, mpir8, &
                          levkntl            , numrecv, displs, mpir8, 0, mpicom)
#endif
         if (masterproc) then
            do k = 1,plev
               do kk = 1,plev
                  levknt(kk,k) = 0.
                  do lat = 1,plat
                     levknt(kk,k) = levknt(kk,k) + levkntl(kk,k,lat)
                  end do
               end do
            end do

            lsum = 0
            do lat = 1,plat
               lsum = lsum + nlon(lat)
            end do

            vflag = .false.
            do k = 1,plev
               vsum = 0.
               do kk = 1,plev
                  vsum = vsum + levknt(kk,k)
               end do
               lvsum = vsum + 0.01
               if( lvsum .ne. lsum*kstep ) vflag = .true.
            end do

            do k = 1,plev
               do kk = 1,plev
                  levknt(kk,k) = levknt(kk,k)/float(lsum*kstep)
                  if (levknt(kk,k) .eq. 0.) levknt(kk,k) = 1.e+36
               end do
            end do

            if(vflag) then
               write(6,*) '********************************************'
               write(6,1000)
               write(6,*) '********************************************'
            else
               write(6,2000)
               k = 1
               write(6,3001) hypm(k)/100.,(levknt(kk,k),kk = 1, 6)
               k = 2
               write(6,3002) hypm(k)/100.,(levknt(kk,k),kk = 1, 7)
               k = 3
               write(6,3003) hypm(k)/100.,(levknt(kk,k),kk = 1, 8)
               k = 4
               write(6,3004) hypm(k)/100.,(levknt(kk,k),kk = 1, 9)
               k = 5
               write(6,3005) hypm(k)/100.,(levknt(kk,k),kk = 1,10)
               do k = 6,plev
                  kk1 = k-5
                  kk2 = k+5
                  if(kk2 .gt. plev) kk2 = plev
                  write(6,3000) hypm(k)/100.,(levknt(kk,k),kk = kk1,kk2)
               end do
            end if
         end if
      end if

      call t_stopf('stepon_single')
!
! Check for end of run
!
      if (nlend) then
         deallocate(phys_state)
         deallocate(phys_tend)
         deallocate(t2)
         deallocate(fu)
         deallocate(fv)
#ifdef COUP_CSM
         call ccsmfin
#endif
         return
      end if

   end do  ! End of timestep loop

1000 format(' ERROR in binning departure points; sums were done'// &
            ' incorrectly. not printing results.')
2000 format(/40x,' BINNING STATISTICS FOR VERTICAL DEPARTURE POINTS'// &
                 '    level          -5        -4        -3        -2     ', &
                 '   -1    0 (arr pt)     1         2         3         4        ', &
                 ' 5'/)
3000 format(' ',f9.4, 4x,11f10.5)
3001 format(' ',f9.4,54x,11f10.5)
3002 format(' ',f9.4,44x,11f10.5)
3003 format(' ',f9.4,34x,11f10.5)
3004 format(' ',f9.4,24x,11f10.5)
3005 format(' ',f9.4,14x,11f10.5)
end subroutine stepon
