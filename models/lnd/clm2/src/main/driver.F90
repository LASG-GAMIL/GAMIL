#include <misc.h>
#include <preproc.h>

subroutine driver (doalb, eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
!
! Purpose:
! clm model driver
!
! Method:
! Calling sequence:
!
! -> histend   Determines if current time step is the end of history interval
!
! -> calendr   Generate the calendar day (1.00 -> 365.99), month (1 -> 12),
!              and day (1 -> 31) used to calculate the surface albedos and
!              leaf and stem areas for the next time step
!
! -> loop over patch points calling for each patch point:
!    -> Hydrology1          canopy interception and precip on ground
!    -> Biogeophysics1      leaf temperature and surface fluxes
!    -> Biogeophysics_Lake  lake temperature and surface fluxes
!    -> Biogeophysics2      soil/snow and ground temp and update surface fluxes
!    -> Hydrology2          surface and soil hydrology
!    -> Hydrology_Lake      lake hydrology
!    -> EcosystemDyn:       ecosystem dynamics: phenology, vegetation, soil carbon
!    -> SurfaceAlbedo:      albedos for next time step
!      -> SnowAlbedo:       snow albedos: direct beam
!      -> SnowAlbedo:       snow albedos: diffuse
!      -> SoilAlbedo:       soil/lake albedos
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (vis dir)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (vis dif)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (nir dir)
!      -> TwoStream:        absorbed, reflected, transmitted solar fluxes (nir dif)
!    -> BalanceCheck        check for errors in energy and water balances
!    -> histUpdate:         accumulate history fields over history time interval
!
!  -> Rtmriverflux          calls RTM river routing model
!
!  -> histHandler           write history and restart files if appropriate
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: driver.F90,v 1.11.2.6.6.1 2002/10/03 20:07:36 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varder
  use clm_varpar    , only : maxpatch
  use clm_varmap    , only : begpatch, endpatch, numpatch, numland, landvec
  use clm_varctl    , only : fsurdat, wrtdia, csm_doflxave
  use histHandlerMod, only : histHandler, histend, do_restwrite
  use restFileMod   , only : restwrt
  use inicFileMod   , only : inicwrt, do_inicwrite
  use mvegFileMod   , only : interpmonthlyveg
  use time_manager  , only : get_step_size, get_curr_calday, get_curr_date, get_nstep
#if (defined RTM)
  use RtmMod        , only : Rtmriverflux
#endif
#if (defined SPMD)
  use spmdMod       , only : masterproc, npes, compute_mpigs_patch
  use mpishorthand  , only : mpir8, mpilog, mpicom
#else
  use spmdMod       , only : masterproc
#endif
#ifdef COUP_CSM
  use clm_csmMod    , only : csm_dosndrcv, csm_recv, csm_send, csm_flxave, &
                             dorecv, dosend, csmstop_now
#endif
  use shr_sys_mod   , only : shr_sys_flush
  implicit none

! ------------------- arguments -----------------------------------
  logical , intent(in) :: doalb   !true if time for surface albedo calculation
  real(r8), intent(in) :: eccen   !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr  !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0  !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp  !Earth's moving vernal equinox long. of perihelion + pi (radians)
! -----------------------------------------------------------------

! ---------------------- local variables --------------------------
  integer  :: i,j,k,l,m           !loop/array indices
  integer  :: yrp1                !year (0, ...) for nstep+1
  integer  :: monp1               !month (1, ..., 12) for nstep+1
  integer  :: dayp1               !day of month (1, ..., 31) for nstep+1
  integer  :: secp1               !seconds into current date for nstep+1
  real(r8) :: caldayp1            !calendar day for nstep+1
  integer  :: dtime               !timestep size [seconds]
  integer  :: nstep               !timestep index
  real(r8) :: buf1d(begpatch:endpatch) !temporary buffer
  real(r8) :: tsxyav              !average ts for diagnostic output
#if (defined SPMD)
  real(r8) :: gather1d(numpatch)  !temporary
  integer  :: numrecvv(0:npes-1)  !vector of items to be received
  integer  :: displsv(0:npes-1)   !displacement vector
  integer  :: numsend             !number of items to be sent
  integer  :: ier                 !error code
#endif
! -----------------------------------------------------------------

  call t_startf('clm_driver')

! determine time step

  nstep = get_nstep()

! ----------------------------------------------------------------------
! Coupler receive
! ----------------------------------------------------------------------

#ifdef COUP_CSM
!
! Determine if information should be sent/received to/from flux coupler
!
  call csm_dosndrcv(doalb)
!
! Get atmospheric state and fluxes from flux coupler

!
  if (dorecv) then
     call csm_recv()
     if (csmstop_now) then
        call t_stopf('clm_driver')
        RETURN
     endif
  endif
#endif

! ----------------------------------------------------------------------
! Determine if end of history interval
! ----------------------------------------------------------------------

  call histend ()

! ----------------------------------------------------------------------
! Calendar information for next time step
! o caldayp1 = calendar day (1.00 -> 365.99) for cosine solar zenith angle
!              calday is based on Greenwich time
! o monp1    = month (1 -> 12) for leaf area index and stem area index
! o dayp1    = day (1 -> 31)   for leaf area index and stem area index
! ----------------------------------------------------------------------

  dtime = get_step_size()
  caldayp1 = get_curr_calday(offset=dtime)
  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=dtime)

! ----------------------------------------------------------------------
! Determine weights for time interpolation of monthly vegetation data.
! This also determines whether it is time to read new monthly vegetation and
! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
! weights obtained here are used in subroutine ecosystemdyn to obtain time
! interpolated values.
! ----------------------------------------------------------------------

  if (doalb) call interpMonthlyVeg (fsurdat, monp1, dayp1)

! ----------------------------------------------------------------------
! LOOP 1
! ----------------------------------------------------------------------

  call t_startf('clm_loop1')
!$OMP PARALLEL DO PRIVATE (K,J)
  do k = begpatch, endpatch    ! begin 1st loop over patches
!
! Initial set of variables
!
     clm(k)%nstep = nstep
     clm(k)%h2osno_old = clm(k)%h2osno  ! snow mass at previous time step
     clm(k)%frac_veg_nosno = clm(k)%frac_veg_nosno_alb
!
! Determine if will cap snow
!
     if (clm(k)%h2osno > 1000.) then
        clm(k)%do_capsnow = .true.
     else
        clm(k)%do_capsnow = .false.
     endif
!
! Energy for non-lake points
!
     if (.not. clm(k)%lakpoi) then
!
! Initial set of previous time step variables
!
        do j = clm(k)%snl+1, 0       ! ice fraction of snow at previous time step
           clm(k)%frac_iceold(j) = clm(k)%h2osoi_ice(j)/(clm(k)%h2osoi_liq(j)+clm(k)%h2osoi_ice(j))
        enddo
!
! Determine beginning water balance (water balance at previous time step)
!
        clm(k)%begwb = clm(k)%h2ocan + clm(k)%h2osno
        do j = 1, nlevsoi
           clm(k)%begwb = clm(k)%begwb + clm(k)%h2osoi_ice(j) + clm(k)%h2osoi_liq(j)
        enddo
!
! Determine canopy interception and precipitation onto ground surface.
! Determine the fraction of foliage covered by water and the fraction
! of foliage that is dry and transpiring. Initialize snow layer if the
! snow accumulation exceeds 10 mm.
!
        call Hydrology1(clm(k))
!
! Determine leaf temperature and surface fluxes based on ground
! temperature from previous time step.
!
        call Biogeophysics1(clm(k))

     else if (clm(k)%lakpoi) then
!
! Determine lake temperature and surface fluxes
!
        call Biogeophysics_Lake (clm(k))

     endif

     if (.not. clm(k)%lakpoi) then
!
! Ecosystem dynamics: phenology, vegetation, soil carbon.
! Also updates snow fraction
!
        call EcosystemDyn (clm(k), doalb, .false.)

     else if (clm(k)%lakpoi) then

! Needed for global history output

        clm(k)%fpsn = 0.

     endif
!
! Albedos for next time step
!
     if (doalb) then
        call SurfaceAlbedo (clm(k), caldayp1, eccen, obliqr, lambm0, mvelpp)
     endif
!
! THIS WILL EVENTUALLY MARK THE END OF THE PATCH LOOP AND
! THE BEGINNING OF THE SINGLE COLUMN SOIL LOOP(S)
!
! Determine soil/snow temperatures including ground temperature and
! update surface fluxes for new ground temperature.
!
     if (.not. clm(k)%lakpoi) call Biogeophysics2(clm(k))

  end do
!$OMP END PARALLEL DO
  call t_stopf('clm_loop1')

! ----------------------------------------------------------------------
! Coupler send
! ----------------------------------------------------------------------

#ifdef COUP_CSM
!
! Average fluxes over interval if appropriate
! Surface states sent to the flux coupler states are not time averaged
!
     if (csm_doflxave) call csm_flxave()
!
! Send fields to flux coupler
! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
! albedos[n+1], and ocnrof_vec[n]
!
     if (dosend) call csm_send()

#endif

! ----------------------------------------------------------------------
! LOOP 2
! ----------------------------------------------------------------------

  call t_startf('clm_loop2')
!$OMP PARALLEL DO PRIVATE (K)
  do k = begpatch, endpatch   ! begin 2nd loop over patches
!
! Vertical (column) soil and surface hydrology
!
     if (.not. clm(k)%lakpoi) call Hydrology2 (clm(k))
!
! Lake hydrology
!
     if (clm(k)%lakpoi) call Hydrology_Lake (clm(k))
!
! Update Snow Age (needed for surface albedo calculation - but is
! really a column type property
!
     call SnowAge (clm(k))
!
! Fraction of soil covered by snow - really a column property
!
     clm(k)%frac_sno = clm(k)%snowdp/(10.*clm(k)%zlnd + clm(k)%snowdp)
!
! Check the energy and water balance
!
     call BalanceCheck (clm(k))
  end do
!$OMP END PARALLEL DO
  call t_stopf('clm_loop2')

! ----------------------------------------------------------------------
! Update history fields and internally accumulated fields
! ----------------------------------------------------------------------

  call t_startf('histup')
  call histUpdate ()
  call t_stopf('histup')

! ----------------------------------------------------------------------
! Write global average diagnostics to standard output
! ----------------------------------------------------------------------

  if (wrtdia) then
     buf1d(begpatch:endpatch) = clm(begpatch:endpatch)%t_rad
#if (defined SPMD)
     call compute_mpigs_patch(1, numsend, numrecvv, displsv)
     call mpi_gatherv (buf1d(begpatch), numsend , mpir8, &
          gather1d, numrecvv, displsv, mpir8, 0, mpicom, ier)
#endif
     if (masterproc) then
        tsxyav = 0._r8
        do m = 1, maxpatch
           do l = 1, numland
              k = landvec%patch(l,m)
#if (defined SPMD)
              tsxyav = tsxyav + gather1d(k)*landvec%wtxy(l,m)
#else
              tsxyav = tsxyav + buf1d(k)*landvec%wtxy(l,m)
#endif
           end do
        end do
        tsxyav = tsxyav / numland
        write (6,1000) nstep,tsxyav
1000    format (1x,'nstep = ',i10,'   TS = ',e21.15)
     end if
  else
#if (!defined COUP_CAM)
     if (masterproc) then
        write(6,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(6)
     endif
#endif
  endif

#if (defined RTM)
! ----------------------------------------------------------------------
! Route surface and subsurface runoff into rivers
! ----------------------------------------------------------------------

  call t_startf('clm_rtm')
  call Rtmriverflux ()
  call t_stopf('clm_rtm')
#endif

! ----------------------------------------------------------------------
! Write history, restart files and initial conditions file if appropriate
! ----------------------------------------------------------------------

  call t_startf('clm_output')

  call histhandler ()

  if (do_restwrite()) call restwrt ()

  if (do_inicwrite()) call inicwrt ()

  call t_stopf('clm_output')

  call t_stopf('clm_driver')

  return
end subroutine driver


