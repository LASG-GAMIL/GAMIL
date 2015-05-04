#include <misc.h>
#include <preproc.h>

module clm_csmMod

!-----------------------------------------------------------------------
!
! Purpose:
! Set of routines that define communication between the land model
! and the flux coupler. The order of sends/receives is as follows:
!  - receive orbital data from coupler (csm_recvorb)
!  - send control data (grids and masks) to coupler (csm_sendcontrol)
!    land grid does not have valid data, runoff grid does
!  - receive valid land grid from flux coupler (csm_recvgrid)
!  - send compressed runoff information to flux coupler (csm_sendrunoff)
!  - send first land model data to flux coupler (csm_send_alb)
!  - start normal send/recv communication pattern
!      => csm_dosndrcv
!      => csm_recv
!      => csm_flxave
!      => csm_send
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
! $Id: clm_csmMod.F90,v 1.12.2.7 2002/06/15 13:50:27 erik Exp $
!-----------------------------------------------------------------------

#ifdef COUP_CSM

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE infnan
  USE mpishorthand
  USE clm_varpar         !parameters
  USE spmdMod            !spmd routines and variables
  USE shr_msg_mod        !csm_share message passing routines
  USE shr_sys_mod, only: shr_sys_irtc !csm_share system utility routines

  implicit none

! Buffer information

  integer, parameter :: nibuff = 100    ! cpl ->atm msg, initial
  integer, parameter :: ncbuff_max=2000 ! Max size of character data from cpl

  integer   :: ncbuff             !Size of character data from cpl
  integer   :: ibuffr(nibuff)     !Integer buffer from cpl
  integer   :: ibuffs(nibuff)     !Integer buffer to   cpl
  real(r8)  :: rbuff(nibuff)      !Floating pt buffer from cpl
  character :: cbuff(ncbuff_max)  !Character data recieved from cpl

! Timing information

  logical :: csm_timing
  integer :: irtc_w               !rc ticks when waiting for msg
  integer :: irtc_r               !rtc tics when msg recved
  integer :: irtc_s               !rtc tics when msg sent

! Send/recv buffers

  integer, parameter :: nsnd = 14
  integer, parameter :: nrcv = 16

  real(r8) :: send2d(lsmlon,lsmlat,nsnd)  !2d send buffer
  real(r8) :: recv2d(lsmlon,lsmlat,nrcv)  !2d recv buffer

  real(r8), private, pointer :: send1d(:,:)    !1d send buffer
  real(r8), private, pointer :: gather1d(:,:)
  real(r8), private, pointer :: recv1d(:,:)    !1d recv buffer
  real(r8), private, pointer :: scatter1d(:,:)

  logical :: debug_flag   !received from coupler

! Flux averaging arrays and counters

  integer  :: icnt  !step counter for flux averager
  integer  :: ncnt  !number of steps over which to average output fluxes
  real(r8) :: rncnt !reciprocal of ncnt

  real(r8), allocatable :: taux_ave(:)     !averaged array
  real(r8), allocatable :: tauy_ave(:)     !averaged array
  real(r8), allocatable :: lhflx_ave(:)    !averaged array
  real(r8), allocatable :: shflx_ave(:)    !averaged array
  real(r8), allocatable :: lwup_ave(:)     !averaged array
  real(r8), allocatable :: qflx_ave(:)     !averaged array
  real(r8), allocatable :: swabs_ave(:)    !averaged array

! When to send/receive messages to coupler and when to make restart and stop

  integer :: ncpday         !number of send/recv calls per day
  logical :: dorecv         !receive data from coupler this step
  logical :: dosend         !send data to coupler this step
  logical :: csmstop_next   !received stop at eod signal and will stop on next ts
  logical :: csmstop_now    !received stop now signal from coupler
  logical :: csmrstrt       !restart write signal received from coupler

! Indices for send/recv fields

  integer, parameter :: irecv_hgt    = 1    !zgcmxy       Atm state m
  integer, parameter :: irecv_u      = 2    !forc_uxy     Atm state m/s
  integer, parameter :: irecv_v      = 3    !forc_vxy     Atm state m/s
  integer, parameter :: irecv_th     = 4    !forc_thxy    Atm state K
  integer, parameter :: irecv_q      = 5    !forc_qxy     Atm state kg/kg
  integer, parameter :: irecv_pbot   = 6    !ptcmxy       Atm state Pa
  integer, parameter :: irecv_t      = 7    !forc_txy     Atm state K
  integer, parameter :: irecv_lwrad  = 8    !flwdsxy      Atm flux  W/m^2
  integer, parameter :: irecv_rainc  = 9    !rainxy       Atm flux  mm/s
  integer, parameter :: irecv_rainl  = 10   !rainxy       Atm flux  mm/s
  integer, parameter :: irecv_snowc  = 11   !snowfxy      Atm flux  mm/s
  integer, parameter :: irecv_snowl  = 12   !snowfxl      Atm flux  mm/s
  integer, parameter :: irecv_soll   = 13   !forc_sollxy  Atm flux  W/m^2
  integer, parameter :: irecv_sols   = 14   !forc_solsxy  Atm flux  W/m^2
  integer, parameter :: irecv_solld  = 15   !forc_solldxy Atm flux  W/m^2
  integer, parameter :: irecv_solsd  = 16   !forc_solsdxy Atm flux  W/m^2

  integer, parameter :: isend_trad   = 1
  integer, parameter :: isend_asdir  = 2
  integer, parameter :: isend_aldir  = 3
  integer, parameter :: isend_asdif  = 4
  integer, parameter :: isend_aldif  = 5
  integer, parameter :: isend_sno    = 6
  integer, parameter :: isend_taux   = 7
  integer, parameter :: isend_tauy   = 8
  integer, parameter :: isend_lhflx  = 9
  integer, parameter :: isend_shflx  = 10
  integer, parameter :: isend_lwup   = 11
  integer, parameter :: isend_qflx   = 12
  integer, parameter :: isend_tref2m = 13
  integer, parameter :: isend_swabs  = 14

! csm timers

  logical  :: timer_lnd_sendrecv = .false. !true => timer is on
  logical  :: timer_lnd_recvsend = .false. !true => timer is on

  SAVE

!===============================================================================
CONTAINS
!===============================================================================

  SUBROUTINE csm_recvorb (eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
!
! Purpose:
! receive the initial integer, real and character control data from
! the flux coupler.  Then parse it out into the variables used by the
! land model.
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

! ---------------------- arguments------- -------------------------
    real(r8), intent(out) :: eccen  !Earth's eccentricity of orbit
    real(r8), intent(out) :: obliqr !Earth's obliquity in radians
    real(r8), intent(out) :: lambm0 !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(out) :: mvelpp !Earth's moving vernal equinox long of perihelion plus pi (radians)
! -----------------------------------------------------------------

! ---------------------- Local variables --------------------------
    integer  :: i,j,k                !indices
    integer  :: ierr                 !error code
    integer  :: info_time            !T => turn on msg-passing timing
    integer  :: maj_vers             !Major version of message passed from the cpl
    integer  :: min_vers             !Minor version of message passed from the cpl
    real(r8) :: spval                !float-pt buffer special value from coupler
! -----------------------------------------------------------------

    if (masterproc) then

! Receive control message from coupler

       ibuffr(:) = 0

       call shr_msg_recv_i (ibuffr, size(ibuffr), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)

       ierr       = ibuffr( 1)  !error code
       info_time  = ibuffr(11)  !T => turn on msg-passing timing
       maj_vers   = ibuffr(40)  !Coupler message major version
       min_vers   = ibuffr(41)  !Coupler message minor version
       ncbuff     = ibuffr(42)  !Size of character data to receive
       write(6,*) '(CSM_RECVORB): recd d->l initial ibuffr msg_id = ',SHR_MSG_TAG_C2LI

! Determine debug flag

       if (ibuffr(12) >= 2) then
          debug_flag = .true.
       else
          debug_flag = .false.
       endif

! Check that the version of the message from the coupler is valid

       call csm_compat(maj_vers, min_vers,SHR_MSG_L_MAJ_V04, SHR_MSG_L_MIN_V00)

! Receive orbital parameters from coupler

       rbuff(:) = 0.0

       call shr_msg_recv_r (rbuff, size(rbuff), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)

       spval  = rbuff(1)      !Coupler float-pt special data flag
       eccen  = rbuff(2)      !Earth's eccentricity of orbit
       obliqr = rbuff(3)      !Earth's Obliquity radians
       lambm0 = rbuff(4)      !Earth's Long. of prehelian at v-equinox
       mvelpp = rbuff(5)      !Earth's Moving vernal equinox of orbit + pi
       write(6,*)'(CSM_RECVORB): recd d->l initial real buff msg_id = ',SHR_MSG_TAG_C2LI

! Check that data received is good data and not the special value

       call compat_check_spval(spval, eccen ,'Eccentricity'     )
       call compat_check_spval(spval, obliqr,'Obliquity'        )
       call compat_check_spval(spval, lambm0,'Long of perhelion')
       call compat_check_spval(spval, mvelpp,'Move long of perh')

       write(6,*)'(CSM_RECVORB): eccen:  ', eccen
       write(6,*)'(CSM_RECVORB): obliqr: ', obliqr
       write(6,*)'(CSM_RECVORB): lambm0: ', lambm0
       write(6,*)'(CSM_RECVORB): mvelpp: ', mvelpp

! Receive character data from coupler and determine if will output csm timing info

       if (ncbuff > 0) then
          call shr_msg_recv_c (cbuff, ncbuff, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
          write(6,*)'(CSM_RECVORB): recd d->a initial char. buf msg_id = ',SHR_MSG_TAG_C2LI
          write(6,*)'(CSM_RECVORB): Char: ',(cbuff(i), i = 1,ncbuff)
       end if

       if (info_time == 0) then
          csm_timing = .false.
       else
          csm_timing = .true.
       endif

    end if
#if ( defined SPMD )
   call mpi_bcast(spval , 1, mpir8, 0, mpicom, ierr)
   call mpi_bcast(eccen , 1, mpir8, 0, mpicom, ierr)
   call mpi_bcast(obliqr, 1, mpir8, 0, mpicom, ierr)
   call mpi_bcast(lambm0, 1, mpir8, 0, mpicom, ierr)
   call mpi_bcast(mvelpp, 1, mpir8, 0, mpicom, ierr)
#endif

 END SUBROUTINE csm_recvorb

!===============================================================================

 SUBROUTINE csm_sendcontrol(irad)

!-----------------------------------------------------------------------
!
! Purpose:
! send first control data to flux coupler and "invalid" grid
! containing special value data
!
! Method:
! The coupler treats points where the mask is nonzero as points where
! you could possibly do a calculation (in the case of the runoff, this
! corresponds to all RTM ocean points). The coupler then defines a "key"
! as points where the model can give you valid data (in the case of runoff,
! this corresponds to points where the land model will give you valid
! compressed data points). The key can be 0 where the mask is 1. However,
! the key cannot be 1 where the mask is 0 unless the data is also zero.
! In the case of runoff, the key the coupler builds is time invariant.
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varctl   , only : csm_doflxave, nsrest
    use RtmMod       , only : area_r, longxy_r, latixy_r, mask_r
    use clm_varcon   , only : re
    use time_manager , only : get_step_size
    use controlMod   , only : csm_dtime
    use shr_const_mod, only : SHR_CONST_CDAY

! ---------------------- arguments--------------------------------------
    integer , intent(in) :: irad  !frequency of radiation computation
!-----------------------------------------------------------------------

! ---------------------- Local variables --------------------------
    real(r8) rtemp_lnd(lsmlon*lsmlat*4) !temporary vector
    integer  itemp_lnd(lsmlon*lsmlat)   !temporary vector
    real(r8) rtemp_rtm(rtmlon*rtmlat*4) !temporary vector
    real(r8) temp_area(rtmlon,rtmlat)   !temporary area
    integer  temp_mask(rtmlon,rtmlat)   !temporary mask
    real(r8) dtime                      !step size
! -----------------------------------------------------------------

    if (masterproc) then

       rtemp_lnd(:) = 1.e30
       itemp_lnd(:) = 999
       rtemp_rtm(:) = 1.e30

! Determine number of send/recv calls steps per day to flux coupler

       if (nsrest == 0) then
          dtime = get_step_size()
       else
          dtime = csm_dtime
       endif
       if (csm_doflxave) then
          ncpday = nint(SHR_CONST_CDAY/dtime)/irad
       else
          ncpday = nint(SHR_CONST_CDAY/dtime)
       endif

! Send integer control information

       ibuffs(:)  = 0                !initialize ibuffs
       ibuffs(7)  = lsmlon           !number of land longitudes
       ibuffs(8)  = lsmlat           !number of land latitudes
       ibuffs(9)  = ncpday           !number of land send/recv calls per day
       ibuffs(34) = 1                !T(or 1) => requests cpl to send valid land domain info back
       ibuffs(36) = 0                !size of compressed runoff vector, if zero then not sending compressed info
       ibuffs(37) = rtmlon           !number of longitudes in uncompressed 2d runoff array
       ibuffs(38) = rtmlat           !number of latitudes  in uncompressed 2d runoff array

       call shr_msg_send_i (ibuffs ,nibuff, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

! Send "invalid" land model grid and mask data

       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat*4), lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat*4), lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i (itemp_lnd(1:lsmlon*lsmlat  ), lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

! Send "valid" RTM grid and mask data

       temp_area(:,:) = area_r(:,:)/(re*re)  !convert from km^2 to radians^2 before sending to coupler
       temp_mask(:,:) = 1 - mask_r(:,:)      !make coupler runoff mask 1 over ocean and 0 over land

       call shr_msg_send_r (longxy_r , rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (latixy_r , rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_rtm, rtmlon*rtmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (rtemp_rtm, rtmlon*rtmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_r (temp_area, rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i (temp_mask, rtmlon*rtmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       write(6,*)'(CSM_SENDCONTROL): there will be ',ncpday, &
            ' send/recv calls per day from the land model to the flux coupler'
       write(6,*)'(CSM_SENDCONTROL):sent l->d control data msg_id = ',SHR_MSG_TAG_L2CI

    endif

  END SUBROUTINE csm_sendcontrol

!===============================================================================

  SUBROUTINE csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

!-----------------------------------------------------------------------
!
! Purpose:
! Receive valid land grid and land mask from coupler
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

! ---------------------- arguments--------------------------------------
    integer , intent(out) :: cam_numlon(lsmlat)           !cam number of longitudes
    real(r8), intent(out) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8), intent(out) :: cam_latixy(lsmlon,lsmlat)    !cam lat values
    real(r8), intent(out) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer , intent(out) :: cam_landmask(lsmlon,lsmlat)  !cam land mask
!-----------------------------------------------------------------------

! ---------------------- Local variables -------------------------------
    integer  i,j                      !loop indices
    real(r8) xe(4,lsmlon,lsmlat)      !coupler land grid edges
    real(r8) ye(4,lsmlon,lsmlat)      !coupler land grid edges
    real(r8) area_a(lsmlon,lsmlat)    !coupler atm grid areas
    integer  mask_a(lsmlon,lsmlat)    !coupler atm valid grid mask
!-----------------------------------------------------------------------

    if (masterproc) then

       ibuffr(:) = 0

       call shr_msg_recv_i (ibuffr       ,nibuff         , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_longxy   ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_latixy   ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (xe           ,lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (ye           ,lsmlon*lsmlat*4, SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (area_a       ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_r (cam_landfrac ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_i (cam_landmask ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)
       call shr_msg_recv_i (mask_a       ,lsmlon*lsmlat  , SHR_MSG_TID_CPL, SHR_MSG_TAG_C2LI)

       write(6,*)'(CSM_SENDGRID):recd d->l land grid, msg_id= ',SHR_MSG_TAG_C2L

! USE mask_a to determine number of valid longitudes for each latitude band
! this is the only use for mask_a

       cam_numlon(:) = 0
       do j = 1,lsmlat
          do i= 1,lsmlon
             if (mask_a(i,j) /= 0) cam_numlon(j) = cam_numlon(j)+1
          end do
       end do

    endif  !end of if-masterproc block

    return
  END SUBROUTINE csm_recvgrid

!===============================================================================

  SUBROUTINE csm_sendrunoff()

!-----------------------------------------------------------------------
!
! Purpose:
! Send valid runoff information back to flux coupler
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use RtmMod, only : ocnrof_iindx, ocnrof_jindx, ocnrof_vec

    if (masterproc) then

! Send integer buffer control info

       ibuffs(7)  = lsmlon           !number of model longitudes
       ibuffs(8)  = lsmlat           !number of model latitudes
       ibuffs(36) = size(ocnrof_vec) !number of data points in compressed runoff data
       ibuffs(37) = rtmlon           !number of longitudes in uncompressed 2d runoff array
       ibuffs(38) = rtmlat           !number of latitudes  in uncompressed 2d runoff array

       call shr_msg_send_i(ibuffs       ,nibuff          , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

! Send runoff vector compression info

       call shr_msg_send_i(ocnrof_iindx, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)
       call shr_msg_send_i(ocnrof_jindx, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2CI)

       write(6,*) '(CSM_SENDROF):sent l->d valid initial runoff info msg_id = ',SHR_MSG_TAG_L2CI

    endif

  END SUBROUTINE csm_sendrunoff

!===============================================================================

  SUBROUTINE csm_sendalb

!-----------------------------------------------------------------------
!
! Purpose:
! Send initial albedos, surface temperature and snow data to the
! flux coupler
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varder
    use clm_varsur
    use clm_varmap
    use clm_varctl  , only : csm_doflxave, nsrest
    use RtmMod      , only : ocnrof_vec
    use time_manager, only : get_curr_date, get_prev_date

! --------------------------- Local variables ---------------------
    integer  :: i,j,k,l,m,n !loop indices
    real(r8) :: wt          !weight
    integer  :: yr          !current year
    integer  :: mon         !current month
    integer  :: day         !current day (0, 1, ...)
    integer  :: ncsec       !current seconds of current date (0, ..., 86400)
    integer  :: ncdate      !current date (yymmdd format) (e.g., 021105)
#if (defined SPMD)
    integer  :: numrecvv(0:npes-1) !vector of items to be received
    integer  :: displsv(0:npes-1)  !displacement vector
    integer  :: numsend            !number of items to be sent
#endif
    integer  :: ier                !return error code
! -----------------------------------------------------------------

! Allocate dynamic memory

    allocate (send1d(nsnd,begpatch:endpatch), STAT=ier)
    if (ier /= 0) then
       write(6,*)'CSM_SENDALB error: send1d allocation error'
       call endrun
    endif
    send1d(:,:)=inf

#if (defined SPMD)
    if (masterproc) then
       allocate (gather1d(nsnd,numpatch), STAT=ier)
       if (ier /= 0) then
          write(6,*)'CSM_SENDALB error: gather1d allocation error'
          call endrun
       endif
       gather1d(:,:) = inf
    endif
#else
    gather1d => send1d
#endif

    if (masterproc) then
       allocate (recv1d(nrcv,numpatch), STAT=ier)
       if (ier /= 0) then
          write(6,*)'CSM_SENDALB error: recv1d allocation error'
          call endrun
       endif
       recv1d(:,:)=inf
    endif

#if (defined SPMD)
    allocate (scatter1d(nrcv,begpatch:endpatch), STAT=ier)
    if (ier /= 0) then
       write(6,*)'CSM_SENDALB error: scatter1d allocation error'
       call endrun
    endif
    scatter1d(:,:) = inf
#else
    scatter1d => recv1d
#endif

! Send first data to coupler

    if (nsrest == 0) then   !initial run

! On initial timestep ONLY: determine 1d vector of states that will be sent
! to coupler and map fields from 1d subgrid vector to 2d [lsmlon]x[lsmlat] grid.

       do k = begpatch,endpatch
          send1d(:,k)            = 1.e30                !don't want to send NaN
          send1d(isend_trad ,k)  = clm(k)%t_grnd        !tsxy
          send1d(isend_asdir,k)  = clm(k)%albd(1)       !asdir
          send1d(isend_aldir,k)  = clm(k)%albd(2)       !aldir
          send1d(isend_asdif,k)  = clm(k)%albi(1)       !asdif
          send1d(isend_aldif,k)  = clm(k)%albi(2)       !aldif
          send1d(isend_sno  ,k)  = clm(k)%h2osno/1000.  !snow (convert from mm to m)
       end do

#if (defined SPMD)
       call compute_mpigs_patch(nsnd, numsend, numrecvv, displsv)
       if (masterproc) then
          call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
               gather1d, numrecvv, displsv, mpir8, 0, mpicom, ier)
       else
          call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
               0._r8, numrecvv, displsv, mpir8, 0, mpicom, ier)
       endif
#else
       gather1d => send1d
#endif

       if (masterproc ) then
          do n=1,nsnd
             where (landmask(:,:) > 0)
                send2d(:,:,n) = 0.
             elsewhere
                send2d(:,:,n) = 1.e30
             end where
          end do
          send2d(:,:,isend_sno) = 0. ! snow initialized to 0 everywhere
          do k = 1,numpatch
             if (patchvec%wtxy(k) /= 0.) then
                i  = patchvec%ixy(k)
                j  = patchvec%jxy(k)
                wt = patchvec%wtxy(k)
                do n = 1,nsnd
                   send2d(i,j,n) = send2d(i,j,n) + gather1d(n,k)*wt
                end do
             end if
          end do
       endif

    else  ! restart run

! On a restart run, no meaningful data is sent to the flux coupler -
! this includes the ocean runoff vector (which should only contain zero values)
! since the runoff code (riverfluxrtm) has not been called yet

       if (masterproc) send2d(:,:,:) = 1.e30  ! this will be sent on a restart timestep

    endif

! Send data to coupler
! Determine time index to send to coupler. Note that for a restart run,
! the next time step is nstep+1. But must send current time step to
! flux coupler here.

    if (masterproc) then

       if (nsrest == 0) then
          call get_curr_date (yr, mon, day, ncsec)
       else
          call get_prev_date (yr, mon, day, ncsec)
       endif
       ncdate = yr*10000 + mon*100 + day

       ibuffs(4) = ncdate  !model date (yyyymmdd)
       ibuffs(5) = ncsec   !elapsed seconds in model date

       call shr_msg_send_i (ibuffs    , size(ibuffs)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (send2d    , size(send2d)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (ocnrof_vec, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)

       if (csm_timing) then
          irtc_s = shr_sys_irtc()
          write(6,9099) irtc_s,'l->d sending'
9099      format('[mp timing]  irtc = ',i20,' ',a)
       end if

    endif  ! end of if_masterproc

    return
  END SUBROUTINE csm_sendalb

!===============================================================================

  SUBROUTINE csm_dosndrcv (doalb)

!-----------------------------------------------------------------------
!
! Purpose:
! Determine when to send and receive messages to/from the
! flux coupler on this time-step.
!
! Method:
! Determine if send/receive information to/from flux coupler
! Send msgs (land state and fluxes) to the flux coupler only when
! doalb is true (i.e. on time steps before the atm does a solar
! radiation computation). Receive msgs (atm state) from the
! flux coupler only when dorad is true (i.e. on time steps
! when the atm does a solar radiation computation).
! The fluxes are then averaged between the send and receive calls.
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varctl   , only : csm_doflxave
    use time_manager , only : get_step_size, get_nstep
    use shr_const_mod, only : SHR_CONST_CDAY

!----------------------------- Arguments --------------------------
    logical, intent(in) :: doalb  !true=>next timestep a radiation time step
!-----------------------------------------------------------------

! ---------------------- Local variables --------------------------
    integer  :: ntspday           !model steps per day
    real(r8) :: dtime             !step size (seconds)
    integer  :: nstep             !time step
!-----------------------------------------------------------------

! -----------------------------------------------------------------

! Determine if send/receive information to/from flux coupler

    nstep = get_nstep()
    if (csm_doflxave) then
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = doalb
       else
          dorecv = dosend
          dosend = doalb
       endif
    else
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = .true.
       else
          dorecv = .true.
          dosend = .true.
       endif
    endif

! If at end of day: check if should write restart file or stop at next time step
! Note these statements must appear here since ibuffr is not received at every time
! step when flux averaging occurs.

    csmstop_next = .false.
    csmrstrt     = .false.
    dtime        = get_step_size()
    ntspday      = nint(SHR_CONST_CDAY/dtime)
    if (mod(nstep,ntspday) == 0) then
       if (ibuffr(2) /= 0) then  !stop at end of day
          csmstop_next = .true.  !will stop on next time step
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
       if (ibuffr(21) /= 0) then !write restart at end of day
          csmrstrt = .true.      !will write restart now
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
    endif

    return
  END SUBROUTINE csm_dosndrcv

!===============================================================================

  SUBROUTINE csm_recv()

!-----------------------------------------------------------------------
!
! Purpose:
! Receive data from the flux coupler
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varder           !derived type definition
    use clm_varcon           !physical constants
    use clm_varmap           !mapping variables

! --------------------------- Local variables ---------------------
    integer :: i,j,k,n             !indices
    real(r8):: forc_rainc          !rainxy Atm flux mm/s
    real(r8):: forc_rainl          !rainxy Atm flux mm/s
    real(r8):: forc_snowc          !snowfxy Atm flux  mm/s
    real(r8):: forc_snowl          !snowfxl Atm flux  mm/s
    integer :: ier                 !return error code
#if (defined SPMD)
    integer  :: numsendv(0:npes-1)    !vector of items to be sent
    integer  :: displsv(0:npes-1)     !displacement vector
    integer  :: numrecv               !number of items to be received
#endif
! -----------------------------------------------------------------

! Start timers

     if (timer_lnd_sendrecv) then
        call t_stopf ('lnd_sendrecv') ; timer_lnd_sendrecv = .false.
     endif

     call t_startf('lnd_recv')

! Receive message from flux coupler

    if (masterproc) then
       ibuffr(:)     = 0
       recv2d(:,:,:) = 1.e30
       if (csm_timing) irtc_w = shr_sys_irtc()
       call shr_msg_recv_i (ibuffr, size(ibuffr), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
       call shr_msg_recv_r (recv2d, size(recv2d), SHR_MSG_TID_CPL, SHR_MSG_TAG_C2L)
       if (csm_timing) then
          irtc_r = shr_sys_irtc()
          write(6,9099) irtc_w,'d->l waiting'
9099      format('[mp timing]  irtc = ',i20,' ',a)
          write(6,9099) irtc_r,'d->l received'
       end if

! Do global integrals of fluxes if flagged

       if (debug_flag) then
          write(6,*)
          write(6,100) 'lnd','recv', irecv_lwrad, global_sum(recv2d(1,1,irecv_lwrad),1.e30), ' lwrad'
          write(6,100) 'lnd','recv', irecv_rainc, global_sum(recv2d(1,1,irecv_rainc),1.e30), ' rainc'
          write(6,100) 'lnd','recv', irecv_rainl, global_sum(recv2d(1,1,irecv_rainl),1.e30), ' rainl'
          write(6,100) 'lnd','recv', irecv_snowc, global_sum(recv2d(1,1,irecv_snowc),1.e30), ' snowc'
          write(6,100) 'lnd','recv', irecv_snowl, global_sum(recv2d(1,1,irecv_snowl),1.e30), ' snowl'
          write(6,100) 'lnd','recv', irecv_soll , global_sum(recv2d(1,1,irecv_soll ),1.e30), ' soll '
          write(6,100) 'lnd','recv', irecv_sols , global_sum(recv2d(1,1,irecv_sols ),1.e30), ' sols '
          write(6,100) 'lnd','recv', irecv_solld, global_sum(recv2d(1,1,irecv_solld),1.e30), ' solld'
          write(6,100) 'lnd','recv', irecv_solsd, global_sum(recv2d(1,1,irecv_solsd),1.e30), ' solsd'
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
          write(6,*)
       endif
    endif  ! end of if-masteproc
#if (defined SPMD)
    call mpi_bcast (ibuffr, size(ibuffr), mpiint, 0, mpicom, ier)
#endif

! Stop timer

     call t_stopf('lnd_recv')

! Check if end of run now, if so stop (each processor does this)

    csmstop_now = .false.
    if (ibuffr(3) /= 0) then
       csmstop_now = .true.
       if (timer_lnd_recvsend) call t_stopf('lnd_recvsend')
       if (timer_lnd_sendrecv) call t_stopf('lnd_sendrecv')
       write(6,*)'(CSM_RECV) stop now signal from flux coupler'
       write(6,*)'(CSM_RECV) ibuffr(3) = ',ibuffr(3)
       if (masterproc) then
          write(6,9001)
          write(6,9002) ibuffr(4)
          write(6,9003)
9001      format(/////' ===========> Terminating CLM Model')
9002      format(     '      Date: ',i8)
9003      format(/////' <=========== CLM Model Terminated')
       endif
       RETURN
    endif

! More timer logic

     if (.not. timer_lnd_recvsend) then
        call t_startf('lnd_recvsend') ; timer_lnd_recvsend = .true.
     endif

! Map 2d received fields on [lsmlon]x[lsmlat] grid to subgrid vectors

    if (masterproc) then
       do k = 1,numpatch
          i = patchvec%ixy(k)
          j = patchvec%jxy(k)
          do n = 1,nrcv
             recv1d(n,k) = recv2d(i,j,n)
          end do
       end do
    end if

#if (defined SPMD)
    call compute_mpigs_patch(nrcv, numrecv, numsendv, displsv)
    if (masterproc) then
       call mpi_scatterv (recv1d, numsendv, displsv, mpir8, &
            scatter1d(1,begpatch), numrecv, mpir8 , 0, mpicom, ier)
    else
       call mpi_scatterv (0._r8, numsendv, displsv, mpir8, &
            scatter1d(1,begpatch), numrecv , mpir8, 0, mpicom, ier)
    endif
#else
    scatter1d => recv1d
#endif

! Split data from coupler into component arrays. Note that the precipitation fluxes received
! from the coupler are in units of kg/s/m^2. To convert these precipitation rates in units of
! mm/sec, one must divide by 1000 kg/m^3 and multiply by 1000 mm/m resulting in an overall
! factor of unity. Below the units are therefore given in mm/s.

    do k = begpatch, endpatch
       clm(k)%forc_hgt      = scatter1d(irecv_hgt  ,k)
       clm(k)%forc_u        = scatter1d(irecv_u    ,k)
       clm(k)%forc_v        = scatter1d(irecv_v    ,k)
       clm(k)%forc_th       = scatter1d(irecv_th   ,k)
       clm(k)%forc_q        = scatter1d(irecv_q    ,k)
       clm(k)%forc_pbot     = scatter1d(irecv_pbot ,k)
       clm(k)%forc_t        = scatter1d(irecv_t    ,k)
       clm(k)%forc_lwrad    = scatter1d(irecv_lwrad,k)
       forc_rainc           = scatter1d(irecv_rainc,k)
       forc_rainl           = scatter1d(irecv_rainl,k)
       forc_snowc           = scatter1d(irecv_snowc,k)
       forc_snowl           = scatter1d(irecv_snowl,k)
       clm(k)%forc_solad(2) = scatter1d(irecv_soll ,k)
       clm(k)%forc_solad(1) = scatter1d(irecv_sols ,k)
       clm(k)%forc_solai(2) = scatter1d(irecv_solld,k)
       clm(k)%forc_solai(1) = scatter1d(irecv_solsd,k)

       ! determine derived quantities

       clm(k)%forc_hgt_u = clm(k)%forc_hgt   !observational height of wind [m]
       clm(k)%forc_hgt_t = clm(k)%forc_hgt   !observational height of temperature [m]
       clm(k)%forc_hgt_q = clm(k)%forc_hgt   !observational height of humidity [m]
       clm(k)%forc_vp    = clm(k)%forc_q*clm(k)%forc_pbot / (0.622+0.378*clm(k)%forc_q)
       clm(k)%forc_rho   = (clm(k)%forc_pbot-0.378*clm(k)%forc_vp) / (rair*clm(k)%forc_t)
       clm(k)%forc_co2   = pco2*clm(k)%forc_pbot
       clm(k)%forc_o2    = po2*clm(k)%forc_pbot

       ! Determine precipitation needed by clm

       clm(k)%forc_rain = forc_rainc + forc_rainl
       clm(k)%forc_snow = forc_snowc + forc_snowl
       if ( clm(k)%forc_snow > 0.0_r8  .and. clm(k)%forc_rain > 0.0_r8 ) then
          write(6,*) 'kpatch= ',k, &
               ' snow= ',clm(k)%forc_snow,' rain= ',clm(k)%forc_rain, &
               ' CLM cannot currently handle both non-zero rain and snow'
          call endrun
       elseif (clm(k)%forc_rain > 0.) then
          clm(k)%itypprc = 1
       elseif (clm(k)%forc_snow > 0.) then
          clm(k)%itypprc = 2
       else
          clm(k)%itypprc = 0
       endif
    end do

    return
  END SUBROUTINE csm_recv

!===============================================================================

  SUBROUTINE csm_send()

!-----------------------------------------------------------------------
!
! Purpose:
! Send data to the flux coupler
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varder
    use clm_varmap            !mapping arrays
    use clm_varctl            !run control variables
    use clm_varsur            !surface variables
    use RtmMod      , only : ocnrof_vec
    use time_manager, only : get_curr_date

! --------------------------- Local variables ---------------------
    integer :: i,j,k,l,m,n !loop indices
    integer :: yr          !current year
    integer :: mon         !current month
    integer :: day         !current day (0, 1, ...)
    integer :: ncsec       !current seconds of current date (0, ..., 86400)
    integer :: ncdate      !current date (yymmdd format) (e.g., 021105)
#if (defined SPMD)
    integer :: numrecvv(0:npes-1)   !vector of items to be received
    integer :: displsv(0:npes-1)    !displacement vector
    integer :: numsend              !number of items to be sent
    integer :: ier                  !error return status
#endif
! -----------------------------------------------------------------

! Send data to the flux coupler

    if (timer_lnd_recvsend) then
       call t_stopf ('lnd_recvsend') ; timer_lnd_recvsend = .false.
    endif

! Start timer

    call t_startf('lnd_send')

! Determine 1d vector of fields that will be sent to coupler.
! Coupler has convention that fluxes are positive downward.

    do k = begpatch,endpatch
       send1d(isend_trad ,k) = clm(k)%t_rad         !tsxy
       send1d(isend_asdir,k) = clm(k)%albd(1)       !asdir
       send1d(isend_aldir,k) = clm(k)%albd(2)       !aldir
       send1d(isend_asdif,k) = clm(k)%albi(1)       !asdif
       send1d(isend_aldif,k) = clm(k)%albi(2)       !aldif
       send1d(isend_sno  ,k) = clm(k)%h2osno/1000.  !snow (convert from mm to m)
       if (csm_doflxave) then
          send1d(isend_taux ,k) = -taux_ave(k)
          send1d(isend_tauy ,k) = -tauy_ave(k)
          send1d(isend_lhflx,k) = -lhflx_ave(k)
          send1d(isend_shflx,k) = -shflx_ave(k)
          send1d(isend_lwup ,k) = -lwup_ave(k)
          send1d(isend_qflx ,k) = -qflx_ave(k)
          send1d(isend_swabs,k) = -swabs_ave(k)
       else
          send1d(isend_taux ,k) = -clm(k)%taux
          send1d(isend_tauy ,k) = -clm(k)%tauy
          send1d(isend_lhflx,k) = -clm(k)%eflx_lh_tot
          send1d(isend_shflx,k) = -clm(k)%eflx_sh_tot
          send1d(isend_lwup ,k) = -clm(k)%eflx_lwrad_out
          send1d(isend_qflx ,k) = -clm(k)%qflx_evap_tot
          send1d(isend_swabs,k) = -clm(k)%fsa
       endif
       send1d(isend_tref2m,k) =  clm(k)%t_ref2m          !tref
    end do

#if (defined SPMD)
       call compute_mpigs_patch(nsnd, numsend, numrecvv, displsv)
       if (masterproc) then
          call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
               gather1d, numrecvv, displsv, mpir8, 0, mpicom, ier)
       else
          call mpi_gatherv (send1d(1,begpatch), numsend , mpir8, &
               0._r8, numrecvv, displsv, mpir8, 0, mpicom, ier)
       endif
#else
       gather1d => send1d
#endif

! Send data to flux coupler
! First, map fields from 1d subgrid vector to 2d [lsmlon]x[lsmlat] grid, weighting
! by subgrid fraction. Use only points with wt > 0 so SPMD code will not use
! uninitialized stack memory values for arrays like taux. NOTE: snow is sent as
! zero over non-land because currently the ocn and sea-ice send no snow cover
! to coupler and so the coupler sends back zero snow over non-land to
! the atm (atm and land grid are currently assumed to be identical)

    if (masterproc) then

       do n = 1,nsnd
          where( landmask(:,:) > 0 )
             send2d(:,:,n) = 0.
          elsewhere
             send2d(:,:,n) = 1.e30
          endwhere
       end do
       send2d(:,:,isend_sno) = 0.     !reset snow to 0 everywhere
       do k = 1, numpatch
          if (patchvec%wtxy(k) /= 0.) then
             i = patchvec%ixy(k)
             j = patchvec%jxy(k)
             do n=1,nsnd
                send2d(i,j,n) = send2d(i,j,n) + gather1d(n,k)*patchvec%wtxy(k)
             end do
          end if
       end do

       call get_curr_date (yr, mon, day, ncsec)
       ncdate = yr*10000 + mon*100 + day

       ibuffs(4)  = ncdate           !model date (yyyymmdd)
       ibuffs(5)  = ncsec            !elapsed seconds in current date

       call shr_msg_send_i (ibuffs    , size(ibuffs)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (send2d    , size(send2d)    , SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)
       call shr_msg_send_r (ocnrof_vec, size(ocnrof_vec), SHR_MSG_TID_CPL, SHR_MSG_TAG_L2C)

       if (csm_timing) then
          irtc_s = shr_sys_irtc()
          write(6,9099) irtc_s,'l->d sending'
9099      format('[mp timing]  irtc = ',i20,' ',a)
       end if

! Do global integrals if flag is set

       if (debug_flag) then
          write(6,*)
          write(6,100) 'lnd','send', isend_taux , global_sum(send2d(1,1,isend_taux ),1.e30), ' taux'
          write(6,100) 'lnd','send', isend_tauy , global_sum(send2d(1,1,isend_tauy ),1.e30), ' tauy'
          write(6,100) 'lnd','send', isend_lhflx, global_sum(send2d(1,1,isend_lhflx),1.e30), ' lhflx'
          write(6,100) 'lnd','send', isend_shflx, global_sum(send2d(1,1,isend_shflx),1.e30), ' shflx'
          write(6,100) 'lnd','send', isend_lwup , global_sum(send2d(1,1,isend_lwup ),1.e30), ' lwup'
          write(6,100) 'lnd','send', isend_qflx , global_sum(send2d(1,1,isend_qflx ),1.e30), ' qflx'
          write(6,100) 'lnd','send', isend_swabs, global_sum(send2d(1,1,isend_swabs),1.e30), ' swabs'
          write(6,*)
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
       endif

    endif  ! end of if_masterproc

! Stop timers

    call t_stopf('lnd_send')

    if (.not. timer_lnd_recvsend) then
       call t_startf('lnd_sendrecv') ; timer_lnd_sendrecv = .true.
    endif

    return
  END SUBROUTINE csm_send

!===============================================================================

  SUBROUTINE csm_flxave()

!-----------------------------------------------------------------------
!
! Purpose:
! Average output fluxes for flux coupler
!
! Method:
! Add land surface model output fluxes to accumulators every time step.
! When icnt==ncnt, compute the average flux over the time interval.
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varder
    use clm_varctl            !run control variables
    use clm_varmap            !mapping arrays
    use RtmMod, only: ocnrof_vec
    use time_manager, only : get_nstep

! ------------------------ local variables -----------------------------
    integer :: i,k,lat,n          !indices
    integer :: nstep              !model time step
! ----------------------------------------------------------------------

! Allocate dynamic memory if necessary

    if (.not. allocated(taux_ave)) then
       allocate (taux_ave(numpatch)) ; taux_ave(:) = inf
    endif
    if (.not. allocated(tauy_ave)) then
       allocate (tauy_ave(numpatch)) ; tauy_ave(:) = inf
    endif
    if (.not. allocated(lhflx_ave)) then
       allocate (lhflx_ave(numpatch)); lhflx_ave(:) = inf
    endif
    if (.not. allocated(shflx_ave)) then
       allocate (shflx_ave(numpatch)); shflx_ave(:) = inf
    endif
    if (.not. allocated(lwup_ave)) then
       allocate (lwup_ave(numpatch)) ; lwup_ave(:) = inf
    endif
    if (.not. allocated(qflx_ave)) then
       allocate (qflx_ave(numpatch)) ; qflx_ave(:) = inf
    endif
    if (.not. allocated(swabs_ave)) then
       allocate (swabs_ave(numpatch)) ; swabs_ave(:) = inf
    endif

! Determine output flux averaging interval

    nstep = get_nstep()
    if (dorecv) then
       icnt = 1
       if ( nstep==0 ) then
          ncnt = irad + 1
       else
          ncnt = irad
       endif
       rncnt = 1./ncnt
    endif

! Initial call of averaging interval, copy data to accumulators

    if (icnt == 1) then
       do k = begpatch,endpatch
          taux_ave(k)  = clm(k)%taux
          tauy_ave(k)  = clm(k)%tauy
          lhflx_ave(k) = clm(k)%eflx_lh_tot
          shflx_ave(k) = clm(k)%eflx_sh_tot
          lwup_ave(k)  = clm(k)%eflx_lwrad_out
          qflx_ave(k)  = clm(k)%qflx_evap_tot
          swabs_ave(k) = clm(k)%fsa
       end do

! Final call of averaging interval, complete averaging

    else if (icnt == ncnt) then
       do k = begpatch,endpatch
          taux_ave (k) = rncnt * (taux_ave(k)  + clm(k)%taux)
          tauy_ave (k) = rncnt * (tauy_ave(k)  + clm(k)%tauy)
          lhflx_ave(k) = rncnt * (lhflx_ave(k) + clm(k)%eflx_lh_tot)
          shflx_ave(k) = rncnt * (shflx_ave(k) + clm(k)%eflx_sh_tot)
          lwup_ave (k) = rncnt * (lwup_ave(k)  + clm(k)%eflx_lwrad_out)
          qflx_ave (k) = rncnt * (qflx_ave(k)  + clm(k)%qflx_evap_tot)
          swabs_ave(k) = rncnt * (swabs_ave(k) + clm(k)%fsa)
       end do

! Intermediate call, add data to accumulators

    else
       do k = begpatch,endpatch
          taux_ave (k) = taux_ave(k)  + clm(k)%taux
          tauy_ave(k)  = tauy_ave(k)  + clm(k)%tauy
          lhflx_ave(k) = lhflx_ave(k) + clm(k)%eflx_lh_tot
          shflx_ave(k) = shflx_ave(k) + clm(k)%eflx_sh_tot
          lwup_ave(k)  = lwup_ave(k)  + clm(k)%eflx_lwrad_out
          qflx_ave(k)  = qflx_ave(k)  + clm(k)%qflx_evap_tot
          swabs_ave(k) = swabs_ave(k) + clm(k)%fsa
       end do
    end if

! Increment counter

    icnt = icnt + 1

    return
  END SUBROUTINE csm_flxave

!===============================================================================

  SUBROUTINE compat_check_spval( spval, data, string )

!-----------------------------------------------------------------------
!
! Purpose:
! Check that the given piece of real data sent from the coupler is valid
! data and not the special data flag set by the coupler. This ensures that
! the expected data is actually being sent by the coupler.
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

!------------------ Arguments ------------------------------------------
    real(r8), intent(in) :: spval
    real(r8), intent(in) :: data
    character(len=*)     :: string
!-----------------------------------------------------------------------

    if ( spval == data )then
       write(6,*)'ERROR:(compat_check_spval) msg incompatibility'
       write(6,*)'ERROR: I expect to recieve the data type: ',string
       write(6,*)'from CPL, but all I got was the special data flag'
       write(6,*)'coupler must not be sending this data, you are'
       write(6,*)'running with an incompatable version of the coupler'
       call endrun
    end if
    return
  END SUBROUTINE compat_check_spval

!===============================================================================

  SUBROUTINE csm_compat(cpl_maj_vers, cpl_min_vers, expect_maj_vers, expect_min_vers)

!-----------------------------------------------------------------------
! Checks that the message recieved from the coupler is compatable
! with the type of message that I expect to recieve.  If the minor
! version numbers differ I print a warning message.  If the major
! numbers differ I abort since that means that the change is
! drastic enough that I can't run with the differences.
! Original Author: Erik Kluzek Dec/97
!-----------------------------------------------------------------------

!----------------------- Arguments -------------------------------------
    integer, intent(in) :: cpl_maj_vers    ! major version from coupler initial ibuffr array
    integer, intent(in) :: cpl_min_vers    ! minor version from coupler initial ibuffr array
    integer, intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
    integer, intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting
!-----------------------------------------------------------------------

    write(6,*)'(cpl_COMPAT): This is revision: $Revision: 1.12.2.7 $'
    write(6,*)'              Tag: $Name: GAMIL1.0 $'
    write(6,*)'              of the message compatability interface:'

    if ( cpl_min_vers /= expect_min_vers )then
       write(6,*) 'WARNING(cpl_compat):: Minor version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_min_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_min_vers
    end if

    if ( cpl_maj_vers /= expect_maj_vers )then
       write(6,*) 'ERROR(cpl_compat):: Major version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_maj_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_maj_vers
       call endrun
    end if

    return
  END SUBROUTINE csm_compat

!===============================================================================

  real(r8) function global_sum(flux, spval)

!-----------------------------------------------------------------------
! Perform global integral
!-----------------------------------------------------------------------

    use clm_varsur, only : area                 !km^2

    real(r8), intent(in) :: flux(lsmlon,lsmlat) !W/m2, Kg/m2-s or N/m2
    real(r8), intent(in) :: spval               !points to not include in global sum

    integer :: i,j                              !indices

    global_sum = 0.
    do j = 1,lsmlat
       do i = 1,lsmlon
          if (flux(i,j) /= spval) then
             global_sum = global_sum + flux(i,j)*area(i,j)*1.e6
          endif
       end do
    end do

    return
  end function global_sum

!===============================================================================

#endif

end module clm_csmMod
