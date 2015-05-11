#include <misc.h>

module time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc
   use esmf_timemgmtmod, only: &
      esmf_errhandlersettype, esmf_err_return, esmf_errprint, esmf_success, &
      esmf_time, esmf_timeinit, esmf_timeget, esmf_timegetdays, &
        esmf_timeincrement, esmf_timedecrement, &
      esmf_date, esmf_dateinit, esmf_gregorian, esmf_no_leap, esmf_dateget, &
        esmf_dateincrementsec, esmf_dateincrementday, esmf_datedecrement, &
        esmf_datediff, esmf_dategetfltdayofyear, &
      esmf_timemgr, esmf_timemgrinit, esmf_timemgradvance, esmf_timemgrgetnstep, &
        esmf_timemgrgetstepsize, esmf_timemgrgetstartdate, esmf_timemgrgetbasedate, &
        esmf_timemgrlaststep, esmf_timemgrgetcurrdate, esmf_timemgrgetprevdate, esmf_dateislater, &
        esmf_timemgrrestartwrite, esmf_timemgrrestartread
   use string_utils, only: to_upper
   use dycore, only: dycore_is
#ifdef SPMD
   use mpishorthand, only: mpicom, mpiint, mpilog
#endif

   implicit none
   private
   save

! Public methods

   public ::&
      timemgr_preset,           &! time manager initialization before namelist input
      timemgr_init,             &! time manager initialization
      advance_timestep,         &! increment timestep number
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_ref_date,             &! return components of the reference date
      get_perp_date,            &! return components of the perpetual date, and current time of day
      get_curr_time,            &! return components of elapsed time since reference date
      get_curr_calday,          &! return calendar day at end of current timestep
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart or branch run
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_perpetual,             &! return true if perpetual calendar is in use
      timemgr_write_restart,    &! write info to file needed to restart the time manager
      timemgr_read_restart,     &! read info from file needed to restart the time manager
      timemgr_restart            ! restart the time manager

! Public data for namelist input

   character(len=32), public ::&
      calendar   = 'NO_LEAP'     ! Calendar to use in date calculations.
                                 ! 'NO_LEAP' or 'GREGORIAN'
   integer, parameter :: uninit_int = -999999999
   integer, public ::&
      dtime         = uninit_int,  &! timestep in seconds
      nestep        = uninit_int,  &! final timestep (or day if negative) number
      nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
      start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
      start_tod     = 0,           &! starting time of day for run in seconds
      stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
      stop_tod      = 0,           &! stopping time of day for run in seconds
      ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
      ref_tod       = 0,           &! reference time of day for time coordinate in seconds
      perpetual_ymd = uninit_int    ! date used by perpetual calendar
   logical, public ::&
      perpetual_run = .false.       ! .true. => use a perpetual calendar
!!
   real(r8),public  ::  dtdy = -9999.9     ! timestep of the dycore in seconds !!(wh 2004.04.14)
   integer   ::  tmp
!!
! Public data for communicating with modules that don't have 'get' methods.

   integer, public ::&
      ic_ymd  = uninit_int,     &! date of initial conditions in yearmmdd format
      ic_tod  = 0                ! time of day of initial conditions in seconds
                                 ! ic_ymd, ic_tod set in control/readinitial.F90
   logical, public ::&
      tm_aqua_planet = .false.   ! flag for aqua-planet simulation

! Private module data

   type(esmf_timemgr) :: tm_id          ! time manager ID
   type(esmf_date)    :: tm_perp_date   ! reference date for perpetual calendar day calc

   integer ::&                      ! Data required to restart time manager:
      rst_type      = uninit_int,  &! calendar type
      rst_nstep     = uninit_int,  &! current step number
      rst_step_days = uninit_int,  &! days component of timestep size
      rst_step_sec  = uninit_int,  &! seconds component of timestep size
      rst_start_ymd = uninit_int,  &! start date
      rst_start_tod = uninit_int,  &! start time of day
      rst_stop_ymd  = uninit_int,  &! stop date
      rst_stop_tod  = uninit_int,  &! stop time of day
      rst_ref_ymd   = uninit_int,  &! reference date
      rst_ref_tod   = uninit_int,  &! reference time of day
      rst_curr_ymd  = uninit_int,  &! current date
      rst_curr_tod  = uninit_int,  &! current time of day
      rst_perp_ymd  = uninit_int    ! perpetual date
   logical ::&
      rst_perp_cal  = .false.       ! true when using perpetual calendar

   logical :: tm_first_restart_step = .false.  ! true for first step of a restart or branch run
   logical :: tm_perp_calendar = .false.       ! true when using perpetual calendar

!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_preset()

! Initialize variables before namelist input.

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'timemgr_preset'
!-----------------------------------------------------------------------------------------

   if ( dtime == uninit_int ) then
      if (dycore_is ('EUL')) then
         dtime  = 1200
         dtdy   = 240.0                 !!(wh 2004.04.14)
      else if (dycore_is ('SLD')) then
         dtime  = 3600
      else if (dycore_is ('LR')) then
         dtime  = 3600
      else
         write(6,*)sub,': need valid dycore to set default DTIME'
         call endrun
      end if
   end if

end subroutine timemgr_preset
!=========================================================================================

subroutine timemgr_init()

! Initialize the ESMF time manager.
!
! NOTE - Assumptions:
!      1) The namelist variables have been set before this routine is called.  (set in control/parse_namelist.F90)
!      2) ic_ymd, ic_tod are set before this routine is called.  (set in control/readinitial.F90)

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'timemgr_init'
   character(len=len(calendar)) :: cal
   integer :: rc                 ! return code
   integer :: cal_type           ! calendar type
   type(esmf_time) :: step_size   ! timestep size
   type(esmf_date) :: start_date  ! start date for run
   type(esmf_date) :: stop_date   ! stop date for run
   type(esmf_date) :: ref_date    ! reference date for time coordinate
! Some backwards compatibility stuff:
   type(esmf_time) :: diff
   integer :: ntspday, ndays, nsecs
   logical :: islater
!-----------------------------------------------------------------------------------------

! Initialize error handling.

   call esmf_errhandlersettype(esmf_err_return)

! Initialize calendar type.

   cal = to_upper(calendar)
   if ( trim(cal) == 'NO_LEAP' ) then
      cal_type = esmf_no_leap
   else if ( trim(cal) == 'GREGORIAN' ) then
      cal_type = esmf_gregorian
   else
      write(6,*)sub,': unrecognized calendar specified: ',calendar
      call endrun
   end if

! Initialize timestep size.

   if ( mod(86400,dtime) /=0 ) then
      write(6,*)sub,': timestep must divide evenly into 1 day'
      call endrun
   end if
!!
   tmp = int(dtdy+0.1)                                      !!
   if ( mod(dtime,tmp) /= 0 ) then                          !!
      write(6,*)sub,': dtdy must divide evenly into dtime'  !! (wh 2004.04.14)
      call endrun                                           !!
   end if                                                   !!
!!
     
   step_size = esmf_timeinit(0, dtime, rc)
   call chkrc(rc, sub//': error return from esmf_timeinit: setting step_size')


! Initialize start date.

   if ( start_ymd == uninit_int ) then
      start_ymd = ic_ymd
      start_tod = ic_tod
   end if
   start_date = esmf_dateinit(cal_type, start_ymd, start_tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateinit: setting start_date')

! Initialize reference date for time coordinate.

   if ( ref_ymd /= uninit_int ) then
      ref_date = esmf_dateinit(cal_type, ref_ymd, ref_tod, rc)
   else
      ref_date = esmf_dateinit(start_date, rc)
   end if
   call chkrc(rc, sub//': error return from esmf_dateinit: setting ref_date')

! Initialize stop date.

   if ( stop_ymd /= uninit_int ) then
      stop_date = esmf_dateinit(cal_type, stop_ymd, stop_tod, rc)
   else if ( nestep /= uninit_int ) then
      if ( nestep >= 0 ) then
         stop_date = esmf_dateincrementsec(start_date, dtime*nestep, rc)
      else
         stop_date = esmf_dateincrementday(start_date, -nestep, rc)
      end if
   else if ( nelapse /= uninit_int ) then
      if ( nelapse >= 0 ) then
         stop_date = esmf_dateincrementsec(start_date, dtime*nelapse, rc)
      else
         stop_date = esmf_dateincrementday(start_date, -nelapse, rc)
      end if
   else
      write(6,*)sub,': Must specify one of stop_ymd, nestep, or nelapse'
      call endrun
   end if
   call chkrc(rc, sub//': error return setting stop_date')

! Initialize a time manager.

   tm_id = esmf_timemgrinit(step_size, start_date, stop_date, ref_date, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrinit')

! Initialize date used for perpetual calendar day calculation.

   if ( tm_aqua_planet ) then
      tm_perp_date = esmf_dateinit(cal_type, 321, 0, rc)
      call chkrc(rc, sub//': error return from esmf_dateinit: setting tm_perp_date')
      tm_perp_calendar = .true.
   else if ( perpetual_run ) then
      if ( perpetual_ymd /= uninit_int ) then
         tm_perp_date = esmf_dateinit(cal_type, perpetual_ymd, 0, rc)
         call chkrc(rc, sub//': error return from esmf_dateinit: setting tm_perp_date')
         tm_perp_calendar = .true.
      else
         tm_perp_date = esmf_dateinit(cal_type, ic_ymd, 0, rc)
         call chkrc(rc, sub//': error return from esmf_dateinit: setting tm_perp_date')
         tm_perp_calendar = .true.
      end if
   end if

! Set variables from "comtim.h interface" for backwards compatibility.

! Calculation of ending timestep number (nestep) assumes a constant stepsize.
   ntspday = 86400/dtime
   diff = esmf_timeinit()
   call esmf_datediff(start_date, stop_date, diff, islater, rc)
   call chkrc(rc, sub//': error return from esmf_datediff calculating nestep')
   call esmf_timeget(diff, ndays, nsecs, rc)
   call chkrc(rc, sub//': error return from esmf_timeget calculating nestep')
   nestep = ntspday*ndays + nsecs/dtime
   if ( mod(nsecs,dtime) /= 0 ) nestep = nestep + 1

! Print configuration summary to log file (stdout).

   if (masterproc) then
      call timemgr_print()
   end if

end subroutine timemgr_init
!=========================================================================================

subroutine timemgr_restart()

! Restart the ESMF time manager.
!
! NOTE - Assumptions:
!      1) The namelist variables have been set before this routine is called.
!         (set in control/parse_namelist.F90)
!         The stop date is the only thing that can be changed by the user on a restart.
!      2) Restart data have been read on the master process before this routine is called.
!         (timemgr_read_restart called from control/restart.F90::read_restart)

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'timemgr_restart'
   integer :: rc                 ! return code
   type(esmf_date) :: start_date  ! start date for run
   type(esmf_date) :: stop_date   ! stop date for run
   type(esmf_date) :: curr_date   ! date of data in restart file
   logical :: islater
   integer :: ymd, tod
! Some backwards compatibility stuff:
   type(esmf_time) :: diff
   integer :: ntspday, ndays, nsecs
!-----------------------------------------------------------------------------------------

#if ( defined SPMD ) 
   call mpibcast(rst_type,      1, mpiint, 0, mpicom)
   call mpibcast(rst_nstep,     1, mpiint, 0, mpicom)
   call mpibcast(rst_step_days, 1, mpiint, 0, mpicom)
   call mpibcast(rst_step_sec,  1, mpiint, 0, mpicom)
   call mpibcast(rst_start_ymd, 1, mpiint, 0, mpicom)
   call mpibcast(rst_start_tod, 1, mpiint, 0, mpicom)
   call mpibcast(rst_stop_ymd,  1, mpiint, 0, mpicom)
   call mpibcast(rst_stop_tod,  1, mpiint, 0, mpicom)
   call mpibcast(rst_ref_ymd,   1, mpiint, 0, mpicom)
   call mpibcast(rst_ref_tod,   1, mpiint, 0, mpicom)
   call mpibcast(rst_curr_ymd,  1, mpiint, 0, mpicom)
   call mpibcast(rst_curr_tod,  1, mpiint, 0, mpicom)
   call mpibcast(rst_perp_ymd,  1, mpiint, 0, mpicom)
   call mpibcast(rst_perp_cal,  1, mpilog, 0, mpicom)
#endif

! Initialize error handling.

   call esmf_errhandlersettype(esmf_err_return)

! Initialize calendar type.

   if ( rst_type == esmf_no_leap ) then
      calendar = 'NO_LEAP'
   else if ( rst_type == esmf_gregorian ) then
      calendar = 'GREGORIAN'
   else
      write(6,*)sub,': unrecognized calendar type in restart file: ',rst_type
      call endrun
   end if

! Initialize the timestep.

   dtime = rst_step_days*86400 + rst_step_sec

! Initialize start date.

   start_date = esmf_dateinit(rst_type, rst_start_ymd, rst_start_tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateinit: setting start_date')

! Initialize current date.

   curr_date = esmf_dateinit(rst_type, rst_curr_ymd, rst_curr_tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateinit: setting curr_date')

! Initialize stop date.

   if ( stop_ymd /= uninit_int ) then
      stop_date = esmf_dateinit(rst_type, stop_ymd, stop_tod, rc)
   else if ( nestep /= uninit_int ) then
      if ( nestep >= 0 ) then
         stop_date = esmf_dateincrementsec(start_date, dtime*nestep, rc)
      else
         stop_date = esmf_dateincrementday(start_date, -nestep, rc)
      end if
   else if ( nelapse /= uninit_int ) then
      if ( nelapse >= 0 ) then
         stop_date = esmf_dateincrementsec(curr_date, dtime*nelapse, rc)
      else
         stop_date = esmf_dateincrementday(curr_date, -nelapse, rc)
      end if
   else
      stop_date = esmf_dateinit(rst_type, rst_stop_ymd, rst_stop_tod, rc)
   end if
   call chkrc(rc, sub//': error return setting stop_date')

! Check that stop date is later than current date.

   call esmf_dateislater(curr_date, stop_date, islater, rc)
   call chkrc(rc, sub//': error return from esmf_dateislater: comparing start and stop dates')
   if ( .not. islater ) then
      write(6,*)sub,': stop date must be specified later than current date: '
      call esmf_dateget(curr_date, ymd, tod)
      write(6,*)' Current date (ymd tod): ', ymd, tod
      call esmf_dateget(stop_date, ymd, tod)
      write(6,*)' Stop date (ymd tod): ', ymd, tod
      call endrun
   end if
   call esmf_dateget(stop_date, rst_stop_ymd, rst_stop_tod)

! Restart a time manager.

   tm_id = esmf_timemgrrestartread(rst_type, rst_nstep, rst_step_days, rst_step_sec, rst_start_ymd, &
                                  rst_start_tod, rst_stop_ymd, rst_stop_tod, rst_ref_ymd, rst_ref_tod, &
                                  rst_curr_ymd, rst_curr_tod, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrrestartread')

! Advance the timestep.  Data from the restart file corresponds to the
! last timestep of the previous run.

   call advance_timestep()

!  Set flag that this is the first timestep of the restart run.

   tm_first_restart_step = .true.

! Initialize date used for perpetual calendar day calculation.

   if ( rst_perp_cal ) then
      tm_perp_date = esmf_dateinit(rst_type, rst_perp_ymd, 0, rc)
      call chkrc(rc, sub//': error return from esmf_dateinit: setting tm_perp_date')
      tm_perp_calendar = .true.
   end if

! Set variables from "comtim.h interface" for backwards compatibility.

! Calculation of ending timestep number (nestep) assumes a constant stepsize.
   ntspday = 86400/dtime
   diff = esmf_timeinit()
   call esmf_datediff(start_date, stop_date, diff, islater, rc)
   call chkrc(rc, sub//': error return from esmf_datediff calculating nestep')
   call esmf_timeget(diff, ndays, nsecs, rc)
   call chkrc(rc, sub//': error return from esmf_timeget calculating nestep')
   nestep = ntspday*ndays + nsecs/dtime
   if ( mod(nsecs,dtime) /= 0 ) nestep = nestep + 1

! Print configuration summary to log file (stdout).

   if (masterproc) then
      call timemgr_print()
   end if

end subroutine timemgr_restart
!=========================================================================================

subroutine timemgr_print()

   implicit none

! Local variables
   character(len=*), parameter :: sub = 'timemgr_print'
   integer :: rc
   integer :: day, sec, ymd, tod
   character(len=32) :: cal     ! Calendar to use in date calculations.
   integer ::&                  ! Data required to restart time manager:
      type      = uninit_int,  &! calendar type
      nstep     = uninit_int,  &! current step number
      step_days = uninit_int,  &! days component of timestep size
      step_sec  = uninit_int,  &! seconds component of timestep size
      start_ymd = uninit_int,  &! start date
      start_tod = uninit_int,  &! start time of day
      stop_ymd  = uninit_int,  &! stop date
      stop_tod  = uninit_int,  &! stop time of day
      ref_ymd   = uninit_int,  &! reference date
      ref_tod   = uninit_int,  &! reference time of day
      curr_ymd  = uninit_int,  &! current date
      curr_tod  = uninit_int    ! current time of day
!-----------------------------------------------------------------------------------------

   call esmf_timemgrrestartwrite(tm_id, type, nstep, step_days, step_sec, &
                                start_ymd, start_tod, stop_ymd, stop_tod, ref_ymd, &
                                ref_tod, curr_ymd, curr_tod, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrrestartwrite')

   write(6,*)' ********** Time Manager Configuration **********'

   if ( type == esmf_no_leap ) then
      cal = 'NO_LEAP'
   else if ( type == esmf_gregorian ) then
      cal = 'GREGORIAN'
   end if

   write(6,*)'  Calendar type:            ',trim(cal)
   write(6,*)'  Timestep size (seconds):  ', (step_days*86400 + step_sec)
   write(6,*)'  .. of dy_core (seconds):  ', int(dtdy+0.1)
   write(6,*)'  Start date (ymd tod):     ', start_ymd, start_tod
   write(6,*)'  Stop date (ymd tod):      ', stop_ymd, stop_tod
   write(6,*)'  Reference date (ymd tod): ', ref_ymd, ref_tod
   write(6,*)'  Current step number:      ', nstep
   write(6,*)'  Ending step number:       ', nestep
   write(6,*)'  Current date (ymd tod):   ', curr_ymd, curr_tod

   if ( tm_perp_calendar ) then
      call esmf_dateget(tm_perp_date, ymd, tod)
      write(6,*)'  Use perpetual diurnal cycle date (ymd): ', ymd
   end if

   write(6,*)' ************************************************'

end subroutine timemgr_print
!=========================================================================================

subroutine advance_timestep()

! Increment the timestep number.

   implicit none
   
! Local variables
   character(len=*), parameter :: sub = 'advance_timestep'
   integer :: rc
!-----------------------------------------------------------------------------------------

   call esmf_timemgradvance(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgradvance')

! Set first step flag off.

   tm_first_restart_step = .false.

end subroutine advance_timestep
!=========================================================================================

function get_step_size()

! Return the step size in seconds.

   implicit none
   
! Return value
   integer :: get_step_size

! Local variables
   character(len=*), parameter :: sub = 'get_step_size'
   integer :: days, seconds
   integer :: rc
!-----------------------------------------------------------------------------------------

   call esmf_timemgrgetstepsize(tm_id, days, seconds, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetstepsize')
   get_step_size = 86400*days + seconds

end function get_step_size
!=========================================================================================

function get_nstep()

! Return the timestep number.

   implicit none
   
! Return value
   integer :: get_nstep

! Local variables
   character(len=*), parameter :: sub = 'get_nstep'
   integer :: rc
!-----------------------------------------------------------------------------------------

   get_nstep = esmf_timemgrgetnstep(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetnstep')

end function get_nstep
!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

! Return date components valid at end of current timestep with an optional
! offset (positive or negative) in seconds.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

! Local variables
   character(len=*), parameter :: sub = 'get_curr_date'
   integer :: rc
   type(esmf_date) :: date
   type(esmf_time) :: off
   integer :: ymd
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetcurrdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetcurrdate')

   if (present(offset)) then
      if (offset > 0) then
         date = esmf_dateincrementsec(date, offset, rc)
         call chkrc(rc, sub//': error incrementing current date')
      else if (offset < 0) then
         off = esmf_timeinit(0, -offset, rc)
         call chkrc(rc, sub//': error setting offset time type')
         date = esmf_datedecrement(date, off, rc)
         call chkrc(rc, sub//': error decrementing current date')
      end if
   end if

   call esmf_dateget(date, ymd, tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   yr = ymd/10000
   mon = mod(ymd, 10000) / 100
   day = mod(ymd, 100)

end subroutine get_curr_date
!=========================================================================================

subroutine get_perp_date(yr, mon, day, tod)

! Return time of day valid at end of current timestep and the components
! of the perpetual date.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_perp_date'
   integer :: rc
   type(esmf_date) :: date
   integer :: ymd, idum
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetcurrdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetcurrdate')

   call esmf_dateget(date, idum, tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   call esmf_dateget(tm_perp_date, ymd, idum, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   yr = ymd/10000
   mon = mod(ymd, 10000) / 100
   day = mod(ymd, 100)

end subroutine get_perp_date
!=========================================================================================

subroutine get_prev_date(yr, mon, day, tod)

! Return date components valid at beginning of current timestep.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_prev_date'
   integer :: rc
   type(esmf_date) :: date
   integer :: ymd
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetprevdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetprevdate')

   call esmf_dateget(date, ymd, tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   yr = ymd/10000
   mon = mod(ymd, 10000) / 100
   day = mod(ymd, 100)

end subroutine get_prev_date
!=========================================================================================

subroutine get_start_date(yr, mon, day, tod)

! Return date components valid at beginning of initial run.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_start_date'
   integer :: rc
   type(esmf_date) :: date
   integer :: ymd
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetstartdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetstartdate')

   call esmf_dateget(date, ymd, tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   yr = ymd/10000
   mon = mod(ymd, 10000) / 100
   day = mod(ymd, 100)

end subroutine get_start_date
!=========================================================================================

subroutine get_ref_date(yr, mon, day, tod)

! Return date components of the reference date.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_ref_date'
   integer :: rc
   type(esmf_date) :: date
   integer :: ymd
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetbasedate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetbasedate')

   call esmf_dateget(date, ymd, tod, rc)
   call chkrc(rc, sub//': error return from esmf_dateget')
   yr = ymd/10000
   mon = mod(ymd, 10000) / 100
   day = mod(ymd, 100)

end subroutine get_ref_date
!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

   implicit none
   
! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'get_curr_time'
   integer :: rc
   type(esmf_date) :: cdate, rdate
   type(esmf_time) :: diff
   logical :: islater
!-----------------------------------------------------------------------------------------

   cdate = esmf_timemgrgetcurrdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetcurrdate')

   rdate = esmf_timemgrgetbasedate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetbasedate')

   call esmf_datediff(rdate, cdate, diff, islater, rc)
   call chkrc(rc, sub//': error return from esmf_datediff')

   call esmf_timeget(diff, days, seconds, rc)
   call chkrc(rc, sub//': error return from esmf_timeget')

end subroutine get_curr_time
!=========================================================================================

function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

   implicit none
   
! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
   real(r8) :: get_curr_calday

! Local variables
   character(len=*), parameter :: sub = 'get_curr_calday'
   integer :: rc
   type(esmf_date) :: date
   type(esmf_time) :: off
   integer :: ymd, tod
!-----------------------------------------------------------------------------------------

   date = esmf_timemgrgetcurrdate(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetcurrdate')

   if (present(offset)) then
      if (offset > 0) then
         date = esmf_dateincrementsec(date, offset, rc)
         call chkrc(rc, sub//': error incrementing current date')
      else if (offset < 0) then
         off = esmf_timeinit(0, -offset, rc)
         call chkrc(rc, sub//': error setting offset time type')
         date = esmf_datedecrement(date, off, rc)
         call chkrc(rc, sub//': error decrementing current date')
      end if
   end if

   if ( tm_perp_calendar ) then
      call esmf_dateget(date, ymd, tod, rc)
      call chkrc(rc, sub//': error return from esmf_dateget')
      date = esmf_dateincrementsec(tm_perp_date, tod, rc)
      call chkrc(rc, sub//': error return from esmf_dateincrementsec')
   end if

   get_curr_calday = esmf_dategetfltdayofyear(date, rc)
   call chkrc(rc, sub//': error return from esmf_dategetfltdayofyear')

end function get_curr_calday
!=========================================================================================

function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

   implicit none
   
! Return value
   logical :: is_end_curr_day

! Local variables
   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
!-----------------------------------------------------------------------------------------

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_day = (tod == 0)

end function is_end_curr_day
!=========================================================================================

function is_end_curr_month()

! Return true if current timestep is last timestep in current month.

   implicit none
   
! Return value
   logical :: is_end_curr_month

! Local variables
   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
!-----------------------------------------------------------------------------------------

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_month = (day == 1  .and.  tod == 0)

end function is_end_curr_month
!=========================================================================================

function is_first_step()

! Return true on first step of initial run only.

   implicit none
   
! Return value
   logical :: is_first_step

! Local variables
   character(len=*), parameter :: sub = 'is_first_step'
   integer :: rc
   integer :: nstep
!-----------------------------------------------------------------------------------------

   nstep = esmf_timemgrgetnstep(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrgetnstep')
   is_first_step = (nstep == 0)

end function is_first_step
!=========================================================================================

function is_first_restart_step()

! Return true on first step of restart run only.

   implicit none
   
! Return value
   logical :: is_first_restart_step
!-----------------------------------------------------------------------------------------

   is_first_restart_step = tm_first_restart_step

end function is_first_restart_step
!=========================================================================================

function is_last_step()

! Return true on last timestep.

   implicit none
   
! Return value
   logical :: is_last_step

! Local variables
   character(len=*), parameter :: sub = 'is_last_step'
   integer :: rc
!-----------------------------------------------------------------------------------------

   is_last_step = esmf_timemgrlaststep(tm_id, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrlaststep')

end function is_last_step
!=========================================================================================

function is_perpetual()

! Return true on last timestep.

   implicit none
   
! Return value
   logical :: is_perpetual
!-----------------------------------------------------------------------------------------

   is_perpetual = tm_perp_calendar

end function is_perpetual
!=========================================================================================

subroutine timemgr_write_restart(ftn_unit)

! Write information needed on restart to a binary Fortran file.
! It is assumed that this routine is called only from the master proc if in SPMD mode.

   implicit none

! Arguments
   integer, intent(in) :: ftn_unit  ! Fortran unit number

! Local variables
   character(len=*), parameter :: sub = 'timemgr_write_restart'
   integer :: rc   ! return code
   integer :: tod
!-----------------------------------------------------------------------------------------

   call esmf_timemgrrestartwrite(tm_id, rst_type, rst_nstep, rst_step_days, rst_step_sec, &
                                rst_start_ymd, rst_start_tod, rst_stop_ymd, rst_stop_tod, rst_ref_ymd, &
                                rst_ref_tod, rst_curr_ymd, rst_curr_tod, rc)
   call chkrc(rc, sub//': error return from esmf_timemgrrestartwrite')

   if ( tm_perp_calendar ) then
      call esmf_dateget(tm_perp_date, rst_perp_ymd, tod, rc)
      call chkrc(rc, sub//': error return from esmf_dateget')
      rst_perp_cal = tm_perp_calendar
   end if

   write(ftn_unit, iostat=rc) rst_type, rst_nstep, rst_step_days, rst_step_sec, &
                              rst_start_ymd, rst_start_tod, rst_stop_ymd, rst_stop_tod, rst_ref_ymd, &
                              rst_ref_tod, rst_curr_ymd, rst_curr_tod, rst_perp_ymd, rst_perp_cal

   if (rc /= 0 ) then
      write (6,*) 'WRITE iostat= ',rc,' on i/o unit = ',ftn_unit
      call endrun
   end if

end subroutine timemgr_write_restart
!=========================================================================================

subroutine timemgr_read_restart(ftn_unit)

! Read information needed to restart from a binary Fortran file.
! It is assumed that this routine is called only from the master proc if in SPMD mode.

   implicit none

! Arguments
   integer, intent(in) :: ftn_unit  ! Fortran unit number

! Local variables
   character(len=*), parameter :: sub = 'timemgr_read_restart'
   integer :: rc   ! return code
!-----------------------------------------------------------------------------------------

   read(ftn_unit, iostat=rc) rst_type, rst_nstep, rst_step_days, rst_step_sec, &
                             rst_start_ymd, rst_start_tod, rst_stop_ymd, rst_stop_tod, rst_ref_ymd, &
                             rst_ref_tod, rst_curr_ymd, rst_curr_tod, rst_perp_ymd, rst_perp_cal

   if (rc /= 0 ) then
      write (6,*) 'READ iostat= ',rc,' on i/o unit = ',ftn_unit
      call endrun
   end if

end subroutine timemgr_read_restart
!=========================================================================================

subroutine chkrc(rc, mes)
   implicit none
   integer, intent(in)          :: rc   ! return code from time management library
   character(len=*), intent(in) :: mes  ! error message
   if ( rc == esmf_success ) return
   write(6,*) mes
   call esmf_errprint(rc)
   call endrun
end subroutine chkrc

end module time_manager
