#include <misc.h>
#include <preproc.h>

module histHandlerMod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : maxhist
  implicit none

! Variables

  logical, private :: lrestwrt           !true => write restart file at this time step
  integer, public  :: mcdate_i(maxhist)  !start of history interval date (yyyymmdd format)
  integer, public  :: mcsec_i(maxhist)   !start of history interval seconds of current date
  integer, public  :: mdcur_i(maxhist)   !start of history interval day
  integer, public  :: mscur_i(maxhist)   !start of history interval seconds of current day

! Methods

   public  :: histhandler         !main history and restart handler
   public  :: histend             !determines if end of history interval
   public  :: do_restwrite        !returns logical setting if time to write restart data
   private :: close_and_disp      !determine if file needs to be closed and archived
   private :: set_hist_filename   !sets history dataset filenames

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine histhandler ()

!-----------------------------------------------------------------------
!
! Purpose:
! History file handler
!
! Method:
! This code does the following for every time step:
!   o increments history field accumulation counters
!   o determines if next time step is beginning of history interval
!
! This code does the following at the end of a history interval:
!   o increments the current time sample counter: ntim <= hist_mfilt
!   o opens a new history file if needed (i.e., when ntim = 1)
!   o writes history data to current history file
!   o resets field accumulation counters to zero
!
! This code does the following when a history file is full
! (i.e., ntim = hist_mfilt) or at the last time step of the simulation:
!   o closes history file and disposes to mass store if appropriate
!   o resets time sample counter to zero if file is full
!
! This code does the following when the primary history file is full
! (i.e., ntim = hist_mfilt) or at the last time step of the simulation:
!   o writes restart file and disposes to mass store if appropriate
!
! Author: Gordon Bonan
!
!-----------------------------------------------------------------------
! $Id: histHandlerMod.F90,v 1.13.2.6.6.1 2002/10/03 20:07:38 erik Exp $
!-----------------------------------------------------------------------

    use clm_varctl, only : locfnh, hist_nhtfrq, hist_mfilt, &
                           archive_dir, mss_wpass, mss_irt
    use histFileMod, only : nhist, nbeghis, ntim, timcom, ncid, ehi, &
                            histcrt, histwrt, histcls , histzero
    use fileutils, only : set_filename, putfil
    use spmdMod, only : masterproc
    use time_manager, only : get_nstep, get_curr_date, get_curr_time
    use shr_const_mod, only: SHR_CONST_CDAY
    include 'netcdf.inc'

! ------------------------ local variables ------------------------
    integer  :: i                          !loop index
    integer  :: m                          !history file do loop counter
    integer  :: ier                        !error code
    integer  :: day                        !current day (1 -> 31)
    integer  :: mon                        !current month (1 -> 12)
    integer  :: yr                         !current year (0 -> ...)
    integer  :: mcsec                      !seconds of current date
    integer  :: mdcur                      !current day
    integer  :: mscur                      !seconds of current day
    integer  :: mcdate                     !current date
    logical  :: lhisdisp(nhist)            !true => save and dispose history file
    logical  :: lremovh(nhist)             !true => remove local history file after dispose
    logical  :: lstop                      !true => last time step of run
    real(r8) :: frac_i                     !fractional day, start of time sample
    real(r8) :: frac                       !current fractional day
    real(r8) :: hour_i                     !fractional hour,start of time sample
    real(r8) :: hour                       !current fractional hour
    character(len=256) :: loc_fn           !local and remote filenames
    character(len=256) :: rem_dir          !remote (archive) directory
    character(len=256) :: rem_fn           !remote (archive) filename
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Loop over history files: Increment history field counters,
! create new history files if necessary and write data to history
! files if end of history interval.
! -----------------------------------------------------------------

    do m = 1, nhist

! Set [nbeghis] to one to indicate next time step is start of time sample

       if (ehi(m)) then
          nbeghis(m) = 1
       else
          nbeghis(m) = 0
       end if

! End of history interval?

       if (ehi(m)) then

! Increment current time sample counter.  If first time sample
! generate unique history file name and open netCDF file.

          ntim(m) = ntim(m) + 1
          if (masterproc) then
             if (ntim(m) == 1) then
                locfnh(m) = set_hist_filename (hist_freq=hist_nhtfrq(m), hist_file=m)
                write(6,*)'(HISTHANDLER): Creating history file ', &
                     trim(locfnh(m)),' at nstep = ',get_nstep()
                call histcrt (m)
             endif
          end if

! Build time comment for current time sample based on start of
! time slice calendar info and  current time calendar info

          if (masterproc) then
             call get_curr_time (mdcur, mscur)
             call get_curr_date (yr, mon, day, mcsec)
             mcdate = yr*10000 + mon*100 + day
             frac_i = float(mscur_i(m))/SHR_CONST_CDAY
             frac   = float(mscur     )/SHR_CONST_CDAY
             hour_i = float(mcsec_i(m))/SHR_CONST_CDAY*24. !/3600. s per hr
             hour   = float(mcsec     )/SHR_CONST_CDAY*24. !/3600. s per hr
             write(timcom(m),200) &
                  mdcur_i(m),frac_i,mdcur,frac,hour_i,mcdate_i(m),hour,mcdate
200          format ('TIME AVG FOR DAYS: ',i6.6,f4.3,'-',i6.6,f4.3, &
                  ' DATES:',f6.3,'Z ',i8.8,'-',f6.3,'Z ',i8.8)
             write(6,*)
             write(6,*)'(HISTHANDLER): Writing current time sample to local history file ',&
                  trim(locfnh(m)),' at nstep = ',get_nstep()
             write(6,*) trim(timcom(m))
             write(6,*)
          endif

! Write history time sample

          call histwrt(m)

! Zero necessary history buffer information

          call histzero(m)

       end if

    end do  ! end loop over history tapes

! -----------------------------------------------------------------
! Loop over history files and if history file is full
! (ntim = hist_mfilt) or if last time step
! o close history file and dispose to mass store if appropriate
! o reset [ntim] time sample counter to zero (only if file is full)
! -----------------------------------------------------------------

! determine if file needs to be closed and disposed

    call close_and_disp (lhisdisp, lstop, lremovh)

! loop over history files

    do m = 1, nhist

       if (lhisdisp(m)) then

! Close open history file and dispose to mass store
! Auxilary files may have been closed and saved off without being full,
! must reopen the files

          if (masterproc) then
             if (ntim(m) /= 0) then
                write(6,*)
                write(6,*) '(HISTHANDLER): Closing local history file ',&
                     trim(locfnh(m)),' at nstep = ', get_nstep()
                write(6,*)
                call histcls(m)
                rem_dir = trim(archive_dir) // '/hist/'
                rem_fn = set_filename(rem_dir, locfnh(m))
                call putfil (locfnh(m), rem_fn, mss_wpass, mss_irt, lremovh(m))
                if (.not.lstop .and. (ntim(m)/=hist_mfilt(m))) then
                   call wrap_open (trim(locfnh(m)), NF_WRITE, ncid(m))
                end if
             else
                write(6,*)'(HISTHANDLER): history file ',m,': no open file to close'
             end if
          endif

! reset number of time samples to zero if file is full

          if (ntim(m) == hist_mfilt(m)) then
             ntim(m) = 0
             locfnh(m) = ' '
          end if

       end if  ! end of if-dispose block

    end do     ! end of loop over tapes

    return
  end subroutine histhandler

!=======================================================================

  subroutine histend ()

!-----------------------------------------------------------------------
!
! Purpose:
! Determine if end of history interval
!
! Method:
!
! Daily-averaged data for the first day in September
! are written on mcdate = 00/09/02 with mscur = 0
!
! In general: daily-averaged data for the first day in
! month mm are written on mcdate = yyyy/mm/02 with mscur = 0
!
! Daily-averaged data for the 30th day (last day in September)
! are written on mcdate = 0000/10/01 mscur = 0
!
! In general: daily-averaged data for the last day in
! month mm are written on mcdate = yyyy/mm+1/01 with mscur = 0
!
! Author: Gordon Bonan
!
!-----------------------------------------------------------------------

    use clm_varctl , only : hist_mfilt, hist_nhtfrq, hist_crtinic
    use histFileMod, only : ehi, nhist, nbeghis
    use time_manager, only : get_nstep, get_curr_date, get_prev_date, get_curr_time

! ------------------------ local variables ------------------------
    integer :: m          !history file (1, ..., maxhist)
    integer :: nstep      !current step
    integer :: day        !nstep day (1 -> 31)
    integer :: mon        !nstep month (1 -> 12)
    integer :: yr         !nstep year (0 -> ...)
    integer :: mcsec      !nstep time of day [seconds]
    integer :: daym1      !nstep-1 day (1 -> 31)
    integer :: monm1      !nstep-1 month (1 -> 12)
    integer :: yrm1       !nstep-1 year (0 -> ...)
    integer :: mcsecm1    !nstep-1 time of day [seconds]
    integer :: mcdate     !nstep date in integer format [yyyymmdd]
    integer :: mcdatem1   !nstep-1 date in integer format [yyyymmdd]
    integer :: mdcur      !day for current nstep
    integer :: mscur      !seconds of current day
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Calendar calculations
! -----------------------------------------------------------------

! get current step

    nstep = get_nstep()

! Set calendar for current time step

    call get_curr_date (yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day

! Set calendar for current for previous time step

    call get_prev_date (yrm1, monm1, daym1, mcsecm1)
    mcdatem1 = yrm1*10000 + monm1*100 + daym1

! Set elapased time since reference date

    call get_curr_time(mdcur, mscur)

! -----------------------------------------------------------------
! Process all active history files
! -----------------------------------------------------------------

! loop over history tapes

    do m = 1, nhist

! Skip nstep=0 if monthly average

       if (nstep==0 .and. hist_nhtfrq(m)==0) cycle

! Determine if end of history interval

       ehi(m) = .false.
       if (hist_nhtfrq(m)==0) then   !monthly average
          if (mon /= monm1) ehi(m) = .true.
       else
          if (mod(nstep,hist_nhtfrq(m)) == 0) ehi(m) = .true.
       end if

! Calendar info for current time sample: start of time interval

       if (nbeghis(m) == 1) then
          mcdate_i(m) = mcdate
          mcsec_i(m) = mcsec
          mdcur_i(m) = mdcur
          mscur_i(m) = mscur
       end if

    end do

    return
  end subroutine histend

!=======================================================================

  subroutine close_and_disp (lhisdisp, lstop, lremovh)

!-----------------------------------------------------------------------
!
! Purpose:
! Determine logic for closeing and/or disposing history file
! Sets values for logical variables lhisdisp, lstop, lremovh (arguments)
! and lrestwrt (module variable)
!
! Method:
!
! Logic to determine if time to dispose history file:
! o if history file is full (ntim = mfilt)
! o if have reached end of run
!
! Logic to determine if time to write a restart file:
! o if run through flux coupler
!   o when coupler gives signal to write restart
!   o when coupler gives signal to stop
! o if not run through flux coupler
!   o when primary history file is closed and disposed
!
! Logic to determine if time to remove history file:
! o remove history files unless this is end of run or
!   history file is not full.
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varctl  , only : hist_mfilt
    use histFileMod , only : nhist, ntim
#ifdef COUP_CSM
    use clm_csmMod
#endif
    use time_manager, only : is_last_step

! ------------------------ arguments ------------------------------
    logical, intent(out) :: lhisdisp(nhist) !true => save and dispose history file
    logical, intent(out) :: lstop           !true => last time step of run
    logical, intent(out) :: lremovh(nhist)  !true => remove local history file after dispose
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer m        ! loop index
! -----------------------------------------------------------------

    lstop = .false.
    lrestwrt = .false.
    lhisdisp(:) = .false.
    lremovh(:) = .false.

#if (defined OFFLINE) || (defined COUP_CAM)

! If end of run dispose all history files and write restart

    if (is_last_step()) then
       lstop = .true.
       lrestwrt = .true.
       lhisdisp(:) = .true.
       lremovh(:) = .false.
       RETURN
    endif

! If time to dispose master history file then dispose all
! history files and write restart

    if (ntim(1)==hist_mfilt(1))  then
       lrestwrt = .true.
       lhisdisp(:) = .true.
       lremovh(1) = .true.
       do m = 2,nhist
          if (ntim(m)==hist_mfilt(m)) then
             lremovh(m) = .true.
          else
             lremovh(m) = .false.
          endif
       end do
       RETURN
    endif

! If not end of run or time to dispose master history file
! then determine if time to dispose individual auxillary files

    lrestwrt = .false.
    do m = 2,nhist
       if (ntim(m)==hist_mfilt(m)) then
          lhisdisp(m) = .true.
          lremovh(m) = .true.
       endif
    end do
    RETURN

#elif (defined COUP_CSM)

! If coupler says that next time step is end of run then
! dispose all history files and write restart

    if (csmstop_next) then
       lstop = .true.
       if (csmrstrt) then
          lrestwrt = .true.
       else
          lrestwrt = .false.
       endif
       lhisdisp(:) = .true.
       lremovh(:) = .false.
       RETURN
    endif

! If coupler says to write restart then dispose all history
! files and write restart

    if (csmrstrt) then
       lhisdisp(:) = .true.
       lrestwrt = .true.
       lremovh(:) = .true.
       RETURN
    endif

! Otherwise check if file is full and dispose if it is

    lrestwrt = .false.
    lremovh(:) = .false.
    do m=1,nhist
       if (ntim(m) == hist_mfilt(m)) then
          lhisdisp(m) = .true.
          lremovh(m) = .true.
       endif
    end do
    lremovh(:) = .true.
    RETURN

#endif

    return
100 continue

  end subroutine close_and_disp

!=======================================================================

  logical function do_restwrite()

!-----------------------------------------------------------------------
!
! Purpose:
! Determine if restart dataset is to be written at this time step
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    do_restwrite = .false.
    if (lrestwrt) do_restwrite = .true.

  end function do_restwrite

!=======================================================================

  character(len=256) function set_hist_filename (hist_freq, hist_file)

!-----------------------------------------------------------------------
!
! Purpose:
! Determine history dataset filenames
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------

    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date, get_prev_date

! ------------------------ arguments ------------------------------
    integer, intent(in)  :: hist_freq   !history file frequency
    integer, intent(in)  :: hist_file   !history file index
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    character(len=256) :: cdate       !date char string
    character(len=  1) :: hist_index  !0, 1 or 2 (currently)
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
! -----------------------------------------------------------------

    if (hist_freq == 0 ) then   !monthly
       call get_prev_date (yr, mon, day, sec)
       write(cdate,'(i4.4,"-",i2.2)') yr,mon
    else                        !other
       call get_curr_date (yr, mon, day, sec)
       write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    endif
    write(hist_index,'(i1.1)') hist_file - 1
    set_hist_filename = "./"//trim(caseid)//".clm2.h"//hist_index//"."//&
         trim(cdate)//".nc"

  end function set_hist_filename

!=======================================================================

end module histHandlerMod


