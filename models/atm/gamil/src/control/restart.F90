#include <misc.h>
#include <params.h>

module restart
!----------------------------------------------------------------------- 
! 
! module to handle reading and writing of the master restart files.
!
!----------------------------------------------------------------------- 
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: masterproc, plev, plevp, plond, plat,plevstd,plond,plat
  use rgrid,        only: nlon
  use ioFileMod,    only: putfil, getfil, opnfil
  use filenames,    only: get_archivedir

#ifdef SPMD
  use mpishorthand, only: mpicom, mpir8, mpiint, mpilog
#endif

  implicit none

  public write_restart
  public read_restart
  public set_restart_filepath

  private

  integer, parameter :: nlen = 256 ! Length of character strings
  character(len=nlen) :: pname     ! Full restart pathname
!
! Filename specifiers for master restart filename
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
!
  character(256) :: rfilename_spec = '%c.gamil.r.%y-%m-%d-%s'

#include <comctl.h>
#include <comlun.h>
#include <comhyb.h>

contains

  subroutine write_restart
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Write the primary, secondary, and history buffer regeneration files.
! 
! Method: 
! The cpp SPMD definition provides for the funnelling of all program i/o
! through the master processor. Processor 0 either reads restart/history
! data from the disk and distributes it to all processors, or collects
! data from all processors and writes it to disk.
! 
! Author: 
! 
!----------------------------------------------------------------------- 
    use history,          only: write_restart_history
    use restart_physics,  only: write_restart_physics
    use restart_dynamics, only: write_restart_dynamics
    use time_manager,     only: timemgr_write_restart
    use filenames,        only: caseid, mss_wpass, mss_irt, get_archivedir, &
                                interpret_filename_spec

    integer ioerr         ! write error status
    character(256) fname  ! Restart filename

    ! Write the primary restart datasets
    !
    ! Everything inside the following "if" is from the obsolete RGNFLS
    !
    if (masterproc) then
      ! Open master restart file.
      fname = interpret_filename_spec(rfilename_spec)
      call opnfil(fname, nrg, 'u')

      ! Write the master restart dataset
      call timemgr_write_restart(nrg)

      write (nrg, iostat=ioerr) caseid, hyai, hybi, hyam, hybm, sig, sigl, dsig, pmtop*100, aeres
      if (ioerr /= 0) then
        write(6, "('[Error]: write_restart: Failed with ioerr ', I,' on i/o unit = ', I)") ioerr, nrg
        call endrun
      end if

      ! Reduced grid stuff
      write (nrg, iostat=ioerr) nlon
      if (ioerr /= 0) then
        write(6, "('[Error]: write_restart: Failed with ioerr ', I,' on i/o unit = ', I)") ioerr, nrg
        call endrun
      end if
    end if

    ! Dynamics, physics, History
    call write_restart_dynamics(nrg)
    call write_restart_physics(nrg, nrg2)
    call write_restart_history(nrg, luhrest)

    if (masterproc) then
      close(nrg)
      pname = trim(get_archivedir('rest')) // fname
      call putfil(fname, pname, mss_wpass, mss_irt, (.not. nlend))
      call write_rest_pfile
    end if

    return
  end subroutine write_restart

  subroutine read_restart
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Acquire and position the restart, master, primary and secondary
! datasets for a continuation run
! 
!-----------------------------------------------------------------------
    use phys_grid,        only: phys_grid_init
    use history,          only: read_restart_history
    use filenames,        only: caseid
    use restart_physics,  only: read_restart_physics
    use restart_dynamics, only: read_restart_dynamics
    use time_manager,     only: timemgr_read_restart, timemgr_restart

    integer i
    integer ioerr
    integer decomp                ! decomposition type

    character(256) locfn          ! Local filename
    character(33)  tcase          ! Read in previous case name

    ! Determine full restart pathname
    !
    ! Test for existence of a restart pointer dataset (nsrest=1 only).
    ! If restart pinter dataset exists, obtain it and read it.
    !
    if (masterproc) then
      if (.not. nlhst) then
        call read_rest_pfile
      end if
      ! Obtain and read the master restart dataset 
      !
      ! Obtain master restart dataset
      !
      call getfil(pname, locfn)
      call opnfil(locfn, nrg, 'u')

      ! Read comtim, comtfc, comhst, comhed and history buffer variables 
      call timemgr_read_restart(nrg)

      read (nrg, iostat=ioerr) tcase, hyai, hybi, hyam, hybm, sig, sigl, dsig, pmtop, aeres
      if (ioerr /= 0) then
        write(6, "('[Error]: read_restart: Failed with ioerr ', I,' on i/o unit = ', I)") ioerr, nrg
        call endrun
      end if

      if (lbrnch .and. tcase == caseid) then
        write(6, *) 'READ_RESTART: Must change case name on branch run'
        write(6, *) 'Prev case = ',tcase,' current case = ',caseid
        call endrun
      end if

      ! Reduced grid stuff
      read (nrg, iostat=ioerr) nlon
      if (ioerr /= 0) then
        write(6, "('[Error]: read_restart: Failed with ioerr ', I,' on i/o unit = ', I)") ioerr, nrg
        call endrun
      end if

      write(6,*) 'Files for restart:'

    end if

#if ( defined SPMD ) 
    call mpibcast(hyai   ,plev+1   ,mpir8  ,0,mpicom) 
    call mpibcast(hybi   ,plev+1   ,mpir8  ,0,mpicom) 
    call mpibcast(hyam   ,plev     ,mpir8  ,0,mpicom)   
    call mpibcast(hybm   ,plev     ,mpir8  ,0,mpicom)   
    call mpibcast(sig    ,plev+1   ,mpir8  ,0,mpicom)   
    call mpibcast(sigl   ,plev     ,mpir8  ,0,mpicom)   
    call mpibcast(dsig   ,plev     ,mpir8  ,0,mpicom)   
    call mpibcast(pmtop  ,1        ,mpir8  ,0,mpicom)
    call mpibcast(aeres  ,1        ,mpilog ,0,mpicom)
    call mpibcast(nlon,   plat,     mpiint, 0,mpicom)
#endif

    ! Restart the time manager.
    call timemgr_restart()

    ! Dynamics, physics, history
    call read_restart_dynamics (nrg)
    call initcom()
    call phys_grid_init
    call read_restart_physics (nrg, nrg2, aeres )
    if (nsrest == 1) then
      !call read_restart_history (nrg, luhrest)
    end if
    if (masterproc) close(nrg)

    return
  end subroutine read_restart

  subroutine write_rest_pfile
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Write out the restart pointer file
!
!----------------------------------------------------------------------- 
    use filenames,       only: rest_pfile
    use restart_physics, only: get_abs_restart_filepath
    use history,         only: get_mtapes, get_hist_restart_filepath, &
                               hstwr, get_hfilepath, nfils, mfilt
!-----------------------------------------------------------------------
! for nsds
#include <comlun.h>
!-----------------------------------------------------------------------
! for aeres
#include <comctl.h>
!-----------------------------------------------------------------------
    integer t      ! Tape number
    integer mtapes ! Number of tapes that are active

    call opnfil(rest_pfile, nsds, 'f')
    rewind nsds
    write (nsds,'(a)') trim(pname)
    write (nsds,'(//a,a)') '# The following lists the other files needed for restarts', &
                           ' (gamil only reads the first line of this file).'
    write (nsds,'(a,a)') '# The files below refer to the files needed for the master restart file:', &
                         trim(pname)
    if (aeres) then
      write (nsds,'(a,a)') '# ', trim(get_abs_restart_filepath())
    end if
!
! History files: Need restart history files when not a time-step to write history info
! Need: history files if they are not full
!
    mtapes = get_mtapes( )
    do t = 1, mtapes
      if ( .not. hstwr(t) ) then
         write (nsds,'(a,a)') '# ', trim(get_hist_restart_filepath( t ))
      end if
      if ( nfils(t) > 0 .and. nfils(t) < mfilt(t) ) then
         write (nsds,'(a,a)') '# ', trim(get_hfilepath( t ))
      end if
    end do
    close (nsds)

    write(6,*)'(WRITE_REST_PFILE): successfully wrote local restart pointer file ', trim(rest_pfile)
    write(6,'("---------------------------------------")')

   end subroutine write_rest_pfile

   subroutine read_rest_pfile
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Read the master restart file from the restart pointer file
!
!----------------------------------------------------------------------- 
    use filenames, only: rest_pfile
!-----------------------------------------------------------------------
! for nsds
#include <comlun.h>

    character(len=256) :: locfn    ! Local pathname for restart pointer file

    call opnfil(rest_pfile, nsds, 'f', status="old")
    read (nsds,'(a)') pname
    close(nsds)

  end subroutine read_rest_pfile

!-----------------------------------------------------------------------
! BOP
!
! !ROUTINE: set_restart_filepath
!
! !DESCRIPTION: Set the filepath of the specific type of restart file.
!
!-----------------------------------------------------------------------
  subroutine set_restart_filepath(rgpath)

    character(len=*), intent(in)  :: rgpath ! Full pathname to restart file

    if (trim(rgpath) == '') then
      write(6,*) 'SET_RESTART_FILEPATH: rgpath sent into subroutine is empty'
      call endrun
    end if
    if (rgpath(1:1) /= '/')then
      write(6,*) 'SET_RESTART_FILEPATH: rgpath sent into subroutine is not an absolute pathname'
      call endrun
    end if
    if (len_trim(rgpath) > nlen)then
      write(6,*) 'SET_RESTART_FILEPATH: rgpath is too long :', rgpath
      call endrun
    end if
    pname = trim(rgpath)

  end subroutine set_restart_filepath

end module restart
