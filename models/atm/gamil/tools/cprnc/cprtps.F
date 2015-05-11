program cprnc
!
! $Id: cprtps.F,v 1.3.6.1 2002/05/13 17:59:24 erik Exp $
!
  use header
  use nldat

  implicit none
      
  include 'netcdf.inc'
!
! Local workspace
!
  character*128 cvsid
  character*180 arg
  character*180 file(2)

  integer, dimension(2) :: nlon  = -1
  integer, dimension(2) :: nlat  = -1 
  integer, dimension(2) :: nlev  = -1
  integer, dimension(2) :: nlevp = -1
  integer n, nn, nargs
  integer numcases
  integer ntimemax
  integer ret(2)

  logical twocases
  logical itexists
  logical lnorm             ! print level by level L2 and L-INF norm stats
  logical constps
  logical verbose

  integer ierr    ! Error code for allocate

  namelist/nlcpr/nsteps, nstepe, iprs, ipre, kprs, kpre, latprs, &
                 latpre, matchts
!
! Externals
!
  integer lenchr,iargc
  external lenchr,iargc

  cvsid = '$Id: cprtps.F,v 1.3.6.1 2002/05/13 17:59:24 erik Exp $'
  write(6,*)'You are using cprnc with cvsid:',cvsid
!
! Default settings before parsing argument list
!
  file(:) = ' '
  nsteps = 0
  nstepe = 99999999
  lnorm = .false.
  constps = .false.
  matchts = .true.
  twocases = .false.
  verbose = .false.

  nargs = iargc()
  n = 1
  do while (n .le. nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1
    if (arg .eq. '-c') then
      constps = .true.
      write (6,*)'Assuming constant PS everywhere'
    else if (arg .eq. '-e') then
      call getarg (n, arg)
      n = n + 1
      if (n.gt.nargs .or. arg(1:1).lt.'0' .or. arg(1:1).gt.'9') then
        call usage_exit ('The -e option requires an input ending timestep')
      end if
      read (arg, '(f15.7)') nstepe
      write (6,*) 'Stopping analysis after NSTEP=', nstepe
    else if (arg .eq. '-m') then
      matchts = .true.
    else if (arg .eq. '-n') then
      lnorm = .true.
    else if (arg .eq. '-s') then
      call getarg (n, arg)
      n = n + 1
      if (n.gt.nargs .or. arg(1:1).lt.'0' .or. arg(1:1).gt.'9') then
        call usage_exit ('The -s option requires an input start timestep')
      end if
      read (arg, '(f15.7)') nsteps
      write (6,*) 'Starting analysis at NSTEP=', nsteps
    else if (arg .eq. '-v') then
      verbose = .true.
    else
      if (file(1).eq.' ') then
        file(1) = arg(1:lenchr(arg))//char(0)
        write (6,*) 'file 1=',file(1)(1:lenchr(file(1)))
      else if (file(2).eq.' ') then
        file(2) = arg(1:lenchr(arg))//char(0)
        write (6,*) 'file 2=',file(2)(1:lenchr(file(2)))
        twocases = .true.
      else
        call usage_exit (' ')
      end if
    end if
  end do
      
  if (file(1).eq.' ') then
    call usage_exit ('You must enter at least 1 input file')
  else if (n.lt.nargs) then
    call usage_exit ('Either -d or -g must be specified')
  end if
!
! Default namelist settings
!
  iprs = 0
  ipre = 0
  kprs = 0
  kpre = 0
  latprs = 0
  latpre = 0

  inquire (file='nl.cpr', exist=itexists)
  if (itexists) then
    open (unit=9, file='nl.cpr', status='OLD')
    read (9,nlcpr)
    close (9)
    if (iprs.gt.0) then

      if (ipre.le.0 .or. ipre.lt.iprs) then
        write(6,*)'cprnc: IPRS set but IPRE invalid or not set'
        stop 99
      end if

      if (latprs.le.0 .or. latpre.lt.latprs) then
        write(6,*)'cprnc: LATPRS and/or LATPRE invalid or not set'
        stop 99
      end if

      if (kprs.le.0 .or. kpre.lt.kprs) then
        write(6,*)'cprnc: KPRS and/or KPRE invalid or not set'
        stop 99
      end if
    end if
  end if

  numcases = 1
  if (twocases) numcases = 2
!
! Open files, then get dimension and variable id info
!
  do n=1,numcases
    inquire (file=file(n), exist=itexists)
    if (.not.itexists) then
      write(6,*)'cprnc: Unable to find input file ',file(n)
      stop 99
    end if

    call wrap_open (file(n)(1:lenchr(file(n))), NF_NOWRITE, ncid(n))

    call wrap_get_att_text (ncid(n), NF_GLOBAL, 'case', case(n))
    call wrap_get_att_text (ncid(n), NF_GLOBAL, 'title', title(n))

    call wrap_inq_dimid (ncid(n), 'lon', londimid(n))
    call wrap_inq_dimid (ncid(n), 'lat', latdimid(n))
    call wrap_inq_dimid (ncid(n), 'lev', levdimid(n))
    call wrap_inq_dimid (ncid(n), 'ilev', ilevdimid(n))
    call wrap_inq_dimid (ncid(n), 'time', unlimdimid(n))

    call wrap_inq_dimlen (ncid(n), londimid(n), nlon(n))
    call wrap_inq_dimlen (ncid(n), latdimid(n), nlat(n))
    call wrap_inq_dimlen (ncid(n), levdimid(n), nlev(n))
    call wrap_inq_dimlen (ncid(n), ilevdimid(n), nlevp(n))
    call wrap_inq_dimlen (ncid(n), unlimdimid(n), ntime(n))

    call wrap_inq_varid (ncid(n), 'lon', lonid(n))
    call wrap_inq_varid (ncid(n), 'lat', latid(n))
    call wrap_inq_varid (ncid(n), 'lev', levid(n))
    call wrap_inq_varid (ncid(n), 'ilev', ilevid(n))
    call wrap_inq_varid (ncid(n), 'time', timeid(n))

    call wrap_inq_varid (ncid(n), 'P0', p0id(n))
    call wrap_inq_varid (ncid(n), 'ntrm', ntrmid(n))
    call wrap_inq_varid (ncid(n), 'ntrn', ntrnid(n))
    call wrap_inq_varid (ncid(n), 'ntrk', ntrkid(n))
    call wrap_inq_varid (ncid(n), 'hyai', hyaiid(n))
    call wrap_inq_varid (ncid(n), 'hybi', hybiid(n))
    call wrap_inq_varid (ncid(n), 'hyam', hyamid(n))
    call wrap_inq_varid (ncid(n), 'hybm', hybmid(n))
    call wrap_inq_varid (ncid(n), 'ndbase', ndbaseid(n))
    call wrap_inq_varid (ncid(n), 'nsbase', nsbaseid(n))
    call wrap_inq_varid (ncid(n), 'nbdate', nbdateid(n))
    call wrap_inq_varid (ncid(n), 'nbsec', nbsecid(n))
    call wrap_inq_varid (ncid(n), 'mdt', mdtid(n))

    call wrap_inq_varid (ncid(n), 'date', ncdateid(n))
    call wrap_inq_varid (ncid(n), 'datesec', ncsecid(n))
    call wrap_inq_varid (ncid(n), 'gw', gwid(n))
!
! Variables that should not result in a fatal error if not found
!
    ret(1) = nf_inq_varid (ncid(n), 'nsteph', nstephid(n))
    if (ret(1).ne.NF_NOERR) nstephid(n) = -1

    ret(1) = nf_inq_varid (ncid(n), 'date_written', date_writtenid(n))
    if (ret(1).ne.NF_NOERR) date_writtenid(n) = -1

    ret(2) = nf_inq_varid (ncid(n), 'time_written', time_writtenid(n))
    if (ret(2).ne.NF_NOERR) time_writtenid(n) = -1
  end do
!
! Allocate space for header variables
!
  ntimemax = ntime(1)
  if (twocases .and. ntime(2).gt.ntime(1)) then
    ntimemax = ntime(2)
  end if

  ierr = 0
  if ( .not. allocated(date_written) ) &
  allocate(date_written(ntimemax,numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(time_written) ) &
  allocate(time_written(ntimemax,numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(nsteph) ) &
  allocate(nsteph(ntimemax,numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(ncdate) ) &
  allocate(ncdate(ntimemax,numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(ncsec) ) &
  allocate(ncsec(ntimemax,numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(time) ) &
  allocate(time(ntimemax), stat=ierr)
  if ( ierr /= 0 ) stop 78

  if ( .not. allocated(lon) ) &
  allocate(lon(nlon(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(lat) ) &
  allocate(lat(nlat(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(lev) ) &
  allocate(lev(nlev(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(ilev) ) &
  allocate(ilev(nlevp(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(gw) ) &
  allocate(gw(nlat(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(hyai) ) &
  allocate(hyai(nlevp(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(hybi) ) &
  allocate(hybi(nlevp(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(hyam) ) &
  allocate(hyam(nlev(1),numcases), stat=ierr)
  if ( ierr /= 0 ) stop 78
  if ( .not. allocated(hybm) ) &
  allocate(hybm(nlev(1),numcases), stat=ierr)
  if (twocases .and. .not. allocated(time2) ) allocate (time2(ntime(2)), stat=ierr)
  if ( ierr /= 0 ) stop 78
!
! Read in header variables
!
  do n=1,numcases
    call wrap_get_var_realx (ncid(n), p0id(n), p0(n))

#ifndef SUN
!
! Sun dies with seg-fault on following statement, don't do this if
! on Sun as a work-around. Really need a better long-term solution.
!
    ret(1) = nf_get_var_text (ncid(n), date_writtenid(n), date_written(:,n))
    if (ret(1).ne.NF_NOERR) then
      do nn=1,ntimemax
        date_written(nn,n) = 'UNKNOWN'
      end do
    end if
#endif

    ret(2) = nf_get_var_text (ncid(n), time_writtenid(n), time_written(:,n))
    if (ret(2).ne.NF_NOERR) then
      do nn=1,ntimemax
        time_written(nn,n) = 'UNKNOWN'
      end do
    end if

    ret(1) = nf_get_var_int (ncid(n), nstephid(n), nsteph(:,n))
    if (ret(1).ne.NF_NOERR) then
      do nn=1,ntimemax
        nsteph(nn,n) = -9999
      end do
    end if

    call wrap_get_var_int (ncid(n), ncdateid(n), ncdate(:,n))
    call wrap_get_var_int (ncid(n), ncsecid(n),  ncsec(:,n))
    call wrap_get_var_int (ncid(n), ntrmid(n), ntrm(n))
    call wrap_get_var_int (ncid(n), ntrnid(n), ntrn(n))
    call wrap_get_var_int (ncid(n), ntrkid(n), ntrk(n))
    call wrap_get_var_int (ncid(n), ndbaseid(n), ndbase(n))
    call wrap_get_var_int (ncid(n), nsbaseid(n), nsbase(n))
    call wrap_get_var_int (ncid(n), nbdateid(n), nbdate(n))
    call wrap_get_var_int (ncid(n), nbsecid(n), nbsec(n))
    call wrap_get_var_int (ncid(n), mdtid(n), mdt(n))

    call wrap_get_var_realx (ncid(n), lonid(n), lon(:,n))
    call wrap_get_var_realx (ncid(n), latid(n), lat(:,n))
    call wrap_get_var_realx (ncid(n), levid(n), lev(:,n))
    call wrap_get_var_realx (ncid(n), ilevid(n), ilev(:,n))
    call wrap_get_var_realx (ncid(n), gwid(n), gw(:,n))
    call wrap_get_var_realx (ncid(n), hyaiid(n), hyai(:,n))
    call wrap_get_var_realx (ncid(n), hybiid(n), hybi(:,n))
    call wrap_get_var_realx (ncid(n), hyamid(n), hyam(:,n))
    call wrap_get_var_realx (ncid(n), hybmid(n), hybm(:,n))
  end do

  call wrap_get_var_realx (ncid(1), timeid(1), time)
  if (twocases) call wrap_get_var_realx (ncid(2), timeid(2), time2)
  
  if (iprs.gt.0) then
    if (ipre.gt.nlon(1)) then
      write(6,*)'Bad ipre'
      write(6,*)'ipre,nlon=  ', ipre, nlon(1)
      stop 99
    end if

    if (iprs.lt.1) then
      write(6,*)'Bad iprs'
      write(6,*)'iprs=  ', iprs
      stop 99
    end if
  end if

  if (twocases) then
!
! Print header diffs
!
    call prhddiff (nlev(1), nlevp(1), nlat(1))
    if (.not.constps) then
      ret(1) = nf_inq_varid (ncid(1), 'PS', psid(1))
      ret(2) = nf_inq_varid (ncid(2), 'PS', psid(2))

      if (ret(1).ne.NF_NOERR .or. ret(2).ne.NF_NOERR) then
        write(6,*) 'WARNING: PS not found on both tapes: assuming constant'
        psid(:) = -1
      end if
    end if
  end if
!
! Print 1-time header
!
  write(6,'(/,a)')' SUMMARY OF FIELD DIFFERENCES:'
  write(6,802)'FIELD','# DIFFS','# POSS','MAX','MIN','DIFFMAX','VALUES','RDIFMAX','VALUES'
!
! Call the routine which loops through the time variables
!
  call cpr (nlon(1), nlat(1), nlev(1), nlevp(1), numcases, lnorm, verbose)
802 format(1x ,a5,2x,a7,1x,a6,3x,a3,19x,a3,20x,a7,1x,a7,16x,a7,1x,a6)
end program cprnc

subroutine usage_exit (arg)
  implicit none
  character*(*) arg
  
  if (arg.ne.' ') write (6,*) arg
  write (6,*) 'Usage: cprnc [-c] [-s nsteps] [-e nstepe] [-m] ', &
              'file1 file2'
  stop 999
end subroutine usage_exit
