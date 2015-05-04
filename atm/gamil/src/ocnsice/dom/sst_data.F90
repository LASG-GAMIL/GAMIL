#include <misc.h>
#include <params.h>
!-----------------------------------------------------------------------
!
! MODULE: sst_data
!
! DESCRIPTION:
!
!   Module to handle dealing with the Sea-Surface Temperature
!   datasets. This module also figures out the location of
!   sea-ice from these datasets where it is assumed that
!   seawater at freezing or below is a flag for the existence of sea-ice.
!   SST datasets that are created for use with the stand-alone CCM should
!   take this into account and set grid-points where sea-ice fraction is
!   greater than 50% to -1.8C and ensure that other grid points where sea-ice
!   is less than 50% have SST's greater than -1.8C.
!
! Public interfaces:
!
!   sstini -- Initialization and reading of dataset.
!   sstint -- Interpolate dataset SST to current time.
!   sstan --- Apply the interpolated SST to the model state.
!
!-----------------------------------------------------------------------

module sst_data

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid,       only: plon, plat, masterproc
    use ppgrid,       only: pcols, begchunk, endchunk
    use phys_grid,    only: scatter_field_to_chunk, get_ncols_p
    use comsrf,       only: plevmx, icefrac
    use physconst,    only: tmelt
    use commap,       only: clat, clon

    implicit none

    public sst ! needed in ice_data

    public sstan    ! Set the surface temperature, oro, and sea-ice fraction
    public sstini   ! Initialization
    public sstint   ! Time interpolation of SST data

    private

    integer, parameter :: totsstsz = 2000
    real(r8), parameter :: daysperyear = 365.0

    real(r8), allocatable :: sstbdy(:,:,:) ! SST values on boundary dataset
    real(r8), allocatable :: sst(:,:) ! interpolated model sst values
    real(r8) cdaysstm   ! Calendar day for prv. month SST values read in
    real(r8) cdaysstp   ! Calendar day for nxt. month SST values read in

    integer nm, np            ! array indices for prv., nxt month sst data
    integer sstid             ! netcdf id for sst variable
    integer lonsiz            ! size of longitude dimension on sst dataset
    integer levsiz            ! size of level dimension on sst dataset
    integer latsiz            ! size of latitude dimension on sst dataset
    integer timesiz           ! size of time dimension on sst dataset
    integer np1               ! current forward time index of sst dataset
    integer date_sst(totsstsz)! date on sst dataset (YYYYMMDD)
    integer sec_sst(totsstsz) ! seconds of date on sst dataset (0-86399)

    real(r8), parameter :: tsice = -1.7999 ! freezing point of sea ice degrees C
                                           ! use this with global sst data

contains

    !-----------------------------------------------------------------------
    !
    ! SUBROUTINE: sstan
    !
    ! DESCRIPTION:
    !
    !   Update sea surface temperatures (sst's) and sea ice distribution
    !
    ! Method:
    !   Assume that the sst data exists in a two dimensional field
    !   encoded as follows:
    !     Land               values where oro field says so ("valid sst's
    !                        are provided globally, the model's land mask
    !                        determines whether the sst is used or not)
    !     Ocean without      values degrees celcius (greater than tsice)
    !      sea ice
    !     Ocean with         values less than tsice
    !      sea ice
    ! New sea ice has a constant 0.5 cm value for snow cover prescribed
    !
    ! Author: CCM1
    !
    !-----------------------------------------------------------------------

    subroutine sstan(lchnk, ncol, ocnfrac, ts)
        integer , intent(in) :: lchnk           ! chunk identifier
        integer , intent(in) :: ncol            ! number of atmospheric columns
        real(r8), intent(in) :: ocnfrac(pcols)  ! Surface type flag array
        real(r8), intent(inout) :: ts(pcols)    ! Surface temperature

        integer i ! Column index

        ! Open ocean
        do i = 1, ncol
            if (ocnfrac(i) > 0.) then
                ts(i) = sst(i,lchnk)+tmelt
            end if
        end do

        return
    end subroutine sstan

!-----------------------------------------------------------------------
!
! BOP
!
! !IROUTINE: sstini
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying sst boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method:
!
! Author: L.Bath
!
!-----------------------------------------------------------------------
!
! !INTERFACE
!
    subroutine sstini

        use rgrid,          only: nlon
        use error_messages, only: alloc_err, handle_ncerr
        use time_manager,   only: get_curr_date, get_curr_calday, &
                                  is_perpetual, get_perp_date
#if ( defined SPMD )
        use mpishorthand,   only: mpicom, mpiint, mpir8
#endif

        !---------------------------Common blocks-------------------------------
#include <comctl.h>
#include <comlun.h>
        !---------------------------Local variables-----------------------------
        integer dateid              ! netcdf id for date variable
        integer secid               ! netcdf id for seconds variable
        integer londimid            ! netcdf id for longitude variable
        integer latdimid            ! netcdf id for latitude variable
        integer lonid               ! netcdf id for longitude variable
        integer latid               ! netcdf id for latitude variable
        integer timeid              ! netcdf id for time variable
        integer nlonid              ! netcdf id for nlon variable (rgrid)
        integer cnt3(3)             ! array of counts for each dimension
        integer strt3(3)            ! array of starting indices
        integer n                   ! indices
        integer nlon_sst(plat)      ! number of lons per lat on bdy dataset
        integer j                   ! latitude index
        integer istat               ! error return
        integer yr, mon, day        ! components of a date
        integer ncdate              ! current date in integer format [yyyymmdd]
        integer ncsec               ! current time of day [seconds]
        real(r8) calday             ! calendar day (includes yr if no cycling)
        real(r8) caldayloc          ! calendar day (includes yr if no cycling)
        real(r8) xvar(plon,plat,2)  ! work space

        ! initialize time indices
        nm = 1
        np = 2

        ! allocate space for data.
        allocate(sst(pcols,begchunk:endchunk), stat=istat)
        call alloc_err(istat, 'sstini', 'sst', pcols*(endchunk-begchunk+1))

        if(aqua_planet) return

        allocate(sstbdy(pcols,begchunk:endchunk,2), stat=istat)
        call alloc_err(istat, 'sstini', 'sstbdy', pcols*(endchunk-begchunk+1)*2)

        ! SPMD: Master does all the work.
        if (masterproc) then
            ! use year information only if not cycling sst dataset
            calday = get_curr_calday()
            if (is_perpetual()) then
                call get_perp_date(yr, mon, day, ncsec)
            else
                call get_curr_date(yr, mon, day, ncsec)
            end if
            ncdate = yr*10000+mon*100+day
            if (sstcyc) then
                caldayloc = calday
            else
                caldayloc = calday + yr*daysperyear
            end if
            ! get and check dimension info
            call wrap_inq_dimid(ncid_sst, 'lon', londimid)
            call wrap_inq_dimid(ncid_sst, 'time', timeid )
            call wrap_inq_dimid(ncid_sst, 'lat', latdimid)

            call wrap_inq_dimlen(ncid_sst, londimid, lonsiz)
            if (lonsiz /= plon) then
                write(6,*) 'SSTINI: lonsiz=',lonsiz,' must = plon=',plon
                call endrun
            end if
            call wrap_inq_dimlen(ncid_sst, latdimid, latsiz)
            if (latsiz /= plat) then
                write(6,*) 'SSTINI: latsiz=',latsiz,' must = plat=',plat
                call endrun
            end if
            call wrap_inq_dimlen(ncid_sst, timeid, timesiz)

            ! check to make sure space allocated for time variables is sufficient
            if (timesiz > totsstsz) then
                write(6, *) 'SSTINI:  Allocated space for sst data is insufficient.'
                write(6, *) 'Please increase parameter totsstsz to',timesiz,' and recompile.'
                call endrun
            end if

            ! check to ensure reduced or not grid of dataset matches that of model
            if (fullgrid) then
                call wrap_inq_varid(ncid_sst, 'lon', lonid)
            else
                call wrap_inq_varid(ncid_sst, 'nlon', nlonid)
                call wrap_get_var_int(ncid_sst, nlonid, nlon_sst)
                do j = 1, plat
                    if (nlon_sst(j) /= nlon(j)) then
                        write(6, *) 'SSTINI: model grid does not match dataset grid'
                        call endrun
                    end if
                end do
            end if

            call wrap_inq_varid(ncid_sst, 'date', dateid)
            call wrap_inq_varid(ncid_sst, 'datesec', secid)
            call wrap_inq_varid(ncid_sst, 'SST_cpl', sstid)
            call wrap_inq_varid(ncid_sst, 'lat', latid)

            ! retrieve entire date and sec variables.

            call wrap_get_var_int (ncid_sst,dateid,date_sst)
            call wrap_get_var_int (ncid_sst,secid,sec_sst)
            if (sstcyc) then
                if (timesiz < 12) then
                    write(6, *) 'SSTINI: ERROR'
                    write(6, *) 'When cycling sst, sst data set must have 12'
                    write(6, *) 'consecutive months of data starting with Jan'
                    write(6, *) 'Current dataset has only ',timesiz,' months'
                    call endrun
                end if
                do n = 1, 12
                    if (mod(date_sst(n), 10000)/100 /= n) then
                        write(6, *) 'SSTINI: ERROR'
                        write(6, *) 'When cycling sst, sst data set must have 12'
                        write(6, *) 'consecutive months of data starting with Jan'
                        write(6, *) 'Month ', n, ' of sst data set is out of order'
                        call endrun
                    end if
                end do
            end if

            strt3(1) = 1
            strt3(2) = 1
            strt3(3) = 1
            cnt3(1)  = lonsiz
            cnt3(2)  = latsiz
            cnt3(3)  = 1

            ! special code for interpolation between December and January
            if (sstcyc) then
                n = 12
                np1 = 1
                call bnddyi(date_sst(n), sec_sst(n), cdaysstm)
                call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)
                if (caldayloc <= cdaysstp .or. caldayloc > cdaysstm) then
                    strt3(3) = n
                    call wrap_get_vara_realx(ncid_sst, sstid, strt3, cnt3, xvar(1,1,nm))
                    strt3(3) = np1
                    call wrap_get_vara_realx(ncid_sst, sstid, strt3, cnt3, xvar(1,1,np))
                    goto 10
                end if
            end if

            ! normal interpolation between consecutive time slices.
            do n = 1, timesiz-1
                np1 = n+1
                call bnddyi(date_sst(n  ), sec_sst(n  ), cdaysstm)
                call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)
                if (.not.sstcyc) then
                    yr = date_sst(n)/10000
                    cdaysstm = cdaysstm + yr*daysperyear
                    yr = date_sst(np1)/10000
                    cdaysstp = cdaysstp + yr*daysperyear
                end if
                if (caldayloc > cdaysstm .and. caldayloc <= cdaysstp) then
                    strt3(3) = n
                    call wrap_get_vara_realx(ncid_sst, sstid, strt3, cnt3, xvar(1,1,nm))
                    strt3(3) = np1
                    call wrap_get_vara_realx(ncid_sst, sstid, strt3, cnt3, xvar(1,1,np))
                    goto 10
                end if
            end do
            write(6, *) 'SSTINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
            call endrun
10          continue
            write(6, *) 'SSTINI: Read sst data for dates ', date_sst(n), sec_sst(n), &
                ' and ', date_sst(np1), sec_sst(np1)
#if (defined SPMD )
            call mpibcast(timesiz, 1, mpiint, 0, mpicom )
            call mpibcast(date_sst, totsstsz, mpiint, 0, mpicom )
            call mpibcast(sec_sst, totsstsz, mpiint, 0, mpicom )
            call mpibcast(cdaysstm, 1, mpir8, 0, mpicom )
            call mpibcast(cdaysstp, 1, mpir8, 0, mpicom )
            call mpibcast(np1, 1, mpiint, 0, mpicom )
        else
            call mpibcast(timesiz, 1, mpiint, 0, mpicom )
            call mpibcast(date_sst, totsstsz, mpiint, 0, mpicom )
            call mpibcast(sec_sst, totsstsz, mpiint, 0, mpicom )
            call mpibcast(cdaysstm, 1, mpir8, 0, mpicom )
            call mpibcast(cdaysstp, 1, mpir8, 0, mpicom )
            call mpibcast(np1, 1, mpiint, 0, mpicom )
#endif
        end if

        call scatter_field_to_chunk(1, 1, 2, plon, xvar, sstbdy)

        return
    end subroutine sstini

    !-----------------------------------------------------------------------
    !
    ! SUBROUTINE: sstint
    !
    ! DESCRIPTION:
    !
    ! if "aqua_planet", specify SST's analytically (Jerry Olson).
    ! Otherwise, time interpolate SST's to current time, reading in new
    ! monthly data if necessary.
    !
    ! Method:
    !
    ! Author: L.Bath
    !
    !-----------------------------------------------------------------------

    subroutine sstint

        use rgrid,          only: nlon
        use comsrf,         only: ocnfrac
        use time_manager,   only: get_curr_date, get_curr_calday, &
                                  is_perpetual, get_perp_date

        !---------------------------Common blocks-------------------------------
#include <comctl.h>
#include <comlun.h>
!---------------------------Local variables-----------------------------
        integer cnt3(3)        ! array of counts for each dimension
        integer strt3(3)       ! array of starting indices
        integer i, j, lchnk    ! indices
        integer ncol           ! number of columns in current chunk
        integer ntmp           ! temporary
        real(r8) fact1, fact2  ! time interpolation factors
        integer yr, mon, day   ! components of a date
        integer ncdate         ! current date in integer format [yyyymmdd]
        integer ncsec          ! current time of day [seconds]
        real(r8) calday        ! current calendar day
        real(r8) caldayloc     ! calendar day (includes yr if no cycling)
        real(r8) deltat        ! time (days) between interpolating sst data

        ! aqua planet variables
        real(r8) pi            ! 3.14159...
        real(r8) pio180        ! pi/180.
        real(r8) tmp           ! temporary
        real(r8) tmp1          ! temporary
        real(r8) t0_max        ! max reference temperature
        real(r8) t0_min        ! min reference temperature
        real(r8) t0_max6       ! max asymmetric reference temperature for option 6
        real(r8) t0_max7       ! max asymmetric reference temperature for option 7
        real(r8) maxlat        ! cutoff latitude poleward of which SST = 0 deg C
        real(r8) shift         ! number of degrees peak SST is shifted off equator
        real(r8) shift9        ! number of degrees peak SST is shifted off equator for opt. 9
        real(r8) shift10       ! number of degrees peak SST is shifted off equator for opt. 10
        real(r8) latcen        ! center of asymmetric SST forcing
        real(r8) latrad6       ! radius of asymmetric SST forcing for option 6
        real(r8) latrad8       ! radius of asymmetric SST forcing for option 8
        real(r8) loncen        ! center of asymmetric SST forcing
        real(r8) lonrad        ! radius of asymmetric SST forcing
        real(r8) xvar(plon,plat,2)    ! work space
        integer  sst_option    ! option of analytical SST algorithm

        ! SPMD: Master does all the work.  Sends needed info to slaves

        if (aqua_planet) then

            if (masterproc) then

                sst_option = 1
                pi         = 4.*atan(1.)
                pio180     = pi/180.
                if(sst_option < 1 .or. sst_option > 10) then
                    write(6,*) 'ERROR SSTINT:  sst_option must be between 1 and 10'
                call endrun
            end if

            ! parameters for zonally symmetric experiments
            t0_max     = 27.
            t0_min     = 0.
            maxlat     = 60.
            shift      = 5.
            shift9     = 10.
            shift10    = 15.

            ! parameters for zonally asymmetric experiments
            t0_max6    = 1.
            t0_max7    = 3.
            latcen     = 0.
            loncen     = 90.
            latrad6    = 15.
            latrad8    = 30.
            lonrad     = 30.

            maxlat     = maxlat *pio180
            shift      = shift  *pio180
            shift9     = shift9 *pio180
            shift10    = shift10*pio180
            latcen     = latcen *pio180
            loncen     = loncen *pio180
            latrad6    = latrad6*pio180
            latrad8    = latrad8*pio180
            lonrad     = lonrad *pio180

            if(sst_option == 1 .or. sst_option == 6 .or. &
                sst_option == 7 .or. sst_option == 8     ) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                else
                    tmp = clat(j)*pi*0.5/maxlat
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if

            if(sst_option == 2) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                else
                    tmp = clat(j)*pi*0.5/maxlat
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp*tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if

            if(sst_option == 3) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                else
                    tmp = clat(j)*pi*0.5/maxlat
                    tmp = sin(tmp)
                    tmp = (2. - tmp*tmp*tmp*tmp - tmp*tmp)*0.5
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if

            if(sst_option == 4) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                else
                    tmp  = (maxlat - abs(clat(j)))/maxlat
                    tmp1 = 1. - tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_max*tmp + t0_min*tmp1
                    end do
                end if
            end do
            end if

            if(sst_option == 5) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                elseif(clat(j) > shift) then
                    tmp = (clat(j)-shift)*pi*0.5/(maxlat-shift)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                else
                    tmp = (clat(j)-shift)*pi*0.5/(maxlat+shift)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if

            if(sst_option == 6) then
            do j = 1,plat
                if(abs(clat(j)-latcen) <= latrad6) then
                    tmp1 = (clat(j)-latcen)*pi*0.5/latrad6
                    tmp1 = cos(tmp1)
                    tmp1 = tmp1*tmp1
                    do i=1,nlon(j)
                        if(abs(clon(i,j)-loncen) <= lonrad) then
                        tmp = (clon(i,j)-loncen)*pi*0.5/lonrad
                        tmp = cos(tmp)
                        tmp = tmp*tmp
                        xvar(i,j,1) = xvar(i,j,1) + t0_max6*tmp*tmp1
                        end if
                    end do
                end if
            end do
            end if

            if(sst_option == 7) then
            do j = 1,plat
                if(abs(clat(j)-latcen) <= latrad6) then
                    tmp1 = (clat(j)-latcen)*pi*0.5/latrad6
                    tmp1 = cos(tmp1)
                    tmp1 = tmp1*tmp1
                    do i=1,nlon(j)
                        if(abs(clon(i,j)-loncen) <= lonrad) then
                        tmp = (clon(i,j)-loncen)*pi*0.5/lonrad
                        tmp = cos(tmp)
                        tmp = tmp*tmp
                        xvar(i,j,1) = xvar(i,j,1) + t0_max7*tmp*tmp1
                        end if
                    end do
                end if
            end do
            end if

            if(sst_option == 8) then
            do j = 1,plat
                if(abs(clat(j)-latcen) <= latrad8) then
                    tmp1 = (clat(j)-latcen)*pi*0.5/latrad8
                    tmp1 = cos(tmp1)
                    tmp1 = tmp1*tmp1
                    do i=1,nlon(j)
                        tmp = cos(clon(i,j)-loncen)
                        xvar(i,j,1) = xvar(i,j,1) + t0_max7*tmp*tmp1
                    end do
                end if
            end do
            end if

            if(sst_option == 9) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                elseif(clat(j) > shift9) then
                    tmp = (clat(j)-shift9)*pi*0.5/(maxlat-shift9)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                else
                    tmp = (clat(j)-shift9)*pi*0.5/(maxlat+shift9)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if

            if(sst_option == 10) then
            do j = 1,plat
                if(abs(clat(j)) > maxlat) then
                    do i=1,nlon(j)
                        xvar(i,j,1) = t0_min
                    end do
                elseif(clat(j) > shift10) then
                    tmp = (clat(j)-shift10)*pi*0.5/(maxlat-shift10)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                else
                    tmp = (clat(j)-shift10)*pi*0.5/(maxlat+shift10)
                    tmp = sin(tmp)
                    tmp = 1. - tmp*tmp
                    do i=1,nlon(j)
                        xvar(i,j,1) = tmp*(t0_max - t0_min) + t0_min
                    end do
                end if
            end do
            end if
            end if

            call scatter_field_to_chunk(1,1,1,plon,xvar(1,1,1),sst(1,begchunk))

        else

            ! use year information only if a multiyear dataset
            calday = get_curr_calday()
            if (is_perpetual()) then
                call get_perp_date(yr, mon, day, ncsec)
            else
                call get_curr_date(yr, mon, day, ncsec)
            end if
            ncdate = yr*10000+mon*100+day
            if (sstcyc) then
                caldayloc = calday
            else
                caldayloc = calday+yr*daysperyear
            end if

            if (masterproc) then
                strt3(1) = 1
                strt3(2) = 1
                strt3(3) = 1
                cnt3(1)  = lonsiz
                cnt3(2)  = latsiz
                cnt3(3)  = 1
            end if
!
! If model time is past current forward sst timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for sstcyc = .true.
! interpolation between December and January (np1==1).  Note that
! np1 is never 1 when sstcyc is .false.
!
            if (caldayloc > cdaysstp .and. .not. (np1 == 1 .and. caldayloc > cdaysstm)) then
                if (sstcyc) then
                    np1 = mod(np1, 12)+1
                else
                    np1 = np1+1
                end if
                if (np1 > timesiz) then
                    if (masterproc) then
                        write(6, *) 'SSTINT: Attempt to read past end of SST dataset'
                    end if
                    call endrun
                end if
                cdaysstm = cdaysstp
                call bnddyi(date_sst(np1), sec_sst(np1), cdaysstp)

                if (.not. sstcyc) then
                    yr = date_sst(np1)/10000
                    cdaysstp = cdaysstp+yr*daysperyear
                end if

                if (np1 == 1 .or. caldayloc <= cdaysstp) then
                    ntmp = nm
                    nm = np
                    np = ntmp
                    if (masterproc) then
                        strt3(3) = np1
                        call wrap_get_vara_realx(ncid_sst, sstid, strt3, cnt3, xvar(1,1,np))
                        write(6, *) 'SSTINT: Read sst for date (yyyymmdd) ',date_sst(np1), &
                            ' sec ',sec_sst(np1)
                    end if
                    call scatter_field_to_chunk(1, 1, 1, plon, xvar(1,1,np), sstbdy(1,begchunk,np))
                else
                    if (masterproc) then
                        write(6, *) 'SSTINT: Input sst for date',date_sst(np1), &
                            ' sec ',sec_sst(np1), 'does not exceed model date',ncdate,&
                            ' sec ',ncsec,' Stopping.'
                    end if
                    call endrun
                end if
            end if
!
! Time interpolation.  Account for December-January interpolation if
! cycling sst dataset.  Again note that np1 is never 1 when sstcyc is false
!
            if (np1 == 1) then                    ! Dec-Jan interpolation
                deltat = cdaysstp+daysperyear-cdaysstm
                if (caldayloc > cdaysstp) then      ! We're in December
                    fact1 = (cdaysstp+daysperyear-caldayloc)/deltat
                    fact2 = (caldayloc-cdaysstm)/deltat
                else                                ! We're in January
                    fact1 = (cdaysstp-caldayloc)/deltat
                    fact2 = (caldayloc+daysperyear-cdaysstm)/deltat
                end if
            else
                deltat = cdaysstp-cdaysstm
                fact1 = (cdaysstp-caldayloc)/deltat
                fact2 = (caldayloc-cdaysstm)/deltat
            end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
            if (abs(fact1+fact2-1.) > 1.e-6 .or. &
                fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
                fact2 > 1.000001 .or. fact2 < -1.e-6) then
                if (masterproc) then
                    write(6, *) 'SSTINT: Bad fact1 and/or fact2=', fact1, fact2
                end if
                call endrun
            end if

            do lchnk = begchunk, endchunk
                ncol = get_ncols_p(lchnk)
                do i = 1, ncol
                    if (ocnfrac(i,lchnk) >= 0.) then
                        sst(i,lchnk) = sstbdy(i,lchnk,nm)*fact1+sstbdy(i,lchnk,np)*fact2
                        ! bound the sst temp by the freezing point of sea water
                        sst(i,lchnk) = max(sst(i,lchnk), tsice)
                    end if
                end do
            end do
        end if

        return
    end subroutine sstint

end module sst_data

