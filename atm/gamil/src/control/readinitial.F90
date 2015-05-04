#include <misc.h>
#include <params.h>

! ***************************************************************************** !
! Purpose: 
!
!   Ensure that requisite netcdf variables are on the initial dataset.
!   Set base day and date info using the "current" values from it.
! 
! Method:
!
!   Issue proper netcdf wrapper calls.  Broadcast to slaves if SPMD
! 
! Author:
!
!   CCM Core Group
!
! Development records:
!
!   1. old unclarified records
!      WAN Hui, 2003/04/30, 2003/07/10, 2003/10/20, 2003/10/27, 2003/11/22
!
! ***************************************************************************** !

subroutine readinitial(ncid)

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use rgrid
    use time_manager, only: ic_ymd, ic_tod
#if ( defined SPMD )
    use mpishorthand
#endif

    implicit none

#include <comctl.h>
#include <comhyb.h>

    include 'netcdf.inc'

    integer, intent(in) :: ncid

    integer :: lonid            !------------------------------------------------------
    integer :: levid            ! 
    integer :: latid            ! 
    integer :: ncdateid         ! Netcdf variable and dimension ids for variable of that
    integer :: ncsecid          ! name with "id" tacked on to the end
    integer :: hyaiid           ! 
    integer :: hybiid           ! 
    integer :: hyamid           ! 
    integer :: hybmid           ! 
    integer :: ilev             ! 
    integer :: ilevid           ! 
    integer :: rlonid           ! 
    integer :: nlonid           ! 

    integer :: phisid
    real(r8):: phis_tmp(plond,plat)  !!(wh 2003.10.20)

    integer strt2d(3)           ! start lon, lat, time indices for netcdf 2-d !!(wh10.20)
    data strt2d/3*1/            ! only index 2 will ever change               !!(wh10.20)
    integer cnt2d(3)            ! lon, lat, time counts for netcdf 2-d        !!(wh10.20)
    data cnt2d/plon,1,1/        ! 2-d arrs: Always grab only a "plon" slice   !!(wh10.20)

    integer :: i, j

    integer :: mlon             ! longitude dimension length from dataset
    integer :: mlev             ! level dimension length from dataset
    integer :: morec            ! latitude dimension length from dataset

#ifdef STAGGERED
    integer :: slonid           ! staggered longitude dimension length from dataset
    integer :: slatid           ! staggered latitude dimension length from dataset
#endif

    if (masterproc) then
        !
        ! Get and check dimension/date info
        !
        call wrap_inq_dimid(ncid, 'lon' , lonid)
        call wrap_inq_dimid(ncid, 'lev' , levid)
        call wrap_inq_dimid(ncid, 'ilev', ilevid)
        call wrap_inq_dimid(ncid, 'lat' , latid)

#ifdef STAGGERED
        call wrap_inq_dimid(ncid, 'slon' , slonid)
        call wrap_inq_dimid(ncid, 'slat' , slatid)
#endif
        call wrap_inq_varid(ncid, 'date'   , ncdateid)
        call wrap_inq_varid(ncid, 'datesec', ncsecid)
        call wrap_inq_varid(ncid, 'hyai'   , hyaiid)
        call wrap_inq_varid(ncid, 'hybi'   , hybiid)
        call wrap_inq_varid(ncid, 'hyam'   , hyamid)
        call wrap_inq_varid(ncid, 'hybm'   , hybmid)
        call wrap_inq_varid(ncid, 'PHIS'   , phisid)

        call wrap_inq_dimlen(ncid, lonid , mlon)
        call wrap_inq_dimlen(ncid, levid , mlev)
        call wrap_inq_dimlen(ncid, ilevid, ilev)
        call wrap_inq_dimlen(ncid, latid , morec)
        !
        ! Check for reduced grid info on initial dataset.  If not present, define
        ! variables for full grid
        !
        if (nf_inq_varid(ncid, 'nlon', nlonid) == nf_noerr) then
            call wrap_get_var_int(ncid, nlonid, nlon)
        else
            nlon(:) = plon
        end if

        if (mlev /= plev .or. mlon /= plon .or. morec /= plat) then
            write(6, "('Error: readinitial: model parameters do not match initial dataset parameters')")
            write(6, "('  Model Parameters:   plev = ', I3, ' plon = ', I3, ' plat  = ', I3)") plev, plon, plat
            write(6, "('  Dataset Parameters: mlev = ', I3, ' mlon = ', I3, ' morec = ', I3)") mlev, mlon, morec
            call endrun
        end if

        call wrap_get_var_int(ncid, ncdateid, ic_ymd)
        call wrap_get_var_int(ncid, ncsecid , ic_tod)

        call wrap_get_var_realx(ncid, hyamid, hyam)
        call wrap_get_var_realx(ncid, hybmid, hybm)
        call wrap_get_var_realx(ncid, hyaiid, hyai)
        call wrap_get_var_realx(ncid, hybiid, hybi)

        !!(wh 2003.10.20)

        do j = 1, plat
            strt2d(2) = j
            call wrap_get_vara_realx(ncid, phisid, strt2d, cnt2d, phis_tmp(1,j))
        end do
        !!(wh)

        !! (wh 2003.10.20)
        !!      do j=1,plat
        !!        do i=1,plon
        !!            ghs(i,j) = phis_tmp(i,plat+1-j)
        !!        enddo
        !!      enddo
        !!
        !!      do j=1,plat
        !!         ghs(plond-1,j) = ghs(1,j)
        !!         ghs(plond,  j) = ghs(2,j)
        !!      enddo
        !! (wh)

    end if

    !!(wh 2003.10.20)

    call copy_phis_ghs(phis_tmp)

#if ( defined SPMD )
    call mpibcast (ic_ymd,  1,    mpiint, 0, mpicom)
    call mpibcast (ic_tod,  1,    mpiint, 0, mpicom)
    call mpibcast (nlon,    plat, mpiint, 0, mpicom)
    !! call mpibcast (wnummax, plat, mpiint, 0, mpicom)

    call mpibcast (hyam  ,plev ,mpir8,  0, mpicom)
    call mpibcast (hybm  ,plev ,mpir8,  0, mpicom)
    call mpibcast (hyai  ,plevp,mpir8,  0, mpicom)
    call mpibcast (hybi  ,plevp,mpir8,  0, mpicom)
#endif

    return
end subroutine readinitial

!!#############################
!!(wanhui 2003.10.20)
!!(wanhui 2003.10.27)

subroutine copy_phis_ghs( phis_tmp )

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use prognostics,  only: phis
    use comfm1,       only: ghs
#if ( defined SPMD )
    use mpishorthand
    use spmd_dyn,     only: npes, compute_gsfactors
#endif

    implicit none

    real(r8), intent(in) :: phis_tmp(plond,plat)

    integer i, jdyn, jcam, begj, endj

    character(50) filename

#if ( defined SPMD )
    integer numperlat         ! number of values per latitude band
    integer numsend(0:npes-1) ! number of items to be sent
    integer numrecv           ! number of items to be received
    integer displs(0:npes-1)  ! displacement array
#endif
    !-----------------------------------------------------------------
    !! scatter phis_tmp to phis
    !!----------------------------------------------------------------

#if ( defined SPMD )

    numperlat = plond
    call compute_gsfactors(numperlat, numrecv, numsend, displs)

    call mpiscatterv (phis_tmp, numsend, displs, mpir8, &
        phis(1,beglat), numrecv, mpir8, 0, mpicom)

#else
    phis(:plon,:) = phis_tmp(:plon,:)
#endif

    !!----------------------------------------------------------------
    !! copy phis to ghs
    !!----------------------------------------------------------------

    begj = beglatexdyn + numbnd
    endj = endlatexdyn - numbnd
    do jdyn = begj,endj
        jcam = beglatexdyn + endlatex - jdyn
        do i=1,plon
            ghs(i,jdyn) = phis(i,jcam )
        enddo
        ghs(plond-1,jdyn) = ghs(1,jdyn)
        ghs(plond,  jdyn) = ghs(2,jdyn)
    enddo

    !- check -------------------------------------------------------
    !
    !#if (defined SPMD)
    !   
    !       write(filename,12) 'ghs-p-',iam,'.out'
    !12     format(a6,i1,a4)
    !       open(10,file=trim(filename))
    !#else
    !       open(10,file='ghs-s.out')
    !#endif
    !
    !        write(10,*) '----------------- ghs -----------------'
    !        do jdyn=begj,endj
    !           write(10,11) jdyn,(ghs(i,jdyn),i=1,3)
    !        enddo
    !11      format(1x,i3,3e30.20)
    !       close(10)
    !!       stop
    !--------------------------------------------------------------

    return
end subroutine copy_phis_ghs
