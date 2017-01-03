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
#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none

#include <comctl.h>
#include <comhyb.h>

  include 'netcdf.inc'

  integer, intent(in) :: ncid

  integer lonid, latid, levid
  integer hyaiid, hybiid, hyamid, hybmid
  integer ilev, ilevid, pmtopid
  integer nlonid
  integer phisid
  real(r8) phis_tmp(plond,plat)

  integer strt2d(3)           ! start lon, lat, time indices for netcdf 2-d !!(wh10.20)
  data strt2d/3*1/            ! only index 2 will ever change               !!(wh10.20)
  integer cnt2d(3)            ! lon, lat, time counts for netcdf 2-d        !!(wh10.20)
  data cnt2d/plon,1,1/        ! 2-d arrs: Always grab only a "plon" slice   !!(wh10.20)

  integer i, j, k

  integer mlon             ! longitude dimension length from dataset
  integer mlev             ! level dimension length from dataset
  integer morec            ! latitude dimension length from dataset

#ifdef STAGGERED
  integer slonid           ! staggered longitude dimension length from dataset
  integer slatid           ! staggered latitude dimension length from dataset
#endif

  if (masterproc) then
    call wrap_inq_dimid(ncid, 'lon' , lonid)
    call wrap_inq_dimid(ncid, 'lev' , levid)
    call wrap_inq_dimid(ncid, 'ilev', ilevid)
    call wrap_inq_dimid(ncid, 'lat' , latid)

#ifdef STAGGERED
    call wrap_inq_dimid(ncid, 'slon' , slonid)
    call wrap_inq_dimid(ncid, 'slat' , slatid)
#endif
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
      write(6, "('[Error]: readinitial: model parameters do not match initial dataset parameters')")
      write(6, "('  Model Parameters:   plev = ', I3, ' plon = ', I3, ' plat  = ', I3)") plev, plon, plat
      write(6, "('  Dataset Parameters: mlev = ', I3, ' mlon = ', I3, ' morec = ', I3)") mlev, mlon, morec
      call endrun
    end if

    call wrap_get_var_realx(ncid, hyamid, hyam)
    call wrap_get_var_realx(ncid, hybmid, hybm)
    call wrap_get_var_realx(ncid, hyaiid, hyai)
    call wrap_get_var_realx(ncid, hybiid, hybi)

    ! Read sigma levels.
    call wrap_inq_varid(ncid, "pmtop", pmtopid)
    call wrap_get_var_realx(ncid, pmtopid, pmtop)
    call wrap_get_var_realx(ncid, levid, sigl)
    call wrap_get_var_realx(ncid, ilevid, sig)

    do k = 1, plev
      dsig(k) = sig(k+1) - sig(k)      
    end do

    do j = 1, plat
      strt2d(2) = j
      call wrap_get_vara_realx(ncid, phisid, strt2d, cnt2d, phis_tmp(1,j))
    end do

  end if

  call copy_phis_ghs(phis_tmp)

#if ( defined SPMD )
  call mpibcast(nlon,  plat,  mpiint, 0, mpicom)
  call mpibcast(hyam,  plev,  mpir8,  0, mpicom)
  call mpibcast(hybm,  plev,  mpir8,  0, mpicom)
  call mpibcast(hyai,  plevp, mpir8,  0, mpicom)
  call mpibcast(hybi,  plevp, mpir8,  0, mpicom)
  call mpibcast(sigl,  plev,  mpir8,  0, mpicom)
  call mpibcast(sig,   plevp, mpir8,  0, mpicom)
  call mpibcast(dsig,  plev,  mpir8,  0, mpicom)
  call mpibcast(pmtop, 1,     mpir8,  0, mpicom)
#endif

end subroutine readinitial

subroutine copy_phis_ghs( phis_tmp )

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use prognostics,  only: phis
  use comfm1,       only: ghs
  use mpi_gamil
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

  ! scatter phis_tmp to phis
  call gamil_scatter_2D_array_phys(phis_tmp, beglonex, endlonex, beglat, endlat, phis(beglonex,beglat)) 

  ! copy phis to ghs
  begj = beglatexdyn + numbnd
  endj = endlatexdyn - numbnd
  do jdyn = begj, endj
    jcam = beglatexdyn + endlatex - jdyn
    do i = beglonex, endlonex
      ghs(i,jdyn) = phis(i,jcam )
    end do
  end do

  call gamil_arrays_comm(COMM_ROTATE_LEFT,2,ghs(:,beglatexdyn))

  return
end subroutine copy_phis_ghs
