#include <misc.h>
#include <preproc.h>

module mvegFileMod

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

! Interpolation data

  integer  :: InterpMonths1 = -999   !saved month index
  real(r8) :: timwt(2)               !time weights for month 1 and month 2

! Monthly Vegetation

  real(r8), allocatable :: mlai1(:) !lai for interpolation (month 1)
  real(r8), allocatable :: mlai2(:) !lai for interpolation (month 2)
  real(r8), allocatable :: msai1(:) !sai for interpolation (month 1)
  real(r8), allocatable :: msai2(:) !sai for interpolation (month 2)
  real(r8), allocatable :: mhvt1(:) !top vegetation height for interpolation (month 1)
  real(r8), allocatable :: mhvt2(:) !top vegetation height for interpolation (month 2)
  real(r8), allocatable :: mhvb1(:) !bottom vegetation height for interpolation(month 1)
  real(r8), allocatable :: mhvb2(:) !bottom vegetation height for interpolation(month 2)

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine interpMonthlyVeg (fveg, kmo, kda)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Determine if 2 new months of data are to be read
! 
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: mvegFileMod.F90,v 1.3.6.5 2002/06/15 13:50:35 erik Exp $
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: fveg  !file with monthly vegetation data
    integer, intent(in) :: kmo             !month (1, ..., 12)
    integer, intent(in) :: kda             !day of month (1, ..., 31)
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------  
    real(r8):: t           !a fraction: kda/ndaypm
    integer :: it(2)       !month 1 and month 2 (step 1)
    integer :: months(2)   !months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
! -----------------------------------------------------------------

    t = (kda-0.5) / ndaypm(kmo)
    it(1) = t + 0.5
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5) - t
    timwt(2) = 1.-timwt(1)

    if (InterpMonths1 /= months(1)) then
       call readMonthlyVegetation (fveg, kmo, kda, months)
       InterpMonths1 = months(1)
    endif

  end subroutine interpMonthlyVeg

!=======================================================================

  subroutine readMonthlyVegetation (fveg, kmo, kda, months)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read monthly vegetation data for two consec. months
! 
! Method: 
! 
! Author: Sam Levis
! 
!-----------------------------------------------------------------------
! $Id: mvegFileMod.F90,v 1.3.6.5 2002/06/15 13:50:35 erik Exp $
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar  , only : lsmlon, lsmlat, maxpatch_pft, maxpatch 
    use clm_varmap  , only : landvec
    use clm_varsur  , only : landmask, numlon
    use fileutils   , only : getfil
    use spmdMod     , only : masterproc
#if (defined SPMD)
    use mpishorthand, only : mpir8, mpicom
#endif
    use time_manager, only : get_nstep
    implicit none

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: fveg    !file with monthly vegetation data
    integer, intent(in) :: kmo              !month (1, ..., 12)
    integer, intent(in) :: kda              !day of month (1, ..., 31)
    integer, intent(in) :: months(2)        !months to be interpolated (1 to 12)
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    character(len=256) :: locfn       !local file name

    integer :: iland                  !number of land points
    integer :: patch_index            !patch index 
    integer :: k                      !month do loop index
    integer :: i                      !longitude do loop index
    integer :: j                      !latitude do loop index
    integer :: m                      !patch do loop index
    integer :: ier                    !error code 

    integer :: ncid,dimid,varid       !input netCDF id's
    integer :: beg4d(4),len4d(4)      !netCDF variable edges
    integer :: ntim                   !number of input data time samples
    integer :: nlon_i                 !number of input data longitudes
    integer :: nlat_i                 !number of input data latitudes
    integer :: numpft_i               !number of input data pft types

    real(r8) :: laixy(lsmlon,lsmlat,maxpatch_pft) !lai read from input files
    real(r8) :: saixy(lsmlon,lsmlat,maxpatch_pft) !sai read from input files
    real(r8) :: hvtxy(lsmlon,lsmlat,maxpatch_pft) !top vegetation height
    real(r8) :: hvbxy(lsmlon,lsmlat,maxpatch_pft) !bottom vegetation height
! -----------------------------------------------------------------

! ----------------------------------------------------------------------
! Open monthly vegetation file
! Read data and convert from [lsmlon] x [lsmlat] grid to patch data 
! ----------------------------------------------------------------------

    if (masterproc) then

       write (6,*) 'Attempting to read monthly vegetation data .....'
       write (6,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda

       call getfil (fveg, locfn, 0)
       call wrap_open(locfn, 0, ncid)

       call wrap_inq_dimid  (ncid, 'lsmlon', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlon_i)
       if (nlon_i /= lsmlon) then
          write(6,*)'ReadMonthlyVegetation: parameter lsmlon= ', &
               lsmlon,'does not equal input nlat_i= ',nlon_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmlat', dimid)
       call wrap_inq_dimlen (ncid, dimid, nlat_i)
       if (nlat_i /= lsmlat) then
          write(6,*)'ReadMonthlyVegetation: parameter lsmlat= ', &
               lsmlat,'does not equal input nlat_i= ',nlat_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'lsmpft', dimid)
       call wrap_inq_dimlen (ncid, dimid, numpft_i)
       if (numpft_i /= maxpatch_pft) then
          write(6,*)'ReadMonthlyVegetation: parameter maxpatch_pft = ', &
               maxpatch_pft,'does not equal input numpft_i= ',numpft_i
          call endrun
       endif

       call wrap_inq_dimid  (ncid, 'time', dimid)
       call wrap_inq_dimlen (ncid, dimid, ntim)

       do k=1,2   !loop over months and read vegetated data

          beg4d(1) = 1         ; len4d(1) = nlon_i
          beg4d(2) = 1         ; len4d(2) = nlat_i
          beg4d(3) = 1         ; len4d(3) = numpft_i
          beg4d(4) = months(k) ; len4d(4) = 1

          call wrap_inq_varid (ncid, 'MONTHLY_LAI', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, laixy )

          call wrap_inq_varid (ncid, 'MONTHLY_SAI', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, saixy )

          call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_TOP', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, hvtxy)

          call wrap_inq_varid (ncid, 'MONTHLY_HEIGHT_BOT', varid)
          call wrap_get_vara_realx (ncid, varid, beg4d, len4d, hvbxy)

! store data directly in patch vectors processing vegetated 
! patches first and then setting non-vegetated patches to zero

          iland = 0
          do j = 1, lsmlat
             do i = 1, numlon(j)
                if (landmask(i,j) == 1) then
                   iland = iland + 1
                   do m = 1, maxpatch_pft
                      if (landvec%wtxy(iland,m) > 0.) then
                         patch_index =  landvec%patch(iland,m)
                         if (k == 1) then
                            mlai1(patch_index) = laixy(i,j,m)
                            msai1(patch_index) = saixy(i,j,m)
                            mhvt1(patch_index) = hvtxy(i,j,m)
                            mhvb1(patch_index) = hvbxy(i,j,m)
                         else !if (k == 2) then
                            mlai2(patch_index) = laixy(i,j,m)
                            msai2(patch_index) = saixy(i,j,m)
                            mhvt2(patch_index) = hvtxy(i,j,m)
                            mhvb2(patch_index) = hvbxy(i,j,m)
                         end if
                      end if
                   end do
                   do m = maxpatch_pft+1, maxpatch
                      if (landvec%wtxy(iland,m) > 0.) then
                         patch_index =  landvec%patch(iland,m)
                         if (k == 1) then
                            mlai1(patch_index) = 0.
                            msai1(patch_index) = 0.
                            mhvt1(patch_index) = 0.
                            mhvb1(patch_index) = 0.
                         else !if (k == 2) then
                            mlai2(patch_index) = 0.
                            msai2(patch_index) = 0.
                            mhvt2(patch_index) = 0.
                            mhvb2(patch_index) = 0.
                         end if
                      end if
                   end do
                end if
             end do
          end do

       end do   ! end of loop over months

       call wrap_close(ncid)

       write (6,*) 'Successfully read monthly vegetation data for'
       write (6,*) 'month ', months(1), ' and month ', months(2)
       write (6,*)

    endif ! end of if-masterproc if block


#if ( defined SPMD )
    call mpi_bcast (mlai1, size(mlai1), mpir8, 0, mpicom, ier)
    call mpi_bcast (mlai2, size(mlai2), mpir8, 0, mpicom, ier)
    call mpi_bcast (msai1, size(msai1), mpir8, 0, mpicom, ier)
    call mpi_bcast (msai2, size(msai2), mpir8, 0, mpicom, ier)
    call mpi_bcast (mhvt1, size(mhvt1), mpir8, 0, mpicom, ier)
    call mpi_bcast (mhvt2, size(mhvt2), mpir8, 0, mpicom, ier)
    call mpi_bcast (mhvb1, size(mhvb1), mpir8, 0, mpicom, ier)
    call mpi_bcast (mhvb2, size(mhvb2), mpir8, 0, mpicom, ier)
#endif

    return
  end subroutine readMonthlyVegetation


!=======================================================================

  subroutine monthveg_ini

    use shr_kind_mod, only: r8 => shr_kind_r8
    use infnan
    use clm_varmap, only : numpatch
    implicit none

! Dynamically allocate memory

    allocate (mlai1(numpatch))
    allocate (mlai2(numpatch))
    allocate (msai1(numpatch))
    allocate (msai2(numpatch))
    allocate (mhvt1(numpatch))
    allocate (mhvt2(numpatch))
    allocate (mhvb1(numpatch))
    allocate (mhvb2(numpatch))

! Set to infinity

    mlai1(:) = inf
    mlai2(:) = inf
    msai1(:) = inf
    msai2(:) = inf
    mhvt1(:) = inf
    mhvt2(:) = inf
    mhvb1(:) = inf
    mhvb2(:) = inf

  end subroutine monthveg_ini

!=======================================================================

end module mvegFileMod
