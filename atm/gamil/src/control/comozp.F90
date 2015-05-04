#include <misc.h>
#include <params.h>

      module comozp
!----------------------------------------------------------------------- 
! 
! Purpose: Variables associated with ozone dataset
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid

      implicit none

      real(r8) cdayozm  ! dataset calendar day previous month
      real(r8) cdayozp  ! dataset calendar day next month
      real(r8) cplos    ! constant for ozone path length integral
      real(r8) cplol    ! constant for ozone path length integral

      integer nm        ! Array indices for previous month ozone data
      integer np        ! Array indices for next month ozone data
      integer oznid     ! netcdf id for ozone variable
      integer lonsiz    ! size of longitude dimension on ozone dataset
      integer levsiz    ! size of level dimension on ozone dataset
      integer latsiz    ! size of latitude dimension on ozone dataset
      integer timesiz   ! size of time dimension on ozone dataset
      integer np1       ! current forward time index of ozone dataset

      real(r8), allocatable :: ozmixm(:,:,:,:)  ! monthly mixing ratios
      real(r8), allocatable :: ozmix(:,:,:)     ! mixing ratio
      real(r8), allocatable :: pin(:)           ! ozone pressure level
      real(r8), allocatable :: ozlon(:)         ! Longitudes of bdy dataset
      real(r8), allocatable :: ozlat(:)         ! Latitudes of bdy dataset

      integer, allocatable :: date_oz(:)        ! Date on ozone dataset (YYYYMMDD)
      integer, allocatable :: sec_oz(:)         ! seconds of date (0-86399)

      end module comozp
