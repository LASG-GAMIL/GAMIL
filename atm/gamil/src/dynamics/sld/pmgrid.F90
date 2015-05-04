#include <misc.h>
#include <params.h>

module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize grid point resolution parameters
! 
! Author: 
! 
!-----------------------------------------------------------------------

! Grid point resolution parameters

  integer, parameter :: plon       = PLON  ! number of longitudes
  integer, parameter :: plev       = PLEV  ! number of vertical levels
  integer, parameter :: plat       = PLAT  ! number of latitudes
  integer, parameter :: plevp      = plev + 1 ! plev + 1
  integer, parameter :: nxpt       = 1     ! no.of pts outside active domain of interpolnt
  integer, parameter :: jintmx     = 2     ! # of extra lats in polar region
  integer, parameter :: plond      = plon + 1 + 2*nxpt ! slt extended domain longitude
  integer, parameter :: platd      = plat + 2*nxpt + 2*jintmx ! slt extended domain lat.
  integer, parameter :: i1         = 1 + nxpt          ! model starting longitude index
  integer, parameter :: j1         = 1 + nxpt + jintmx ! model starting latitude index
  integer, parameter :: numbnd     = nxpt + jintmx ! no.of lats passed N/S of forecast lat
  integer, parameter :: plnlv      = plon*plev     ! Length of multilevel field slice
  integer, parameter :: plndlv     = plond*plev    ! Length of multilevel 3-d field slice

  integer iam
  integer begirow    ! beg. index for lat pairs owned by a proc
  integer endirow    ! end. index for lat pairs owned by a proc
  integer beglat     ! beg. index for lats owned by a given proc
  integer endlat     ! end. index for lats owned by a given proc
  integer numlats    ! number of lats owned by a given proc
  integer numlatsex  ! number of ext lats owned by a given proc
  integer beglatex   ! extended grid beglat
  integer endlatex   ! extended grid endlat 
  logical masterproc 
  logical :: dyngrid_set = .false. ! flag indicates dynamics grid has been set
#if ( ! defined SPMD )
  parameter (iam        = 0)
  parameter (begirow    = 1)
  parameter (endirow    = plat/2)
  parameter (beglat     = 1)
  parameter (endlat     = plat)
  parameter (numlats    = plat)
  parameter (numlatsex  = platd)
  parameter (beglatex   = 1)
  parameter (endlatex   = platd)
  parameter (masterproc = .true.)
#endif
end module pmgrid


