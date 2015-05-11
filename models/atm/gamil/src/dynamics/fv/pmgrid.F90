#include <misc.h>
#include <params.h>

module pmgrid
!BOP
!
! !MODULE: pmgrid --- Initialize grid point resolution parameters
!
! !USES:
   use decompmodule, only : decomptype

   implicit none

!
! !DESCRIPTION:  Initialize grid point resolution parameters
! 
! !REVISION HISTORY:
!
!   ??.??.??   ??????     Creation
!   01.03.26   Sawyer     Added ProTeX documentation
!   01.06.27   Mirin      Added 2D decomposition material
!   01.10.16   Sawyer     Added Y-at-Z subdomain decompositions
!   02.05.02   Sawyer     Added XY variable definitions for non-SPMD
!
! !PUBLIC DATA MEMBERS:
   integer, parameter :: plon   = PLON       ! number of longitudes
   integer, parameter :: plev   = PLEV       ! number of vertical levels
   integer, parameter :: plat   = PLAT       ! number of latitudes

   integer, parameter :: plevp  = plev+1     ! plev + 1
   integer, parameter :: plond  = plon       ! no ghost pts in LR so plond=plon
   integer, parameter :: platd  = plat       ! no ghost pts in LR so platd=plat
   integer, parameter :: numbnd = 0          ! no.of latitudes passed N and S of forecast lat
   integer, parameter :: plnlv = plon*plev   ! Length of multilevel field slice
   integer, parameter :: plndlv = plond*plev ! Length of multilevel 3-d field slice

   integer beglat     ! beg. index for latitudes owned by a given proc
   integer endlat     ! end. index for latitudes owned by a given proc
   integer numlats    ! number of latitudes owned by a given proc

   integer beglev     ! beg. index for levels owned by a given task
   integer endlev     ! end. index for levels owned by a given task
   integer endlevp1   ! end. index + 1 for levels owned by a given task
   integer endlevp    ! equals endlev, except in last subdomain where equals endlevp1

   integer myid_y     ! subdomain index (0-based) in latitude (y)
   integer myid_z     ! subdomain index (0 based) in level (z)
   integer npr_y      ! number of subdomains in y
   integer npr_z      ! number of subdomains in z

! Install secondary xy decomposition - to be consistent with segmented decomposition
! This is temporary; will generalize at later date

   integer twod_decomp  ! 1 for multi-2D decompositions, 0 otherwise
#if defined ( TWOD_YZ )
   parameter(twod_decomp=1)
#else
   parameter(twod_decomp=0)
#endif

   integer myidxy_x     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   integer myidxy_y     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   integer nprxy_x      ! number of subdomains in x (second. decomp.)
   integer nprxy_y      ! number of subdomains in y (second. decomp.)
   integer beglonxy     ! beg. index for longitudes (second. decomp.)
   integer endlonxy     ! end. index for longitudes (second. decomp.)
   integer beglatxy     ! beg. index for latitudes (second. decomp.)
   integer endlatxy     ! end. index for latitudes (second. decomp.)

   integer iam
   logical masterproc ! Flag for (iam eq 0)
   logical :: dyngrid_set = .false. ! flag indicates dynamics grid has been set

   type(decomptype) :: strip2d, strip3dxyz, strip3dxzy, strip3dq3,         &
                       strip3dxyzp, strip3zaty, strip3dxzyp,               &
                       strip3yatz, strip3yatzp, strip3dq3old,              &
                       strip3kxyz, strip3kxzy, strip3kxyzp, strip3kxzyp,   &
                       strip3kq3, strip3zatypt, strip3zatyj2,              &
                       strip3zatypj1, strip3zatypj2, strip3zatyt4

!
#if ( ! defined SPMD )
   parameter (iam      = 0)
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (numlats  = plat)
   parameter (masterproc = .true.)
   parameter (beglev = 1)
   parameter (endlev = plev)
   parameter (endlevp1 = plev+1)
   parameter (endlevp = plev+1)
   parameter (myid_y = 0)
   parameter (myid_z = 0)
   parameter (npr_y = 1)
   parameter (npr_z = 1)
!
! These are needed to pass strict run-time error checking
!
   parameter (myidxy_x=0)     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   parameter (myidxy_y=0)     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   parameter (nprxy_x=1)      ! number of subdomains in x (second. decomp.)
   parameter (nprxy_y=1)      ! number of subdomains in y (second. decomp.)
   parameter (beglonxy=1)     ! beg. index for longitudes (second. decomp.)
   parameter (endlonxy=1)     ! end. index for longitudes (second. decomp.)
   parameter (beglatxy=1)     ! beg. index for latitudes (second. decomp.)
   parameter (endlatxy=1)     ! end. index for latitudes (second. decomp.)
#endif

! Staggered grid parameters
! splon and splat may eventually need to become new parameters
! in params.h to define the size of the staggered grid arrays. - gg

   integer, parameter :: splon = plon     ! Number of longitudes on the staggered grid
   integer, parameter :: splat = plat     ! Number of latitudes on the staggered grid

! Note: In reality, the staggered latitude array for Lin-Rood dynamics only
! uses PLAT-1 latitudes, the first one being ignored.  So ideally the line
! above should read:
!   parameter (splat = plat-1)
! to define the staggered latitude winds with the correct dimension.
! However, the assumption that the staggered latitude grid has one extra
! latitude (making it the same dimension as the non-staggered grid) is
! pervasive throughout the Lin-Rood dynamical core, necessitating the
! extra latitude.
!
! A temporary parameter for I/O purposes is defined below. When the Lin-Rood
! dynamical core has been changed to use a staggered-U wind arrays with the
! correct dimension (plat-1), replace all occurances of platm1 with splat,
! set splat to plat-1, and delete the following definition for platm1:

!  integer platm1
!  parameter (platm1 = plat-1)

!EOP
end module pmgrid
