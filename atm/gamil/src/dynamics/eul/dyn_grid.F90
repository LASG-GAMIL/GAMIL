#include <misc.h>
module dyn_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of dynamics computational grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!
!      get_block_bounds_d  get first and last indices in global 
!                          block ordering
!      get_block_owner_d   get process "owning" given block
!      get_block_col_cnt_d get number of columns in given block
!      get_block_lvl_cnt_d get number of vertical levels in column
!      get_block_levels_d  get vertical levels in column
!
!      get_block_coord_cnt_d   get number of blocks containing
!                              data from a given vertical column
!      get_block_coord_d   get global block indices and local columns 
!                          index for given coordinates
!
!      get_xxx_d           get global indices or coordinates for a single
!                          column
!      where xxx is
!       lat                for global latitude index
!       lon                for global longitude index
!
! Author: John Drake and Patrick Worley
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plev

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore
   ! sxj--- add ptimelevels 2008-11-10

contains
!========================================================================
!
   integer function get_block_coord_cnt_d(glon,glat)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of blocks contain data for the vertical column
!          at location (glon,glat)
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: glon     ! global longitude index
   integer, intent(in) :: glat     ! global latitude index
!-----------------------------------------------------------------------
!  latitude slice block
   get_block_coord_cnt_d = 1

   return
   end function get_block_coord_cnt_d
!
!========================================================================
!
   subroutine get_block_coord_d(glon,glat,cnt,blockid,bcid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return global block index and local column index
!          for column at location (glon,glat)
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: glon     ! global longitude index
   integer, intent(in) :: glat     ! global latitude index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
!---------------------------Local workspace-----------------------------
!
    integer jb                     ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (cnt < 1) then
      write(6,*)'GET_BLOCK_COORD_D: arrays not large enough (', &
                          cnt,' < ',1,' ) '
      call endrun
   else
      blockid(1) = glat
      bcid(1)    = glon
      do jb=2,cnt
         blockid(jb) = -1
         bcid(jb)    = -1
      enddo
   endif
!
   return
   end subroutine get_block_coord_d
!
!========================================================================
!
   Subroutine get_block_bounds_d(block_first,block_last)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return first and last indices used in global block ordering
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plat

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
   block_first = 1
   block_last  = plat

   return
   end subroutine get_block_bounds_d

!
!========================================================================
!
   integer function get_block_col_cnt_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of lon/lat coordinates in indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use rgrid, only: nlon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_col_cnt_d = nlon(blockid)

   return
   end function get_block_col_cnt_d
!
!========================================================================
!
   integer function get_block_lvl_cnt_d(blockid,bcid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of levels in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_lvl_cnt_d = plev + 1

   return
   end function get_block_lvl_cnt_d
!
!========================================================================
!
   subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return level indices in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block
   integer, intent(in) :: lvlsiz   ! dimension of levels array

   integer, intent(out) :: levels(lvlsiz) ! levels indices for block

!---------------------------Local workspace-----------------------------
!
    integer k                      ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (lvlsiz < plev + 1) then
      write(6,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
                          lvlsiz,' < ',plev + 1,' ) '
      call endrun
   else
      do k=0,plev
         levels(k+1) = k
      enddo
      do k=plev+2,lvlsiz
         levels(k) = -1
      enddo
   endif

   return
   end subroutine get_block_levels_d
!
!========================================================================
!
   integer function get_lon_d(blockid,bcid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return global longitude index for the column at location bcid
!          in block blockid.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! block column index

!-----------------------------------------------------------------------
!  latitude slice block
   get_lon_d = bcid

   return
   end function get_lon_d
!
!========================================================================
!
   integer function get_lat_d(blockid,bcid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return global latitude index for the column at location bcid
!          in block blockid.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! block column index

!-----------------------------------------------------------------------
!  latitude slice block
   get_lat_d = blockid

   return
   end function get_lat_d

!
!========================================================================
!
   integer function get_block_owner_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return id of processor that "owns" the indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
#if ( defined SPMD )
   use spmd_dyn, only: proc
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
   get_block_owner_d = proc(blockid)
#else
   get_block_owner_d = 0
#endif

   return
   end function get_block_owner_d
!#######################################################################
end module dyn_grid







