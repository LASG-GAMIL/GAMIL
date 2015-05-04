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

contains
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
   use pmgrid, only: plat, twod_decomp, nprxy_x, nprxy_y

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
   block_first = 1
   if (twod_decomp .eq. 1) then
! Assume 1 block per subdomain
      block_last  = nprxy_x*nprxy_y
   else
      block_last  = plat
   endif

   return
   end subroutine get_block_bounds_d
!
!========================================================================
!
   integer function get_block_coord_cnt_d(glon,glat)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of blocks that contain data for the vertical 
!          column at location (glon,glat)
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
   use pmgrid, only: twod_decomp, nprxy_x, nprxy_y
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: glon     ! global longitude index
   integer, intent(in) :: glat     ! global latitude index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
!---------------------------Local workspace-----------------------------
!
    integer i,j,ii,jj                   ! loop indices
    integer ddlon                       ! number of longitudes in block
!-----------------------------------------------------------------------
!  lon/lat block
   if (cnt < 1) then
      write(6,*)'GET_BLOCK_COORD_D: arrays not large enough (', &
                          cnt,' < ',1,' ) '
      call endrun
   else
      if (twod_decomp .eq. 1) then
! Determine block coordinates (ii,jj), where ii ranges from 1 to
! nprxy_x and jj ranges from 1 to nprxy_y.
#if ( defined SPMD )
         ii=0
         do i=1,nprxy_x
            if (lonrangexy(1,i) .le. glon .and. glon .le. lonrangexy(2,i)) ii=i
         enddo
         jj=0
         do j=1,nprxy_y
            if (latrangexy(1,j) .le. glat .and. glat .le. latrangexy(2,j)) jj=j
         enddo
         if (ii .eq. 0 .or. jj .eq. 0) then
            write(6,*)'GET_BLOCK_COORD_D: could not find block indices for (', &
                      glon,',',glat,' ) '
            call endrun
         endif

! Global block index
         blockid(1) = (jj-1)*nprxy_x+ii

! Local coordinates in block
         j = glat-latrangexy(1,jj)+1
         i = glon-lonrangexy(1,ii)+1
         ddlon = lonrangexy(2,ii)-lonrangexy(1,ii)+1

! Local column index in block
         bcid(1) = (j-1)*ddlon+i
!
#endif
      else
         blockid(1) = glat
         bcid(1)    = glon
      endif
!
      do j=2,cnt
         blockid(j) = -1
         bcid(j)    = -1
      enddo
!
   endif
!
   return
   end subroutine get_block_coord_d
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
   use shr_kind_mod, only: r8 => shr_kind_r8
   use rgrid, only: nlon
   use pmgrid, only: twod_decomp, nprxy_x
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none
!
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id


!---------------------------Local workspace-----------------------------
   integer i, j

!-----------------------------------------------------------------------
   if (twod_decomp .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      get_block_col_cnt_d = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
         (latrangexy(2,j)-latrangexy(1,j)+1)
#endif
   else
      get_block_col_cnt_d = nlon(blockid)
   endif

   return
   end function get_block_col_cnt_d
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
   use pmgrid, only: twod_decomp, nprxy_x
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! block column index

!-----------------------------------------------------------------------
   integer i, j, ddlon, ddlat, ii, jj
!  latitude slice block
   if (twod_decomp .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      ddlon = lonrangexy(2,i)-lonrangexy(1,i)+1
      ddlat = latrangexy(2,j)-latrangexy(1,j)+1
      jj = (bcid-1) / ddlon + 1
      ii = bcid - (jj-1) * ddlon
      get_lon_d = lonrangexy(1,i) - 1 + ii
#endif
   else
      get_lon_d = bcid
   endif

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
   use pmgrid, only: twod_decomp, nprxy_x
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! block column index

!-----------------------------------------------------------------------
   integer i, j, ddlon, ddlat, ii, jj
!  latitude slice block
   if (twod_decomp .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      ddlon = lonrangexy(2,i)-lonrangexy(1,i)+1
      ddlat = latrangexy(2,j)-latrangexy(1,j)+1
      jj = (bcid-1) / ddlon + 1
      ii = bcid - (jj-1) * ddlon
      get_lat_d = latrangexy(1,j) - 1 + jj
#endif
   else
      get_lat_d = blockid
   endif

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
   use pmgrid, only: twod_decomp
#if ( defined SPMD )
   use spmd_dyn, only: proc
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
   if (twod_decomp .eq. 1) then
      get_block_owner_d = blockid - 1
   else
      get_block_owner_d = proc(blockid)
   endif
#else
   get_block_owner_d = 0
#endif

   return
   end function get_block_owner_d
!#######################################################################
end module dyn_grid







