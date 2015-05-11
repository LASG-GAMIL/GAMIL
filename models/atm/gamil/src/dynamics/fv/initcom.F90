#include <misc.h>
#include <params.h>

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  initcom --- Initialize model commons
!
! !INTERFACE:
subroutine initcom

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use commap
   use physconst,    only: rair
   use dynconst,     only: dynconsti
   use constituents, only: ppcnst, qmin, qmincg

   implicit none

!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------

! !DESCRIPTION:
!
!   Initialize Model commons.
! 
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.02.01      Buja          Modifications
!   01.01.19      Lin           incorporated Rasch's bug fix for the 
!                               "Gaussian" weights
!   02.04.04      Sawyer        Removed comspe
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer m           ! Index for legendre array
#ifdef PVP
   integer n           ! Index for legendre array
#endif
   integer irow        ! Latitude pair index
   integer lat         ! Latitude index

   real(r8) dp             ! Spacing between latitudes
   real(r8) pi             ! Mathematical pi (3.14...)
   real(r8) sum
!
!-----------------------------------------------------------------------
   call dynconsti
!
!
! Initialize commap.  Set hybrid level dependent arrays
!
   call hycoef

! L-R dynamics uses a regular latitude distribution (not gausian).
! The algorithm below is a bastardized version of LSM: map.F.

   pi = 4.0*atan(1.0)
   dp = 180./(plat-1)
   do lat = 1, plat
      latdeg(lat) = -90.0 + (lat-1)*dp
      clat(lat) = latdeg(lat)*pi/180.
   end do

! Calculate latitudes for the staggered grid

   do lat = 1, plat-1
      clat_staggered(lat) = (clat(lat) + clat(lat+1)) / 2.0
   end do

! Weights are defined as cos(phi)*(delta-phi)
! For a sanity check, the sum of w across all lats should be 2, or 1 across
! half of the latitudes.

   do lat = 2, plat-1
      w(lat) = sin(clat_staggered(lat)) - sin(clat_staggered(lat-1))
   end do
   w(1) = sin(clat_staggered(1)) + 1.
   w(plat) = w(1)

   sum = 0.
   do lat=1,plat
      if (iam .eq. 0) write (6,*) 'initcom: lat, clat, w ', lat, clat(lat), w(lat)
      sum = sum + w(lat)
   end do

   if (abs(sum - 2._r8) > 1.e-8) then
      write(6,*) 'INITCOM 1: weights do not sum to 2. sum=',sum
      call endrun
   end if

   dp = pi / float(plat-1)
   do lat = 1, plat-1
      w_staggered(lat) = sin(clat(lat+1)) - sin(clat(lat))
   end do

   sum = 0.
   do lat=1,plat-1
      sum = sum + w_staggered(lat)
   end do

   if (abs(sum - 2._r8) > 1.e-8) then
      write(6,*) 'INITCOM 2: weights do not sum to 2. sum=',sum
      call endrun
   end if
!
! Set minimum mixing ratio for moisture and advected tracers
!
   qmin(1) = 1.e-12          ! Minimum mixing ratio for moisture
   do m=2,ppcnst
      qmin(m) = 0.0
   end do
!
! Set the minimum mixing ratio for the counter-gradient term.  
! Normally this should be the same as qmin, but in order to 
! match control case 414 use zero for water vapor.
!
   qmincg(1) = 0.
   do m=2,ppcnst
      qmincg(m) = qmin(m)
   end do
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   if (iam .eq. 0) write(6,*) 'Number of longitudes per latitude = ', plon
!
! Longitude array
!
   do lat=1,plat
      do i=1,plon
         londeg(i,lat) = (i-1)*360./plon
         clon(i,lat)   = (i-1)*2.0*pi/plon
      end do
   end do
!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!
   dyngrid_set = .true.

   return
!EOC
end subroutine initcom
!-----------------------------------------------------------------------
