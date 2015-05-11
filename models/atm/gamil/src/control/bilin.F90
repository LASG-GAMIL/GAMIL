#include <misc.h>
#include <params.h>

subroutine bilin (arrin, xin, yin, nlondin, nlonin, &
                  nlevdin, nlev, nlatin, arrout, xout, &
                  yout, nlondout, nlonout, nlevdout, nlatout)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! Do a bilinear interpolation from input mesh defined by xin, yin to output
! mesh defined by xout, yout.  Circularity is assumed in the x-direction so
! input x-grid must be in degrees east and must start from Greenwich.  When
! extrapolation is necessary in the N-S direction, values will be copied 
! from the extreme edge of the input grid.  Vectorization is over the
! longitude dimension.  Input grid is assumed rectangular. Output grid
! is assumed ragged, i.e. xout is a 2-d mesh.
! 
! Method: Interpolate first in longitude, then in latitude.
! 
! Author: Jim Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: nlondin                        ! longitude dimension of input grid
   integer, intent(in) :: nlonin                         ! number of real longitudes (input)
   integer, intent(in) :: nlevdin                        ! vertical dimension of input grid
   integer, intent(in) :: nlev                           ! number of vertical levels
   integer, intent(in) :: nlatin                         ! number of input latitudes
   integer, intent(in) :: nlatout                        ! number of output latitudes
   integer, intent(in) :: nlondout                       ! longitude dimension of output grid
   integer, intent(in) :: nlonout(nlatout)               ! number of output longitudes per lat
   integer, intent(in) :: nlevdout                       ! vertical dimension of output grid

   real(r8), intent(in) :: arrin(nlondin,nlevdin,nlatin) ! input array of values to interpolate
   real(r8), intent(in) :: xin(nlondin)                  ! input x mesh
   real(r8), intent(in) :: yin(nlatin)                   ! input y mesh
   real(r8), intent(in) :: xout(nlondout,nlatout)        ! output x mesh
   real(r8), intent(in) :: yout(nlatout)                 ! output y mesh
!
! Output arguments
!
   real(r8), intent(out) :: arrout(nlondout,nlevdout,nlatout) ! interpolated array
!
! Local workspace
!
   integer :: i, ii, iw, ie, iiprev ! longitude indices
   integer :: j, jj, js, jn, jjprev ! latitude indices
   integer :: k                     ! level index
   integer :: icount                ! number of bad values

   real(r8) :: extrap               ! percent grid non-overlap
   real(r8) :: dxinwrap             ! delta-x on input grid for 2-pi
   real(r8) :: avgdxin              ! avg input delta-x
   real(r8) :: ratio                ! compare dxinwrap to avgdxin
   real(r8) :: sum                  ! sum of weights (used for testing)
!
! Dynamic
!
   integer :: iim(nlondout)         ! interpolation index to the left
   integer :: iip(nlondout)         ! interpolation index to the right
   integer :: jjm(nlatout)          ! interpolation index to the south
   integer :: jjp(nlatout)          ! interpolation index to the north

   real(r8) :: wgts(nlatout)        ! interpolation weight to the north
   real(r8) :: wgtn(nlatout)        ! interpolation weight to the north
   real(r8) :: wgte(nlondout)       ! interpolation weight to the north
   real(r8) :: wgtw(nlondout)       ! interpolation weight to the north
   real(r8) :: igrid(nlonin)        ! interpolation weight to the north
!
! Check validity of input coordinate arrays: must be monotonically increasing,
! and have a total of at least 2 elements
!
   if (nlonin < 2 .or. nlatin < 2) then
      write(6,*)'BILIN: Must have at least 2 input points for ', &
           'interpolation'
      call endrun
   end if

   if (xin(1) < 0. .or. xin(nlonin) > 360.) then
      write(6,*)'BILIN: Input x-grid must be between 0 and 360'
      call endrun
   end if

   icount = 0
   do i=1,nlonin-1
      if (xin(i) >= xin(i+1)) icount = icount + 1
   end do

   do j=1,nlatin-1
      if (yin(j) >= yin(j+1)) icount = icount + 1
   end do
  
   write(*,*) nlatout
   do j=1,nlatout-1
      if (yout(j) >= yout(j+1)) icount = icount + 1
   end do

   do j=1,nlatout
      do i=1,nlonout(j)-1
         if (xout(i,j) >= xout(i+1,j)) icount = icount + 1
      end do
   end do

   if (icount > 0) then
      write(6,*)'BILIN: Non-monotonic coordinate array(s) found'
      call endrun
   end if

   if (yout(nlatout) <= yin(1) .or. yout(1) >= yin(nlatin)) then
      write(6,*)'BILIN: No overlap in y-coordinates'
      call endrun
   end if

   do j=1,nlatout
      if (xout(1,j) < 0. .or. xout(nlonout(j),j) > 360.) then
         write(6,*)'BILIN: Output x-grid must be between 0 and 360'
         call endrun
      end if

      if (xout(nlonout(j),j) <= xin(1) .or.  &
          xout(1,j)          >= xin(nlonin)) then
         write(6,*)'BILIN: No overlap in x-coordinates'
         call endrun
      end if
   end do
!
! Initialize index arrays for later checking
!
   do j=1,nlatout
      jjm(j) = 0
      jjp(j) = 0
   end do
!
! For values which extend beyond N and S boundaries, set weights
! such that values will just be copied.
!
   do js=1,nlatout
      if (yout(js) > yin(1)) exit
      jjm(js) = 1
      jjp(js) = 1
      wgts(js) = 1.
      wgtn(js) = 0.
   end do

   do jn=nlatout,1,-1
      if (yout(jn) <= yin(nlatin)) exit
      jjm(jn) = nlatin
      jjp(jn) = nlatin
      wgts(jn) = 1.
      wgtn(jn) = 0.
   end do
!
! Loop though output indices finding input indices and weights
!
   jjprev = 1
   do j=js,jn
      do jj=jjprev,nlatin-1
         if (yout(j) > yin(jj) .and. yout(j) <= yin(jj+1)) then
            jjm(j) = jj
            jjp(j) = jj + 1
            wgts(j) = (yin(jj+1) - yout(j)) / (yin(jj+1) - yin(jj))
            wgtn(j) = (yout(j)   - yin(jj)) / (yin(jj+1) - yin(jj))
            goto 30
         end if
      end do
      write(6,*)'BILIN: Failed to find interp values'
      call endrun
30    jjprev = jj
   end do

   dxinwrap = xin(1) + 360. - xin(nlonin)
!
! Check for sane dxinwrap values.  Allow to differ no more than 10% from avg
!
   avgdxin = (xin(nlonin)-xin(1))/(nlonin-1.)
   ratio = dxinwrap/avgdxin
   if (ratio < 0.9 .or. ratio > 1.1) then
      write(6,*)'BILIN: Insane dxinwrap value =',dxinwrap,' avg=', &
           avgdxin
      call endrun
   end if
!
! Check grid overlap
!
   extrap = 100.*((js - 1.) + float(nlatout - jn))/nlatout
   if (extrap > 20.) then
      write(6,*)'BILIN:',extrap,' % of N/S output', &
           ' grid will have to be extrapolated'
   end if
!
! Check that interp/extrap points have been found for all outputs, and that
! they are reasonable (i.e. within 32-bit roundoff)
!
   icount = 0
   do j=1,nlatout
      if (jjm(j) == 0 .or. jjp(j) == 0) icount = icount + 1
      sum = wgts(j) + wgtn(j)
      if (sum < 0.99999 .or. sum > 1.00001) icount = icount + 1
      if (wgts(j) < 0. .or. wgts(j) > 1.) icount = icount + 1
      if (wgtn(j) < 0. .or. wgtn(j) > 1.) icount = icount + 1
   end do

   if (icount > 0) then
      write(6,*)'BILIN: something bad in latitude indices or weights'
      call endrun
   end if
!
! Do the bilinear interpolation
!
   do j=1,nlatout
!
! Initialize index arrays for later checking
!
      do i=1,nlondout
         iim(i) = 0
         iip(i) = 0
      end do
!
! For values which extend beyond E and W boundaries, set weights such that
! values will be interpolated between E and W edges of input grid.
!
      do iw=1,nlonout(j)
         if (xout(iw,j) > xin(1)) exit
         iim(iw) = nlonin
         iip(iw) = 1
         wgtw(iw) = (xin(1)        - xout(iw,j))   /dxinwrap
         wgte(iw) = (xout(iw,j)+360. - xin(nlonin))/dxinwrap
      end do

      do ie=nlonout(j),1,-1
         if (xout(ie,j) <= xin(nlonin)) exit
         iim(ie) = nlonin
         iip(ie) = 1
         wgtw(ie) = (xin(1)+360. - xout(ie,j))   /dxinwrap
         wgte(ie) = (xout(ie,j)    - xin(nlonin))/dxinwrap
      end do
!
! Loop though output indices finding input indices and weights
!
      iiprev = 1
      do i=iw,ie
         do ii=iiprev,nlonin-1
            if (xout(i,j) > xin(ii) .and. xout(i,j) <= xin(ii+1)) then
               iim(i) = ii
               iip(i) = ii + 1
               wgtw(i) = (xin(ii+1) - xout(i,j)) / (xin(ii+1) - xin(ii))
               wgte(i) = (xout(i,j) - xin(ii))   / (xin(ii+1) - xin(ii))
               goto 60
            end if
         end do
         write(6,*)'BILIN: Failed to find interp values'
         call endrun
60       iiprev = ii
      end do

      icount = 0
      do i=1,nlonout(j)
         if (iim(i) == 0 .or. iip(i) == 0) icount = icount + 1
         sum = wgtw(i) + wgte(i)
         if (sum < 0.99999 .or. sum > 1.00001) icount = icount + 1
         if (wgtw(i) < 0. .or. wgtw(i) > 1.) icount = icount + 1
         if (wgte(i) < 0. .or. wgte(i) > 1.) icount = icount + 1
      end do

      if (icount > 0) then
         write(6,*)'BILIN: j=',j,' Something bad in longitude indices or weights'
         call endrun
      end if
!
! Do the interpolation, 1st in longitude then latitude
!
      do k=1,nlev
         do i=1,nlonin
            igrid(i) = arrin(i,k,jjm(j))*wgts(j) + arrin(i,k,jjp(j))*wgtn(j)
         end do

         do i=1,nlonout(j)
            arrout(i,k,j) = igrid(iim(i))*wgtw(i) + igrid(iip(i))*wgte(i)
         end do
      end do
   end do

   return
end subroutine bilin
