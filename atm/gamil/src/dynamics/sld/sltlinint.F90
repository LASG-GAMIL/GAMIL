#include <misc.h>
#include <params.h>
subroutine sltlinint(kdim    ,ikcnst  ,fb      ,xl      ,xr      , &
                     ys      ,yn      ,zt      ,zb      ,idp     , &
                     jdp     ,kdp     ,fdp     ,nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Interpolate field to departure points using tri-linear interpolation
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id: sltlinint.F90,v 1.3.42.1 2002/06/15 13:48:31 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  implicit none

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: kdim             ! vertical dimension
  integer , intent(in)   :: ikcnst           ! constituent index

  real(r8), intent(in)   :: fb (plond,kdim*ikcnst,beglatex:endlatex) ! input field
  real(r8), intent(in)   :: xl (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(in)   :: xr (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(in)   :: ys (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(in)   :: yn (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(in)   :: zt (plon,plev)   ! top vertical interpolation weight 
  real(r8), intent(in)   :: zb (plon,plev)   ! bot vertical interpolation weight 
  integer , intent(in)   :: idp(plon,plev,4) ! index of x-coordinate of dep pt
  integer , intent(in)   :: jdp(plon,plev)   ! index of y-coordinate of dep pt
  integer , intent(in)   :: kdp(plon,plev)   ! index of z-coordinate of dep pt
  real(r8), intent(out)  :: fdp(plon,plev)   ! interpolant
  integer , intent(in)   :: nlon             ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i, k, ii2, ii3, jj, kk             ! indices
!
!-----------------------------------------------------------------------
!
  do k=1,plev
     do i = 1,nlon
        ii2 = idp(i,k,2)
        ii3 = idp(i,k,3)
        jj  = jdp(i,k)
        kk  = kdp(i,k)
        fdp(i,k) = (( fb (ii2  ,kk  ,jj  )*xl (i,k,2) &
                 +    fb (ii2+1,kk  ,jj  )*xr (i,k,2) )*ys(i,k) &
                 +  ( fb (ii3  ,kk  ,jj+1)*xl (i,k,3) &
                 +    fb (ii3+1,kk  ,jj+1)*xr (i,k,3) )*yn(i,k) )*zt(i,k) &
                 + (( fb (ii2  ,kk+1,jj  )*xl (i,k,2) &
                 +    fb (ii2+1,kk+1,jj  )*xr (i,k,2) )*ys(i,k) &
                 +  ( fb (ii3  ,kk+1,jj+1)*xl (i,k,3) &
                 +    fb (ii3+1,kk+1,jj+1)*xr (i,k,3) )*yn(i,k) )*zb(i,k)
     end do
  end do

  return
end subroutine sltlinint

