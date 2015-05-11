#include <misc.h>
#include <params.h>

module prognostics

!! (wanhui 2003.05.14)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prognostic variables held in-core for convenient access.
! q3 is specific humidity (water vapor) and other constituents.
! pcnst is advected constituents, pnats is non-advected.
! 
! Author: G. Grant
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid 
   use infnan
   use constituents, only: pcnst, pnats

   implicit none

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore
   integer :: n3   = 3
   integer :: n3m1 = 2
   integer :: n3m2 = 1

   real(r8), allocatable :: ps(:,:,:)
   real(r8), allocatable :: u3(:,:,:,:)
   real(r8), allocatable :: v3(:,:,:,:)
   real(r8), allocatable :: t3(:,:,:,:)
   real(r8), allocatable :: t31(:,:,:)
   real(r8), allocatable :: t32(:,:,:)
   real(r8), allocatable :: q3(:,:,:,:,:)
   real(r8), allocatable :: q31(:,:,:)
   real(r8), allocatable :: q32(:,:,:)
   real(r8), allocatable :: qminus(:,:,:,:)

!! real(r8), allocatable :: vort(:,:,:,:)   ! vorticity
!! real(r8), allocatable :: div(:,:,:,:)    ! divergence

!! real(r8), allocatable :: dpsl(:,:)       ! longitudinal pressure gradient
!! real(r8), allocatable :: dpsm(:,:)       ! meridional pressure gradient
!! real(r8), allocatable :: dps(:,:)        ! pressure gradient
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity

CONTAINS

   subroutine initialize_prognostics
!
! Purpose:  Allocate and initialize the prognostic arrays.
!
      allocate (ps    (plond                 ,beglat:endlat    ,ptimelevels))
      allocate (u3    (plond,plev            ,beglatex:endlatex,ptimelevels))
      allocate (v3    (plond,plev            ,beglatex:endlatex,ptimelevels))
      allocate (t3    (plond,plev            ,beglatex:endlatex,ptimelevels))
      allocate (t31   (plond,plev            ,beglatex:endlatex))
      allocate (t32   (plond,plev            ,beglatex:endlatex))
      allocate (q3    (plond,plev,pcnst+pnats,beglatex:endlatex,ptimelevels))
      allocate (q31   (plond,plev            ,beglatex:endlatex))
      allocate (q32   (plond,plev            ,beglatex:endlatex))
      allocate (qminus(plond,plev,pcnst      ,beglatex:endlatex  ))

!!    allocate (vort  (plond,plev,beglat:endlat,ptimelevels))   
!!    allocate (div   (plond,plev,beglat:endlat,ptimelevels))    

!!    allocate (dpsl  (plond,beglat:endlat))        
!!    allocate (dpsm  (plond,beglat:endlat))        
!!    allocate (dps   (plond,beglat:endlat))         
      allocate (phis  (plond,beglat:endlat))        
      allocate (omga  (plond,plev,beglat:endlat))    

      ps(:,:,:)       = inf
      u3(:,:,:,:)     = inf
      v3(:,:,:,:)     = inf
      t3(:,:,:,:)     = inf
      t31(:,:,:)      = inf
      t32(:,:,:)      = inf
      q3(:,:,:,:,:)   = inf
      q31(:,:,:)      = inf
      q32(:,:,:)      = inf
      qminus(:,:,:,:) = inf

!!    vort(:,:,:,:)   = inf
!!    div (:,:,:,:)   = inf

!!    dpsl  (:,:) = inf
!!    dpsm  (:,:) = inf
!!    dps   (:,:) = inf
      phis  (:,:) = inf
      omga  (:,:,:) = inf

      return
   end subroutine initialize_prognostics

   subroutine shift_time_indices
!
! Purpose: 
! Shift the indices that keep track of which index stores
! the relative times (current time, previous, time before previous etc).
!
      integer :: itmp

!!    write(6,*) 'before shifting-----'
!!    write(6,*) 'n3m2 = ',n3m2,';  n3 = ',n3

      itmp = n3m2

!!  (wanhui 2003.06.13)
!
!      n3m2 = n3m1
!      n3m1 = n3

      n3m2 = n3
      n3   = itmp

!!    write(6,*) 'after shifting-----'
!!    write(6,*) 'n3m2 = ',n3m2,';  n3 = ',n3
!!
   end subroutine shift_time_indices

end module prognostics
