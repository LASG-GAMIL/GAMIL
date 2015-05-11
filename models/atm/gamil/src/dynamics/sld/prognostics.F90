#include <misc.h>
#include <params.h>

module prognostics

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

   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore
   integer :: n3   = 2
   integer :: n3m1 = 1

   real(r8), allocatable :: ps(:,:,:)
   real(r8), allocatable :: u3(:,:,:,:)
   real(r8), allocatable :: v3(:,:,:,:)
   real(r8), allocatable :: t3(:,:,:,:)
   real(r8), allocatable :: q3(:,:,:,:,:)
   real(r8), allocatable :: u3sld(:,:,:)
   real(r8), allocatable :: v3sld(:,:,:)
   real(r8), allocatable :: lnpssld(:,:,:)
   real(r8), allocatable :: prhssld(:,:,:)
   real(r8), allocatable :: tarrsld(:,:,:)
   real(r8), allocatable :: parrsld(:,:,:)
   real(r8), allocatable :: etadot(:,:,:,:)
   
   real(r8), allocatable :: div(:,:,:,:)    ! divergence

   real(r8), allocatable :: urhs(:,:,:)
   real(r8), allocatable :: vrhs(:,:,:)
   real(r8), allocatable :: trhs(:,:,:)
   real(r8), allocatable :: prhs(:,:,:)
   
   real(r8), allocatable :: ql(:,:,:)
   real(r8), allocatable :: qm(:,:,:)
   real(r8), allocatable :: ed1(:,:,:)
   real(r8), allocatable :: tlm1(:,:,:)
   real(r8), allocatable :: tl(:,:,:)
   real(r8), allocatable :: tmm1(:,:,:)
   real(r8), allocatable :: tm(:,:,:)
   real(r8), allocatable :: omga(:,:,:)     ! vertical velocity
   real(r8), allocatable :: dpsl(:,:)       ! longitudinal pressure gradient
   real(r8), allocatable :: dpslm1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable :: dpslp1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable :: dpsm(:,:)       ! meridional pressure gradient
   real(r8), allocatable :: dpsmm1(:,:)     ! meridional pressure gradient
   real(r8), allocatable :: dpsmp1(:,:)     ! meridional pressure gradient
   real(r8), allocatable :: dps(:,:)        ! pressure gradient
   real(r8), allocatable :: phis(:,:)       ! surface geopotential
   real(r8), allocatable :: phisl(:,:)      ! surface geopotential
   real(r8), allocatable :: phism(:,:)      ! surface geopotential
   
CONTAINS

  subroutine initialize_prognostics
!
! Purpose: Allocate and initialize the prgnostic variables
!
    allocate (ps(plond                 ,beglat:endlat    ,ptimelevels))
    allocate (u3(plond,plev            ,beglatex:endlatex,ptimelevels))
    allocate (v3(plond,plev            ,beglatex:endlatex,ptimelevels))
    allocate (t3(plond,plev            ,beglatex:endlatex,ptimelevels))
    allocate (q3(plond,plev,pcnst+pnats,beglatex:endlatex,ptimelevels))
    allocate (u3sld(plond,plev,beglatex:endlatex))
    allocate (v3sld(plond,plev,beglatex:endlatex))
    allocate (lnpssld(plond,plev,beglatex:endlatex))
    allocate (prhssld(plond,plev,beglatex:endlatex))
    allocate (tarrsld(plond,plev,beglatex:endlatex))
    allocate (parrsld(plond,plev,beglatex:endlatex))
    allocate (etadot(plond,plevp,beglatex:endlatex,ptimelevels))

    allocate (div   (plond,plev,beglat:endlat,ptimelevels))

    allocate (urhs(plond,plev,beglat:endlat))
    allocate (vrhs(plond,plev,beglat:endlat))
    allocate (trhs(plond,plev,beglat:endlat))
    allocate (prhs(plond,plev,beglat:endlat))

    allocate (ql     (plond,plev,beglat:endlat))    
    allocate (qm     (plond,plev,beglat:endlat))    
    allocate (ed1    (plond,plev,beglat:endlat))  
    allocate (tlm1   (plond,plev,beglat:endlat))   
    allocate (tl     (plond,plev,beglat:endlat))   
    allocate (tmm1   (plond,plev,beglat:endlat))   
    allocate (tm     (plond,plev,beglat:endlat))   
    allocate (omga   (plond,plev,beglat:endlat))    
    allocate (dpsl   (plond,beglat:endlat))        
    allocate (dpslm1 (plond,beglat:endlat))        
    allocate (dpslp1 (plond,beglat:endlat))        
    allocate (dpsm   (plond,beglat:endlat))        
    allocate (dpsmm1 (plond,beglat:endlat))        
    allocate (dpsmp1 (plond,beglat:endlat))        
    allocate (dps    (plond,beglat:endlat))         
    allocate (phis   (plond,beglat:endlat))        
    allocate (phisl  (plond,beglat:endlat))        
    allocate (phism  (plond,beglat:endlat))        

    ps(:,:,:)     = inf
    u3(:,:,:,:)   = inf
    v3(:,:,:,:)   = inf
    t3(:,:,:,:)   = inf
    q3(:,:,:,:,:) = inf
    u3sld(:,:,:) = inf
    v3sld(:,:,:) = inf
    lnpssld(:,:,:) = inf
    prhssld(:,:,:) = inf
    tarrsld(:,:,:) = inf
    parrsld(:,:,:) = inf
    etadot(:,:,:,:) = inf

    div(:,:,:,:) = inf

    urhs(:,:,:) = inf
    vrhs(:,:,:) = inf
    trhs(:,:,:) = inf
    prhs(:,:,:) = inf

    ql     (:,:,:) = inf
    qm     (:,:,:) = inf
    ed1    (:,:,:) = inf
    tlm1   (:,:,:) = inf
    tl     (:,:,:) = inf
    tmm1   (:,:,:) = inf
    tm     (:,:,:) = inf
    omga   (:,:,:) = inf
    dpsl   (:,:) = inf
    dpslm1 (:,:) = inf
    dpslp1 (:,:) = inf
    dpsm   (:,:) = inf
    dpsmm1 (:,:) = inf
    dpsmp1 (:,:) = inf
    dps    (:,:) = inf
    phis   (:,:) = inf
    phisl  (:,:) = inf
    phism  (:,:) = inf

    return
  end subroutine initialize_prognostics

  subroutine shift_time_indices
!
! Purpose: Shift the time indices that track the current and previous times.
!
     integer :: itmp

     itmp = n3m1
     n3m1 = n3
     n3   = itmp
   end subroutine shift_time_indices

end module prognostics




