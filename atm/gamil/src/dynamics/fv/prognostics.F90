#include <misc.h>
#include <params.h>

module prognostics
!BOP
!
! !MODULE: prognostics --- dynamics-physics coupling module
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid 
   use infnan
   use constituents, only: ppcnst
   use dynamics_vars, only: ng_d, ng_s

   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:

   public initialize_prognostics, shift_time_indices

   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore

!-----------------------------------------------------------------------
! Time index pointers: only used for prognostic cloud water
!---------------------------------------------------------------------------

   integer :: n3   = 2
   integer :: n3m1 = 1

!-----------------------------------------------------------------------
! The next 5 arrays are the REQUIRED state variables for the f-v dynamics
!-----------------------------------------------------------------------
! dyn_state:
   real(r8), allocatable :: u3s(:,:,:)  ! Staggered grid winds, latitude
   real(r8), allocatable :: v3s(:,:,:)  ! Satggered grid winds, longitude
   real(r8), allocatable :: delp(:,:,:) ! delta pressure
   real(r8), allocatable :: pt(:,:,:)   ! virtual potential temperature (ghosted)
   real(r8), allocatable :: q3(:,:,:,:) ! Moisture and constituents

!-----------------------------------------------------------------------
! Auxilliary pressure (for definition of vertical coordinate) arrays
! These can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_grid:
   real(r8), allocatable :: pe  (:,:,:)     ! edge pressure
   real(r8), allocatable :: pk  (:,:,:)     ! pe**cappa
   real(r8), allocatable :: piln(:,:,:)     ! ln(pe)
   real(r8), allocatable :: pkz (:,:,:)     ! finite-volume mean pk

!-----------------------------------------------------------------------
! The following arrays also existed in phys_state
!-----------------------------------------------------------------------
! phys_state:

   real(r8), allocatable :: phis(:,:)      ! Surface geopotential
   real(r8), allocatable ::   ps(:,:)      ! Surface pressure
   real(r8), allocatable ::   t3(:,:,:)    ! (virtual) Temperature
   real(r8), allocatable :: omga(:,:,:)    ! vertical pressure velocity

! !DESCRIPTION:
!
!   {\bf Purpose:} Prognostic variables held in-core for convenient 
!   access. q3 is specific humidity (water vapor) and other 
!   constituents. pcnst is advected constituents, pnats is non-advected, ppcnst is total.
! 
! !REVISION HISTORY:
!
!   00.08.15     Lin        Modifications
!   00.12.14     Sawyer     SPMD bug: do j=1,plat => beglat,endlat
!   01.03.26     Sawyer     Added ProTeX documentation
!   01.12.10     Sawyer     Ghosted PT
!   01.12.20     Sawyer     Changed index order of Q3
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: initialize_prognostics --- initialize prognostic variables
!
! !INTERFACE:
subroutine initialize_prognostics

! !DESCRIPTION:
!
!   Initialize the prognostic variables
!
! !REVISION HISTORY:
!   00.10.15  Lin      modified to declare only one-time-level; 
!                      may want to declare them as 3D arrays eventually
!   01.03.26  Sawyer   Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer :: i, j, k, m

      allocate (phis(plon      , beglat:endlat))
      allocate (ps  (plon      , beglat:endlat))
      allocate (t3  (plon, beglev:endlev, beglat:endlat))
      allocate (omga(plon, beglev:endlev, beglat:endlat))    

      allocate ( u3s(plon, beglat-ng_d:endlat+ng_s, beglev:endlev))
      allocate ( v3s(plon, beglat-ng_s:endlat+ng_d, beglev:endlev))
      allocate (delp(plon, beglat:endlat, beglev:endlev))
      allocate ( pt(plon, beglat-ng_d:endlat+ng_d, beglev:endlev))
      allocate ( q3(plon, beglat-ng_d:endlat+ng_d, beglev:endlev, ppcnst) )

      allocate (pkz (plon, beglat:endlat, beglev:endlev))
      allocate (pk  (plon, beglat:endlat, beglev:endlevp1))
      allocate (pe  (plon, beglev:endlevp1, beglat:endlat)) 
      allocate (piln(plon, beglev:endlevp1, beglat:endlat)) 

      phis (:,:) = inf         ! phis is used in SMP k-loop

! The order of the "first touch" initialization may need to be r-arranged
! in order to maximize the parallel efficiency

!$omp parallel do private (i, j, k, m)

  do j=beglat, endlat

        do i=1,plon
           ps(i,j) = inf
        enddo

      do k=beglev,endlev
        do i=1,plon
             t3(i,k,j) = inf
           omga(i,k,j) = inf
        enddo
      enddo

      do m=1,ppcnst
         do k=beglev,endlev
           do i=1,plon
              q3(i,j,k,m) = inf
           enddo
         enddo
      enddo

      do k=beglev,endlevp1
         do i=1,plon
            pe  (i,k,j) = inf
            piln(i,k,j) = inf
         enddo
      enddo
  enddo

!$omp parallel do private (i, j, k)

  do k=beglev,endlev
      do j=beglat,endlat
         do i=1,plon
            u3s(i,j,k) = inf
            v3s(i,j,k) = inf
             pt(i,j,k) = inf
           delp(i,j,k) = inf
            pkz(i,j,k) = inf
         enddo
      enddo
  enddo

!$omp parallel do private (i, j, k)

  do k=beglev,endlevp1
      do j=beglat,endlat
        do i=1,plon
           pk(i,j,k) = inf
        enddo
      enddo
  enddo

   return
!EOC
 end subroutine initialize_prognostics
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: shift_time_indices --- Shift time indices
!
! !INTERFACE:
 subroutine shift_time_indices

! !DESCRIPTION:
!
!   Shift time indices
!
! !REVISION HISTORY:
!   01.03.26  Sawyer   Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer :: itmp

      itmp = n3m1
      n3m1 = n3
      n3   = itmp
      return
!EOC
   end subroutine shift_time_indices
!-----------------------------------------------------------------------
 end module prognostics

