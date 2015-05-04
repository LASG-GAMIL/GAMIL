module comhd

!! (wanhui 2003.10.23)
!! (wanhui 2003.10.24)
!! (wanhui 2004.04.14)
!----------------------------------------------------------------------- 
! Purpose: horizontal diffusion module
!-----------------------------------------------------------------------

!! radically changed by wanhui for lasg framework, 2003.04.28
!! changed into a module by wanhui, 2003.10.23
  
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond, beglatexdyn, endlatexdyn
   use infnan

   implicit none

      real(r8) :: dfs0 
      real(r8) :: dthdfs 

      real(r8),allocatable :: frdt (:,:)
      real(r8),allocatable :: frds (:,:)
      real(r8),allocatable :: frdu (:,:)
      real(r8),allocatable :: frdv (:,:)
      real(r8),allocatable :: frdp (:,:)

      real(r8),allocatable :: dxvpn(:)
      real(r8),allocatable :: dxvps(:)

CONTAINS

  subroutine initialize_comhd
!
! Purpose:  Allocate and initialize the arrays of the horizontal diffusion
!
      allocate (frdt  (beglatexdyn:endlatexdyn,3))
      allocate (frds  (beglatexdyn:endlatexdyn,3))
      allocate (frdu  (beglatexdyn:endlatexdyn,3))
      allocate (frdv  (beglatexdyn:endlatexdyn,3))
      allocate (frdp  (beglatexdyn:endlatexdyn,3))

      allocate (dxvpn (beglatexdyn:endlatexdyn))
      allocate (dxvps (beglatexdyn:endlatexdyn))

      dthdfs = inf

      frdt(:,:) = inf
      frds(:,:) = inf
      frdu(:,:) = inf
      frdv(:,:) = inf
      frdp(:,:) = inf

      dxvpn(:)  = inf
      dxvps(:)  = inf

      return
   end subroutine initialize_comhd

end module comhd
