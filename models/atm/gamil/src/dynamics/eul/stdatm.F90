module stdatm

!! (wanhui 2003.10.23)
!-----------------------------------------------------------------------
! Purpose:  variables of the LASG std.atm.
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plond, beglatexdyn, endlatexdyn, plevstd
   use infnan

   implicit none

      real(r8) :: tbb (plevstd)         
      real(r8) :: hbb (plevstd)         
      real(r8) :: cbb (plevstd)         
      real(r8) :: dcbb(plevstd)         
      real(r8) :: p00, t00              

      real(r8),allocatable :: psb (:,:)      
      real(r8),allocatable :: tsb (:,:)      

CONTAINS

    subroutine initialize_stdatm
!
! Purpose:  Allocate and initialize the arrays of the standard atmosphere.
!
      allocate (psb (plond, beglatexdyn:endlatexdyn))
      allocate (tsb (plond, beglatexdyn:endlatexdyn))

      tbb (:) = inf
      hbb (:) = inf
      cbb (:) = inf
      dcbb(:) = inf
      
      p00     = inf
      t00     = inf

      psb (:,:) = inf
      tsb (:,:) = inf
 
    end subroutine initialize_stdatm

end module stdatm


