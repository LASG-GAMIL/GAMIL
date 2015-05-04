#include <misc.h>
#include <params.h>

subroutine vertinterp(ncol, ncold, nlev, pmid, pout, arrin, arrout)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Vertically interpolate input array to output pressure level
! Copy values at boundaries.
! 
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: ncol              ! column dimension
  integer , intent(in)  :: ncold             ! declared column dimension
  integer , intent(in)  :: nlev              ! vertical dimension
  real(r8), intent(in)  :: pmid(ncold,nlev)  ! input level pressure levels 
  real(r8), intent(in)  :: pout              ! output pressure level 
  real(r8), intent(in)  :: arrin(ncold,nlev) ! input  array
  real(r8), intent(out) :: arrout(ncold)     ! output array (interpolated)
!--------------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k               ! indices
  integer kupper(ncold)     ! Level indices for interpolation
  real(r8) dpu              ! upper level pressure difference
  real(r8) dpl              ! lower level pressure difference
  logical found(ncold)      ! true if input levels found
  logical error             ! error flag 
!-----------------------------------------------------------------
!
! Initialize index array and logical flags
!
  do i=1,ncol
     found(i)  = .false.
     kupper(i) = 1
  end do
  error = .false.
!     
! Store level indices for interpolation. 
! If all indices for this level have been found, 
! do the interpolation 
!     
  do k=1,nlev-1
     do i=1,ncol
        if ((.not. found(i)) .and. pmid(i,k)<pout .and. pout<=pmid(i,k+1)) then
           found(i) = .true.
           kupper(i) = k
        end if
     end do
  end do
!
! If we've fallen through the k=1,nlev-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top data level for at least some
! of the longitude points.
!
  do i=1,ncol
     if (pout <= pmid(i,1)) then
        arrout(i) = arrin(i,1)
     else if (pout >= pmid(i,nlev)) then
        arrout(i) = arrin(i,nlev)
     else if (found(i)) then
        dpu = pout - pmid(i,kupper(i))
        dpl = pmid(i,kupper(i)+1) - pout
        arrout(i) = (arrin(i,kupper(i)  )*dpl + arrin(i,kupper(i)+1)*dpu)/(dpl + dpu)
     else
        error = .true.
     end if
  end do
!     
! Error check
!
  if (error) then
     write(6,*)'VERTINTERP: ERROR FLAG'
     call endrun
  end if

  return
end subroutine vertinterp
