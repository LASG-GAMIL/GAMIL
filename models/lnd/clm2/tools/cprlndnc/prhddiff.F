subroutine prhddiff (nlev, nlat)
                     
  use shr_kind_mod, only: r8 => shr_kind_r8
  use header

  implicit none
!
! Input arguments
!
  integer, intent(in) :: nlev, nlat
!
! Local workspace
!
  integer k, j, n

  write(6,*)'SUMMARY OF IMPORTANT HEADER DIFFERENCES:'

  do n=1,ntime(1)
    if (nsteph(n,1) .ne. nsteph(n,2)) write(6,100)'nsteph:', nsteph(n,1), nsteph(n,2)
    if (ncdate(n,1).ne.ncdate(n,2)) write(6,100)'ncdate:', ncdate(n,1), ncdate(n,2)
    if (ncsec(n,1).ne.ncsec(n,2)) write(6,100)'ncsec:', ncsec(n,1) , ncsec(n,2) 
  end do

  write(6,*)
  do j=1,nlat
    if (gw(j,1) .ne. gw(j,2)) write(6,200)'gw(',j,'):',gw(j,1), gw(j,2)
  end do

  return
100 format(a,2i9)
200 format(a,i3,a,1p,2e23.15)
  end subroutine prhddiff
