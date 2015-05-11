subroutine prhddiff (nlev, nlevp, nlat)
                     
  use shr_kind_mod, only: r8 => shr_kind_r8
  use header

  implicit none
!
! Input arguments
!
  integer, intent(in) :: nlev, nlevp, nlat
!
! Local workspace
!
  integer k, j, n

  write(6,*)'SUMMARY OF IMPORTANT HEADER DIFFERENCES:'
  if (ntrm(1)   .ne. ntrm(2)  ) write(6,100)'ntrm  :', ntrm(1)  , ntrm(2)  
  if (ntrn(1)   .ne. ntrn(2)  ) write(6,100)'ntrn  :', ntrn(1)  , ntrn(2)  
  if (ntrk(1)   .ne. ntrk(2)  ) write(6,100)'ntrk  :', ntrk(1)  , ntrk(2)  
  if (ndbase(1) .ne. ndbase(2)) write(6,100)'ndbase:', ndbase(1), ndbase(2)
  if (nsbase(1) .ne. nsbase(2)) write(6,100)'nsbase:', nsbase(1), nsbase(2)
  if (nbdate(1) .ne. nbdate(2)) write(6,100)'nbdate:', nbdate(1), nbdate(2)
  if (nbsec(1)  .ne. nbsec(2) ) write(6,100)'nbsec :', nbsec(1) , nbsec(2) 
  if (mdt(1)    .ne. mdt(2)   ) write(6,100)'mdt   :', mdt(1)   , mdt(2)   
  if (p0(1)     .ne. p0(2))     write(6,100)'p0    :', p0(1), p0(2)

  do n=1,ntime(1)
    if (nsteph(n,1) .ne. nsteph(n,2)) write(6,100)'nsteph:', nsteph(n,1), nsteph(n,2)
    if (ncdate(n,1).ne.ncdate(n,2)) write(6,100)'ncdate:', ncdate(n,1), ncdate(n,2)
    if (ncsec(n,1).ne.ncsec(n,2)) write(6,100)'ncsec:', ncsec(n,1) , ncsec(n,2) 
  end do

  write(6,*)
  do k=1,nlevp
    if (hyai(k,1) .ne. hyai(k,2)) write(6,200)'hyai(',k,'):',hyai(k,1),hyai(k,2)
    if (hybi(k,1) .ne. hybi(k,2)) write(6,200)'hybi(',k,'):',hybi(k,1),hybi(k,2)
  end do

  do k=1,nlev
    if (hyam(k,1) .ne. hyam(k,2)) write(6,200)'hyam(',k,'):',hyam(k,1),hyam(k,2)
    if (hybm(k,1) .ne. hybm(k,2)) write(6,200)'hybm(',k,'):',hybm(k,1),hybm(k,2)
  end do

  write(6,*)
  do j=1,nlat
    if (gw(j,1) .ne. gw(j,2)) write(6,200)'gw(',j,'):',gw(j,1), gw(j,2)
  end do

  return
100 format(a,2i9)
200 format(a,i3,a,1p,2e23.15)
  end subroutine prhddiff
