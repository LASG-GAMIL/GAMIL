subroutine snowdp2lev

!----------------------------------------------------------------------- 
! 
! Purpose: 
! create snow layers and interfaces given snow depth
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: snowdp2lev.F90,v 1.2.10.2 2002/06/15 13:50:37 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : nlevsoi, nlevsno, nlevlak
  use clm_varmap, only : begpatch, endpatch
  use clm_varder, only : clm
  implicit none

! ------------------- local variables -----------------------------
  integer i,k    !indices
! -----------------------------------------------------------------

! note that clm%zi(0) is set in routine iniTimeConst

  do k = begpatch,endpatch
     clm(k)%dz(-nlevsno+1:0) = 1.e36 
     clm(k)%z (-nlevsno+1:0) = 1.e36 
     clm(k)%zi(-nlevsno:-1)  = 1.e36 
     if (.not. clm(k)%lakpoi) then  !not lake
        if (clm(k)%snowdp < 0.01) then
           clm(k)%snl = 0
           clm(k)%dz(-nlevsno+1:0) = 0.
           clm(k)%z (-nlevsno+1:0) = 0.
           clm(k)%zi(-nlevsno+0:0) = 0.
        else
           if ((clm(k)%snowdp >= 0.01) .AND. (clm(k)%snowdp <= 0.03)) then
              clm(k)%snl = -1
              clm(k)%dz(0)  = clm(k)%snowdp
           else if ((clm(k)%snowdp > 0.03) .AND. (clm(k)%snowdp <= 0.04)) then
              clm(k)%snl = -2
              clm(k)%dz(-1) = clm(k)%snowdp/2.
              clm(k)%dz( 0) = clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.04) .AND. (clm(k)%snowdp <= 0.07)) then
              clm(k)%snl = -2
              clm(k)%dz(-1) = 0.02
              clm(k)%dz( 0) = clm(k)%snowdp - clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.07) .AND. (clm(k)%snowdp <= 0.12)) then
              clm(k)%snl = -3
              clm(k)%dz(-2) = 0.02
              clm(k)%dz(-1) = (clm(k)%snowdp - 0.02)/2.
              clm(k)%dz( 0) = clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.12) .AND. (clm(k)%snowdp <= 0.18)) then
              clm(k)%snl = -3
              clm(k)%dz(-2) = 0.02
              clm(k)%dz(-1) = 0.05
              clm(k)%dz( 0) = clm(k)%snowdp - clm(k)%dz(-2) - clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.18) .AND. (clm(k)%snowdp <= 0.29)) then
              clm(k)%snl = -4
              clm(k)%dz(-3) = 0.02
              clm(k)%dz(-2) = 0.05
              clm(k)%dz(-1) = (clm(k)%snowdp - clm(k)%dz(-3) - clm(k)%dz(-2))/2.
              clm(k)%dz( 0) = clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.29) .AND. (clm(k)%snowdp <= 0.41)) then
              clm(k)%snl = -4
              clm(k)%dz(-3) = 0.02
              clm(k)%dz(-2) = 0.05
              clm(k)%dz(-1) = 0.11
              clm(k)%dz( 0) = clm(k)%snowdp - clm(k)%dz(-3) - clm(k)%dz(-2) - clm(k)%dz(-1)
           else if ((clm(k)%snowdp > 0.41) .AND. (clm(k)%snowdp <= 0.64)) then
              clm(k)%snl = -5
              clm(k)%dz(-4) = 0.02
              clm(k)%dz(-3) = 0.05
              clm(k)%dz(-2) = 0.11
              clm(k)%dz(-1) = (clm(k)%snowdp - clm(k)%dz(-4) - clm(k)%dz(-3) - clm(k)%dz(-2))/2.
              clm(k)%dz( 0) = clm(k)%dz(-1)
           else if (clm(k)%snowdp > 0.64) then 
              clm(k)%snl = -5
              clm(k)%dz(-4) = 0.02
              clm(k)%dz(-3) = 0.05
              clm(k)%dz(-2) = 0.11
              clm(k)%dz(-1) = 0.23
              clm(k)%dz( 0)=clm(k)%snowdp-clm(k)%dz(-4)-clm(k)%dz(-3)-clm(k)%dz(-2)-clm(k)%dz(-1)
           endif
           do i = 0, clm(k)%snl+1, -1
              clm(k)%z(i)    = clm(k)%zi(i) - 0.5*clm(k)%dz(i)
              clm(k)%zi(i-1) = clm(k)%zi(i) - clm(k)%dz(i)
           enddo
        endif
     else   !lake points
        clm(k)%snl = 0
        clm(k)%dz(-nlevsno+1:0) = 0.
        clm(k)%z (-nlevsno+1:0) = 0.
        clm(k)%zi(-nlevsno+0:0) = 0.
     endif
  end do

end subroutine snowdp2lev
