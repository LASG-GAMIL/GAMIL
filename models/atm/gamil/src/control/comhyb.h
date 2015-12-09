!----------------------------------------------------------------------- 
! 
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
! 
!-----------------------------------------------------------------------
      real(r8) hyai(plevp)       ! ps0 component of hybrid coordinate - interfaces
      real(r8) hybi(plevp)       ! ps component of hybrid coordinate - interfaces
      real(r8) hyam(plev)        ! ps0 component of hybrid coordinate - midpoints
      real(r8) hybm(plev)        ! ps component of hybrid coordinate - midpoints

      real(r8) hypi(plevp)       ! reference pressures at interfaces
      real(r8) hypm(plev)        ! reference pressures at midpoints

      real(r8) ps0               ! base state sfc pressure for level definitions

!     Currently, GAMIL uses sigma coordinate.
      real(r8) pmtop
      real(r8) sig (plevp)
      real(r8) sigl(plev)
      real(r8) dsig(plev)

      common /comhyb/ hyai, hybi , hyam, hybm
      common /comhyb/ hypi, hypm
      common /comhyb/ ps0 , pmtop
      common /comhyb/ sig , sigl , dsig
 
