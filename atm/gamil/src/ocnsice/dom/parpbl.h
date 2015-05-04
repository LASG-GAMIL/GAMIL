!
! $Id: parpbl.h,v 1.4 2000/06/02 16:19:18 jet Exp $
! $Author: jet $
!
!
! Constants for the surface temperature calculation
!
      real(r8), parameter :: zzocen = 0.0001   ! Ocean aerodynamic roughness length
      real(r8), parameter :: zzsice = 0.0400   ! Sea ice aerodynamic roughness length
      real(r8), parameter :: xkar   = 0.4      ! Von Karman constant
      real(r8), parameter :: hmixmn = 500.0    ! Minimum boundary layer height
      real(r8), parameter :: ric    = 3.05     ! Critical Richardson number
      real(r8), parameter :: rimax  = 0.50*ric ! Maximum Richardson number
      real(r8), parameter :: epsi   = 1.0e-12  ! Small factor to prevent exact zeros
      real(r8), parameter :: zref   = 10.0     ! 10m reference height 
      real(r8), parameter :: umin   = 1.0      ! Minimum wind speed at bottom level
!
 
