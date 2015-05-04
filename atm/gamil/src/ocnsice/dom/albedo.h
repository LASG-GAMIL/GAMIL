!
! $Id: albedo.h,v 1.4 2000/06/02 16:19:16 jet Exp $
! $Author: jet $
!
      real(r8), parameter :: snws  = 0.95  ! Snow albedo for 0.2-0.7 micro-meters
      real(r8), parameter :: snwl  = 0.70  ! Snow albedo for 0.7-5.0 micro-meters    
      real(r8), parameter :: sices = 0.70  ! Sea ice albedo for 0.2-0.7 micro-meters
      real(r8), parameter :: sicel = 0.50  ! Sea ice albedo for 0.7-5.0 micro-meters
!
! Slab ocean model mods
!
      real(r8), parameter :: sicsns = 0.84 ! Sea-ice snow albedo for 0.2-0.7 micro-meters
      real(r8), parameter :: sicsnl = 0.60 ! Sea-ice snow albedo for 0.7-5.0 micro-meters
!
      real(r8), parameter :: sicsmn = 0.50 ! min Sea-ice albedo for 0.2-0.7 micro-meters
      real(r8), parameter :: siclmn = 0.26 ! min Sea-ice albedo for 0.7-5.0 micro-meters
!
 
