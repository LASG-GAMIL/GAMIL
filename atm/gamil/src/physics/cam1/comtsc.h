!
! $Id: comtsc.h,v 1.4 2000/06/02 16:20:40 jet Exp $
! $Author: jet $
!
!
! Constants for surface temperature/energy exchange calculations
!
      common/comtsc/latice  ,tmelt   ,latvap  ,rair    ,stebol  
      common/comtsc/snwedp  
!
      real(r8) latice     ! Latent heat of fusion
      real(r8) tmelt      ! Melting temperature of snow and ice
      real(r8) latvap     ! Latent heat of vaporization
      real(r8) rair       ! Gas constant for dry air
      real(r8) stebol     ! Stefan-Boltzmann constant
      real(r8) snwedp     ! Snow equivalent depth factor 
!
 
