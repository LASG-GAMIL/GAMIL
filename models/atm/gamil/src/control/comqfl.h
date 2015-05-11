!----------------------------------------------------------------------- 
! 
! Purpose: Global integrals for moisture and mass conservation and geopotential height
! 
!-----------------------------------------------------------------------
!
   real(r8) tmass(plat)  ! Mass integral for each latitude pair
   real(r8) tmass0       ! Specified dry mass of atmosphere
   real(r8) tmassf       ! Global mass integral
   real(r8) qmassf       ! Global moisture integral
   real(r8) fixmas       ! Proportionality factor for ps in dry mass fixer
   real(r8) qmass1       ! Contribution to global moisture integral (mass
                         !  weighting is based upon the "A" part of the hybrid grid)
   real(r8) qmass2       ! Contribution to global moisture integral (mass
                         !  weighting is based upon the "B" part of the hybrid grid)
   real(r8) pdela(plond,plev)   ! pressure difference between interfaces (pressure
                         !  defined using the "A" part of hybrid grid only)
   real(r8) zgsint       ! global integral of geopotential height
 
   common /comqfl/ tmass, tmass0, tmassf
   common /comqfl/ qmassf ,fixmas ,qmass1, qmass2, pdela
   common /comqfl/ zgsint
