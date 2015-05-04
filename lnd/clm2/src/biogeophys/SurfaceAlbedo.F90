#include <misc.h>
#include <preproc.h>

subroutine SurfaceAlbedo (clm, caldayp1, eccen, obliqr, lambm0, mvelpp)

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely
!  L                           M  available land surface process model.
!  M --COMMUNITY LAND MODEL--  C
!  C                           L
!  LMCLMCLMCLMCLMCLMCLMCLMCLMCLM
!
!-----------------------------------------------------------------------
! Purpose: 
! Surface albedo and two-stream fluxes
! 
! Method: 
! Surface albedos. Also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation. 
! Also sunlit fraction of the canopy. 
!
! The calling sequence is:
!   -> SurfaceAlbedo:     albedos for next time step
!        -> SnowAge:      snow age
!        -> SnowAlbedo:   snow albedos: direct beam
!        -> SnowAlbedo:   snow albedos: diffuse
!        -> SoilAlbedo:   soil/lake/glacier/wetland albedos
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (vis dif)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dir)
!        -> TwoStream:    absorbed, reflected, transmitted solar fluxes (nir dif)
! 
! Author: 
! Gordon Bonan
! April 2002: Vertenstein/Oleson/Levis; Final form
! 
!-----------------------------------------------------------------------
! $Id: SurfaceAlbedo.F90,v 1.1.10.3 2002/06/15 13:50:20 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use shr_orb_mod
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm  !CLM 1-D Module

  real(r8), intent(in) :: caldayp1 ! calendar day at Greenwich (1.00, ..., 365.99)
  real(r8), intent(in) :: eccen    ! Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   ! Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   ! Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   ! Earth's moving vernal equinox long. of perihelion + pi (radians)

!----Local Variables----------------------------------------------------

  integer  :: ib              ! band index
  integer  :: ic              ! direct beam: ic=0; diffuse: ic=1
  integer  :: nband = numrad  ! number of solar radiation wave bands
  real(r8) :: wl              ! fraction of LAI+SAI that is LAI
  real(r8) :: ws              ! fraction of LAI+SAI that is SAI
  real(r8) :: mpe = 1.e-06    ! prevents overflow for division by zero 
  real(r8) :: vai             ! elai+esai
  real(r8) :: rho(numrad)     ! leaf/stem refl weighted by fraction LAI and SAI
  real(r8) :: tau(numrad)     ! leaf/stem tran weighted by fraction LAI and SAI
  real(r8) :: ftdi(numrad)    ! down direct flux below veg per unit dif flux = 0
  real(r8) :: albsnd(numrad)  ! snow albedo (direct)
  real(r8) :: albsni(numrad)  ! snow albedo (diffuse)
  real(r8) :: gdir            ! aver projected leaf/stem area in solar direction
  real(r8) :: ext             ! optical depth direct beam per unit LAI+SAI
  real(r8) :: delta           ! solar declination angle in radians
  real(r8) :: eccf            ! earth orbit eccentricity factor
  real(r8) :: coszen          ! cosine solar zenith angle for next time step

!----End Variable List--------------------------------------------------

!
! Solar declination  for next time step
!

  call shr_orb_decl (caldayp1, eccen, mvelpp, lambm0, obliqr, &
                     delta   , eccf )

!
! Cosine solar zenith angle for next time step
!

  coszen = shr_orb_cosz(caldayp1, clm%lat, clm%lon, delta)

!
! Initialize output because solar radiation only done if coszen > 0
!

  do ib = 1, nband
     clm%albd(ib)   = 1.
     clm%albi(ib)   = 1.
     clm%albgrd(ib) = 0._r8
     clm%albgri(ib) = 0._r8
     clm%fabd(ib)   = 0._r8
     clm%fabi(ib)   = 0._r8
     clm%ftdd(ib)   = 0._r8
     clm%ftid(ib)   = 0._r8
     clm%ftii(ib)   = 0._r8
     if (ib==1) clm%fsun = 0.
  end do

!
! Return if coszen is not positive 
!

  if (coszen <= 0._r8) RETURN 

!
! Weight reflectance/transmittance by lai and sai
!

  do ib = 1, nband
     vai = clm%elai + clm%esai
     wl = clm%elai / max( vai,mpe )
     ws = clm%esai / max( vai,mpe )
     rho(ib) = max( clm%rhol(ib)*wl + clm%rhos(ib)*ws, mpe )
     tau(ib) = max( clm%taul(ib)*wl + clm%taus(ib)*ws, mpe )
  end do

!
! Snow albedos: only if h2osno > 0
!

  if ( clm%h2osno > 0._r8 ) then
     ic=0; call SnowAlbedo (clm, coszen, nband, ic, albsnd)
     ic=1; call SnowAlbedo (clm, coszen, nband, ic, albsni)  
  else
     albsnd(:) = 0._r8
     albsni(:) = 0._r8
  endif
     
!
! Ground surface albedos
!

  call SoilAlbedo (clm, coszen, nband, albsnd, albsni)      

  if (vai /= 0.) then  ! vegetated patch

!
! Loop over nband wavebands to calculate surface albedos and solar 
! fluxes for vegetated patch for unit incoming direct 
! (ic=0) and diffuse flux (ic=1)
!

     do ib = 1, nband
        ic = 0
        call TwoStream (clm,      ib,  ic,       coszen,   vai,      &
                        rho,      tau, clm%fabd, clm%albd, clm%ftdd, &
                        clm%ftid, gdir )
        ic = 1
        call TwoStream (clm,      ib,  ic,       coszen,   vai,  &
                        rho,      tau, clm%fabi, clm%albi, ftdi, &
                        clm%ftii, gdir )
     end do
     
!
! Sunlit fraction of canopy. Set fsun = 0 if fsun < 0.01.
!
     
     ext = gdir/coszen * sqrt(1.-rho(1)-tau(1))
     clm%fsun = (1.-exp(-ext*vai)) / max(ext*vai,mpe)
     ext = clm%fsun                                       !temporary fsun
     if (ext < 0.01) then 
        wl = 0._r8                                        !temporary fsun
     else
        wl = ext                                          !temporary fsun
     end if
     clm%fsun = wl

  else     ! non-vegetated patch

     do ib = 1,numrad
        clm%fabd(ib) = 0.
        clm%fabi(ib) = 0.
        clm%ftdd(ib) = 1.
        clm%ftid(ib) = 0.
        clm%ftii(ib) = 1.
        clm%albd(ib) = clm%albgrd(ib)
        clm%albi(ib) = clm%albgri(ib)
        clm%fsun     = 0.
     end do

  endif

  return
end subroutine SurfaceAlbedo
