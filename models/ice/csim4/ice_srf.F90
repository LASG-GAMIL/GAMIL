! Questions: Should ice model reduce snow thickness from evap?
! I am including it
! I am also adding snowfall to hsnow and limiting the snowfall to 0.5m

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute sea ice to atmosphere surface fluxes; then compute
! sea ice temperature change.
!
! Method: 
! Temperatures over sea-ice surfaces are specified in 'plevmx' layers of
! fixed thickness and thermal properties.  The forecast temperatures are
! determined from a backward/implicit diffusion calculation using
! linearized sensible/latent heat fluxes. The bottom ocean temperature
! is fixed at -2C, allowing heat flux exchange with underlying ocean.
! Temperature over sea ice is not allowed to exceed melting temperature.
! 
! The spectral components of shortwave and albedos must be sent to 
! this sea ice model because the sea ice extinction coefficient depends
! on the wavelength.
!
! Author: C.M. Bitz
! 
!-----------------------------------------------------------------------
  subroutine seaice (c, ncol, dtime, icefrac, Tsice,       &
                     hi, snowh, ubot, vbot, tbot,          &
                     qbot, thbot, zbot, pbot, flwds,       &
                     swvdr, swidr, swvdf, swidf, alvdr,    &
                     alidr, alvdf, alidf, snowfall, tssub, &
                     qflx, taux, tauy, ts, shflx,          &
                     lhflx, lwup, tref)

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols
  use constituents, only: pcnst, pnats
  use ice_constants
  use ice_sfc_flux
  use ice_tstm
  use ice_dh

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: c               ! chunk index
  integer , intent(in) :: ncol           ! number of columns this chunk

  real(r8), intent(in) :: dtime              ! Land/ocean/seaice flag
  real(r8), intent(in) :: icefrac(pcols)    ! Land/ocean/seaice flag
  real(r8), intent(inout) :: tsice(pcols)   ! ice/snow surface temperature (K)
  real(r8), intent(inout) :: snowh(pcols)   ! Snow depth (liquid water equivalent)
  real(r8), intent(inout) :: hi(pcols)      ! Ice thickness
  real(r8), intent(in) :: ubot(pcols)          ! Bottom level u wind
  real(r8), intent(in) :: vbot(pcols)          ! Bottom level v wind
  real(r8), intent(in) :: tbot(pcols)          ! Bottom level temperature
  real(r8), intent(in) :: qbot(pcols)          ! Bottom level specific humidity
  real(r8), intent(in) :: thbot(pcols)         ! Bottom level potential temperature
  real(r8), intent(in) :: zbot(pcols)          ! Bottom level height above surface
  real(r8), intent(in) :: pbot(pcols)        ! Bottom level pressure
  real(r8), intent(in) :: flwds(pcols)   ! net down longwave radiation at surface
  real(r8), intent(in) :: swvdr(pcols)   ! direct beam solar radiation onto srf (sw)
  real(r8), intent(in) :: swidr(pcols)   ! direct beam solar radiation onto srf (lw)
  real(r8), intent(in) :: swvdf(pcols)   ! diffuse solar radiation onto srf (sw)
  real(r8), intent(in) :: swidf(pcols)   ! diffuse solar radiation onto srf (lw)
  real(r8), intent(in) :: snowfall(pcols)  ! total snow rate (m h2o/s) 
  real(r8), intent(inout) :: alvdr(pcols)   ! ocean + ice albedo: shortwave, direct
  real(r8), intent(inout) :: alvdf(pcols)   ! ocean + ice albedo: shortwave, diffuse
  real(r8), intent(inout) :: alidr(pcols)   ! ocean + ice albedo: longwave, direct
  real(r8), intent(inout) :: alidf(pcols)   ! ocean + ice albedo: longwave, diffuse

  real(r8), intent(inout):: tssub(pcols,plevmx)  ! Surface/sub-surface temperatures

! fluxes/quantities summed over surface types
  real(r8), intent(out):: qflx(pcols,pcnst+pnats)    ! Constituent flux (kg/m2/s)
  real(r8), intent(out):: taux(pcols)          ! X surface stress (N/m2)
  real(r8), intent(out):: tauy(pcols)          ! Y surface stress (N/m2)
  real(r8), intent(out):: ts(pcols)            ! surface temperature (K)
  real(r8), intent(out):: shflx(pcols)         ! Surface sensible heat flux (J/m2/s)
  real(r8), intent(out):: lhflx(pcols)         ! Surface latent   heat flux (J/m2/s)
  real(r8), intent(out):: lwup(pcols)          ! surface longwave up flux (W/m2)
  real(r8), intent(out):: tref(pcols)          ! 2m reference temperature

!---------------------------Local variables-----------------------------
! fluxes/quantities over sea ice only
  real(r8) :: tauxice(pcols)       ! X surface stress (N/m2)
  real(r8) :: tauyice(pcols)       ! Y surface stress (N/m2)
  real(r8) :: shflxice(pcols)      ! Surface sensible heat flux (J/m2/s)
  real(r8) :: lhflxice(pcols)      ! Surface latent   heat flux (J/m2/s)
  real(r8) :: trefice(pcols)       ! 2m reference temperature
  real(r8) :: flwup(pcols)
  real(r8) :: evap(pcols)

  real(r8) :: Fnet
  real(r8) :: condb, swbot
  real(r8) :: Flwdabs,Fswabs,Fswabsv,Fswabsi
  real(r8) :: dflhdT,dfshdT,dflwdT
  real(r8) :: Tiz(0:plevmx) ! local 1D ice temperature profile (C)
  real(r8) :: Tsfc          ! local ice surface temperature (C)
  real(r8) :: Tbasal        ! ice bottom temp (C)
  real(r8) :: asnow ! snow fractional coverage
  real(r8) :: hs   ! snow thickness 
  real(r8) :: dhs  ! change in snow thickness from melting
  real(r8) :: subs ! change in snow thickness from sublimation/condensation

  integer :: npts ! number of gridcells with sea ice
  integer :: linpts ! counter for number of ice points reset to linear profile
  integer :: indx(pcols)

! Sea ice thickness is fixed, so
! bottom melt/growth is not computed and the ice-ocean flux 
! is ignored. A flag is set in ice_dh.F accordingly
  real(r8), parameter :: Fbot = 0.

! dummies sent to melt/grow routine but ignored when ice is fixed thickness
  real(r8)  :: Focn
  real(r8)  :: dhib, dhit, subi, dhif, dhsf
  real(r8)  :: qi(plevmx)

  integer  :: i,k,m,ii
!-----------------------------------------------------------------------
!  call t_startf ('seaice_other')
  npts = 0
  do i=1,ncol
     if (icefrac(i)>0) then
!        write(6,*) '(ice_srf)',c,i,icefrac(i),hi(i)
!        write(6,*) '(ice_srf) tssub:', i,c,(tssub(i,k),k=1,plevmx)
        npts = npts + 1
        indx(npts) = i
     else
        Tsice(i) = TfrezK
     end if
  end do

  flwup(:)=0.
  shflxice(:)=0.
  lhflxice(:)=0.
  tauxice(:)=0.
  tauyice(:)=0.
!  call t_stopf ('seaice_other')

  if (npts.gt.0) then
  linpts=0
  do ii=1,npts
!     call t_startf ('seaice_other')
     i = indx(ii)
     ! Convert temperatures to C, use 1D array
     Tsfc = Tsice(i) - Tffresh
     do k=1,plevmx
        Tiz(k) = tssub(i,k) - Tffresh
        Tiz(k) = min(Tmelz(k),Tiz(k))
     end do
     ! snow temperature is diagnostic so init to surf
     Tiz(0) = Tsfc
     
     ! snow lands on ice no matter what its temperature
     if (prognostic_icesnow) then
        hs =  (snowh(i) + snowfall(i)*dtime)*rhofresh/rhos 
        if (fixice) hs = min(hs,0.5)  ! do not let the snow get out of hand
     else
        hs =  snowh(i)*rhofresh/rhos 
     endif
     
     ! 1 - snow covered area fraction
     asnow = c1-hs/(hs + snowpatch)
     
     !-----------------------------------------------------------------
     ! compute air to ice heat, momentum, radiative and water fluxes 
     !-----------------------------------------------------------------
!        write(6,*) '(ice_srf) T in C',i,c,Tsfc,(Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) b4 ice_sfc_flux',Tsfc, ubot(i), vbot(i), tbot(i), &
!             qbot(i)    ,thbot(i)   ,zbot(i)  ,pbot(i)  , flwds(i) ,&
!             swvdr(i)   ,swidr(i)   ,swvdf(i) ,swidf(i) , &
!             alvdr(i)   ,alidr(i)   ,alvdf(i) ,alidf(i)
!        write(6,*) '(ice_srf) ',c,i,alvdr(i),alidr(i)   ,alvdf(i)   ,alidf(i)

!     call t_stopf ('seaice_other')
!     call t_startf ('ice_atm_flux')
     call ice_atm_flux(Tsfc, ubot(i), vbot(i), tbot(i), &
             qbot(i)    ,thbot(i)   ,zbot(i)    ,pbot(i)  , flwds(i) ,&
             swvdr(i)   ,swidr(i)   ,swvdf(i)   ,swidf(i) , &
             alvdr(i)   ,alidr(i)   ,alvdf(i)   ,alidf(i) , &
             tauxice(i) ,tauyice(i) ,flwup(i), &
             shflxice(i),lhflxice(i),trefice(i), &
             Flwdabs,Fswabs,Fswabsv,Fswabsi, &
             dflhdT,dfshdT,dflwdT)
!     call t_stopf ('ice_atm_flux')

!        write(6,*) '(ice_srf) after ice_sfc_flux', &
!             tauxice(i) ,tauyice(i) ,flwup(i), &
!             shflxice(i),lhflxice(i),trefice(i), &
!             Flwdabs,Fswabs,Fswabsv,Fswabsi, &
!             dflhdT,dfshdT,dflwdT

      !-----------------------------------------------------------------
      ! solve heat equation
      !-----------------------------------------------------------------
!        write(6,*) '(ice_srf) Tiz', (Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) sensible heat flux',shflxice(i)
!        write(6,*) '(ice_srf) latent heat flux',lhflxice(i)

!     call t_startf ('tstm')
     call tstm( dtime, Tmelz, saltz, Tfrez &
                     , icefrac(i), hi(i), hs &
                     , fswabs, fswabsv, fswabsi &
                     , flwdabs, dflwdT, dflhdT, dfshdT &
                     , asnow,  Tbasal &
                     , swbot, Fnet, condb, Tsfc, Tiz &
                     , flwup(i), lhflxice(i), shflxice(i),linpts)
!     call t_stopf ('tstm')

!        write(6,*) '(ice_srf) Tiz', (Tiz(k),k=0,plevmx)
!        write(6,*) '(ice_srf) sensible heat flux',shflxice(i)
!        write(6,*) '(ice_srf) latent heat flux',lhflxice(i)
      !-----------------------------------------------------------------
      ! compute snow melt and sublimation
      !-----------------------------------------------------------------

!     call t_startf ('dh')
     call dh  ( dtime, saltz, Tiz     &
                     , Tbasal, hi(i), hs, Fbot &
                     , Fnet, condb, lhflxice(i) &
                     , dhib, dhit, dhs, subi &
                     , subs, dhif, dhsf, qi, Focn, i,c)
!     call t_stopf ('dh')

!
! If we are not fixing ice thickness then adjust height 
!
!     call t_startf ('seaice_other')
     if (.not. fixice) then
        hi(i)=hi(i)+subi+dhit+dhib+dhif
     end if
!
! If we are prognosing snow then adjust the height due to snow melt and
! sublimation
!
     if (prognostic_icesnow) then
        if (.not. fixice) then
           hs = hs + subs + dhs + dhsf
        else 
           hs = hs  + subs + dhs 
        endif
     endif

     evap(i) = rhos*subs/dtime 
     snowh(i) = hs*rhos/rhofresh

     ! Convert temperatures to K, filling 2D arrays
     Tsice(i) = Tsfc + Tffresh
     do k=1,plevmx
        tssub(i,k) = Tiz(k) + Tffresh
     end do
!     call t_stopf ('seaice_other')
     
  end do
  if (linpts.gt.0) then
     write(6,*)'WARNING: ice_tstm ::profile reset ',linpts,&
          ' points at chunck ',c,' see NOTE in ice_tstm.F for more info'
  end if
!  call t_startf ('seaice_other')
! Update fluxes, sum ice percentages into total flux values 
! ALL FLUXES IN ICE MODEL ARE POSITIVE DOWN
! CCM DOES NOT USE THIS STANDARD 
  do ii=1,npts
     i = indx(ii)
     Ts(i) = Tsice(i)
     tref(i)=trefice(i)
     lwup(i) = -1.*(flwup(i)-(1-emissivity)*flwds(i))
     shflx(i)=-1.*shflxice(i)
     lhflx(i)=-1.*lhflxice(i)
     taux(i)=-1.*tauxice(i)
     tauy(i)=-1.*tauyice(i)
     qflx(i,1) =-1.*evap(i)
  end do
!  call t_stopf ('seaice_other')

  endif

! Set non-water constituent fluxes to zero

!  call t_startf ('seaice_other')
  do m=2,pcnst+pnats
     do ii=1,npts
        i = indx(ii)
        qflx(i,m) = 0.
     end do
  end do
!  call t_stopf ('seaice_other')

  return

end subroutine seaice


