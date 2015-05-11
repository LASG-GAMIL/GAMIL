c=======================================================================
!---! atmospheric boundary interface (stability based flux calculations)
!---!
!---! based on code by Bill Large
c=======================================================================

      module ice_sfc_flux

      use ice_constants
      use shr_const_mod, only: SHR_CONST_ZVIR, SHR_CONST_RDAIR

      implicit none

c=======================================================================

      contains

c=======================================================================

      subroutine ice_atm_flux(Tsf, uatm, vatm, tair,
     &    qa    ,potT   ,zlvl    ,pbot  , Flw ,
     &    swvdr   ,swidr   ,swvdf   ,swidf , 
     &    alvdr   ,alidr   ,alvdf   ,alidf , 
     &    strx, stry, flwup, fsh, flh, Tref,
     &    flwdabs,fswabs,fswabsv,fswabsi,
     &    dflhdT,dfshdT,dflwdT)

!---!-------------------------------------------------------------------
!---! Compute ice-atm surface fluxes, stress, and reference 
!---! temperature
!---! NOTE: 
!---!  o all fluxes are positive downward
!---!  o net heat flux = fswabs + flwup + (1-emissivity)flwdn + fsh + flh
!---!  o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!---!  o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!---!
!---! ASSUME:
!---!  o The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!---!
!---! code originally based on CSM1
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   Tsf    ! surface temperature of ice 
     &,  uatm  ! surface u wind
     &,  vatm  ! surface v wind
     &,  tair        ! Bottom level temperature
     &,  qa          ! Bottom level specific humidity
     &,  potT         ! Bottom level potential temperature
     &,  zlvl          ! Bottom level height above surface
     &,  pbot        ! Bottom level pressure
     &,  Flw     ! net down longwave radiation at surface
     &,  swvdr   ! direct beam solar radiation onto srf (sw)
     &,  swidr   ! direct beam solar radiation onto srf (lw)
     &,  swvdf   ! diffuse solar radiation onto srf (sw)
     &,  swidf   ! diffuse solar radiation onto srf (lw)
     &,  alvdr   ! ocean + ice albedo: shortwave, direct
     &,  alvdf   ! ocean + ice albedo: shortwave, diffuse
     &,  alidr   ! ocean + ice albedo: longwave, direct
     &,  alidf   ! ocean + ice albedo: longwave, diffuse

      real (kind=dbl_kind), intent(out) ::
     &   strx   ! x surface stress (N)
     &,  stry   ! y surface stress (N)
     &,  Tref    ! reference height temperature  (K)
     &,    flwdabs   ! down long-wave  absorbed heat flx   (W/m**2)
     &,    flwup     ! emitted long-wave upward heat flux  (W/m**2)
     &,    fswabs     ! fswabs sum
     &,    fswabsv    ! fswabs in vis (wvlngth < 700nm)  (W/m**2)
     &,    fswabsi    ! fswabs in nir (wvlngth > 700nm)  (W/m**2)
     &,    fsh      ! sensible         heat flux  (W/m**2)
     &,    flh      ! latent           heat flux  (W/m**2)
     &,    dflhdT     ! d(flh)/d(T)      (W/m**2/K)
     &,    dfshdT     ! d(fsh)/d(T)      (W/m**2/K)
     &,    dflwdT    ! d(flwup)/d(T)     (W/m**2/K)


! local variables
      real (kind=dbl_kind) :: 
     &   dssqdt ! derivative of ssq wrt Ti (kg/kg/K)
     &,  delt   ! potential T difference   (K)
     &,  delq   ! humidity difference      (kg/kg)
 
      integer (kind=int_kind) :: k      ! iteration index
      real (kind=dbl_kind) :: 
     &   TsfK   ! surface temperature in Kelvin (K)
     &,  thva   ! virtual temperature      (K)
     &,  stable ! stability factor
     &,  rdn    ! sqrt of neutral exchange coefficient (momentum)
     &,  rhn    ! sqrt of neutral exchange coefficient (heat)
     &,  ren    ! sqrt of neutral exchange coefficient (water)
     &,  hol    ! H (at zlvl  ) over L
     &,  xsq    ! temporary variable
     &,  xqq    ! temporary variable
     &,  psimh  ! stability function at zlvl   (momentum)
     &,  psixh  ! stability function at zlvl   (heat and water)
     &,  alz    ! ln(zlvl  /z10)
     &,  tau    ! stress at zlvl
     &,  bn     ! exchange coef funct for interpolation
     &,  bh     ! exchange coef funct for interpolation
     &,  fac    ! interpolation factor
     &,  ln0    ! log factor for interpolation
     &,  ln3    ! log factor for interpolation
     &,  ustar  ! ustar (m/s)
     &,  tstar  ! tstar
     &,  qstar  ! qstar
     &,  rd     ! sqrt of exchange coefficient (momentum)
     &,  re     ! sqrt of exchange coefficient (water)            
     &,  rh     ! sqrt of exchange coefficient (heat)
     &,  vmag   ! surface wind magnitude   (m/s)
     &,  ssq    ! sat surface humidity     (kg/kg)
     &,  cp     ! specific heat of moist air
     &,  rhoa   ! air density (kg/m**3)
     &,  lhcoef
     &,  shcoef
     &,  wind  ! surface wind speed 

      real (kind=dbl_kind), parameter ::
     &   cpvir = cpwv/cp_air - c1  ! Defined as cpwv/cp_air - 1.
     &,  zTref  = c2          ! reference height for air temperature (m)
     &,  umin  = c1          ! minimum wind speed (m/s)
     &,  zvir  = SHR_CONST_ZVIR    ! rh2o/rair - 1.0
!     &,  qqq   = 640380._dbl_kind  ! for qsat, dqsatdt
!     &,  TTT   = 5107.4_dbl_kind   ! for qsat, dqsatdt
     &,  qqq   = 11637800._dbl_kind  ! for qsat, dqsatdt
     &,  TTT   = 5897.8_dbl_kind   ! for qsat, dqsatdt

      ! local functions
      real (kind=dbl_kind) :: 
     &   Tk      ! temperature (K)
     &,  qsat    ! the saturation humididty of air (kg/m**3)
     &,  dqsatdt ! derivative of qsat wrt surface temperature
     &,  xd      ! dummy argument  
     &,  psimhu  ! unstable part of psimh
     &,  psixhu  ! unstable part of psimx

      qsat(Tk)    = qqq / exp(TTT/Tk)

      dqsatdt(Tk) = (TTT / Tk**2) * qqq / exp(TTT/Tk)

      psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8)
     $              - c2*atan(xd) + 1.571_dbl_kind

      psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

      ! define some needed variables
      TsfK   = Tsf +Tffresh                 !  surface temp (K)
      wind   = sqrt(uatm*uatm+vatm*vatm)
      vmag   = max(umin, wind)
      rhoa   = pbot/(SHR_CONST_RDAIR*tair)
      thva   = potT * (c1 + zvir * Qa)      ! virtual pot temp (K)
      ssq    = qsat   (TsfK) / rhoa         ! sat surf hum (kg/kg)
      dssqdt = dqsatdt(TsfK) / rhoa         ! deriv of ssq wrt Ti 
      delt   = potT - TsfK                  ! pot temp diff (K)
      delq   = Qa - ssq                     ! spec hum dif (kg/kg)
      alz    = log(zlvl/zref) 
      cp     = cp_air*(c1 + cpvir*ssq)
ccc        write(6,*) 'cp is the problem',cp,ssq,qsat(TsfK),TsfK,rhoa, 
ccc       & Tsf,Tffresh
ccc  
ccc        write(6,*) 'IN ice_sfc_flux', Tsf, uatm, vatm, tair,
ccc       &    qa    ,potT   ,zlvl    ,pbot  , Flw ,
ccc       &    swvdr   ,swidr   ,swvdf   ,swidf , 
ccc       &    alvdr   ,alidr   ,alvdf   ,alidf 
      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      ! neutral coefficients, z/L = 0.0 
      rdn = vonkar/log(zref/iceruf)
      rhn = rdn
      ren = rdn

      ! ustar,tstar,qstar
      ustar = rdn * vmag
      tstar = rhn * delt  
      qstar = ren * delq  

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k=1,5

        ! compute stability & evaluate all stability functions 
        hol    = vonkar * gravit * zlvl
     $           * (tstar/thva+qstar/(c1/zvir+Qa)) / ustar**2
        hol    = sign( min(abs(hol),c10), hol )
        stable = p5 + sign(p5 , hol)
        xsq    = max(sqrt(abs(c1 - c16*hol)) , c1)
        xqq    = sqrt(xsq)
        psimh  = -c5*hol*stable + (c1-stable)*psimhu(xqq)
        psixh  = -c5*hol*stable + (c1-stable)*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
        rd = rdn / (c1+rdn/vonkar*(alz-psimh))
        rh = rhn / (c1+rhn/vonkar*(alz-psixh))
        re = ren / (c1+ren/vonkar*(alz-psixh))

        ! update ustar, tstar, qstar using updated, shifted coeffs 
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

      enddo    ! end iteration

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

      shcoef = rhoa*ustar*cp  *rh
      lhcoef = rhoa*ustar*Lsub*re

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa * ustar * ustar 
      ! strx = tau * uatm / vmag 
      ! stry = tau * vatm / vmag 
      !------------------------------------------------------------

      tau = rhoa * ustar * rd    ! not the stress at zlvl
      strx = tau * uatm 
      stry = tau * vatm

      !------------------------------------------------------------
      ! reference temperature interpolation
      !------------------------------------------------------------
      ! Assume that 
      ! cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
      ! 1/sqrt(cn)=1/rdn and sqrt(cm)/ch=1/rh 
      !------------------------------------------------------------
      bn = vonkar/rdn
      bh = vonkar/rh

      ! Interpolation factor for stable and unstable cases
      ln0 = log(c1 + (zTref/zlvl)*(exp(bn) - c1))
      ln3 = log(c1 + (zTref/zlvl)*(exp(bn - bh) - c1))
      fac = (ln0 - zTref/zlvl*(bn - bh))/bh 
     &         * stable
     $    + (ln0 - ln3)/bh * (c1-stable)
      fac = min(max(fac,c0),c1)

      Tref = TsfK + (Tair - TsfK)*fac

      ! shortwave radiative flux
      fswabsv  = swvdr*(c1-alvdr) + swvdf*(c1-alvdf)
      fswabsi  = swidr*(c1-alidr) + swidf*(c1-alidf)
      fswabs   = fswabsv + fswabsi

      ! longwave radiative flux
      flwdabs = emissivity*Flw
      flwup   = -emissivity*stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      flh = lhcoef * delq
      fsh = shcoef * delt

!      write(6,*) '(ice_sfc_flux)',fsh,shcoef,delt

      ! derivatives wrt surface temp
      dflwdT = - emissivity*stefan_boltzmann * c4*TsfK**3 
      dflhdT = - lhcoef * dssqdt
      dfshdT = - shcoef

!! NOTE
! if the temperature dependence of flh, fsh and flwup are
! included in the iteration for the ice sfc temperature, move
! fsw* out of this routine and use the following expressions instead:
! (NONE of the intent(in) variables will then be needed)

c**        Qsfc = Qcoef * exp(22.47*(c1-Tffresh/TsfK))
c**        fsh  = shcoef*(potT - TsfK)
c**        flh  = lhcoef*(Qa - Qsfc)

c**        dfshdT = - shcoef
c**        dflhdT = - lhcoef*lvrrv*Qsfc/TsfK**2

      end subroutine ice_atm_flux

c=======================================================================

      end module ice_sfc_flux

c=======================================================================

