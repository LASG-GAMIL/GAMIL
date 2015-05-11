subroutine raddedmx(coszrs  ,ndayc   ,idayc   ,abh2o   , &
                    abo3    ,abco2   ,abo2    ,uh2o    ,uo3     , &
                    uco2    ,uo2     ,trayoslp,pflx    ,ns      , &
                    tauxcl  ,wcl     ,gcl     ,fcl     ,tauxci  , &
                    wci     ,gci     ,fci     ,tauxar  ,wa      , &
                    ga      ,fa      ,rdir    ,rdif    ,tdir    , &
                    tdif    ,explay  ,rdirc   ,rdifc   ,tdirc   , &
                    tdifc   ,explayc )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes layer reflectivities and transmissivities, from the top down
! to the surface using the delta-Eddington solutions for each layer
! 
! Method: 
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! Modified for maximum/random cloud overlap by Bill Collins and John
!    Truesdale
! 
! Author: Bill Collins
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid

   implicit none

   integer nspint           ! Num of spctrl intervals across solar spectrum

   parameter ( nspint = 19 )
!
! Minimum total transmission below which no layer computation are done:
!
   real(r8) trmin                ! Minimum total transmission allowed
   real(r8) wray                 ! Rayleigh single scatter albedo
   real(r8) gray                 ! Rayleigh asymetry parameter
   real(r8) fray                 ! Rayleigh forward scattered fraction

   parameter (trmin = 1.e-3)
   parameter (wray = 0.999999)
   parameter (gray = 0.0)
   parameter (fray = 0.1)
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: coszrs(pcols)        ! Cosine zenith angle
   real(r8), intent(in) :: trayoslp             ! Tray/sslp
   real(r8), intent(in) :: pflx(pcols,0:pverp)  ! Interface pressure
   real(r8), intent(in) :: abh2o                ! Absorption coefficiant for h2o
   real(r8), intent(in) :: abo3                 ! Absorption coefficiant for o3
   real(r8), intent(in) :: abco2                ! Absorption coefficiant for co2
   real(r8), intent(in) :: abo2                 ! Absorption coefficiant for o2
   real(r8), intent(in) :: uh2o(pcols,0:pver)   ! Layer absorber amount of h2o
   real(r8), intent(in) :: uo3(pcols,0:pver)    ! Layer absorber amount of  o3
   real(r8), intent(in) :: uco2(pcols,0:pver)   ! Layer absorber amount of co2
   real(r8), intent(in) :: uo2(pcols,0:pver)    ! Layer absorber amount of  o2
   real(r8), intent(in) :: tauxcl(pcols,0:pver) ! Cloud extinction optical depth (liquid)
   real(r8), intent(in) :: wcl(pcols,0:pver)    ! Cloud single scattering albedo (liquid)
   real(r8), intent(in) :: gcl(pcols,0:pver)    ! Cloud asymmetry parameter (liquid)
   real(r8), intent(in) :: fcl(pcols,0:pver)    ! Cloud forward scattered fraction (liquid)
   real(r8), intent(in) :: tauxci(pcols,0:pver) ! Cloud extinction optical depth (ice)
   real(r8), intent(in) :: wci(pcols,0:pver)    ! Cloud single scattering albedo (ice)
   real(r8), intent(in) :: gci(pcols,0:pver)    ! Cloud asymmetry parameter (ice)
   real(r8), intent(in) :: fci(pcols,0:pver)    ! Cloud forward scattered fraction (ice)
   real(r8), intent(in) :: tauxar(pcols,0:pver) ! Aerosol extinction optical depth
   real(r8), intent(in) :: wa(pcols,0:pver)     ! Aerosol single scattering albedo
   real(r8), intent(in) :: ga(pcols,0:pver)     ! Aerosol asymmetry parameter
   real(r8), intent(in) :: fa(pcols,0:pver)     ! Aerosol forward scattered fraction

   integer, intent(in) :: ndayc                 ! Number of daylight columns
   integer, intent(in) :: idayc(pcols)          ! Daylight column indices
   integer, intent(in) :: ns                    ! Index of spectral interval
!
! Input/Output arguments
!
! Following variables are defined for each layer; 0 refers to extra
! layer above top of model:
!
   real(r8), intent(inout) :: rdir(nspint,pcols,0:pver)   ! Layer reflectivity to direct rad
   real(r8), intent(inout) :: rdif(nspint,pcols,0:pver)   ! Layer reflectivity to diffuse rad
   real(r8), intent(inout) :: tdir(nspint,pcols,0:pver)   ! Layer transmission to direct rad
   real(r8), intent(inout) :: tdif(nspint,pcols,0:pver)   ! Layer transmission to diffuse rad
   real(r8), intent(inout) :: explay(nspint,pcols,0:pver) ! Solar beam exp transm for layer
!
! Corresponding quantities for clear-skies
!
   real(r8), intent(inout) :: rdirc(nspint,pcols,0:pver)  ! Clear layer reflec. to direct rad
   real(r8), intent(inout) :: rdifc(nspint,pcols,0:pver)  ! Clear layer reflec. to diffuse rad
   real(r8), intent(inout) :: tdirc(nspint,pcols,0:pver)  ! Clear layer trans. to direct rad
   real(r8), intent(inout) :: tdifc(nspint,pcols,0:pver)  ! Clear layer trans. to diffuse rad
   real(r8), intent(inout) :: explayc(nspint,pcols,0:pver)! Solar beam exp transm clear layer
!
!---------------------------Local variables-----------------------------
!
   integer i                 ! Column indices
   integer k                 ! Level index
   integer nn                ! Index of column loops (max=ndayc)

   real(r8) taugab(pcols)        ! Layer total gas absorption optical depth
   real(r8) tauray(pcols)        ! Layer rayleigh optical depth
   real(r8) taucsc               ! Layer cloud scattering optical depth
   real(r8) tautot               ! Total layer optical depth
   real(r8) wtot                 ! Total layer single scatter albedo
   real(r8) gtot                 ! Total layer asymmetry parameter
   real(r8) ftot                 ! Total layer forward scatter fraction
   real(r8) wtau                 !  rayleigh layer scattering optical depth
   real(r8) wt                   !  layer total single scattering albedo
   real(r8) ts                   !  layer scaled extinction optical depth
   real(r8) ws                   !  layer scaled single scattering albedo
   real(r8) gs                   !  layer scaled asymmetry parameter
!
!---------------------------Statement functions-------------------------
!
! Statement functions and other local variables
!
   real(r8) alpha                ! Term in direct reflect and transmissivity
   real(r8) gamma                ! Term in direct reflect and transmissivity
   real(r8) el                   ! Term in alpha,gamma,n,u
   real(r8) taus                 ! Scaled extinction optical depth
   real(r8) omgs                 ! Scaled single particle scattering albedo
   real(r8) asys                 ! Scaled asymmetry parameter
   real(r8) u                    ! Term in diffuse reflect and
!    transmissivity
   real(r8) n                    ! Term in diffuse reflect and
!    transmissivity
   real(r8) lm                   ! Temporary for el
   real(r8) ne                   ! Temporary for n
   real(r8) w                    ! Dummy argument for statement function
   real(r8) uu                   ! Dummy argument for statement function
   real(r8) g                    ! Dummy argument for statement function
   real(r8) e                    ! Dummy argument for statement function
   real(r8) f                    ! Dummy argument for statement function
   real(r8) t                    ! Dummy argument for statement function
   real(r8) et                   ! Dummy argument for statement function
!
! Intermediate terms for delta-eddington solution
!
   real(r8) alp                  ! Temporary for alpha
   real(r8) gam                  ! Temporary for gamma
   real(r8) ue                   ! Temporary for u
   real(r8) arg                  ! Exponential argument
   real(r8) extins               ! Extinction
   real(r8) amg                  ! Alp - gam
   real(r8) apg                  ! Alp + gam
!
   alpha(w,uu,g,e) = .75_r8*w*uu*((1._r8 + g*(1._r8-w))/(1._r8 - e*e*uu*uu))
   gamma(w,uu,g,e) = .50_r8*w*((3._r8*g*(1._r8-w)*uu*uu + 1._r8)/(1._r8-e*e*uu*uu))
   el(w,g)         = sqrt(3._r8*(1._r8-w)*(1._r8 - w*g))
   taus(w,f,t)     = (1._r8 - w*f)*t
   omgs(w,f)       = (1._r8 - f)*w/(1._r8 - w*f)
   asys(g,f)       = (g - f)/(1._r8 - f)
   u(w,g,e)        = 1.5_r8*(1._r8 - w*g)/e
   n(uu,et)        = ((uu+1._r8)*(uu+1._r8)/et ) - ((uu-1._r8)*(uu-1._r8)*et)
!
!-----------------------------------------------------------------------
!
! Compute layer radiative properties
!
! Compute radiative properties (reflectivity and transmissivity for
!    direct and diffuse radiation incident from above, under clear
!    and cloudy conditions) and transmission of direct radiation
!    (under clear and cloudy conditions) for each layer.
!
   do k=0,pver
      do nn=1,ndayc
         i=idayc(nn)
            tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
            taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) + abco2*uco2(i,k) + abo2*uo2(i,k)
            tautot = tauxcl(i,k) + tauxci(i,k) + tauray(i) + taugab(i) + tauxar(i,k)
            taucsc = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k) + tauxar(i,k)*wa(i,k)
            wtau   = wray*tauray(i)
            wt     = wtau + taucsc
            wtot   = wt/tautot
            gtot   = (wtau*gray + gcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + gci(i,k)*wci(i,k)*tauxci(i,k) + ga(i,k) *wa(i,k) *tauxar(i,k))/wt
            ftot   = (wtau*fray + fcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + fci(i,k)*wci(i,k)*tauxci(i,k) + fa(i,k) *wa(i,k) *tauxar(i,k))/wt
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
!
!     Limit argument of exponential to 25, in case lm very large:
!
            arg  = min(lm*ts,25._r8)
            extins = exp(-arg)
            ne = n(ue,extins)
            rdif(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
            tdif(ns,i,k)   =   4._r8*ue/ne
!
!     Limit argument of exponential to 25, in case coszrs is very small:
!
            arg       = min(ts/coszrs(i),25._r8)
            explay(ns,i,k) = exp(-arg)
            apg = alp + gam
            amg = alp - gam
            rdir(ns,i,k) = amg*(tdif(ns,i,k)*explay(ns,i,k)-1._r8) + apg*rdif(ns,i,k)
            tdir(ns,i,k) = apg*tdif(ns,i,k) + (amg*rdif(ns,i,k)-(apg-1._r8))*explay(ns,i,k)
!
!     Under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
            rdir(ns,i,k) = max(rdir(ns,i,k),0.0_r8)
            tdir(ns,i,k) = max(tdir(ns,i,k),0.0_r8)
            rdif(ns,i,k) = max(rdif(ns,i,k),0.0_r8)
            tdif(ns,i,k) = max(tdif(ns,i,k),0.0_r8)
!
!     Clear-sky calculation
!
            if (tauxcl(i,k) == 0.0_r8 .and. tauxci(i,k) == 0.0_r8) then

               rdirc(ns,i,k) = rdir(ns,i,k)
               tdirc(ns,i,k) = tdir(ns,i,k)
               rdifc(ns,i,k) = rdif(ns,i,k)
               tdifc(ns,i,k) = tdif(ns,i,k)
               explayc(ns,i,k) = explay(ns,i,k)
            else
               tautot = tauray(i) + taugab(i) + tauxar(i,k)
               taucsc = tauxar(i,k)*wa(i,k)
!
! wtau already computed for all-sky
!
               wt     = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + ga(i,k)*wa(i,k)*tauxar(i,k))/wt
               ftot   = (wtau*fray + fa(i,k)*wa(i,k)*tauxar(i,k))/wt
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
!
!     Limit argument of exponential to 25, in case lm very large:
!
               arg  = min(lm*ts,25._r8)
               extins = exp(-arg)
               ne = n(ue,extins)
               rdifc(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
               tdifc(ns,i,k)   =   4._r8*ue/ne
!
!     Limit argument of exponential to 25, in case coszrs is very small:
!
               arg       = min(ts/coszrs(i),25._r8)
               explayc(ns,i,k) = exp(-arg)
               apg = alp + gam
               amg = alp - gam
               rdirc(ns,i,k) = amg*(tdifc(ns,i,k)*explayc(ns,i,k)-1._r8)+ &
                               apg*rdifc(ns,i,k)
               tdirc(ns,i,k) = apg*tdifc(ns,i,k) + (amg*rdifc(ns,i,k) - (apg-1._r8))* &
                               explayc(ns,i,k)
!
!     Under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
               rdirc(ns,i,k) = max(rdirc(ns,i,k),0.0_r8)
               tdirc(ns,i,k) = max(tdirc(ns,i,k),0.0_r8)
               rdifc(ns,i,k) = max(rdifc(ns,i,k),0.0_r8)
               tdifc(ns,i,k) = max(tdifc(ns,i,k),0.0_r8)
            end if
         end do
   end do

   return
end subroutine raddedmx
