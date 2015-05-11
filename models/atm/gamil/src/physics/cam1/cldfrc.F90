#include <misc.h>
#include <params.h>
subroutine cldfrc(lchnk   ,ncol    , &
                  pmid    ,temp    ,q       ,omga    , &
                  cldtop  ,cldbot  ,cloud   ,clc     ,pdel    , &
                  cmfmc   ,landfrac,snowh   ,concld  ,cldst   , &
                  ts      ,ps      ,zdu     ,ocnfrac ,&
                  rhu00   ,relhum  ,dindex )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud fraction using scheme of J.M.Slingo,
! as modified by J.J.Hack and J.T.Kiehl
! 
! Method: 
! This scheme is based on the operational scheme used in the ECMWF model
! A full description of its development can be found in Slingo (1987),
! which appears in the QJRMS July issue.  A number of modifications have
! been introduced to the original scheme in the following implementation
! 
! Author: J. Hack
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use physconst, only: cappa, gravit, rair
   use cldconst
   use wv_saturation, only: aqsat
   use dycore

   implicit none

   real(r8), parameter :: pnot = 1.e5                 ! reference pressure
!
! Arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier
   integer, intent(in) :: ncol                   ! number of atmospheric columns
   integer, intent(in) :: dindex                 ! 0 or 1 to perturb rh

   real(r8), intent(in) :: pmid(pcols,pver)      ! midpoint pressures
   real(r8), intent(in) :: temp(pcols,pver)      ! temperature
   real(r8), intent(in) :: q(pcols,pver)         ! specific humidity
   real(r8), intent(in) :: omga(pcols,pver)      ! vertical pressure velocity
   real(r8), intent(in) :: cldtop(pcols)         ! top level of convection
   real(r8), intent(in) :: cldbot(pcols)         ! bottom level of convection
   real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
   real(r8), intent(in) :: snowh(pcols)          ! snow depth (liquid water equivalent)
   real(r8), intent(in) :: pdel(pcols,pver)      ! pressure depth of layer
   real(r8), intent(in) :: landfrac(pcols)       ! Land fraction
   real(r8), intent(in) :: ocnfrac(pcols)        ! Ocean fraction
   real(r8), intent(in) :: ts(pcols)             ! surface temperature
   real(r8), intent(in) :: ps(pcols)             ! surface pressure
   real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection

!
! Output arguments
!
   real(r8), intent(out) :: cloud(pcols,pver)     ! cloud fraction
   real(r8), intent(out) :: clc(pcols)            ! column convective cloud amount
   real(r8), intent(out) :: cldst(pcols,pver)     ! cloud fraction
   real(r8), intent(out) :: rhu00(pcols,pver)     ! RH threshold for cloud
   real(r8), intent(out) :: relhum(pcols,pver)    ! RH 
!      real(r8) dmudp                 ! measure of mass detraining in a layer
!
!---------------------------Local workspace-----------------------------
!
   real(r8) concld(pcols,pver)    ! convective cloud cover
   real(r8) cld                   ! intermediate scratch variable (low cld)
   real(r8) cld8(pcols)           ! low cloud fraction estimate
   real(r8) cld9(pcols)           ! mid and high cloud fraction estimate
#ifdef STDCONCLD
   real(r8) cck(pcols)            ! convective cloud per level (assuming
!                                  random overlap in convective layer)
   real(r8) zrth                  ! reciprocal of no. of convective layers
   real(r8) ccldt(pcols)          ! estimate of total convective cloud
#endif
   real(r8) dthtdp(pcols,pver)    ! lapse rate (d theta/dp) below 750 mb
   real(r8) dtdpmn(pcols)         ! most stable lapse rate below 750 mb
   real(r8) dthdp                 ! lapse rate (intermediate variable)
   real(r8) es(pcols,pver)        ! saturation vapor pressure
   real(r8) qs(pcols,pver)        ! saturation specific humidity
   real(r8) premib                ! bottom pressure bound of middle cloud
   real(r8) pretop                ! pressure bounding high cloud
   real(r8) rh(pcols,pver)        ! relative humidity
#ifdef OLDLOWCLD
   real(r8) rhb                   ! intermediate scratch variable
   real(r8) pdepth                ! intermediate scratch variable
   real(r8) stratfac              ! intermediate scratch variable
#endif
   real(r8) rhdif                 ! intermediate scratch variable
   real(r8) strat                 ! intermediate scratch variable
   real(r8) theta(pcols,pver)     ! potential temperature
   real(r8) bvf                   ! brunt-vaisalla frequency
   real(r8) rbvflim               ! bound on inverse of bvf
   real(r8) rho                   ! local density (used to calculate bvf)
   real(r8) rhlim                 ! local rel. humidity threshold estimate
   real(r8) rhden                 ! intermediate scratch variable
   real(r8) rhdif2                ! intermediate scratch variable
   real(r8) rhminl                ! minimum rh for low stable clouds
   real(r8) rhminh                ! minimum rh for high stable clouds
   real(r8) mcbar(pcols)          ! mean convective scale motion in column
   real(r8) dpsum(pcols)          ! vertical sum of delta-p (k-1 levels)
   real(r8) coef1                 ! coefficient to convert mass flux to mb/d
   real(r8) clrsky(pcols)         ! temporary used in random overlap calc
   real(r8) rpdeli(pcols,pver-1) ! 1./(pmid(k+1)-pmid(k))
   real(r8) rhpert                !the specified perturbation to rh

   logical lol(pcols)             ! region of low level cloud
   logical cldbnd(pcols)          ! region below high cloud boundary

   integer i,k                    ! longitude, level indices
   integer kp1
   integer kdthdp(pcols)
   integer numkcld                ! number of levels in which to allow clouds

   real(r8) thetas(pcols)
!
! Statement functions
!
   logical land
   logical ocean

   land(i) = nint(landfrac(i)) == 1
   ocean(i) = nint(ocnfrac(i)) == 1
!
! Set bound for inverse of brunt-vaisalla frequency and minimum relative
! humidity thresholds for stable clouds.  These are the principal
! "disposable" parameters for the cloud fraction scheme
!
   rbvflim = 1./0.00035

! set defaults for rhu00

   rhu00(:,:) = 2.0

   if ( dycore_is ('LR') ) then
        rhminl = .90
   else
        rhminl = .85
   endif

!!(ljli0731)   rhminh = .90
    rhminh = .95   !!(ljli0722)
!
! define rh perturbation in order to estimate rhdfda
!
   rhpert = 0.01 
!
! Evaluate potential temperature and relative humidity
!
   call aqsat(temp    ,pmid    ,es      ,qs      ,pcols   , &
              ncol    ,pver    ,1       ,pver    )
   do k=1,pver
      do i=1,ncol
         theta(i,k)  = temp(i,k)*(pnot/pmid(i,k))**cappa
         rh(i,k)     = q(i,k)/qs(i,k)*(1.0+float(dindex)*rhpert)
!
!  record relhum, rh itself will later be modified related with concld
!
         relhum(i,k) = rh(i,k)
         cloud(i,k)  = 0.
         cldst(i,k)  = 0.
         concld(i,k) = 0.
      end do
   end do
!
! Initialize other temporary variables
!
   do i=1,ncol
      thetas(i)  = ts(i)*(pnot/ps(i))**cappa
      clc(i) = 0.0
   end do
   coef1 = gravit*864.0    ! conversion to millibars/day
   do i=1,ncol
      mcbar(i) = 0.0
      dpsum(i) = 0.0
   end do

   do k=1,pver-1
      do i=1,ncol
         rpdeli(i,k) = 1./(pmid(i,k+1) - pmid(i,k))
      end do
   end do
!
! Calculate mean convective motion throughout column (in units of mb/day)
!
   do k=1,pver-1
      do i=1,ncol
         mcbar(i) = mcbar(i) + max(cmfmc(i,k+1)*coef1,0._r8)*pdel(i,k)
         dpsum(i) = dpsum(i) + pdel(i,k)
      end do
   end do
!
! Estimate of total convective cloud cover based on mean convective motion
!
#ifdef STDCONCLD
   do i=1,ncol
      cck(i) = 0.0
      mcbar(i) = max(mcbar(i)/dpsum(i),1.0e-15_r8)
      ccldt(i) = min(0.035*log(1.0+mcbar(i)),0.80_r8)
      if ((cldbot(i) - cldtop(i)) >= 1.0) then
!
! Inverse of depth of convection (depth is expressed in model levels)
!
         zrth = 1.0/(cldbot(i) - cldtop(i))
!
! Compute amount of convective cloud at each level so that
! after random overlap, the total convective cloud cover is ccldt
!
         cck(i) = 1.0 - (1.0 - ccldt(i))**zrth
      end if
   end do
!
! Vertically distribute cloud in convective layer
!
   do k=1,pver-1
      do i=1,ncol
         if (k <= cldbot(i) .and. k >= cldtop(i)) then
            concld(i,k) = cck(i)
            rh(i,k) = (rh(i,k) - concld(i,k))/(1.0 - concld(i,k))
         end if
      end do
   end do
#else
!     make the convective cloud depend on the conv. mass detraining
!     for upper levels only (above 500mb), since Xu and Kreuger showed
!     rh is a very poor predictor of those clouds
   do k = 1,pver-1
      do i = 1,ncol
         if (pmid(i,k) < 5.e4) then
!               dmudp = (cmfmc(i,k+1)-cmfmc(i,k))/pdel(i,k)
            concld(i,k) = min(rh(i,k),min(1._r8,max(0._r8,zdu(i,k)*5.e4)))
!!(ljli0801)            concld(i,k) = min(rh(i,k)*rh(i,k),min(1._r8,max(0._r8,zdu(i,k)*1.e4)))     !!(ljli0722)
         endif
      end do
   end do
#endif
!
! Evaluate effective column-integrated convective cloud cover using
! random overlap assumption (for diagnostic purposes only)
!
   do i=1,ncol
      clrsky(i) = 1.0
   end do
   do k=pver,1,-1
      do i=1,ncol
         clrsky(i) = clrsky(i)*(1. - concld(i,k))
      end do
   end do
   do i=1,ncol
      clc(i) = 1. - clrsky(i)
   end do
!
!          ****** Compute layer cloudiness ******
!
! There is effecively no top for high cloud formation (can for all the way
! up to 1mb)
! The bottom of middle level cloud (or the top of low level cloud) is
! arbitrarily define to be 750 mb (premib)
!
   premib = 750.e2
   pretop = 1.0e2                 ! top of cloud layer is at 1 mb
!
! Find most stable level below 750 mb for evaluating stratus regimes
!
   do i=1,ncol
      dtdpmn(i) = 0.0
      kdthdp(i) = 0
      dthtdp(i,1) = 0.0
   end do
   do k=2,pver-2
      do i=1,ncol
         if (pmid(i,k) >= premib) then
            dthdp = 100.0*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
         else
            dthdp = 0.0
         end if
         if (dthdp < dtdpmn(i)) then
            dtdpmn(i) = dthdp
            kdthdp(i) = k     ! index of interface of max inversion
         end if
         dthtdp(i,k) = dthdp
      end do
   end do
   do k=pver-1,pver
      do i=1,ncol
         if (0.0 < dtdpmn(i)) then
            dtdpmn(i) = 0.0
         end if
         dthtdp(i,k) = 0.0
      end do
   end do
!
! For some reason, allowing clouds in the bottom model layer causes bad
! error growth characteristics
!
#ifdef PERGRO
   numkcld = pver - 1
#else
   numkcld = pver
#endif
!
! bvf => brunt-vaisalla frequency (approx. 1-sided diff.)
! this stability measure is used to set a local relative humidity
! threshold when evaluating the fractional area of layered cloud
!
   do 10 k=2,numkcld
      kp1 = min(k + 1,pver)
      do i=1,ncol
         if (dthtdp(i,k) > dtdpmn(i)) then
            dthtdp(i,k) = 0.
         end if
         cldbnd(i) = pmid(i,k).ge.pretop
         lol(i) = pmid(i,k).ge.premib
         rho = pmid(i,k)/(rair*temp(i,k))
         bvf = -rho*gravit*gravit*((theta(i,k)-theta(i,k-1))* &
                  rpdeli(i,k-1))/theta(i,k)
         if (cldbnd(i)) then
            rhlim = 0.999 - (1.0-rhminh)*(1.0-min(1.0_r8,max(0.0_r8,bvf*rbvflim)))
            rhden = 1.0 - rhlim
         else
            rhlim = 0.999
            rhden = 0.001
         end if
         rhdif = (rh(i,k) - rhlim)/rhden
#ifdef PERGRO
         cld9(i) = 0.1
#else
         cld9(i) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
#endif
!
! Ignore brunt-vaisalla stability estimate of local relative humidity
! threshold when evaluating low cloud where local vertical motion is
! less than some prescribed value (see low cloud section below)
! Relative humidity threshold is fixed at rhminl for this case, except
! over snow-free land, where it is reduced by 10%.  This distinction is
! made to account for enhanced cloud drop nucleation ({\it i.e.,} at
! lower relative humidities) that can occur over CCN rich land areas.
!
         if (lol(i)) then
            if (land(i) .and. (snowh(i) <= 0.000001)) then
               rhlim = rhminl - 0.10
            else
               rhlim = rhminl
            endif
            rhdif2 = (rh(i,k) - rhlim)/(1.0-rhlim)
            cld8(i) = min(0.999_r8,(max(rhdif2,0.0_r8))**2)
         else
            cld8(i) = cld9(i)
         end if
!
! save rhlim to rhu00, it handles well by itself for low/high cloud
!
         rhu00(i,k)=rhlim
      end do
!
! Final evaluation of layered cloud fraction
!
      do i=1,ncol
!
! Low cloud: non-zero only if vertical velocity requirements are satisfied
! Current vertical velocity threshold is omega < +50 mb/day with a 50 mb/day
! linear ramp (other quantities in the class of "disposable" parameters)
!
         if (lol(i)) then
            if (omga(i,k) < 0.05787) then
               cld = cld8(i)* min(1.0_r8,max(0.0_r8,(0.05787-omga(i,k))/0.05787))
            else
               cld = 0.0
!
! give a fake value of rhlim, 2.0, which would never be used
!            
               rhu00(i,k)=2.0
            end if
            cloud(i,k) = cld
#ifdef OLDLOWCLD
!
! Compute cloud associated with low level inversions.
!
            strat = max(0.,min(0.95_r8,-6.67*dthtdp(i,k) - 0.667))
            rhb   = 1.0 - (0.9 - rh(i,k+1))/0.3
            if (rh(i,k+1) < 0.6) then
               strat = 0.0
            end if
            if (rh(i,k+1) >= 0.6 .and. rh(i,k+1) <= 0.9) then
               strat = strat*rhb
            end if
!
! Linear transition from stratus to trade cu as inversion rises.
! Transition starts at 900mb and completes at 750mb (premib)
!
            pdepth = max(pmid(i,k) - 750.e2,0.0)
            stratfac = min(pdepth,150.0e2_r8)/150.e2
            if (dthtdp(i,k) <= -0.125 ) then
               cloud(i,k) = strat*stratfac
            else
               cloud(i,k) = cld
            end if
#endif
         else                  ! Middle and high level cloud
            if ( cldbnd(i) ) then
               cloud(i,k) = cld9(i)
            else
               cloud(i,k) = 0.0
!
! set a fake value for rhlim, 2.0
!
               rhu00(i,k)=2.0
            end if
         end if
      end do
10    continue                  ! k=2,pver-1
#ifdef PERGRO
      rhu00(:ncol,pver)=0.0
#endif
!
! Add in the marine strat
!
      do i=1,ncol
!
!jrbee bugfix?
!
         if (kdthdp(i) /= 0) then
            k = kdthdp(i)
            kp1 = min(k+1,pver)
            strat = min(1._r8,max(0._r8,(theta(i,k700)-thetas(i))*.057-.5573))
!!(ljli0725)strat = min(1._r8,0.8*max(0._r8,(theta(i,k700)-thetas(i))*.057-.5573))
!
! assign the stratus to the layer just below max inversion
! the relative humidity changes so rapidly across the inversion
! that it is not safe to just look immediately below the inversion
! so limit the stratus cloud by rh in both layers below the inversion
!
            if (ocean(i) .and. dthtdp(i,k) <= -0.125 ) then
               cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
               cloud(i,k) = max(cloud(i,k),cldst(i,k))
            else
               cldst(i,k) = 0.
            end if
         end if
      end do
!
! Merge convective and layered cloud fraction for total cloud
!
      do k=1,pver
         do i=1,ncol
!
!          cloud(i,k) = max(0.0,min(0.999_r8,
!     $                 (1.0 - concld(i,k))*cloud(i,k) + concld(i,k)))
! change to a max overlap assumption between convective and strat clouds
!
            cloud(i,k) = max(0.0_r8,min(0.999_r8,max(concld(i,k),cloud(i,k))))
#ifndef PERGRO
            if (rh(i,k) > 0.99) then
               cloud(i,k) = max(0.01_r8,cloud(i,k))
            endif
#endif
         end do
      end do
!
      return
end subroutine cldfrc
