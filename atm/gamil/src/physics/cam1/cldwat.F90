#undef DEBUG
#include <misc.h>
#include <params.h>
module cldwat
!----------------------------------------------------------------------- 
! 
! Purpose: Prognostic cloud water data and methods.
! 
! Public interfaces:
!
! inimc -- Initialize constants
! pcond -- Calculate prognostic condensate
!
! Author: P. Rasch, with Modifications by Minghua Zhang
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,        only: masterproc, plevstd,plond,plat
   use ppgrid,        only: pcols, pver, pverp
   use wv_saturation, only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
                            cp, epsqs, ttrice
   use physconst,     only: tmelt

   implicit none

!-----------------------------------------------------------------------
! PUBLIC: Make default data and interfaces private
!-----------------------------------------------------------------------
   private
   public inimc, pcond          ! Public interfaces
   public cldwat_fice           !! sxj add 1124
   integer, public::  ktop      ! Level above 10 hPa
!-----------------------------------------------------------------------
! PRIVATE: Everything else is private to this module
!-----------------------------------------------------------------------

 !! sxj 1124
   real(r8), private, parameter :: tmax_fice = tmelt - 10._r8       ! max temperature for cloud ice formation
!!$   real(r8), private, parameter :: tmax_fice = tmelt          ! max temperature for cloud ice formation
!!$   real(r8), private, parameter :: tmin_fice = tmax_fice - 20.! min temperature for cloud ice formation
   real(r8), private, parameter :: tmin_fice = tmax_fice - 30._r8   ! min temperature for cloud ice formation
!  pjr
   real(r8), private, parameter :: tmax_fsnow = tmelt            ! max temperature for transition to convective snow
   real(r8), private, parameter :: tmin_fsnow = tmelt-5._r8         ! min temperature for transition to convective snow

   real(r8), private:: rhonot   ! air density at surface
   real(r8), private:: t0       ! Freezing temperature
   real(r8), private:: cldmin   ! assumed minimum cloud amount
   real(r8), private:: small    ! small number compared to unity
   real(r8), private:: c        ! constant for graupel like snow cm**(1-d)/s
   real(r8), private:: d        ! constant for graupel like snow
   real(r8), private:: esi      ! collection efficient for ice by snow
   real(r8), private:: esw      ! collection efficient for water by snow
   real(r8), private:: nos      ! particles snow / cm**4
   real(r8), private:: pi       ! Mathematical constant
   real(r8), private:: gravit   ! Gravitational acceleration at surface
   real(r8), private:: rh2o
   real(r8), private:: prhonos
   real(r8), private:: thrpd    ! numerical three added to d
   real(r8), private:: gam3pd   ! gamma function on (3+d)
   real(r8), private:: gam4pd   ! gamma function on (4+d)
   real(r8), private:: rhoi     ! ice density
   real(r8), private:: rhos     ! snow density
   real(r8), private:: rhow     ! water density
   real(r8), private:: mcon01   ! constants used in cloud microphysics
   real(r8), private:: mcon02   ! constants used in cloud microphysics
   real(r8), private:: mcon03   ! constants used in cloud microphysics
   real(r8), private:: mcon04   ! constants used in cloud microphysics
   real(r8), private:: mcon05   ! constants used in cloud microphysics
   real(r8), private:: mcon06   ! constants used in cloud microphysics
   real(r8), private:: mcon07   ! constants used in cloud microphysics
   real(r8), private:: mcon08   ! constants used in cloud microphysics

   integer, private ::  k1mb    ! index of the eta level near 1 mb
#ifdef DEBUG
   integer, private,parameter ::  nlook = 2  ! Number of points to examine
   integer, private ::  ilook(nlook)         ! Longitude index to examine
   integer, private ::  latlook(nlook)       ! Latitude index to examine
   integer, private ::  lchnklook(nlook)     ! Chunk index to examine
   integer, private ::  icollook(nlook)      ! Column index to examine
#endif

contains

!===============================================================================
subroutine cldwat_fice(ncol, t, fice, fsnow)  !! sxj add 1124
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only. 
! This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------
    implicit none
! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature
    real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    real(r8), intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

! Local variables
    integer :: i,k                                   ! loop indexes
!-----------------------------------------------------------------------
! Define fractional amount of cloud that is ice
    do k=1,pver
       do i=1,ncol
! If warmer than tmax then water phase
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_r8
! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_r8
! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if
! snow fraction partitioning
! If warmer than tmax then water phase
          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_r8
! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do
    return
end subroutine cldwat_fice

subroutine inimc( tmeltx, rhonotx, gravitx, rh2ox )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize constants for the prognostic condensate
! 
! Author: P. Rasch, April 1997
! 
!-----------------------------------------------------------------------
   use phys_grid, only: get_chunk_coord_p
   use pmgrid, only: plev, plevp
   integer k
   real(r8), intent(in) :: tmeltx
   real(r8), intent(in) :: rhonotx
   real(r8), intent(in) :: gravitx
   real(r8), intent(in) :: rh2ox

#include <comhyb.h>

#ifdef CRAY
   real(r8) signgam              ! variable required by cray gamma function
   external gamma
#endif
   rhonot = rhonotx          ! air density at surface (gm/cm3)
   gravit = gravitx
   rh2o   = rh2ox
   rhos = .1                 ! assumed snow density (gm/cm3)
   rhow = 1.                 ! water density
   rhoi = 1.                 ! ice density
   esi = 1.0                 ! collection efficient for ice by snow
   esw = 0.1                 ! collection efficient for water by snow
   t0 = tmeltx               ! approximate freezing temp
   cldmin = 0.02             ! assumed minimum cloud amount
   small = 1.e-22            ! a small number compared to unity
   c = 152.93                ! constant for graupel like snow cm**(1-d)/s
   d = 0.25                  ! constant for graupel like snow
   nos = 3.e-2               ! particles snow / cm**4
   pi = 4.*atan(1.0)
   prhonos = pi*rhos*nos
   thrpd = 3. + d
#ifdef CRAY
   call gamma(3.+d, signgam, gam3pd)
   gam3pd = sign(exp(gam3pd),signgam)
   call gamma(4.+d, signgam, gam4pd)
   gam4pd = sign(exp(gam4pd),signgam)
   write (6,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
#else
   if (d==0.25) then
      gam3pd = 2.549256966718531 ! only right for d = 0.25
      gam4pd = 8.285085141835282
   else
      write (6,*) ' can only use d ne 0.25 on a cray '
      stop
   endif
#endif
   mcon01 = pi*nos*c*gam3pd/4.
   mcon02 = 1./(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4.)))
   mcon03 = -(0.5+d/4.)
   mcon04 = 4./(4.+d)
   mcon05 = (3+d)/(4+d)
   mcon06 = (3+d)/4.
   mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
   mcon08 = -0.5/(4.+d)

!  find the level about 1mb, we wont do the microphysics above this level
   k1mb = 1
   do k=1,pver-1
      if (hypm(k) < 1.e2 .and. hypm(k+1) >= 1.e2) then
         if (1.e2-hypm(k) < hypm(k+1)-1.e2) then
            k1mb = k
         else
            k1mb = k + 1
         end if
         goto 20
      end if
   end do
   if (masterproc) then
      write(6,*)'inimc: model levels bracketing 1 mb not found'
   end if
!  call endrun
   k1mb = 1
20 if( masterproc ) write(6,*)'inimc: model level nearest 1 mb is',k1mb,'which is',hypm(k1mb),'pascals'

#ifdef DEBUG
!
! Set indicies of the point to examine for debugging
!
   latlook(:) = (/64, 32/)   ! Latitude indices to examine
   ilook(:)   = (/1,   1/)   ! Longitude indicex to examine
   call get_chunk_coord_p( nlook, ilook, latlook, icollook, lchnklook )
#endif
   if( masterproc ) write (6,*) 'cloud water initialization by inimc complete '

   return
end subroutine inimc

subroutine pcond (lchnk   ,ncol    , &
                  tn      ,ttend   ,qn      ,qtend   ,omega   , &
                  cwat    ,p       ,pdel    ,cldn    , &
                  cme     ,evapr   ,prain   ,rmelt   , &     
                  deltat  ,pcflx   ,fwaut   ,fsaut   ,fracw   , &
                  fsacw   ,fsaci   ,lctend  ,rhdfda  ,rhu00, icefrac)  
!----------------------------------------------------------------------- 
! 
! Purpose: 
! The public interface to the cloud water parameterization
! returns tendencies to water vapor, temperature and cloud water variables
! 
! For basic method 
!  See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
!  model climate using diagnosed and 
!  predicted condensate parameterizations, 1998, J. Clim., 11,
!  pp1587---1614.
! 
! For important modifications to improve the method of determining
! condensation/evaporation see Zhang et al (2001, in preparation)
!
! Authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
!-----------------------------------------------------------------------
   use wv_saturation, only: vqsatd
!
!---------------------------------------------------------------------
!
! Input Arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: cldn(pcols,pver)     ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: cwat(pcols,pver)     ! cloud water (kg/kg)
   real(r8), intent(in) :: omega(pcols,pver)    ! vert pressure vel (Pa/s)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure          (K)
   real(r8), intent(in) :: pcflx(pcols,pverp)   ! convective precip level by level (kg/m2/s)  (DISABLED)
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure thickness (Pa)
   real(r8), intent(in) :: qn(pcols,pver)       ! new water vapor    (kg/kg)
   real(r8), intent(in) :: qtend(pcols,pver)    ! mixing ratio tend  (kg/kg/s)
   real(r8), intent(in) :: tn(pcols,pver)       ! new temperature    (K)
   real(r8), intent(in) :: ttend(pcols,pver)    ! temp tendencies    (K/s)
   real(r8), intent(in) :: deltat               ! time step to advance solution over
   real(r8), intent(in) :: lctend(pcols,pver)   ! cloud liquid water tendencies   ====wlin
   real(r8), intent(in) :: rhdfda(pcols,pver)   ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
   real(r8), intent(in) :: rhu00 (pcols,pver)   ! Rhlim for cloud                 ====wlin
   real(r8), intent(in) :: icefrac(pcols)       ! sea ice fraction  (fraction)
!
! Output Arguments
!
   real(r8), intent(out) :: cme(pcols,pver)      ! rate of cond-evap within the cloud
   real(r8), intent(out) :: evapr(pcols,pver)    ! rate of evaporation of falling precip (1/s)
   real(r8), intent(out) :: prain(pcols,pver)    ! rate of conversion of condensate to precip (1/s)
   real(r8), intent(out) :: rmelt(pcols,pver)    ! heating rate due to precip phase change (K/s) (DISABLED)

!
! Local workspace
!
   integer i                 ! work variable
   integer iter              ! #iterations for precipitation calculation
   integer k                 ! work variable
   integer l                 ! work variable

   real(r8) cldm(pcols)          ! mean cloud fraction over the time step
   real(r8) cldmax(pcols)        ! max cloud fraction above
   real(r8) coef(pcols)          ! conversion time scale for condensate to rain
   real(r8) conke                ! rate of evaporation of precipitation:
   real(r8) cwm(pcols)           ! cwat mixing ratio at midpoint of time step
   real(r8) cwn(pcols)           ! cwat mixing ratio at end
   real(r8) denom                ! work variable
   real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
   real(r8) es(pcols)            ! sat. vapor pressure
   real(r8) fice(pcols)          ! fraction of cwat that is ice
   real(r8) fracw(pcols,pver)    ! relative importance of collection of liquid by rain
   real(r8) fsaci(pcols,pver)    ! relative importance of collection of ice by snow
   real(r8) fsacw(pcols,pver)    ! relative importance of collection of liquid by snow
   real(r8) fsaut(pcols,pver)    ! relative importance of ice auto conversion
   real(r8) fwaut(pcols,pver)    ! relative importance of warm cloud autoconversion
   real(r8) gamma(pcols)         ! d qs / dT
   real(r8) iceab(pcols)         ! rate of ice only from above
   real(r8) icwc(pcols)          ! in-cloud water content (kg/kg)
   real(r8) mincld               ! a small cloud fraction to avoid / zero
   real(r8) omeps                ! 1 minus epsilon
   real(r8) omsm                 ! a number just less than unity (for rounding)
   real(r8) precab(pcols)        ! rate of precipitation (kg / (m**2 * s))
   real(r8) prect(pcols)         ! rate of precipitation including convection (kg / (m**2 * s))
   real(r8) prprov(pcols)        ! provisional value of precip at btm of layer
   real(r8) prtmp                ! work variable
   real(r8) q(pcols,pver)        ! mixing ratio before time step ignoring condensate
   real(r8) qs(pcols)            ! spec. hum. of water vapor
   real(r8) qsn, esn             ! work variable
   real(r8) qsp(pcols,pver)      ! sat pt mixing ratio
   real(r8) qtl(pcols)           ! tendency which would saturate the grid box in deltat
   real(r8) qtmp, ttmp           ! work variable
   real(r8) relhum1(pcols)        ! relative humidity
   real(r8) relhum(pcols)        ! relative humidity
   real(r8) tc                   ! crit temp of transition to ice
   real(r8) t(pcols,pver)        ! temp before time step ignoring condensate
   real(r8) tsp(pcols,pver)      ! sat pt temperature
   real(r8) pol                  ! work variable
   real(r8) cdt                  ! work variable

! Extra local work space for cloud scheme modification       

   real(r8) cpohl                !Cp/Hlatv
   real(r8) hlocp                !Hlatv/Cp
   real(r8) dto2                 !0.5*deltat (delta=2.0*dt)
   real(r8) calpha(pcols)        !alpha of new C - E scheme formulation
   real(r8) cbeta (pcols)        !beta  of new C - E scheme formulation
   real(r8) cbetah(pcols)        !beta_hat at saturation portion 
   real(r8) cgamma(pcols)        !gamma of new C - E scheme formulation
   real(r8) cgamah(pcols)        !gamma_hat at saturation portion
   real(r8) rcgama(pcols)        !gamma/gamma_hat
   real(r8) csigma(pcols)        !sigma of new C - E scheme formulation
   real(r8) cmec1 (pcols)        !c1    of new C - E scheme formulation
   real(r8) cmec2 (pcols)        !c2    of new C - E scheme formulation
   real(r8) cmec3 (pcols)        !c3    of new C - E scheme formulation
   real(r8) cmec4 (pcols)        !c4    of new C - E scheme formulation
   real(r8) cmeres(pcols)        !residual cond of over-sat after cme and evapr
   real(r8) ctmp                 !a scalar representation of cmeres
   real(r8) clrh2o               ! Ratio of latvap to water vapor gas const
!
!------------------------------------------------------------
#include <comadj.h>              
!------------------------------------------------------------
!
   clrh2o = hlatv/rh2o   ! Ratio of latvap to water vapor gas const
   omeps = 1.0 - epsqs
#ifdef PERGRO
   mincld = 1.e-4
   iter = 1   ! number of times to iterate the precipitation calculation
#else
   mincld = 1.e-4
   iter = 2
#endif
   omsm = 0.99999
   cpohl = cp/hlatv
   hlocp = hlatv/cp
   dto2=0.5*deltat
!
! Constant for computing rate of evaporation of precipitation:
!
   conke = 1.e-5
!
! initialize a few single level fields
!
   do i = 1,ncol
      precab(i) = 0.0
      prect(i) = 0.0
      iceab(i) = 0.0                ! latent heat of precip above
      cldmax(i) = 0.0
   end do
!
! initialize multi-level fields 
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qn(i,k) 
         t(i,k) = tn(i,k)
      end do
   end do
   cme  (:ncol,:) = 0._r8
   evapr(:ncol,:) = 0._r8
   prain(:ncol,:) = 0._r8
   rmelt(:ncol,:) = 0._r8
   fwaut(:ncol,:) = 0._r8
   fsaut(:ncol,:) = 0._r8
   fracw(:ncol,:) = 0._r8
   fsacw(:ncol,:) = 0._r8
   fsaci(:ncol,:) = 0._r8
!
! find the wet bulb temp and saturation value
! for the provisional t and q without condensation
!
   call findsp (lchnk, ncol, qn, tn, p, tsp, qsp)
   do 800 k = k1mb,pver
      call vqsatd (t(1,k), p(1,k), es, qs, gamma, ncol)
      do i = 1,ncol
         relhum(i) = q(i,k)/qs(i)
!
         cldm(i) = max(cldn(i,k),mincld)
!
! the max cloud fraction above this level
!
         cldmax(i) = max(cldmax(i), cldm(i))

! define the coefficients for C - E calculation

         calpha(i) = 1.0/qs(i)
         cbeta (i) = q(i,k)/qs(i)**2*gamma(i)*cpohl
         cbetah(i) = 1.0/qs(i)*gamma(i)*cpohl
         cgamma(i) = calpha(i)+hlatv*cbeta(i)/cp
         cgamah(i) = calpha(i)+hlatv*cbetah(i)/cp
         rcgama(i) = cgamma(i)/cgamah(i)

         if(cldm(i) > mincld) then
            icwc(i) = max(0._r8,cwat(i,k)/cldm(i))
         else
            icwc(i) = 0.0
         endif

!
! initial guess of evaporation, will be updated within iteration
!
         evapr(i,k) = conke*(1. - cldm(i))*sqrt(precab(i)) &
                        *(1. - min(relhum(i),1._r8))

!
! zero cmeres before iteration for each level
!
         cmeres(i)=0.0

      end do
      do i = 1,ncol
!
! fractions of ice at this level
!
         tc = t(i,k) - t0
         fice(i) = max(0._r8,min(-tc*0.05,1.0_r8))
!
! calculate the cooling due to a phase change of the rainwater
! from above
!
         if (t(i,k) >= t0) then
!              rmelt(i,k) =  -hlatf/cp*iceab(i)*gravit/pdel(i,k)
            rmelt(i,k) = 0.
            iceab(i) = 0.
         else
            rmelt(i,k) = 0.
         endif
      end do

!
! calculate cme and formation of precip. 
!
! The cloud microphysics is highly nonlinear and coupled with cme
! Both rain processes and cme are calculated iteratively.
! 
      do 100 l = 1,iter

        do i = 1,ncol

!
! calculation of cme has 4 scenarios
! ==================================
!
           if(relhum(i) > rhu00(i,k)) then
    
           ! 1. whole grid saturation
           ! ========================
              if(relhum(i) >= 0.999_r8 .or. cldm(i) >= 0.999_r8 ) then
                 cme(i,k)=(calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)

           ! 2. fractional saturation
           ! ========================
              else
                  csigma(i) = 1.0/(rhdfda(i,k)+cgamma(i)*icwc(i))
                  cmec1(i) = (1.0-cldm(i))*csigma(i)*rhdfda(i,k)
                  cmec2(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0-rcgama(i)*cldm(i))*   &
                             csigma(i)*calpha(i)*icwc(i)
                  cmec3(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
                           (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigma(i)*icwc(i)
                  cmec4(i) = csigma(i)*cgamma(i)*icwc(i)

                  ! Q=C-E=-C1*Al + C2*Aq - C3* At + C4*Er
  
                  cme(i,k) = -cmec1(i)*lctend(i,k) + cmec2(i)*qtend(i,k)  &
                             -cmec3(i)*ttend(i,k) + cmec4(i)*evapr(i,k)
               endif

           ! 3. when rh < rhu00, evaporate existing cloud water
           ! ================================================== 
           else if(cwat(i,k) > 0.0)then
              ! liquid water should be evaporated but not to exceed 
              ! saturation point. if qn > qsp, not to evaporate cwat
              cme(i,k)=-min(max(0._r8,qsp(i,k)-qn(i,k)),cwat(i,k))/deltat 

           ! 4. no condensation nor evaporation
           ! ==================================
           else
              cme(i,k)=0.0
           endif

  
        end do    !end loop for cme update

! Because of the finite time step, 
! place a bound here not to exceed wet bulb point
! and not to evaporate more than available water
!
         do i = 1, ncol
            qtmp = qn(i,k) - cme(i,k)*deltat

! possibilities to have qtmp > qsp
!
!   1. if qn > qs(tn), it condenses; 
!      if after applying cme,  qtmp > qsp,  more condensation is applied. 
!      
!   2. if qn < qs, evaporation should not exceed qsp,
    
            if(qtmp > qsp(i,k)) then
              cme(i,k) = cme(i,k) + (qtmp-qsp(i,k))/deltat
            endif

!
! if net evaporation, it should not exceed available cwat
!
            if(cme(i,k) < -cwat(i,k)/deltat)  &
               cme(i,k) = -cwat(i,k)/deltat
!
! addition of residual condensation from previous step of iteration
!
            cme(i,k) = cme(i,k) + cmeres(i)

         end do

         do i = 1,ncol
!
! as a safe limit, condensation should not reduce grid mean rh below rhu00
! 
           if(cme(i,k) > 0.0 .and. relhum(i) > rhu00(i,k) )  &
              cme(i,k) = min(cme(i,k), (qn(i,k)-qs(i)*rhu00(i,k))/deltat)
!
! initial guess for cwm (mean cloud water over time step) if 1st iteration
!
           if(l < 2) then
             cwm(i) = max(cwat(i,k)+cme(i,k)*dto2,  0._r8)
           endif

         enddo

! provisional precipitation falling through model layer
         do i = 1,ncol
            prprov(i) = prect(i) + prain(i,k)*pdel(i,k)/gravit
         end do

! calculate conversion of condensate to precipitation by cloud microphysics 
         call findmcnew (lchnk   ,ncol    , &
                         k       ,prprov  ,t       ,p       , &
                         cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
                         fwaut(1,k),fsaut(1,k),fracw(1,k),fsacw(1,k),fsaci(1,k), &
                         icefrac)
!
! calculate the precip rate
!
         do i = 1,ncol
            if (cldm(i) > 0) then  
!
! first predict the cloud water
!
               cdt = coef(i)*deltat
               if (cdt > 0.01) then
                  pol = cme(i,k)/coef(i) ! production over loss
                  cwn(i) = max(0._r8,(cwat(i,k)-pol)*exp(-cdt)+ pol)
               else
                  cwn(i) = max(0._r8,(cwat(i,k) + cme(i,k)*deltat)/(1+cdt))
               endif
!
! now back out the tendency of net rain production
!
               prain(i,k) =  max(0._r8,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
            else
               prain(i,k) = 0.0
               cwn(i) = 0.
            endif
!
! update any remaining  provisional values
!
            cwm(i) = (cwn(i) + cwat(i,k))*0.5
!
! update in cloud water
!
            if(cldm(i) > mincld) then
               icwc(i) = cwm(i)/cldm(i)
            else
               icwc(i) = 0.0
            endif

         end do              ! end of do i = 1,ncol

!
! calculate provisional value of cloud water for
! evaporation of rain (evapr) calculation
!
      do i = 1,ncol
         qtmp = qn(i,k) - cme(i,k)*deltat
         ttmp = tn(i,k) + rmelt(i,k)*deltat + hlocp*deltat*cme(i,k)
         esn = estblf(ttmp)
         qsn = min(epsqs*esn/(p(i,k) - omeps*esn),1._r8)
         qtl(i) = max((qsn - qtmp)/deltat,0._r8)
         relhum1(i) = qtmp/qsn
      end do
!
      do i = 1,ncol
#ifdef PERGRO
         evapr(i,k) = conke*(1. - max(cldm(i),mincld))* &
                      sqrt(precab(i))*(1. - min(relhum1(i),1._r8))
#else
         evapr(i,k) = conke*(1. - cldm(i))*sqrt(precab(i)) &
                      *(1. - min(relhum1(i),1._r8))
#endif
!
! limit the evaporation to the amount which is entering the box
! or saturates the box
!
         prtmp = precab(i)*gravit/pdel(i,k)
         evapr(i,k) = min(evapr(i,k), prtmp, qtl(i))*omsm
#ifdef PERGRO
!           zeroing needed for pert growth
         evapr(i,k) = 0.
#endif
      end do

! now remove the residual of any over-saturation. Normally,
! the oversaturated water vapor should have been removed by 
! cme formulation plus constraints by wet bulb tsp/qsp
! as computed above. However, because of non-linearity,
! addition of (cme-evapr) to update t and q may still cause
! a very small amount of over saturation. It is called a
! residual of over-saturation because theoretically, cme
! should have taken care of all of large scale condensation.
! 

       do i = 1,ncol
          qtmp = qn(i,k)-(cme(i,k)-evapr(i,k))*deltat
          ttmp = tn(i,k)+(rmelt(i,k)+hlocp*(cme(i,k)-evapr(i,k)) )  &
                      *deltat
          esn = estblf(ttmp)
          qsn = min(epsqs*esn/(p(i,k) - omeps*esn),1._r8)
          !
          !Upper stratosphere and mesosphere, qsn calculated
          !above may be negative. Here just to skip it instead
          !of resetting it to 1 as in aqsat
          !
          if(qtmp > qsn .and. qsn > 0) then
             !calculate dqsdt, a more precise calculation
             !which taking into account different range of T 
             !can be found in aqsatd.F. Here follows
             !cond.F to calculate it.
             !
             denom = (p(i,k)-omeps*esn)*ttmp*ttmp
             dqsdt = clrh2o*qsn*p(i,k)/denom
             !
             !now extra condensation to bring air to just saturation
             !
             ctmp = (qtmp-qsn)/(1.+hlocp*dqsdt)/deltat
             cme(i,k) = cme(i,k)+ctmp
!
! save residual on cmeres to addtion to cme on entering next iteration
! cme exit here contain the residual but overrided if back to iteration
!
             cmeres(i) = ctmp
          else
             cmeres(i) = 0.0
          endif
       end do
              
 100 continue              ! end of do l = 1,iter
!
! precipitation
!
      do i = 1,ncol
         prtmp = pdel(i,k) / gravit *(prain(i,k) - evapr(i,k))
         iceab(i) = iceab(i) + fice(i)*prtmp
         precab(i) = precab(i) + prtmp
         prect(i) = prect(i) + prtmp + pcflx(i,k+1)
         if ((precab(i)) < 1.e-10) then      
            precab(i) = 0.
         endif
         if ((prect(i)) < 1.e-10) then      
            prect(i) = 0.
         endif
      end do
 800 continue                ! level loop (k=1,pver)

   return
end subroutine pcond

!##############################################################################

subroutine findmcnew (lchnk   ,ncol    , &
                      k       ,precab  ,t       ,p       , &
                      cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
                      fwaut   ,fsaut   ,fracw   ,fsacw   ,fsaci   , &
                      icefrac )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! calculate the conversion of condensate to precipitate
! 
! Method: 
! See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
!  model climate using diagnosed and 
!  predicted condensate parameterizations, 1998, J. Clim., 11,
!  pp1587---1614.
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
   use phys_grid, only: get_rlat_all_p
   use comsrf,        only: landm
!
! input args
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: k                     ! level index

   real(r8), intent(in) :: precab(pcols)        ! rate of precipitation from above (kg / (m**2 * s))
   real(r8), intent(in) :: t(pcols,pver)        ! temperature       (K)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure          (Pa)
   real(r8), intent(in) :: cldm(pcols)          ! cloud fraction
   real(r8), intent(in) :: cldmax(pcols)        ! max cloud fraction above this level
   real(r8), intent(in) :: cwm(pcols)           ! condensate mixing ratio (kg/kg)
   real(r8), intent(in) :: fice(pcols)          ! fraction of cwat that is ice
   real(r8), intent(in) :: icefrac(pcols)       ! sea ice fraction 

! output arguments
   real(r8), intent(out) :: coef(pcols)          ! conversion rate (1/s)
   real(r8), intent(out) :: fwaut(pcols)         ! relative importance of liquid autoconversion (a diagnostic)
   real(r8), intent(out) :: fsaut(pcols)         ! relative importance of ice autoconversion    (a diagnostic)
   real(r8), intent(out) :: fracw(pcols)         ! relative importance of rain accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsacw(pcols)         ! relative importance of snow accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsaci(pcols)         ! relative importance of snow accreting ice    (a diagnostic)

! work variables

   integer i
   integer ii
   integer ind(pcols)
   integer ncols

   real(r8), parameter :: degrad = 57.296 ! divide by this to convert degrees to radians
   real(r8) alpha                ! ratio of 3rd moment radius to 2nd
   real(r8) capc                 ! constant for autoconversion
   real(r8) capn                 ! local cloud particles / cm3
   real(r8) capnoice             ! local cloud particles when not over sea ice / cm3
   real(r8) capnsi               ! sea ice cloud particles / cm3
   real(r8) capnc                ! cold and oceanic cloud particles / cm3
   real(r8) capnw                ! warm continental cloud particles / cm3
   real(r8) ciaut                ! coefficient of autoconversion of ice (1/s)
   real(r8) ciautb               ! coefficient of autoconversion of ice (1/s)
   real(r8) cldloc(pcols)        ! non-zero amount of cloud
   real(r8) cldpr(pcols)         ! assumed cloudy volume occupied by rain and cloud
   real(r8) con1                 ! work constant
   real(r8) con2                 ! work constant
   real(r8) convfw               ! constant used for fall velocity calculation
   real(r8) cracw                ! constant used for rain accreting water
   real(r8) critpr               ! critical precip rate collection efficiency changes
   real(r8) csacx                ! constant used for snow accreting liquid or ice
   real(r8) dtice                ! interval for transition from liquid to ice
   real(r8) effc                 ! collection efficiency
   real(r8) icemr(pcols)         ! in-cloud ice mixing ratio
   real(r8) icrit                ! threshold for autoconversion of ice
   real(r8) icritc               ! threshold for autoconversion of cold ice
   real(r8) icritw               ! threshold for autoconversion of warm ice
   real(r8) kconst               ! const for terminal velocity (stokes regime)
   real(r8) liqmr(pcols)         ! in-cloud liquid water mixing ratio
   real(r8) pracw                ! rate of rain accreting water
   real(r8) prlloc(pcols)        ! local rain flux in mm/day
   real(r8) prscgs(pcols)        ! local snow amount in cgs units
   real(r8) psaci                ! rate of collection of ice by snow (lin et al 1983)
   real(r8) psacw                ! rate of collection of liquid by snow (lin et al 1983)
   real(r8) psaut                ! rate of autoconversion of ice condensate
   real(r8) ptot                 ! total rate of conversion
   real(r8) pwaut                ! rate of autoconversion of liquid condensate
   real(r8) r3l                  ! volume radius
   real(r8) r3lcrit              ! critical radius at which autoconversion become efficient
   real(r8) rainmr(pcols)        ! in-cloud rain mixing ratio
   real(r8) rat1                 ! work constant
   real(r8) rat2                 ! work constant
   real(r8) rdtice               ! recipricol of dtice
   real(r8) rho(pcols)           ! density (mks units)
   real(r8) rhocgs               ! density (cgs units)
   real(r8) rlat(pcols)          ! latitude (radians)
   real(r8) snowfr               ! fraction of precipate existing as snow
   real(r8) totmr(pcols)         ! in-cloud total condensate mixing ratio
   real(r8) vfallw               ! fall speed of precipitate as liquid
   real(r8) wp                   ! weight factor used in calculating pressure dep of autoconversion
   real(r8) wsi                  ! weight factor for sea ice
   real(r8) wt                   ! fraction of ice
   real(r8) wland                ! fraction of land

!      real(r8) csaci
!      real(r8) csacw
!      real(r8) cwaut
!      real(r8) efact
!      real(r8) lamdas
!      real(r8) lcrit
!      real(r8) rcwm
!      real(r8) r3lc2
!      real(r8) snowmr(pcols)
!      real(r8) vfalls


!     inline statement functions
   real(r8) heavy, heavym, a1, a2, heavyp, heavymp
   heavy(a1,a2) = max(0._r8,sign(1._r8,a1-a2))  ! heavyside function
   heavym(a1,a2) = max(0.01_r8,sign(1._r8,a1-a2))  ! modified heavyside function
!
! New heavyside functions to perhaps address error growth problems
!
   heavyp(a1,a2) = a1/(a2+a1+1.e-36)
   heavymp(a1,a2) = (a1+0.01*a2)/(a2+a1+1.e-36)

! critical precip rate at which we assume the collector drops can change the
! drop size enough to enhance the auto-conversion process (mm/day)
   critpr = 0.5

   convfw = 1.94*2.13*sqrt(rhow*1000.*9.81*2.7e-4)

! liquid microphysics
!      cracw = 6                 ! beheng
   cracw = .884*sqrt(9.81/(rhow*1000.*2.7e-4)) ! tripoli and cotton

! ice microphysics
!      ciautb = 6.e-4
!      ciautb = 1.e-3
   ciautb = 5.e-4
!      icritw = 1.e-5
!      icritw = 5.e-5
   icritw = 4.e-4
!      icritc = 4.e-6
!      icritc = 6.e-6
   icritc = 5.e-6

   dtice = 20.
   rdtice = 1./dtice

   capnw = 400.              ! warm continental cloud particles / cm3
   capnc = 150.              ! cold and oceanic cloud particles / cm3
!  capnsi = 40.              ! sea ice cloud particles density  / cm3
   capnsi =  5.              ! sea ice cloud particles density  / cm3

   kconst = 1.18e6           ! const for terminal velocity

!      effc = 1.                 ! autoconv collection efficiency following boucher 96
!      effc = .55*0.05           ! autoconv collection efficiency following baker 93
   effc = 0.55                ! autoconv collection efficiency following tripoli and cotton
!   effc = 0.    ! turn off warm-cloud autoconv
   alpha = 1.1**4
   capc = pi**(-.333)*kconst*effc *(0.75)**(1.333)*alpha  ! constant for autoconversion

   r3lcrit = 15.0e-6         ! 15.0u  crit radius where liq conversion begins
!
! find all the points where we need to do the microphysics
! and set the output variables to zero
!
   ncols = 0
   do i = 1,ncol
      coef(i) = 0.
      fwaut(i) = 0.
      fsaut(i) = 0.
      fracw(i) = 0.
      fsacw(i) = 0.
      fsaci(i) = 0.
      liqmr(i) = 0.
      rainmr(i) = 0.
      if (cwm(i) > 1.e-20) then
         ncols = ncols + 1
         ind(ncols) = i
      endif
   end do

!cdir$ ivdep
   do ii = 1,ncols
      i = ind(ii)
!
! the local cloudiness at this level
!
      cldloc(i) = max(cldmin,cldm(i))
!
! a weighted mean between max cloudiness above, and this layer
!
      cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5)
!
! decompose the suspended condensate into
! an incloud liquid and ice phase component
!
      totmr(i) = cwm(i)/cldloc(i)
      icemr(i) = totmr(i)*fice(i)
      liqmr(i) = totmr(i)*(1.-fice(i))
!
! density
!
      rho(i) = p(i,k)/(287.*t(i,k))
      rhocgs = rho(i)*1.e-3     ! density in cgs units
!
! decompose the precipitate into a liquid and ice phase
!
      if (t(i,k) > t0) then
         vfallw = convfw/sqrt(rho(i))
         rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
         snowfr = 0
      else
         snowfr = 1
         rainmr(i) = 0.
      endif
!
! local snow amount in cgs units
!
      prscgs(i) = precab(i)/cldpr(i)*0.1*snowfr
!
! local rain amount in mm/day
!
      prlloc(i) = precab(i)*86400./cldpr(i)
   end do

   con1 = 1./(1.333*pi)**0.333 * 0.01 ! meters
!
! calculate the conversion terms
!
   call get_rlat_all_p(lchnk, ncol, rlat)

!cdir$ ivdep
   do ii = 1,ncols
      i = ind(ii)
      rhocgs = rho(i)*1.e-3     ! density in cgs units
!
! exponential temperature factor
!
!        efact = exp(0.025*(t(i,k)-t0))
!
! some temperature dependent constants
!
      wt = min(1._r8,max(0._r8,(t0-t(i,k))*rdtice))
      icrit = icritc*wt + icritw*(1-wt)
!
! linear weight factor in pressure (1 near sfc, 0 at .8 of sfc)
!
      wp = min(1._r8,max(0._r8,(p(i,k)-0.8*p(i,pver))/(0.2*p(i,pver))))
!
! near land near sfc raise the number concentration
! except south of 60S where land is so clean we treat it as ocean
!
      if (rlat(i) < -60./degrad) then
         wland = 0.
      else
         wland = landm(i,lchnk)
      endif
      
      capnoice =  wland*(capnw*wp + capnc*(1-wp))+ &
                  (1.-wland)*capnc
!
!     modify the estimated value to acknowledge cloud water number changes over sea ice
!
      wsi = icefrac(i) ! (1 all sea ice, 0, all ocean or land)
      capn = (capnsi*wsi + capnoice*(1-wsi))*wp + capnoice*(1-wp)

!!!!  enable the following line to get old cloudwater behavior
!!!!  capn = capnoice

      if (icefrac(i) > 0.001) then
         capn = capnsi
      else
         capn = capnoice
      endif

#ifdef DEBUG
      if ( (lat(i) == latlook(1)) .or. (lat(i) == latlook(2)) ) then
         if (i == ilook(1)) then
            write (6,*) ' findmcnew: lat, k, icefrac, landm, wp, capnoice, capn ', &
                 lat(i), k, icefrac(i), landm(i,lat(i)), wp, capnoice, capn
         endif
      endif
#endif

!
! useful terms in following calculations
!
      rat1 = rhocgs/rhow
      rat2 = liqmr(i)/capn
      con2 = (rat1*rat2)**0.333
!
! volume radius
!
!        r3l = (rhocgs*liqmr(i)/(1.333*pi*capn*rhow))**0.333 * 0.01 ! meters
      r3l = con1*con2
!
! critical threshold for autoconversion if modified for mixed phase
! clouds to mimic a bergeron findeisen process
! r3lc2 = r3lcrit*(1.-0.5*fice(i)*(1-fice(i)))
!
! autoconversion of liquid
!
!        cwaut = 2.e-4
!        cwaut = 1.e-3
!        lcrit = 2.e-4
!        lcrit = 5.e-4
!        pwaut = max(0._r8,liqmr(i)-lcrit)*cwaut
!
! pwaut is following tripoli and cotton (and many others)
! we reduce the autoconversion below critpr, because these are regions where
! the drop size distribution is likely to imply much smaller collector drops than
! those relevant for a cloud distribution corresponding to the value of effc = 0.55
! suggested by cotton (see austin 1995 JAS, baker 1993)

! easy to follow form
!        pwaut = capc*liqmr(i)**2*rhocgs/rhow
!    $           *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
!    $           *heavy(r3l,r3lcrit)
!    $           *max(0.10_r8,min(1._r8,prlloc(i)/critpr))
! somewhat faster form
#define HEAVYNEW
#ifdef HEAVYNEW
!#ifdef PERGRO
      pwaut = capc*liqmr(i)**2*rat1*con2*heavymp(r3l,r3lcrit) * &
              max(0.10_r8,min(1._r8,prlloc(i)/critpr))
#else
      pwaut = capc*liqmr(i)**2*rat1*con2*heavym(r3l,r3lcrit)* &
              max(0.10_r8,min(1._r8,prlloc(i)/critpr))
#endif
!
! autoconversion of ice
!
!        ciaut = ciautb*efact
      ciaut = ciautb
!        psaut = capc*totmr(i)**2*rhocgs/rhoi
!     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
!
! autoconversion of ice condensate
!
#ifdef PERGRO
      psaut = heavyp(icemr(i),icrit)*icemr(i)*ciaut
#else
      psaut = max(0._r8,icemr(i)-icrit)*ciaut
#endif
!
! collection of liquid by rain
!
!        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
      pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
!!      pracw = 0.
!
! the following lines calculate the slope parameter and snow mixing ratio
! from the precip rate using the equations found in lin et al 83
! in the most natural form, but it is expensive, so after some tedious
! algebraic manipulation you can use the cheaper form found below
!            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
!     $               *0.01   ! convert from cm/s to m/s
!            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
!            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
!            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25
!            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
!
! coefficient for collection by snow independent of phase
!
      csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05
!
! collection of liquid by snow (lin et al 1983)
!
      psacw = csacx*liqmr(i)*esw
#ifdef PERGRO
! this is necessary for pergro
      psacw = 0.
#endif
!
! collection of ice by snow (lin et al 1983)
!
      psaci = csacx*icemr(i)*esi
!
! total conversion of condensate to precipitate
!
      ptot = pwaut + psaut + pracw + psacw + psaci
!
! the recipricol of cloud water amnt (or zero if no cloud water)
!
!         rcwm =  totmr(i)/(max(totmr(i),small)**2)
!
! turn the tendency back into a loss rate (1/seconds)
!
      if (totmr(i) > 0.) then
         coef(i) = ptot/totmr(i)
      else
         coef(i) = 0.
      endif

      fwaut(i) = pwaut/max(ptot,small)
      fsaut(i) = psaut/max(ptot,small)
      fracw(i) = pracw/max(ptot,small)
      fsacw(i) = psacw/max(ptot,small)
      fsaci(i) = psaci/max(ptot,small)

   end do
#ifdef DEBUG
   i = icollook(nlook)
   if (lchnk == lchnklook(nlook) ) then
      write (6,*)
      write (6,*) '------', k
      write (6,*) ' liqmr, rainmr,precab ', liqmr(i), rainmr(i), precab(i)*8.64e4
      write (6,*) ' frac: waut,saut,racw,sacw,saci ', &
           fwaut(i), fsaut(i), fracw(i), fsacw(i), fsaci(i)
   endif
#endif

   return
end subroutine findmcnew

!##############################################################################

subroutine findsp (lchnk, ncol, q, t, p, tsp, qsp)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! water vapor (kg/kg)
   real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure    (Pa)
!
! output arguments
!
   real(r8), intent(out) :: tsp(pcols,pver)      ! saturation temp (K)
   real(r8), intent(out) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
!
! local variables
!
   integer i                 ! work variable
   integer k                 ! work variable
   logical lflg              ! work variable
   integer iter              ! work variable
   integer l                 ! work variable

   real(r8) omeps                ! 1 minus epsilon
   real(r8) trinv                ! work variable
   real(r8) es                   ! sat. vapor pressure
   real(r8) desdt                ! change in sat vap pressure wrt temperature
!     real(r8) desdp                ! change in sat vap pressure wrt pressure
   real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
   real(r8) dgdt                 ! work variable
   real(r8) g                    ! work variable
   real(r8) weight(pcols)        ! work variable
   real(r8) hlatsb               ! (sublimation)
   real(r8) hlatvp               ! (vaporization)
   real(r8) hltalt(pcols,pver)   ! lat. heat. of vap.
   real(r8) tterm                ! work var.
   real(r8) qs                   ! spec. hum. of water vapor
   real(r8) tc                   ! crit temp of transition to ice

! work variables
   real(r8) t1, q1, dt, dq
   real(r8) dtm, dqm
   real(r8) qvd, a1, tmp
   real(r8) rair
   real(r8) r1b, c1, c2, c3
   real(r8) denom
   real(r8) dttol
   real(r8) dqtol
   integer doit(pcols) 
   real(r8) enin(pcols), enout(pcols)
   real(r8) tlim(pcols)

   omeps = 1.0 - epsqs
   trinv = 1.0/ttrice
   a1 = 7.5*log(10._r8)
   rair =  287.04
   c3 = rair*a1/cp
   dtm = 0.    ! needed for iter=0 blowup with f90 -ei
   dqm = 0.    ! needed for iter=0 blowup with f90 -ei
   dttol = 1.e-4 ! the relative temp error tolerance required to quit the iteration
   dqtol = 1.e-4 ! the relative moisture error tolerance required to quit the iteration
!  tmin = 173.16 ! the coldest temperature we can deal with
!
! max number of times to iterate the calculation
   iter = 8
!
   do k = k1mb,pver

!
! first guess on the wet bulb temperature
!
      do i = 1,ncol

#ifdef DEBUG
         if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
            write (6,*) ' '
            write (6,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
         endif
#endif
! limit the temperature range to that relevant to the sat vap pres tables
         tlim(i) = min(max(t(i,k),173._r8),373._r8)
         es = estblf(tlim(i))
         denom = p(i,k) - omeps*es
         qs = epsqs*es/denom
         doit(i) = 0
         enout(i) = 1.
! make sure a meaningful calculation is possible
         if (p(i,k) > 5.*es .and. qs > 0. .and. qs < 0.5) then
!
! Saturation specific humidity
!
             qs = min(epsqs*es/denom,1._r8)
!
! "generalized" analytic expression for t derivative of es
!  accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
             tc     = tlim(i) - t0
             lflg   = (tc >= -ttrice .and. tc < 0.0)
             weight(i) = min(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0*tc
             if (tlim(i) < t0) then
                hltalt(i,k) = hlatsb
             else
                hltalt(i,k) = hlatvp
             end if
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36.)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
#ifdef DEBUG
             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                write (6,*) ' relative humidity ', q(i,k)/qs
                write (6,*) ' first guess ', tsp(i,k)
             endif
#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
          else
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1.
          endif
       end do   ! end do i
!
! now iterate on first guess
!
      do l = 1, iter
         dtm = 0
         dqm = 0
         do i = 1,ncol
            if (doit(i) == 0) then
               es = estblf(tsp(i,k))
!
! Saturation specific humidity
!
               qs = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
               tc     = tsp(i,k) - t0
               lflg   = (tc >= -ttrice .and. tc < 0.0)
               weight(i) = min(-tc*trinv,1.0_r8)
               hlatsb = hlatv + weight(i)*hlatf
               hlatvp = hlatv - 2369.0*tc
               if (tsp(i,k) < t0) then
                  hltalt(i,k) = hlatsb
               else
                  hltalt(i,k) = hlatvp
               end if
               if (lflg) then
                  tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
               else
                  tterm = 0.0
               end if
               desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
               dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
!              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
               g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
               dgdt = -(cp + hltalt(i,k)*dqsdt)
               t1 = tsp(i,k) - g/dgdt
               dt = abs(t1 - tsp(i,k))/t1
               tsp(i,k) = max(t1,tmin)
               es = estblf(tsp(i,k))
               q1 = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
               dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
               qsp(i,k) = q1
#ifdef DEBUG
               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                  write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
               endif
#endif
               dtm = max(dtm,dt)
               dqm = max(dqm,dq)
! if converged at this point, exclude it from more iterations
               if (dt < dttol .and. dq < dqtol) then
                  doit(i) = 2
               endif
               enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
! bail out if we are too near the end of temp range
               if (tsp(i,k) < 174.16) then
                  doit(i) = 4
               endif
            else
            endif
         end do              ! do i = 1,ncol

         if (dtm < dttol .and. dqm < dqtol) then
            go to 10
         endif

      end do                 ! do l = 1,iter
10    continue

      if (dtm > dttol .or. dqm > dqtol) then
         do i = 1,ncol
            if (doit(i) == 0) then
               write (6,*) ' findsp not converging at point i, k ', i, k
               write (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
               write (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
               call endrun ()
            endif
         end do
      endif
      do i = 1,ncol
         if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4) then
            write (6,*) ' the enthalpy is not conserved for point ', &
                 i, k, enin(i), enout(i)
            write (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
            write (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
            call endrun ()
         endif
      end do
      
   end do                    ! level loop (k=1,pver)

   return
end subroutine findsp

end module cldwat
