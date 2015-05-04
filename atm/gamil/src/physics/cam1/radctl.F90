#include <misc.h>
#include <params.h>

subroutine radctl(lchnk   ,ncol    ,lwup    ,emis    ,          &
                  pmid    ,pint    ,pmln    ,piln    ,t       , &
                  qm1     ,cld     ,clwp    ,coszrs  ,          &
                  asdir   ,asdif   ,aldir   ,aldif   ,pmxrgn  , &
                  nmxrgn  ,fsns    ,fsnt    ,flns    ,flnt    , &
                  qrs     ,qrl     ,flwds   ,rel     ,rei     , &
                  fice    ,sols    ,soll    ,solsd   ,solld   , &
                  landfrac,zm      )
!-----------------------------------------------------------------------
!
! Purpose:
! Driver for radiation computation.
!
! Method:
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Author: CCM1,  CMS Contact: J. Truesdale
!
!-----------------------------------------------------------------------
!!(wh 2003.12.27)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use pspect
   use so4bnd
   use so4bnd_IPCC    !!(wh)
   use commap
   use history, only: outfld
   use tracers,      only: ixcldw
   use constituents, only: ppcnst, cnst_get_ind
   use physconst, only: cpair

    ! *** added by DONG Li *** !
#if (defined TRIAL_RUN_FORCING)
    use SolarForcing, only: currentSolarConstant
#endif
    ! ************************ !

   implicit none

#include <ptrrgrid.h>
#include <comctl.h>
#include <comsol.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: lwup(pcols)          ! Longwave up flux at surface
   real(r8), intent(in) :: emis(pcols,pver)     ! Cloud emissivity
   real(r8), intent(in) :: pmid(pcols,pver)     ! Model level pressures
   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressures
   real(r8), intent(in) :: pmln(pcols,pver)     ! Natural log of pmid
   real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)      ! ice effective drop size (microns)
   real(r8), intent(in) :: fice(pcols,pver)     ! fractional ice content within cloud
   real(r8), intent(in) :: piln(pcols,pverp)    ! Natural log of pint
   real(r8), intent(in) :: t(pcols,pver)        ! Model level temperatures
   real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! Specific humidity and tracers
   real(r8), intent(in) :: cld(pcols,pver)      ! Fractional cloud cover
   real(r8), intent(in) :: clwp(pcols,pver)     ! Cloud liquid water path
   real(r8), intent(in) :: coszrs(pcols)        ! Cosine solar zenith angle
   real(r8), intent(in) :: asdir(pcols)         ! albedo shortwave direct
   real(r8), intent(in) :: asdif(pcols)         ! albedo shortwave diffuse
   real(r8), intent(in) :: aldir(pcols)         ! albedo longwave direct
   real(r8), intent(in) :: aldif(pcols)         ! albedo longwave diffuse
   real(r8), intent(in) :: landfrac(pcols)      ! land fraction
   real(r8), intent(in) :: zm(pcols,pver)       ! Height of midpoints (above surface)
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pmid for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pmid for
!    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer, intent(inout) :: nmxrgn(pcols)     ! Number of maximally overlapped regions
!
! Output solar arguments
!
   real(r8), intent(out) :: fsns(pcols)          ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)          ! Net column abs solar flux at model top
   real(r8), intent(out) :: flns(pcols)          ! Srf longwave cooling (up-down) flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing lw flux at model top
   real(r8), intent(out) :: sols(pcols)          ! Downward solar rad onto surface (sw direct)
   real(r8), intent(out) :: soll(pcols)          ! Downward solar rad onto surface (lw direct)
   real(r8), intent(out) :: solsd(pcols)         ! Downward solar rad onto surface (sw diffuse)
   real(r8), intent(out) :: solld(pcols)         ! Downward solar rad onto surface (lw diffuse)
   real(r8), intent(out) :: qrs(pcols,pver)      ! Solar heating rate
!
! Output longwave arguments
!
   real(r8), intent(out) :: qrl(pcols,pver)      ! Longwave cooling rate
   real(r8), intent(out) :: flwds(pcols)         ! Surface down longwave flux


!
!---------------------------Local variables-----------------------------
!
   integer i, k              ! index

   integer :: in2o, ich4, if11, if12 ! indexes of gases in constituent array

   real(r8) solin(pcols)         ! Solar incident flux
   real(r8) fsds(pcols)          ! Flux Shortwave Downwelling Surface
   real(r8) fsus(pcols)          ! short wave upward flux at surface ! added by DONG Li
   real(r8) fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) fsutoac(pcols)       ! short wave clearsky upward flux at TOA ! added by DONG Li
   real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
   real(r8) fsusc(pcols)         ! short wave clearsky upward flux at surface ! added by DONG Li
   real(r8) flut(pcols)          ! Upward flux at top of model
   real(r8) lwcf(pcols)          ! longwave cloud forcing
   real(r8) swcf(pcols)          ! shortwave cloud forcing
   real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) pbr(pcols,pverr)     ! Model mid-level pressures (dynes/cm2)
   real(r8) pnm(pcols,pverrp)    ! Model interface pressures (dynes/cm2)
   real(r8) o3vmr(pcols,pverr)   ! Ozone volume mixing ratio
   real(r8) o3mmr(pcols,pverr)   ! Ozone mass mixing ratio
   real(r8) eccf                 ! Earth/sun distance factor
   real(r8) n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8) ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8) cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8) cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8) aermmr(pcols,pverr)  ! level aerosol mass mixing ratio
   real(r8) rh(pcols,pverr)      ! level relative humidity (fraction)
   real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

!
! Declare local arrays to which model input arrays are interpolated here.
! Current default is none since radiation grid = model grid.
!
! Declare variables used for indirect forcing calculations:
!
! ++ tls --------------------------------------------------------------2
   real(r8) locrhoair(pcols,pver)  ! dry air density            [kg/m^3 ]
   real(r8) lwcwat(pcols,pver)     ! in-cloud liquid water path [kg/m^3 ]
   real(r8) sulfbio(pcols,pver)    ! biogenic sulfate mmr       [kg/kg  ]
   real(r8) sulfant(pcols,pver)    ! anthropogenic sulfate mmr  [kg/kg  ]
   real(r8) sulf   (pcols,pver)    ! sulfate mmr  [kg/kg  ]  (for IPCC)!!(wh)
   real(r8) sulfscalef             ! sulfate scale factor
   real(r8) sulfmix(pcols,pver)    ! sulfate mass mixing ratio  [kg/kg  ]
   real(r8) so4mass(pcols,pver)    ! sulfate mass concentration [g/cm^3 ]
   real(r8) Aso4(pcols,pver)       ! sulfate # concentration    [#/cm^3 ]
   real(r8) Ntot(pcols,pver)       ! ccn # concentration        [#/cm^3 ]
   real(r8) relmod(pcols,pver)     ! effective radius           [microns]

   real(r8) wrel(pcols,pver)       ! weighted effective radius    [microns]
   real(r8) wlwc(pcols,pver)       ! weighted liq. water content  [kg/m^3 ]
   real(r8) cldfrq(pcols,pver)     ! frequency of occurance of...
!                                  ! clouds (cld => 0.01)         [fraction]
   real(r8) ftem(pcols,pver)       ! temporary array for outfld

   real(r8) locPi                  ! my piece of the pi
   real(r8) Rdryair                ! gas constant of dry air   [J/deg/kg]
   real(r8) rhowat                 ! density of water          [kg/m^3  ]
   real(r8) Acoef                  ! m->A conversion factor; assumes
!                                  ! Dbar=0.10, sigma=2.0      [g^-1    ]
   real(r8) rekappa                ! kappa in evaluation of re(lmod)
   real(r8) recoef                 ! temp. coeficient for calc of re(lmod)
   real(r8) reexp                  ! 1.0/3.0
   real(r8) Ntotb                  ! temp var to hold below cloud ccn
! -- Parameters for background CDNC (from `ambient' non-sulfate aerosols)...
   real(r8) Cmarn                  ! Coef for CDNC_marine         [cm^-3]
   real(r8) Cland                  ! Coef for CDNC_land           [cm^-3]
   real(r8) Hmarn                  ! Scale height for CDNC_marine [m]
   real(r8) Hland                  ! Scale height for CDNC_land   [m]
   parameter ( Cmarn = 50.0, Cland = 100.0 )
   parameter ( Hmarn = 1000.0, Hland = 2000.0 )
   real(r8) bgaer                  ! temp var to hold background CDNC
!
! Statement functions
!
   logical land
   land(i) = nint(landfrac(i)).gt.0.5_r8
!
! -- tls --------------------------------------------------------------2
!
!--------------------------------------------------------------------------
!
! Interpolate ozone volume mixing ratio to model levels
!
   call radozn(lchnk   ,ncol    ,pmid    ,o3vmr   )
   call outfld('O3VMR   ',o3vmr ,pcols, lchnk)
!
! Set chunk dependent radiation input
!
   call radinp(lchnk   ,ncol    ,                                &
               pmid    ,pint    ,o3vmr   , pbr     ,&
               pnm     ,eccf    ,o3mmr   )
!
! Solar radiation computation
!
   if (dosw) then
! ++ tls ---------------------------------------------------------------2
!     write(6,*) 'Sulfate Scale Factor = ', sulfscalef

      locPi = 3.141592654
      Rdryair = 287.04
      rhowat = 1000.0
      Acoef = 1.2930E14
      recoef = 3.0/(4.0*locPi*rhowat)
      reexp = 1.0/3.0
!
      if ( doRamp_so4 ) then
         call getso4bnd( lchnk, ncol, sulfbio, sulfant )
         sulfscalef = so4ramp()
         do k = 1, pver
            do i = 1, ncol
               sulfmix(i,k) = sulfbio(i,k) + sulfscalef*sulfant(i,k)
            end do
         end do
         call outfld('SULFBIO ',sulfbio,pcols,lchnk)
         call outfld('SULFANT ',sulfant,pcols,lchnk)
         call outfld('SULFMMR ',sulfmix,pcols,lchnk)
                                                             !!(wh)
      else if ( doIPCC_so4 ) then
         call getso4bnd_IPCC( lchnk, ncol, sulfmix )
         call outfld('SULFMMR ',sulfmix,pcols,lchnk)
                                                             !!(wh)
      else
         do k = 1, pver
            do i = 1, ncol
               sulfmix(i,k) = 0.
            end do
         end do
      endif

      if ( indirect ) then ! Method of Martin et. al.
         do k=pver,1,-1
            do i = 1,ncol
               locrhoair(i,k) = pmid(i,k)/( Rdryair*t(i,k) )
               lwcwat(i,k) = ( qm1(i,k,ixcldw)*(1.-fice(i,k))/max(0.01_r8,cld(i,k)) )* &
                             locrhoair(i,k)
!                 NOTE: 0.001 converts kg/m3 -> g/cm3
               so4mass(i,k) = sulfmix(i,k)*locrhoair(i,k)*0.001
               Aso4(i,k) = so4mass(i,k)*Acoef

               if (Aso4(i,k) <= 280.0) then
                  Aso4(i,k) = max(36.0_r8,Aso4(i,k))
                  Ntot(i,k) = -1.15E-3*Aso4(i,k)**2 + 0.963*Aso4(i,k)+5.30
                  rekappa = 0.80
               else
                  Aso4(i,k) = min(1500.0_r8,Aso4(i,k))
                  Ntot(i,k) = -2.10E-4*Aso4(i,k)**2 + 0.568*Aso4(i,k)-27.9
                  rekappa = 0.67
               end if
               if (land(i)) then ! Account for local background aerosol;
                  bgaer = Cland*exp(-(zm(i,k)/Hland))
                  Ntot(i,k) = max(bgaer,Ntot(i,k))
               else
                  bgaer = Cmarn*exp(-(zm(i,k)/Hmarn))
                  Ntot(i,k) = max(bgaer,Ntot(i,k))
               end if
!
               if (k == pver) then
                  Ntotb = Ntot(i,k)
               else
                  Ntotb = Ntot(i,k+1)
               end if
!
               relmod(i,k) = (( (recoef*lwcwat(i,k))/(rekappa*Ntotb))**reexp)*10000.0
               relmod(i,k) = max(4.0_r8,relmod(i,k))
               relmod(i,k) = min(20.0_r8,relmod(i,k))
               if (cld(i,k) >= 0.01) then
                  cldfrq(i,k) = 1.0
               else
                  cldfrq(i,k) = 0.0
               end if
               wrel(i,k) = relmod(i,k)*cldfrq(i,k)
               wlwc(i,k) = lwcwat(i,k)*cldfrq(i,k)
            end do
         end do
      else
         do k = 1, pver
            do i = 1, ncol
               relmod(i,k) = rel(i,k)
            end do
         end do
      end if
!
! Specify aerosol mass mixing ratio
!
      call aermix(lchnk   ,ncol    ,pnm     ,sulfmix ,aermmr  ,rh      )

      call t_startf('radcswmx')
#if (defined TRIAL_RUN_FORCING)
      call radcswmx(lchnk   ,ncol    ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aermmr  ,cld     ,clwp    ,rel     ,rei     , &
                    fice    ,eccf    ,coszrs  ,&
                    currentSolarConstant, &
                    solin, &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsutoac ,fsnirt  ,fsnrtc  ,fsnirtsq, &
                    fsns    ,fsnsc   ,fsdsc   ,fsusc   ,fsds    , &
                    sols    ,soll    ,solsd   ,solld   ,fsus    )
#else
      call radcswmx(lchnk   ,ncol    ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aermmr  ,cld     ,clwp    ,rel     ,rei     , &
                    fice    ,eccf    ,coszrs  ,scon    ,solin   , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsutoac ,fsnirt  ,fsnrtc  ,fsnirtsq, &
                    fsns    ,fsnsc   ,fsdsc   ,fsusc   ,fsds    , &
                    sols    ,soll    ,solsd   ,solld   ,fsus    )
#endif
      call t_stopf('radcswmx')

      call outfld('AERMMR  ',aermmr, pcols,lchnk)
      call outfld('REL     ',relmod ,pcols,lchnk)

      if ( indirect ) then
         call outfld('MSO4    ',so4mass,pcols,lchnk)
         call outfld('LWC     ',lwcwat ,pcols,lchnk)
         call outfld('CLDFRQ  ',cldfrq ,pcols,lchnk)
         call outfld('WREL    ',wrel   ,pcols,lchnk)
         call outfld('WLWC    ',wlwc   ,pcols,lchnk)
      end if
! -- tls ---------------------------------------------------------------2
!
! Convert units of shortwave fields needed by rest of model from CGS to MKS
!
      do i=1,ncol
         solin(i) = solin(i)*1.e-3
         fsds(i)  = fsds(i)*1.e-3
         fsus(i)  = fsus(i)*1.0e-3 ! added by DONG Li
         fsnirt(i)= fsnirt(i)*1.e-3
         fsnrtc(i)= fsnrtc(i)*1.e-3
         fsnirtsq(i)= fsnirtsq(i)*1.e-3
         fsnt(i)  = fsnt(i) *1.e-3
         fsns(i)  = fsns(i) *1.e-3
         fsntc(i) = fsntc(i)*1.e-3
         fsnsc(i) = fsnsc(i)*1.e-3
         fsdsc(i) = fsdsc(i)*1.e-3
         fsusc(i) = fsusc(i)*1.e-3 ! added by DONG Li
         fsntoa(i)=fsntoa(i)*1.e-3
         fsntoac(i)=fsntoac(i)*1.e-3
         fsutoac(i)=fsutoac(i)*1.e-3 ! added by DONG Li
      end do
!
! Dump shortwave radiation information to history tape buffer (diagnostics)
!
      ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
      call outfld('QRS     ',ftem  ,pcols,lchnk)
      call outfld('SOLIN   ',solin ,pcols,lchnk)
      call outfld('FSDS    ',fsds  ,pcols,lchnk)
      call outfld('FSUS    ', fsus, pcols, lchnk) ! added by DONG Li
      call outfld('FSNIRTOA',fsnirt,pcols,lchnk)
      call outfld('FSNRTOAC',fsnrtc,pcols,lchnk)
      call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
      call outfld('FSNT    ',fsnt  ,pcols,lchnk)
      call outfld('FSNS    ',fsns  ,pcols,lchnk)
      call outfld('FSNTC   ',fsntc ,pcols,lchnk)
      call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
      call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
      call outfld('FSUSC   ',fsusc ,pcols,lchnk) ! added by DONG Li
      call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
      call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
      call outfld('FSUTOAC ',fsutoac,pcols,lchnk) ! added by DONG Li
      call outfld('SOLS    ',sols  ,pcols,lchnk)
      call outfld('SOLL    ',soll  ,pcols,lchnk)
      call outfld('SOLSD   ',solsd ,pcols,lchnk)
      call outfld('SOLLD   ',solld ,pcols,lchnk)
!
   end if
!
! Longwave radiation computation
!
   if (dolw) then
!
! Convert upward longwave flux units to CGS
!
      do i=1,ncol
         lwupcgs(i) = lwup(i)*1000.
      end do
!
! Do longwave computation. If not implementing greenhouse gas code then
! first specify trace gas mixing ratios. If greenhouse gas code then:
!  o ixtrcg   => indx of advected n2o tracer
!  o ixtrcg+1 => indx of advected ch4 tracer
!  o ixtrcg+2 => indx of advected cfc11 tracer
!  o ixtrcg+3 => indx of advected cfc12 tracer
!
      if (trace_gas) then
         call cnst_get_ind('N2O'  , in2o)
         call cnst_get_ind('CH4'  , ich4)
         call cnst_get_ind('CFC11', if11)
         call cnst_get_ind('CFC12', if12)
         call t_startf("radclwmx")
         call radclwmx(lchnk   ,ncol    ,                            &
                       lwupcgs ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr     ,pnm     ,pmln    ,piln    ,          &
                       qm1(1,1,in2o)    ,qm1(1,1,ich4)    ,          &
                       qm1(1,1,if11)    ,qm1(1,1,if12)    ,          &
                       cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut    ,flutc   )

         call t_stopf("radclwmx")
      else
         call trcmix(lchnk   ,ncol    , &
                     pmid    ,n2o     ,ch4     ,                     &
                     cfc11   ,cfc12   )

         call t_startf("radclwmx")
         call radclwmx(lchnk     ,ncol    ,                            &
                       lwupcgs   ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr       ,pnm     ,pmln    ,piln    ,          &
                       n2o       ,ch4     ,cfc11   ,cfc12   ,          &
                       cld       ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       flns      ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut      ,flutc   )
         call t_stopf("radclwmx")
      endif
!
! Convert units of longwave fields needed by rest of model from CGS to MKS
!
      do i=1,ncol
         flnt(i)  = flnt(i)*1.e-3
         flut(i)  = flut(i)*1.e-3
         flutc(i) = flutc(i)*1.e-3
         flns(i)  = flns(i)*1.e-3
         flntc(i) = flntc(i)*1.e-3
         flnsc(i) = flnsc(i)*1.e-3
         flwds(i) = flwds(i)*1.e-3
         lwcf(i)=flutc(i) - flut(i)
         swcf(i)=fsntoa(i) - fsntoac(i)
      end do
!
! Dump longwave radiation information to history tape buffer (diagnostics)
!
      call outfld('QRL     ',qrl/cpair ,pcols,lchnk)
      call outfld('FLNT    ',flnt  ,pcols,lchnk)
      call outfld('FLUT    ',flut  ,pcols,lchnk)
      call outfld('FLUTC   ',flutc ,pcols,lchnk)
      call outfld('FLNTC   ',flntc ,pcols,lchnk)
      call outfld('FLNS    ',flns  ,pcols,lchnk)
      call outfld('FLNSC   ',flnsc ,pcols,lchnk)
      call outfld('LWCF    ',lwcf  ,pcols,lchnk)
      call outfld('SWCF    ',swcf  ,pcols,lchnk)
!
   end if
!
   return
end subroutine radctl
