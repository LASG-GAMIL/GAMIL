#include <misc.h>
#include <preproc.h>

subroutine iniTimeConst (vegxy)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize time invariant model variables
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: iniTimeConst.F90,v 1.4.6.6.6.1 2002/10/03 20:07:39 erik Exp $
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use infnan
  use clm_varder
  use clm_varpar   , only : nlevsoi, nlevlak, lsmlon, lsmlat, maxpatch, &
                            npatch_urban, npatch_lake, npatch_wet, npatch_gla 
  use clm_varmap   , only : begpatch, endpatch, patchvec
  use clm_varcon   , only : istsoil, istice, istdlak, istslak, istwet, spval
  use pft_varcon   , only : ncorn, nwheat, roota_par, rootb_par,  &
                            z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                            qe25, vcmx25, mp, c3psn
  use clm_varsur   , only : soic2d, sand3d, clay3d, latixy, longxy, &
                            zlak, dzlak, zsoi, dzsoi, zisoi
  use clm_varctl   , only : nsrest
  use time_manager , only : get_step_size
  use shr_const_mod, only : SHR_CONST_PI
  use spmdMod      , only : masterproc
  implicit none

! ------------------------ arguments ---------------------------------
  integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) !vegetation type
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  integer  :: i,j,k,l,m,ib    !indices
  real(r8) :: pi              !3.159...
  integer  :: ivt             !vegetation type index
  real(r8) :: bd              !bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm             !mineral conductivity
  real(r8) :: xksat           !maximum hydraulic conductivity of soil [mm/s]
  real(r8) :: scalez = 0.025  !Soil layer thickness discretization (m)     
  real(r8) :: hkdepth = 0.5   !Length scale for Ksat decrease (m)
  real(r8) :: sand(1:nlevsoi,begpatch:endpatch) !temporary sand
  real(r8) :: clay(1:nlevsoi,begpatch:endpatch) !temporary clay
! --------------------------------------------------------------------

! Initialize local variables for error checking

  sand(:,:) = inf
  clay(:,:) = inf

! Determine itypveg, isoicol, itypwat, sand, and clay from 2-d surface type

  do k = begpatch, endpatch

     i = patchvec%ixy(k)                 !longitude index
     j = patchvec%jxy(k)                 !latitude index
     m = patchvec%mxy(k)                 !patch index

     if (patchvec%wtxy(k) > 0.) then     !valid subgrid patch

        clm(k)%kpatch = k
        clm(k)%itypveg    = vegxy(i,j,m)
        clm(k)%isoicol    = soic2d(i,j)

        if      (m == npatch_urban) then !urban, from pcturb
           clm(k)%itypwat = istsoil
           do l = 1,nlevsoi
              sand(l,k) = sand3d(i,j,l)
              clay(l,k) = clay3d(i,j,l)
           end do
        else if (m == npatch_lake) then  !deep lake, from pctlak
           clm(k)%itypwat = istdlak
           do l = 1,nlevsoi
              sand(l,k) = 0._r8
              clay(l,k) = 0._r8
           end do
        else if (m == npatch_wet) then   !wetland, from pctwet
           clm(k)%itypwat = istwet
           do l = 1,nlevsoi
              sand(l,k) = 0._r8
              clay(l,k) = 0._r8
           end do
        else if (m == npatch_gla) then   !glacier, from pctgla
           clm(k)%itypwat = istice
           do l = 1,nlevsoi
              sand(l,k) = sand3d(i,j,l)
              clay(l,k) = clay3d(i,j,l)
           end do
        else                             !soil 
           clm(k)%itypwat = istsoil
           do l = 1, nlevsoi
              sand(l,k) = sand3d(i,j,l)
              clay(l,k) = clay3d(i,j,l)
           end do
        end if

     end if

  end do

! Tag lake points

  do k = begpatch, endpatch
     if (clm(k)%itypwat==istdlak .or. clm(k)%itypwat==istslak) then
        clm(k)%lakpoi = .true.
     else
        clm(k)%lakpoi = .false.
     end if
  end do

! Determine latitudes and longitudes

  pi = SHR_CONST_PI
  do k = begpatch, endpatch
     i = patchvec%ixy(k)     
     j = patchvec%jxy(k)     
     clm(k)%lat    = latixy(i,j) * pi/180.
     clm(k)%lon    = longxy(i,j) * pi/180.
     clm(k)%latdeg = latixy(i,j) 
     clm(k)%londeg = longxy(i,j) 
  end do

! Initialize TUNABLE constants

  do k = begpatch,endpatch
     clm(k)%zlnd   = 0.01    !Roughness length for soil [m]
     clm(k)%zsno   = 0.0024  !Roughness length for snow [m]
     clm(k)%csoilc = 0.004   !Drag coefficient for soil under canopy [-]
     clm(k)%dewmx  = 0.1     !maximum dew
     clm(k)%wtfact = 0.3     !Fraction of model area with high water table 
     clm(k)%capr   = 0.34    !Tuning factor to turn first layer T into surface T  
     clm(k)%cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
     clm(k)%ssi    = 0.033   !Irreducible water saturation of snow
     clm(k)%wimp   = 0.05    !Water impremeable if porosity less than wimp
     clm(k)%pondmx = 10.0    !Ponding depth (mm)
     clm(k)%smpmax = -1.5e5  !Wilting point potential in mm
     clm(k)%smpmin = -1.e8   !Restriction for min of soil poten. (mm)
     clm(k)%trsmx0 =  2.e-4  !Max transpiration for moist soil+100% veg. [mm/s]
  end do

! Define layer structure for soil and lakes 
! Vertical profile of snow is initialized in routine iniTimeVar

! Check that lake and soil levels are the same for now

  if (nlevlak /= nlevsoi) then
     write(6,*)'number of soil levels and number of lake levels must be the same'
     write(6,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
     call endrun
  endif

! Lake layers (assumed same for all lake patches)
  
  dzlak(1) = 1.               
  dzlak(2) = 2.               
  dzlak(3) = 3.               
  dzlak(4) = 4.               
  dzlak(5) = 5.
  dzlak(6) = 7.
  dzlak(7) = 7.
  dzlak(8) = 7.
  dzlak(9) = 7.
  dzlak(10)= 7.

  zlak(1) =  0.5
  zlak(2) =  1.5
  zlak(3) =  4.5
  zlak(4) =  8.0
  zlak(5) = 12.5
  zlak(6) = 18.5
  zlak(7) = 25.5
  zlak(8) = 32.5
  zlak(9) = 39.5
  zlak(10)= 46.5

! Soil layers and interfaces (assumed same for all non-lake patches)
! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil

  do j = 1, nlevsoi
     zsoi(j) = scalez*(exp(0.5*(j-0.5))-1.)    !node depths
  enddo
  
  dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
  do j = 2,nlevsoi-1
     dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1)) 
  enddo
  dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)
  
  zisoi(0) = 0.
  do j = 1, nlevsoi-1
     zisoi(j) = 0.5*(zsoi(j)+zsoi(j+1))         !interface depths
  enddo
  zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5*dzsoi(nlevsoi)

! Put soil and lake layers into derived type components

  do k = begpatch,endpatch
     if (clm(k)%itypwat == istdlak) then        !assume all lakes are deep lakes
        clm(k)%z(1:nlevlak)  = zlak(1:nlevlak)
        clm(k)%dz(1:nlevlak) = dzlak(1:nlevlak)
     else if (clm(k)%itypwat == istslak) then   !shallow lake (not used)
        clm(k)%dz(1:nlevlak) = NaN
        clm(k)%z(1:nlevlak)  = NaN
     else                                       !soil, ice, wetland
        clm(k)%z(1:nlevsoi)  = zsoi(1:nlevsoi)
        clm(k)%dz(1:nlevsoi) = dzsoi(1:nlevsoi)
        clm(k)%zi(0:nlevsoi) = zisoi(0:nlevsoi) 
     endif
  end do
  
! Initialize root fraction (computing from surface, d is depth in meter):
! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with beta & d_obs
! given in Zeng et al. (1998).

  do k = begpatch, endpatch
     if (.not. clm(k)%lakpoi) then
        ivt = clm(k)%itypveg
        do j = 1, nlevsoi-1
           clm(k)%rootfr(j) = .5*( exp(-roota_par(ivt)*clm(k)%zi(j-1))  &
                                 + exp(-rootb_par(ivt)*clm(k)%zi(j-1))  &
                                 - exp(-roota_par(ivt)*clm(k)%zi(j  ))  &
                                 - exp(-rootb_par(ivt)*clm(k)%zi(j  )) )
        end do
        clm(k)%rootfr(nlevsoi) = .5*( exp(-roota_par(ivt)*clm(k)%zi(nlevsoi-1))  &
                                    + exp(-rootb_par(ivt)*clm(k)%zi(nlevsoi-1)) )
     else
        clm(k)%rootfr(1:nlevsoi) = spval
     end if
  end do

! Initialize soil thermal and hydraulic properties 

  do k = begpatch, endpatch
     if (clm(k)%itypwat == istsoil) then                !soil
        do j = 1, nlevsoi
           clm(k)%bsw(j)    = 2.91 + 0.159*clay(j,k)
           clm(k)%watsat(j) = 0.489 - 0.00126*sand(j,k) 
           xksat            = 0.0070556 *( 10.**(-0.884+0.0153*sand(j,k)) ) ! mm/s
           clm(k)%hksat(j)  = xksat * exp(-clm(k)%zi(j)/hkdepth)
           clm(k)%sucsat(j) = 10. * ( 10.**(1.88-0.0131*sand(j,k)) )
           tkm              = (8.80*sand(j,k)+2.92*clay(j,k))/(sand(j,k)+clay(j,k))          ! W/(m K)
           bd               = (1.-clm(k)%watsat(j))*2.7e3
           clm(k)%tkmg(j)   = tkm ** (1.- clm(k)%watsat(j))           
           clm(k)%tksatu(j) = clm(k)%tkmg(j)*0.57**clm(k)%watsat(j)   
           clm(k)%tkdry(j)  = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)  
           clm(k)%csol(j)   = (2.128*sand(j,k)+2.385*clay(j,k))/ (sand(j,k)+clay(j,k))*1.e6  ! J/(m3 K)
        end do
     else                                               !ice, lakes, wetlands
        do j = 1, nlevsoi
           clm(k)%bsw(j)    =  spval
           clm(k)%watsat(j) =  spval
           clm(k)%hksat(j)  =  spval
           clm(k)%sucsat(j) =  spval
           clm(k)%tkmg(j)   =  spval
           clm(k)%tksatu(j) =  spval
           clm(k)%tkdry(j)  =  spval
           clm(k)%csol(j)   =  spval
        end do
     end if
  end do

! Initialize clm derived type components from pft_varcon 

  do k = begpatch, endpatch
     ivt = clm(k)%itypveg
     clm(k)%z0mr    = z0mr(ivt)
     clm(k)%displar = displar(ivt)
     clm(k)%dleaf  = dleaf(ivt)
     clm(k)%xl     = xl(ivt)
     do ib = 1,numrad 
        clm(k)%rhol(ib) = rhol(ivt,ib)
        clm(k)%rhos(ib) = rhos(ivt,ib)
        clm(k)%taul(ib) = taul(ivt,ib)
        clm(k)%taus(ib) = taus(ivt,ib)
     end do
     clm(k)%qe25   = qe25(ivt)      ! quantum efficiency at 25c (umol co2 / umol photon)
     clm(k)%vcmx25 = vcmx25(ivt)    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
     clm(k)%mp     = mp(ivt)        ! slope for conductance-to-photosynthesis relationship
     clm(k)%c3psn  = c3psn(ivt)     ! photosynthetic pathway: 0. = c4, 1. = c3
  end do

! Initialize other misc derived type components - note that for a 
! restart run, dtime is set in routine restrd()

  if (nsrest == 0) then
     do k = begpatch, endpatch
        clm(k)%dtime = get_step_size()
     end do
  endif
     
! Write out level info 

  if (masterproc) then
     write(6,*)
     write(6,30)
     do j = 1,nlevlak
        write(6,40)zlak(j),dzlak(j)
     end do
     write(6,*)
     write(6,35)
     do j = 1,nlevsoi	
	write(6,45)zsoi(j),dzsoi(j),zisoi(j)
     end do
     write(6,50)
     write(6,*)
  endif
     
30 format(' ',' lake levels ',' lake thickness(m)')
35 format(' ',' soil levels ',' soil thickness(m)',' soil interfaces(m)')
40 format(' ',2(f7.3,8x))
45 format(' ',3(f7.3,8x))
50 format(' ','Note: top level soil interface is set to 0')

  return
end subroutine iniTimeConst 
