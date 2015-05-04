#include <misc.h>
#include <params.h>

module inidat
!BOP
!
! !MODULE: inidat --- dynamics-physics coupling module
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use chemistry, only: chem_init_mix

! !PUBLIC MEMBER FUNCTIONS:
   public read_inidat, copy_inidat

! !PUBLIC DATA MEMBERS:
   real(r8), allocatable :: ps_tmp(:,:)
   real(r8), allocatable :: u3s_tmp(:,:,:)
   real(r8), allocatable :: v3s_tmp(:,:,:)
   real(r8), allocatable :: uv_local(:,:,:)
   real(r8), allocatable :: t3_tmp(:,:,:)
   real(r8), allocatable :: q3_tmp(:,:,:,:)
   real(r8), allocatable :: q3_local(:,:,:,:)
   real(r8), allocatable :: qcwat_tmp(:,:,:)
   real(r8), allocatable :: lcwat_tmp(:,:,:)
   real(r8), allocatable :: phis_tmp(:,:)
   real(r8), allocatable :: landfrac_tmp(:,:)
   real(r8), allocatable :: landm_tmp(:,:)
   real(r8), allocatable :: sgh_tmp(:,:)
   real(r8), allocatable :: ts_tmp(:,:)
   real(r8), allocatable :: tsice_tmp(:,:)
   real(r8), allocatable :: tssub_tmp(:,:,:)
   real(r8), allocatable :: sicthk_tmp(:,:)
   real(r8), allocatable :: snowhice_tmp(:,:)

   real(r8) zgsint_tmp

!
! !DESCRIPTION:
!
!      This module provides
!
!      \begin{tabular}{|l|l|} \hline \hline
!        read\_inidat    &   \\ \hline
!        copy\_inidat    &   \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   YY.MM.DD   ?????      Creation
!   00.06.01   Grant      First attempt at modifying for LRDC
!   01.10.01   Lin        Various revisions
!   01.01.15   Sawyer     Bug fixes for SPMD mode
!   01.03.26   Sawyer     Added ProTeX documentation
!   02.04.04   Sawyer     Removed comspe
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: read_inidat --- read initial dataset
!
! !INTERFACE:
   subroutine read_inidat

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use pspect
      use rgrid
      use comsrf,    only: plevmx,srfflx_state
      use commap
      use physconst, only: gravit
      use history, only: fillvalue
      use constituents, only: pcnst, pnats, cnst_name, qmin
      use tracers, only: nusr_adv, nusr_nad, ixuadv, ixunad, ixcldw

      implicit none

      include 'netcdf.inc'

!------------------------------Commons----------------------------------

#include <comctl.h>
#include <comqfl.h>
#include <comlun.h>
#include <perturb.h>

! !DESCRIPTION:
!
!   Read initial dataset and spectrally truncate as appropriate.
!
! !REVISION HISTORY:
!
!   00.06.01   Grant      First attempt at modifying for LRDC
!   00.10.01   Lin        Various revisions
!   01.01.15   Sawyer     Bug fixes for SPMD mode
!   01.03.09   Eaton      Modifications
!   01.03.26   Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer i,j,k,m,lat       ! grid and constituent indices
      integer ihem              ! hemisphere index
      real(r8) pdelb(plond,plev)! pressure diff between interfaces
      real(r8) pertval          ! perturbation value
      real(r8) zgssum           ! partial sums of phis
      integer ii, ic
!
! Netcdf related variables
!
      integer lonsiz, latsiz, levsiz ! Dimension sizes
      integer londimid, levdimid, latdimid ! Dimension ID's
      integer tid, qid  ! Variable ID's
      integer tracid(pcnst+pnats) ! Variable ID's
      integer phisid, sghid, psid ! Variable ID's
      integer landmid
#ifndef COUP_CSM
      integer tsid, ts1id, ts2id, ts3id, ts4id,tsiceid ! Variable ID's
#endif
      integer sicid,  snowhiceid           ! Variable ID's
      integer landfracid        ! Variable ID's
      integer usid, vsid
      integer strt2d(3)         ! start lon, lat, time indices for netcdf 2-d
      integer strt3d(4)         ! start lon, lev, lat, time for netcdf 3-d
      data strt2d/3*1/          ! Only index 2 will ever change
      data strt3d/4*1/          ! Only indices 2,3 will ever change

      integer cnt2d(3)          ! lon, lat, time counts for netcdf 2-d
      integer cnt3d(4)          ! lon, lat, lev, time counts for netcdf 2-d
      data cnt2d/plon,1,1/      ! 2-d arrs: Always grab only a "plon" slice
      data cnt3d/plon,plev,plat,1/ ! 3-d arrs: Always grab a full time slice

      integer ndims2d           ! number of dimensions
      integer dims2d(NF_MAX_VAR_DIMS) ! variable shape
      integer ndims3d           ! number of dimensions
      integer dims3d(NF_MAX_VAR_DIMS) ! variable shape
      integer tmptype
      integer natt, ret, attlen ! netcdf return values
      logical phis_hires        ! true => PHIS came from hi res topo
      real(r8) arr3d(plon,plev,plat)

      character*(NF_MAX_NAME) tmpname
      character*256 text
      character*80 trunits      ! tracer untis

      real(r8) splon_arr3d(plon,plev,plat)
!     real(r8) splat_arr3d(plon,plev,splat)   ! Glenn Grant's original code
      real(r8) splat_arr3d(plon,plev,plat-1)  ! temporary patch until splat = plat-1

      integer slatid, slatdimid, slatsiz
      integer slonid, slondimid, slonsiz
      integer cnt3dus(4)        ! index counts for netcdf U staggered grid
      integer cnt3dvs(4)        ! index counts for netcdf V staggered grid

!      data cnt3dus/plon,plev,splat,1/ ! 3-d arrs: Always grab a full time slice
! SJL
      integer platm1
      parameter (platm1=plat-1)
      data cnt3dus/plon,plev,platm1,1/ ! temporary patch
      data cnt3dvs/plon,plev,plat,1/ ! 3-d arrs: Always grab a full time slice

!
!-----------------------------------------------------------------------
! Allocate memory for temporary arrays
!-----------------------------------------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter
!
      allocate ( ps_tmp(plond,plat) )
      allocate ( u3s_tmp(plon,plat,plev) )
      allocate ( v3s_tmp(plon,plat,plev) )
      allocate ( t3_tmp(plond,plev,plat) )
      allocate ( q3_tmp(plond,plev,pcnst+pnats,plat) )
      allocate ( qcwat_tmp(plond,plev,plat) )
      allocate ( lcwat_tmp(plond,plev,plat) )
      allocate ( phis_tmp(plond,plat) )
      allocate ( landm_tmp(plond,plat) )
      allocate ( sgh_tmp(plond,plat) )
      allocate ( ts_tmp(plond,plat) )
      allocate ( tsice_tmp(plond,plat) )
      allocate ( tssub_tmp(plond,plevmx,plat) )
      allocate ( sicthk_tmp(plond,plat) )
      allocate ( snowhice_tmp(plond,plat) )
      allocate ( landfrac_tmp(plond,plat) )
!
!-----------------------------------------------------------------------
! Read in input variables
!-----------------------------------------------------------------------

      if (masterproc) then
!
! Get dimension IDs and lengths
!
         call wrap_inq_dimid  (ncid_ini, 'lat', latdimid)
         call wrap_inq_dimlen (ncid_ini, latdimid, latsiz)
         call wrap_inq_dimid  (ncid_ini, 'lev', levdimid)
         call wrap_inq_dimlen (ncid_ini, levdimid, levsiz)
         call wrap_inq_dimid  (ncid_ini, 'lon', londimid)
         call wrap_inq_dimlen (ncid_ini, londimid, lonsiz)
         call wrap_inq_dimid  (ncid_ini, 'slat', slatdimid)
         call wrap_inq_dimlen (ncid_ini, slatdimid, slatsiz)
         call wrap_inq_dimid  (ncid_ini, 'slon', slondimid)
         call wrap_inq_dimlen (ncid_ini, slondimid, slonsiz)

!
! Get variable id's
! Check that all tracer units are in mass mixing ratios
!
!         call wrap_inq_varid (ncid_ini, 'U'   , uid)
!         call wrap_inq_varid (ncid_ini, 'V'   , vid)

         call wrap_inq_varid (ncid_ini, 'slat', slatid)
         call wrap_inq_varid (ncid_ini, 'slon', slonid)
         call wrap_inq_varid (ncid_ini, 'US'  , usid)
         call wrap_inq_varid (ncid_ini, 'VS'  , vsid)

         call wrap_inq_varid (ncid_ini, 'T'   , tid)
         call wrap_inq_varid (ncid_ini, 'Q'   , qid)
         call wrap_inq_varid (ncid_ini, 'PS'  , psid)
         call wrap_inq_varid (ncid_ini, 'PHIS', phisid)
         call wrap_inq_varid (ncid_ini, 'SGH' , sghid)
         call wrap_inq_varid (ncid_ini, 'LANDM', landmid)
#ifndef COUP_CSM
!
! For land-fraction check if the variable name LANDFRAC is on the dataset if not assume FLAND
!
         if ( nf_inq_varid(ncid_ini, 'LANDFRAC', landfracid ) == NF_NOERR ) then
            call wrap_inq_varid (ncid_ini, 'LANDFRAC', landfracid)
         else
            call wrap_inq_varid (ncid_ini, 'FLAND', landfracid)
         end if
         call wrap_inq_varid (ncid_ini, 'TS', tsid)
         call wrap_inq_varid (ncid_ini, 'TSICE', tsiceid)
         call wrap_inq_varid (ncid_ini, 'TS1', ts1id)
         call wrap_inq_varid (ncid_ini, 'TS2', ts2id)
         call wrap_inq_varid (ncid_ini, 'TS3', ts3id)
         call wrap_inq_varid (ncid_ini, 'TS4', ts4id)
         call wrap_inq_varid (ncid_ini, 'SNOWHICE', snowhiceid)
#ifdef COUP_SOM
         call wrap_inq_varid (ncid_ini, 'SICTHK', sicid)
#endif
#endif
         if (readtrace) then
            do m=2,pcnst+pnats
               call wrap_inq_varid (NCID_INI,cnst_name(m), tracid(m))
               call wrap_get_att_text (NCID_INI,tracid(m),'units', trunits)
               if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
                  write(6,*)'INIDAT: tracer units for tracer = ', &
                            cnst_name(m),' must be in KG/KG'
                  call endrun
               endif
            end do
         end if
!
! Check dimension ordering for one 2-d and one 3-d field.
! Assume other arrays of like rank will have dimensions ordered the same.
!
         call wrap_inq_var (ncid_ini, psid, tmpname, tmptype, &
                            ndims2d, dims2d, natt)
         if (dims2d(1).ne.londimid .or. dims2d(2).ne.latdimid .or. &
             ndims2d.gt.3) then
            write(6,*)'INIDAT: Bad number of dims or ordering on 2d fld'
            call endrun
         end if
!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
         ret = nf_inq_attlen (ncid_ini, phisid, 'from_hires', attlen)
         if (ret.eq.NF_NOERR .and. attlen.gt.256) then
            write(6,*)'INIDAT: from_hires attribute length is too long'
            call endrun
         end if
         ret = nf_get_att_text (ncid_ini, phisid, 'from_hires', text)
         if (ret.eq.NF_NOERR .and. text(1:4).eq.'true') then
            phis_hires = .true.
!            write(6,*)'INIDAT: Will filter input PHIS: attribute ', &
!                      'from_hires is true'
         else
            phis_hires = .false.
!            write(6,*)'INIDAT: Will not filter input PHIS: attribute ', &
!                      'from_hires is either false or not present'
         end if
!
! Read in 2d fields.
! For stand alone run: get surface temp and 4 (sub)surface temp fields
! For stand alone run with slab-ocean: get sea-ice thickness and snow cover
!
         do j=1,plat
            strt2d(2) = j
            if (ideal_phys .or. aqua_planet) then
               do i=1,nlon(j)
                  phis_tmp(i,j) = 0.
                  sgh_tmp (i,j) = 0.
               end do
            else
               call wrap_get_vara_realx (ncid_ini, phisid, strt2d, cnt2d, &
                                         phis_tmp(1,j))
               call wrap_get_vara_realx (ncid_ini, sghid , strt2d, cnt2d, &
                                         sgh_tmp(1,j))
            endif
            call wrap_get_vara_realx (ncid_ini, landmid, strt2d, cnt2d, &
                                      landm_tmp(1,j))
            call wrap_get_vara_realx (ncid_ini, psid, strt2d, cnt2d, &
                                      ps_tmp(1,j))
#ifndef COUP_CSM
            if (aqua_planet) then
               do i=1,nlon(j)
                  landfrac_tmp(i,j) = 0.
               end do
            else
               call wrap_get_vara_realx (ncid_ini, landfracid, strt2d, cnt2d, &
                                         landfrac_tmp(1,j))
            endif
            call wrap_get_vara_realx (ncid_ini, tsid, strt2d, cnt2d, &
                                           ts_tmp(1,j))
            call wrap_get_vara_realx (ncid_ini, tsiceid, strt2d, cnt2d, &
                                           tsice_tmp(1,j))
            call wrap_get_vara_realx (ncid_ini, ts1id, strt2d, cnt2d, &
                                           tssub_tmp(1,1,j))
            call wrap_get_vara_realx (ncid_ini, ts2id, strt2d, cnt2d, &
                                           tssub_tmp(1,2,j))
            call wrap_get_vara_realx (ncid_ini, ts3id, strt2d, cnt2d, &
                                           tssub_tmp(1,3,j))
            call wrap_get_vara_realx (ncid_ini, ts4id, strt2d, cnt2d, &
                                           tssub_tmp(1,4,j))
!
! Set sea-ice thickness and snow cover:
!
#ifdef COUP_SOM
            call wrap_get_vara_realx(ncid_ini, sicid, strt2d, cnt2d, sicthk_tmp(1,j))
#endif
            call wrap_get_vara_realx(ncid_ini, snowhiceid, strt2d, cnt2d, snowhice_tmp(1,j))
#endif
         end do
!
! Read in 3d fields.

! Staggered grid variables and transpose

         call wrap_get_vara_realx(ncid_ini,usid, strt3d, cnt3dus, splat_arr3d)

         do k = 1, plev
!
! SJL: initialize j=1 because later on u3s_tmp will be copied to u3s using f90 array syntax
               do i = 1, plon
                  u3s_tmp(i,1,k) = fillvalue
               enddo
!
            do j = 1, plat-1
               do i = 1, plon
                  u3s_tmp(i,j+1,k) = splat_arr3d(i,k,j)
               enddo
            enddo
         enddo

         call wrap_get_vara_realx(ncid_ini,vsid, strt3d, cnt3dvs, splon_arr3d)

         do k = 1, plev
            do j = 1, plat
               do i = 1, plon
                  v3s_tmp(i,j,k) = splon_arr3d(i,k,j)
               enddo
            enddo
         enddo


         call wrap_get_vara_realx(ncid_ini, tid, strt3d, cnt3d, arr3d)
         t3_tmp(:plon,:plev,:plat) = arr3d(:,:,:)

         call wrap_get_vara_realx(ncid_ini, qid, strt3d, cnt3d, arr3d)
         q3_tmp(:plon,:plev,1,:plat) = arr3d(:,:,:)

! Initialize tracers if not read in from input data.
! Initialize all user tracers (advected and non-advectec to 0.)

         if (readtrace) then
            do m=2,pcnst+pnats
               call wrap_get_vara_realx(ncid_ini, tracid(m), strt3d, cnt3d, arr3d)
               q3_tmp(:plon,:plev,m,:plat) = arr3d(:,:,:)
            end do
         endif
!
! Add random perturbation to temperature if required
!
         if (pertlim.ne.0.0) then
            write(6,*)'INIDAT: Adding random perturbation bounded by +/-', &
                      pertlim,' to initial temperature field'
            do lat=1,plat
               do k=1,plev
                  do i=1,nlon(lat)
                     call random_number (pertval)
                     pertval = 2.*pertlim*(0.5 - pertval)
                     t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1. + pertval)
                  end do
               end do
            end do
         endif
!
! Initialize tracers if not read in from input data.
! Initialize all user tracers (advected and non-advectec to 0.)
! Ensure sufficient constituent concentration at all gridpoints
!
         if (.not. readtrace) then
            do lat=1,plat
               q3_tmp(:plon,:plev,ixcldw,lat) = 0.
               if (nusr_adv .gt. 0) then
                  q3_tmp(:,:plev,ixuadv:ixuadv+nusr_adv-1,lat) = 0.
               endif
               if (nusr_nad .gt. 0) then
                  q3_tmp(:,:plev,ixunad:ixunad+nusr_nad-1,lat) = 0.
               end if
               if (trace_gas) then
                  if (doRamp_ghg ) call ramp_ghg
                  call chem_init_mix(lat, ps_tmp(1,lat), q3_tmp(1,1,1,lat), nlon(lat))
               endif
               if (trace_test1 .or. trace_test2 .or. trace_test3) then
                  call initesttr( q3_tmp(1,1,1,lat) ,nlon(lat))
               endif
            end do
         endif

         do lat=1,plat
            call qneg3('INIDAT  ',lat   ,nlon(lat),plond   ,plev    , &
                       pcnst+pnats,qmin ,q3_tmp(1,1,1,lat))
         end do

         zgsint_tmp = 0.
!
! Compute integrals of mass, moisture, and geopotential height
!
!gg  Integrals of mass and moisture should be unnecessary in Lin-Rood dynamics
!gg  because they are conserved. What's left is the global geopotential...
!gg  Dunno if that's necessary or not, so I left it in.

         do lat = 1, plat
!
! Accumulate average mass of atmosphere
!
            zgssum = 0.
            do i=1,nlon(lat)
               zgssum = zgssum + phis_tmp(i,lat)
            end do
            zgsint_tmp = zgsint_tmp + w(lat)*zgssum/nlon(lat)
         end do                  ! end of latitude loop
!
! Normalize average height
!
         zgsint_tmp = zgsint_tmp*.5/gravit
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
!
! SJL:
         tmass0 = 98222./gravit
!
! WS: Average pole information moved here: treat the global arrays
!
!-----------------------------------------------------------
! Average T, PS, PHIS and Q at the poles.       The initial
! conditions *should* already have these variables averaged,
! but do it here for safety -- no harm if it's already done.
!-----------------------------------------------------------

         call xpavg(phis_tmp(1, 1), plon)
         call xpavg(phis_tmp(1,plat), plon)
         call xpavg(ps_tmp(1,   1), plon)
         call xpavg(ps_tmp(1,plat), plon)

!$omp parallel do private(i, j, k, ic)
         do k = 1, plev

            call xpavg(t3_tmp(1,k,   1), plon)
            call xpavg(t3_tmp(1,k,plat), plon)

            do ic = 1, pcnst+pnats
               call xpavg(q3_tmp(1,k,ic,   1),plon)
               call xpavg(q3_tmp(1,k,ic,plat),plon)
            enddo

         enddo

         if (ideal_phys) tmass0 = 100000./gravit
!        write(6,800) tmassf_tmp,tmass0,qmassf_tmp
!        write(6,810) zgsint_tmp
800      format('INIDAT: MASS OF INITIAL DATA BEFORE CORRECTION = ' &
                ,1p,e20.10,/,' DRY MASS WILL BE HELD = ',e20.10,/, &
                ' MASS OF MOISTURE AFTER REMOVAL OF NEGATIVES = ',e20.10)
810      format(/69('*')/'INIDAT: Globally averaged geopotential ', &
                'height = ',f16.10,' meters'/69('*')/)
      endif                     ! end of if-masterproc
!
!-----------------------------------------------------------------------
! Copy temporary arrays to model arrays
!-----------------------------------------------------------------------
!
      call copy_inidat
!
!-----------------------------------------------------------------------
! Deallocate memory for temporary arrays
!-----------------------------------------------------------------------
!
      deallocate ( ps_tmp )
      deallocate ( u3s_tmp )
      deallocate ( v3s_tmp )
      deallocate ( t3_tmp )
      deallocate ( q3_tmp )
      deallocate ( qcwat_tmp )
      deallocate ( lcwat_tmp )
      deallocate ( phis_tmp )
      deallocate ( landfrac_tmp )
      deallocate ( landm_tmp )
      deallocate ( sgh_tmp )
      deallocate ( ts_tmp )
      deallocate ( tsice_tmp )
      deallocate ( tssub_tmp )
      deallocate ( sicthk_tmp )
      deallocate ( snowhice_tmp )
!
      return
!EOC
   end subroutine read_inidat
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: copy_inidat --- Copy temporary arrays to model arrays
!
! !INTERFACE:
   subroutine copy_inidat

! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use prognostics
      use buffer
      use comsrf
      use phys_grid
      use tracers, only: ixcldw
      use phys_grid, only: get_ncols_p
#if ( defined SPMD )
      use mpishorthand
      use spmd_dyn, only : comm_y, comm_z
      use parutilitiesmodule, only: parcollective2d, BCSTOP
#endif

      implicit none
!------------------------------Commons----------------------------------
#include <comqfl.h>

! !DESCRIPTION:
!
! Copy temporary arrays to model arrays
! note that the use statements below contain the definitions
! of the model arrays
!
! !REVISION HISTORY:
!
!   00.06.01   Grant      First attempt at modifying for LRDC
!   00.10.01   Lin        Various revisions
!   00.12.02   Sawyer     Use PILGRIM to scatter data sets
!   01.03.26   Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8), allocatable :: tmpchunk3d(:,:,:)
      real(r8), allocatable :: tmpchunk(:,:)
      integer i, j, ic, k, m
      integer n
      integer ncol
      real(r8) :: pmx, pmn

!-----------------------------------------------------------------------

#if ( defined SPMD )
! dynamics variables
      if (myid_z .eq. 0) then
         call scatter( ps_tmp, strip2d, ps, comm_y )
         call scatter( phis_tmp, strip2d, phis, comm_y )
      endif
#if defined ( TWOD_YZ )
      call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, ps )
      call parcollective2d( comm_z, BCSTOP, plon, endlat-beglat+1, phis )
#endif

! Warning: beglat may be invalid for different sized staggered grids. - gg

! dynamics variables
      allocate( uv_local(plon,beglat:endlat,beglev:endlev) )
      call scatter( u3s_tmp, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               u3s(i,j,k) = uv_local(i,j,k)
            enddo
         enddo
      enddo

      call scatter( v3s_tmp, strip3dxyz, uv_local, mpicom )
!$omp parallel do private(i,j,k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               v3s(i,j,k) = uv_local(i,j,k)
            enddo
         enddo
      enddo
      deallocate(uv_local)

      call scatter( t3_tmp, strip3dxzy, t3, mpicom )
      allocate( q3_local(plon,beglev:endlev,pcnst+pnats,beglat:endlat) )

      call scatter( q3_tmp, strip3dq3old, q3_local, mpicom )
      do m=1,pcnst+pnats
!$omp parallel do private(i,j,k)
         do k=beglev,endlev
            do j=beglat,endlat
               do i=1,plon
                  q3(i,j,k,m) = q3_local(i,k,m,j)
               enddo
            enddo
         enddo
      enddo
      deallocate( q3_local )
#else

! dynamics variables
      ps(:,:) = ps_tmp(:,:)
      phis(:,:) = phis_tmp(:,:)

! dynamics variables

!$omp parallel do private(i, j, k)
      do k=beglev,endlev
         do j=beglat,endlat
            do i=1,plon
               u3s(i,j,k) = u3s_tmp(i,j,k)
               v3s(i,j,k) = v3s_tmp(i,j,k)
            enddo
         enddo
      enddo

!$omp parallel do private(i, j, k, ic)
      do j=beglat,endlat
         do k=beglev,endlev
            do i=1,plon
               t3(i,k,j) = t3_tmp(i,k,j)
            enddo
         enddo

         do ic=1,pcnst+pnats
            do k=beglev,endlev
               do i=1,plon
                  q3(i,j,k,ic) = q3_tmp(i,k,ic,j)
               enddo
            enddo
         enddo
      enddo

#endif
      allocate ( tmpchunk(pcols,begchunk:endchunk) )
      allocate ( tmpchunk3d(pcols,plevmx,begchunk:endchunk) )
! physics variables
      call scatter_field_to_chunk(1,1,1,plond,landfrac_tmp,landfrac(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,landm_tmp,landm(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,sgh_tmp,sgh(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,tsice_tmp,tsice(1,begchunk))
      call scatter_field_to_chunk(1,1,1,plond,ts_tmp,tmpchunk(1,begchunk))
      do i =begchunk,endchunk
         ncol = get_ncols_p(i)
         srfflx_state2d(i)%ts(:ncol) = tmpchunk(:ncol,i)
      end do
#ifdef COUP_SOM
      call scatter_field_to_chunk(1,1,1,plond,sicthk_tmp,sicthk(1,begchunk))
#endif
      call scatter_field_to_chunk(1,1,1,plond,snowhice_tmp,snowhice(1,begchunk))

      call scatter_field_to_chunk(1,plevmx,1,plond,tssub_tmp,tmpchunk3d)
      do i =begchunk,endchunk
         ncol = get_ncols_p(i)
         surface_state2d(i)%tssub(:ncol,:) = tmpchunk3d(:ncol,:,i)
      end do
!
!JR cloud and cloud water initialization.  Does this belong somewhere else?
!
      if (masterproc) then
         qcwat_tmp(:,:,:) = q3_tmp(:,:,1,:)
         lcwat_tmp(:,:,:) = q3_tmp(:,:,ixcldw,:)
      endif
      call scatter_field_to_chunk(1,plev,1,plond,qcwat_tmp,qcwat(1,1,begchunk,1))
      call scatter_field_to_chunk(1,plev,1,plond,lcwat_tmp,lcwat(1,1,begchunk,1))
      call scatter_field_to_chunk(1,plev,1,plond,t3_tmp,tcwat(1,1,begchunk,1))
      cld(:,:,:,1) = 0.
      do n=2,2
         do i =begchunk,endchunk
            ncol = get_ncols_p(i)
            cld(:ncol,:pver,i,n) = 0.
            qcwat(:ncol,:pver,i,n) = qcwat(:ncol,:pver,i,1)
            tcwat(:ncol,:pver,i,n) = tcwat(:ncol,:pver,i,1)
            lcwat(:ncol,:pver,i,n) = lcwat(:ncol,:pver,i,1)
         end do
      end do
!
! Global integerals
!
      if (masterproc) then
         zgsint = zgsint_tmp
      endif
!
#if ( defined SPMD )
      call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif
      deallocate ( tmpchunk )
      deallocate ( tmpchunk3d)
!EOC
   end subroutine copy_inidat
!-----------------------------------------------------------------------
end module inidat

