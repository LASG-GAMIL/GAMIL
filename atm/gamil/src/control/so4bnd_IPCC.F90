#include <misc.h>
#include <params.h>

module so4bnd_IPCC

!----------------------------------------------------------------------- 
! Purpose: 
! SO4 boundary module.  Deals with interpolating SO4 datasets.
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!
!JR Stuck this "only" business in because Compaq compiler barfed on pcnst, pnats having dual
!JR declarations when radctl.F90 gets compiled.
!
   use pmgrid,    only: plon, plat,  masterproc
   use ppgrid,    only: pcols, pver, begchunk, endchunk
   use phys_grid, only: scatter_field_to_chunk, get_ncols_p
   use shr_sys_mod, only : shr_sys_flush
   implicit none
   save
!
! Floating point data
!
   real(r8),private,allocatable,dimension(:,:,:,:) ::  sulfi
                      ! input sulfate bio mixing ratios (pcols,pver,begchunk:endchunk,2)
   real(r8),private,allocatable,dimension(:,:,:)   ::  sulf
                      ! time interpolated sulfate bio mixing ratios (pcols,pver,begchunk:endchunk)

   real(r8), private :: cdaysulfm        ! calendar day for prv. month sulfate values read in
   real(r8), private :: cdaysulfp        ! calendar day for nxt. month sulfate values read in

   real(r8),private :: bndry_m          ! calendar day for the lower time-boundary of inputdata
   real(r8),private :: bndry_p          ! calendar day for the upper time-boundary of inputdata
   logical, private :: outofbndry       ! whether current model time is out of the time-boundary

   integer,private,parameter :: nfile = 19    ! number of input data files
   integer,private,parameter :: yr_sulf(nfile) = (/ &
                      1850, 1900, 1920, 1930, 1940, &
                      1950, 1960, 1970, 1980, 2010, &
                      2020, 2030, 2040, 2050, 2060, &
                      2070, 2080, 2090, 2100     /)

   character*256,private :: file_sulf(nfile)

   integer, private,parameter :: nsample = 12   ! number of samples in each input data file
   real(r8),private           :: sample_sulf(nsample) =      (/&
                1296000.0,  3888000.0,  6480000.0,  9072000.0, &
               11664000.0, 14256000.0, 16848000.0, 19440000.0, &  
               22032000.0, 24624000.0, 27216000.0, 29808000.0 /)
!
! Integer data
!
   integer, private :: nm,np        ! Array indices for prv., nxt month sulfate data
   integer, private :: sp_m, sp_p
   integer, private :: yrm,  yrp
   integer, private :: yr_am, yr_bm, yr_ap, yr_bp, indxbp

   integer, private :: ncid_sulf_am, ncid_sulf_bm  ! sulfate dataset id
   integer, private :: ncid_sulf_ap, ncid_sulf_bp  ! sulfate dataset id
   integer, private :: sulf_id_am,   sulf_id_bm    ! netcdf id for sulfate mmr bio variable
   integer, private :: sulf_id_ap,   sulf_id_bp    ! netcdf id for sulfate mmr bio variable

   integer, private :: lonsiz     ! size of longitude dimension on sulfate dataset
   integer, private :: levsiz     ! size of level dimension on sulfate dataset
   integer, private :: latsiz     ! size of latitude dimension on sulfate dataset
   integer, private :: timsiz     ! size of time dimension on sulfate dataset
 
   character*16,  private ::  IPCC_scenario   ! name of the IPCC scenario
   character*80,  private ::  sulfdata        ! full pathname for sulfate dataset
   character*256, private ::  sulfdatafile_am ! full pathname for sulfate dataset
   character*256, private ::  sulfdatafile_bm ! full pathname for sulfate dataset
   character*256, private ::  sulfdatafile_ap ! full pathname for sulfate dataset
   character*256, private ::  sulfdatafile_bp ! full pathname for sulfate dataset

   character*80,  private,parameter ::  p1='so4.'   ! part of the name of sulfate data file
   character*80,  private,parameter ::  p2='.nc'   ! part of the name of sulfate data file 

contains

subroutine so4bndnl_IPCC( xsulfdata, xIPCC_scenario )
!----------------------------------------------------------------------- 
! 
! Purpose: Set variables from namelist input.
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------

   character*80, intent(in):: xsulfdata       ! full pathname for sulfate dataset
   character*16, intent(in):: xIPCC_scenario  ! name of the IPCC scenario

!-----------------------------------------------------------------------
   sulfdata = xsulfdata
   IPCC_scenario = xIPCC_scenario

   if (masterproc) then 
      write(6,*)'Path of time-variant sulfate dataset is: ',trim(sulfdata)
      write(6,*)'Name of IPCC_scenario is: ',               trim(IPCC_scenario)
   endif

   return
end subroutine so4bndnl_IPCC

!###############################################################################

subroutine sulfini_IPCC
!----------------------------------------------------------------------- 
! 
! Purpose: Do initial read of time-variant sulfate dataset, containing
!          sulfate mixing ratios as a function of time.  It is currently
!          required that the sulfate dataset have the *SAME* horizontal
!          and vertical resolution as the model. Therefore, ONLY a time
!          interpolation of the dataset is currently performed.
! 
!-----------------------------------------------------------------------
   use ioFileMod
   use error_messages, only: alloc_err, handle_ncerr
   use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                           is_perpetual
#if ( defined SPMD )
   use mpishorthand
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!
! Local workspace
!
   character(len=256) locfn_am, locfn_bm, locfn_ap, locfn_bp     ! local filenames

   integer londimid              ! netcdf id for longitude dimension
   integer latdimid              ! netcdf id for latitude dimension
   integer levdimid              ! netcdf id for level dimension
   integer lonid                 ! netcdf id for longitude variable
   integer latid                 ! netcdf id for latitude variable
   integer levid                 ! netcdf id for level variable
   integer timeid                ! netcdf id for time variable
   integer cnt4(4)               ! array of counts for each dimension
   integer strt4(4)              ! array of starting indices
   integer i, k, lat, n          ! longitude, level, latitude, time indices
   integer istat                 ! error return
   integer dimids(nf_max_var_dims) ! netcdf variable shape

   integer yr, mon, day, ncsec ! components of a date
   real(r8) calday                   ! current calendar day
   real(r8) caldayloc                ! calendar day (includes yr if no cycling)
   real(r8) fact1, fact2, deltat

   character*4 yrtmp

   real(r8) xsulfi_am(plon,pver,plat)  ! input sulfate bio mixing ratios
   real(r8) xsulfi_bm(plon,pver,plat)  ! input sulfate ant mixing ratios
   real(r8) xsulfi_ap(plon,pver,plat)  ! input sulfate bio mixing ratios
   real(r8) xsulfi_bp(plon,pver,plat)  ! input sulfate ant mixing ratios

   real(r8) xsulfi(plon,pver,plat,2)   ! input sulfate ant mixing ratios
!
!-----------------------------------------------------------------------
!
! Initialize time counters
!
   nm = 1
   np = 2
!
! Allocate space for data.
!
  allocate( sulfi(pcols,pver,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulfi', pcols*pver*(endchunk-begchunk+1)*2 )
  allocate( sulf(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulf',  pcols*pver*(endchunk-begchunk+1) )
!
! unit shift from second to day
!
   do n=1,nsample
      sample_sulf(n) = sample_sulf(n)/86400.0+1.0
   enddo

!!   if (masterproc) then
!!      do n=1,nsample
!!        write(*,*) n,sample_sulf(n)
!!      enddo
!!       write(*,*) '-----------------------'
!!   endif
!
! the names of rawdata files
!
   do n=1,9
      write(yrtmp,'(i4.4)') yr_sulf(n)
      file_sulf(n) = trim(p1)//yrtmp//trim(p2)
   enddo

   do n=10,19
      write(yrtmp,'(i4.4)') yr_sulf(n)
      file_sulf(n) = trim(p1)//trim(IPCC_scenario)//'.'//yrtmp//trim(p2)
   enddo

!!   if (masterproc) then
!!      do n=1,nfile
!!        write(*,*) n,trim(file_sulf(n))
!!      enddo
!!       write(*,*) '-----------------------'
!!   endif
!
! set the time_boundaries of input dataset
!
   bndry_m = yr_sulf(1    )*365.0 + sample_sulf(1) 
   bndry_p = yr_sulf(nfile)*365.0 + sample_sulf(nsample) 
  
!!   if (masterproc) then
!!       write(*,*) 'bndry_m=',bndry_m
!!       write(*,*) 'bndry_p=',bndry_p
!!   endif
!!    call endrun
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
   if (masterproc) then

     calday = get_curr_calday()
     call get_curr_date(yr, mon, day, ncsec)
     caldayloc = calday + yr*365.

!----------------------------------------------------------------
! locate the model time among the samples of sulfate data file      
!---------------------------------------------------------------
    if ( caldayloc.lt.bndry_m ) then

        outofbndry = .true.
        cdaysulfp  = bndry_m
        sp_p   = 1
        yrp    = yr_sulf(1)
        yr_ap  = yr_sulf(1)
        yr_bp  = yr_sulf(2)
        indxbp = 2
        sulfdatafile_ap = trim(sulfdata)//file_sulf(1)
        sulfdatafile_bp = trim(sulfdata)//file_sulf(2)

    elseif ( caldayloc.ge.bndry_p ) then

        outofbndry = .true.
        sp_p   = nsample 
        yrp    = yr_sulf(nfile)
        yr_ap  = yr_sulf(nfile)
        yr_bp  = yr_sulf(nfile)
        sulfdatafile_ap = trim(sulfdata)//file_sulf(nfile)

    else
        outofbndry = .false.
!
! month
!
       if ( calday.lt.sample_sulf(1) ) then
          yrm  = yr-1
          yrp  = yr
          sp_m = nsample 
          sp_p = 1
          goto 100
       elseif (calday.ge.sample_sulf(nsample) ) then
          yrm = yr
          yrp = yr+1
          sp_m = nsample
          sp_p = 1 
          goto 100
       else 
         do n=1,nsample-1
            if ( calday.ge.sample_sulf(n) .and. calday.lt.sample_sulf(n+1) ) then
               yrm  = yr
               yrp  = yr
               sp_m = n
               sp_p = n+1
               goto 100
            endif
         enddo
       endif
       write(*,*) 'SULFINI: Failed to find seconds bracketing model_time =',calday
       call endrun
100    continue
       cdaysulfm = yrm*365.0+sample_sulf(sp_m)
       cdaysulfp = yrp*365.0+sample_sulf(sp_p)
!
! year ( both yrm & yrp )
!
       do n=1,nfile-1
          yr_am = yr_sulf(n)
          yr_bm = yr_sulf(n+1)
          if ( yrm.ge.yr_am .and. yrm.lt.yr_bm ) goto 200
       enddo
       write(*,*) 'SULFINI: Failed to find years bracketing year  yrm=',yrm
       call endrun
200    continue
       sulfdatafile_am = trim(sulfdata)//trim(file_sulf(n))
       sulfdatafile_bm = trim(sulfdata)//trim(file_sulf(n+1))
 
       do n=1,nfile-1
          yr_ap = yr_sulf(n)
          yr_bp = yr_sulf(n+1)
          if ( yrp.ge.yr_ap .and. yrp.lt.yr_bp ) goto 300
       enddo
       write(*,*) 'SULFINI: Failed to find years bracketing year  yrp=',yrp
       call endrun
300    indxbp = n+1 
       sulfdatafile_ap = trim(sulfdata)//trim(file_sulf(n))
       sulfdatafile_bp = trim(sulfdata)//trim(file_sulf(n+1))

    endif

!--------------------------------------------------------------------------
! open the right data files and read in raw data
!--------------------------------------------------------------------------
    if ( outofbndry ) then

      call getfil(sulfdatafile_ap, locfn_ap)
      call wrap_open(locfn_ap, 0, ncid_sulf_ap)
!
! Obtain dimension id's
!
      call wrap_inq_dimid( ncid_sulf_ap, 'lon', londimid)
      call wrap_inq_dimid( ncid_sulf_ap, 'lat', latdimid)
      call wrap_inq_dimid( ncid_sulf_ap, 'lev', levdimid)
      call wrap_inq_dimid( ncid_sulf_ap, 'time',timeid  )
!
! Obtain size of dimensions.
! Check that horizontal and vertical dimensions are same as model's
! (check one file only and assume that all the other files are of the same setting)
!
      call wrap_inq_dimlen( ncid_sulf_ap, londimid, lonsiz   )
      if (lonsiz /= plon) then
         write(6,*)'SULFINI: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf_ap, latdimid, latsiz   )
      if (latsiz /= plat) then
         write(6,*)'SULFINI: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf_ap, levdimid, levsiz   )
      if (levsiz /= pver) then
         write(6,*)'SULFINI: levsiz=',levsiz,' must = ',pver
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf_ap, timeid, timsiz   )
      if (timsiz /= nsample) then
         write(6,*)'SULFINI: timsiz=',timsiz,' must = ', nsample 
         call endrun
      end if
!
! the variable
!
      call wrap_inq_varid( ncid_sulf_ap, 'sulfmmrbio' , sulf_id_ap )

      call wrap_inq_vardimid (ncid_sulf_ap, sulf_id_ap, dimids)
      if (dimids(1) /= londimid .and. dimids(2) /= levdimid .and. dimids(3) /= latdimid) then
         write(6,*)'SULFINI: Data must be ordered lon, lev, lat, time'
         call endrun
      end if
!
! Set up hyperslab corners
!
      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1
!
! need to read in only 1 sample of rawdata
!
      strt4(4) = sp_p 
      call wrap_get_vara_realx (ncid_sulf_ap,sulf_id_ap,strt4,cnt4,xsulfi_ap)

      do n    = 1,2
       do lat = 1,plat
        do k  = 1,pver
         do i = 1,plon
           xsulfi(i,k,lat,n)= xsulfi_ap(i,k,lat)
         enddo
        enddo
       enddo
      enddo

!---------------------------------------------------------------
! if not out of boundary
!---------------------------------------------------------------
    else   
  
      call getfil(sulfdatafile_am, locfn_am)
      call wrap_open(locfn_am, 0, ncid_sulf_am)

      call getfil(sulfdatafile_bm, locfn_bm)
      call wrap_open(locfn_bm, 0, ncid_sulf_bm)

      call getfil(sulfdatafile_ap, locfn_ap)
      call wrap_open(locfn_ap, 0, ncid_sulf_ap)

      call getfil(sulfdatafile_bp, locfn_bp)
      call wrap_open(locfn_bp, 0, ncid_sulf_bp)
! 
! Obtain sulfate mixing ratio id
!
      call wrap_inq_varid( ncid_sulf_am, 'sulfmmrbio' , sulf_id_am )
      call wrap_inq_varid( ncid_sulf_bm, 'sulfmmrbio' , sulf_id_bm )
      call wrap_inq_varid( ncid_sulf_ap, 'sulfmmrbio' , sulf_id_ap )
      call wrap_inq_varid( ncid_sulf_bp, 'sulfmmrbio' , sulf_id_bp )

      call wrap_inq_vardimid (ncid_sulf_am, sulf_id_am, dimids)
      if (dimids(1) /= londimid .and. dimids(2) /= levdimid .and. dimids(3) /= latdimid) then
         write(6,*)'SULFINI: Data must be ordered lon, lev, lat, time'
         call endrun
      end if
!
! Set up hyperslab corners
!
      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1
!
! read in the 4 sample of rawdata
!
      strt4(4) = sp_m 
      call wrap_get_vara_realx (ncid_sulf_am,sulf_id_am,strt4,cnt4,xsulfi_am)
      call wrap_get_vara_realx (ncid_sulf_bm,sulf_id_bm,strt4,cnt4,xsulfi_bm)
      strt4(4) = sp_p 
      call wrap_get_vara_realx (ncid_sulf_ap,sulf_id_ap,strt4,cnt4,xsulfi_ap)
      call wrap_get_vara_realx (ncid_sulf_bp,sulf_id_bp,strt4,cnt4,xsulfi_bp)

!
! Normal interpolation between consecutive time slices.
!
      deltat = dble( yr_bm - yr_am )
      fact1  = dble( yr_bm - yrm  )/deltat
      fact1  = dble( yrm   - yr_am)/deltat
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
        if (abs(fact1+fact2-1.).gt.1.e-6 .or.  &
            fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or. &
            fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
            write(6,*)'SULFINI--1: Bad fact1 and/or fact2=',fact1,fact2
            call endrun
        end if
!
! interpolation
!      
      do lat = 1,plat
       do k  = 1,pver
        do i = 1,plon
           xsulfi(i,k,lat,nm)= fact1*xsulfi_am(i,k,lat)+fact2*xsulfi_bm(i,k,lat)
        enddo
       enddo
      enddo
!
      deltat = dble( yr_bp - yr_ap )
      fact1  = dble( yr_bp - yrp  )/deltat
      fact1  = dble( yrp   - yr_ap)/deltat
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
        if (abs(fact1+fact2-1.).gt.1.e-6 .or.  &
            fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or. &
            fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
            write(6,*)'SULFINI--2: Bad fact1 and/or fact2=',fact1,fact2
            call endrun
        end if
!
! interpolation
!      
      do lat = 1,plat
       do k  = 1,pver
        do i = 1,plon
           xsulfi(i,k,lat,np)= fact1*xsulfi_ap(i,k,lat)+fact2*xsulfi_bp(i,k,lat)
        enddo
       enddo
      enddo

    endif    ! end if (.not.outofbndry)
!------------------------------------------------------------------------------
! master does the work above
!------------------------------------------------------------------------------
   endif
!------------------------------------------------------------------------------

#if (defined SPMD )
   call mpibcast( cdaysulfm, 1, mpir8,  0, mpicom )
   call mpibcast( cdaysulfp, 1, mpir8,  0, mpicom )
   call mpibcast( sp_m,      1, mpiint, 0, mpicom )
   call mpibcast( sp_p,      1, mpiint, 0, mpicom )
   call mpibcast( yrp,       1, mpiint, 0, mpicom )
#endif

   call scatter_field_to_chunk(1,pver,2,plon,xsulfi,sulfi)

   if (outofbndry) call scatter_field_to_chunk(1,pver,1,plon,xsulfi_ap,sulf)

   return
end subroutine sulfini_IPCC

!###############################################################################

subroutine sulfint_IPCC
!----------------------------------------------------------------------- 
!
! Purpose: Interpolate sulfate mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
!-----------------------------------------------------------------------
   use ioFileMod
   use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                           is_perpetual
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

#include <comctl.h>
#include <comlun.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Local workspace
!
   integer cnt4(4)            ! array of counts for each dimension
   integer strt4(4)           ! array of starting indices
   integer i, k, lchnk,lat    ! column, level, chunk, latitude indices
   integer ncol               ! number of columns in current chunk
   integer ntmp               ! temporary
   integer yr, mon, day       ! components of a date
   integer ncsec              ! current time of day [seconds]

   character(len=256)  locfn_bp     ! local filename

   real(r8) calday            ! current calendar day
   real(r8) caldayloc         ! calendar day (includes yr if no cycling)
   real(r8) deltat            ! time (days) between interpolating sulfate data
   real(r8) fact1, fact2      ! time interpolation factors

   real(r8) xsulfi   (plon,pver,plat)  ! input sulfate bio mixing ratios
   real(r8) xsulfi_ap(plon,pver,plat)  ! input sulfate bio mixing ratios
   real(r8) xsulfi_bp(plon,pver,plat)  ! input sulfate ant mixing ratios

!-----------------------------------------------------------------------
!
! Use year information only if a multiyear dataset
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   caldayloc = calday + yr*365.

   if ((caldayloc.lt.bndry_m).or.(caldayloc.ge.bndry_p)) then
        outofbndry = .true.
   else
        outofbndry = .false.
   endif
!-------------------------------------------------------------------------------
!  there is no need to do anything when current model time is
!   out of the time-boundary of input data files
!-------------------------------------------------------------------------------
 if (.not.outofbndry) then
!
! Initialize hyperslab corners
!
   if (masterproc) then

      strt4(1) = 1
      strt4(2) = 1
      strt4(3) = 1
      cnt4(1)  = lonsiz
      cnt4(2)  = levsiz
      cnt4(3)  = latsiz
      cnt4(4)  = 1

   endif
!
! If model time is past current forward sulfate timeslice, read in the next
! timeslice for time interpolation. 
!
   if ( caldayloc .gt. cdaysulfp ) then
      sp_m = sp_p
      sp_p = sp_p+1

      if ( sp_p.gt.nsample ) then     ! a new year
          sp_p = 1
          yrp  = yrp+1

          if (masterproc) then

            if ( yrp.gt.yr_bp ) then    ! need to locate the new year
             indxbp = indxbp+1
             yr_ap = yr_bp
             yr_bp = yr_sulf(indxbp)

             ncid_sulf_ap = ncid_sulf_bp
             sulf_id_ap = sulf_id_bp

             sulfdatafile_bp = trim(sulfdata)//trim(file_sulf(indxbp))
             call getfil(sulfdatafile_bp, locfn_bp)
             call wrap_open(locfn_bp, 0, ncid_sulf_bp)
             call wrap_inq_varid( ncid_sulf_bp, 'sulfmmrbio' , sulf_id_bp )
            endif

          endif

      endif
      
      cdaysulfm = cdaysulfp
      cdaysulfp = yrp*365.0+sample_sulf(sp_p)

      ntmp = nm
      nm   = np
      np   = ntmp

      if (masterproc) then
        strt4(4) = sp_p 
        call wrap_get_vara_realx (ncid_sulf_ap,sulf_id_ap,strt4,cnt4,xsulfi_ap)
        call wrap_get_vara_realx (ncid_sulf_bp,sulf_id_bp,strt4,cnt4,xsulfi_bp)
!
! Normal interpolation between consecutive time slices.
!
        deltat = dble( yr_bp - yr_ap )
        fact1  = dble( yr_bp - yrp  )/deltat
        fact1  = dble( yrp   - yr_ap)/deltat
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
        if (abs(fact1+fact2-1.).gt.1.e-6 .or.  &
            fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or. &
            fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
            write(6,*)'SULFINT--1: Bad fact1 and/or fact2=',fact1,fact2
            call endrun
        end if
!
! interpolation
!      
        do lat = 1,plat
         do k  = 1,pver
          do i = 1,plon
             xsulfi(i,k,lat)= fact1*xsulfi_ap(i,k,lat)+fact2*xsulfi_bp(i,k,lat)
          enddo
         enddo
        enddo
!
! master does the work above
!
      endif  

      call scatter_field_to_chunk(1,pver,1,plon,xsulfi,sulfi(1,1,begchunk,np))

   end if
!
! Determine time interpolation factor.   
!
   deltat = cdaysulfp - cdaysulfm
   fact1 = (cdaysulfp - caldayloc)/deltat
   fact2 = (caldayloc - cdaysulfm)/deltat
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
   if (abs(fact1+fact2-1.).gt.1.e-6 .or.  &
        fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or. &
        fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
      write(6,*)'SULFINT--2: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! Time interpolation.
!
   do lchnk=begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do k=1,pver
         do i=1,ncol
            sulf(i,k,lchnk) = sulfi(i,k,lchnk,nm)*fact1 +  &
                              sulfi(i,k,lchnk,np)*fact2
         end do
      end do
   end do

!----------------------------------------------------------------------------
 endif
!----------------------------------------------------------------------------

   return
end subroutine sulfint_IPCC

!###############################################################################

subroutine getso4bnd_IPCC( lchnk, ncol, sulfdata )

!----------------------------------------------------------------------- 
! 
! Purpose: Return slice of time interpolated sulfate aerosol.
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns

   real(r8), intent(out) :: sulfdata(pcols,pver)  ! biogenic sulfate
!
! Local variables.
!
   integer :: i, k                 ! longitude, level indices
!-----------------------------------------------------------------------
   do k = 1, pver
      do i = 1, ncol
         sulfdata(i,k)  = sulf(i,k,lchnk)
      end do
   end do
   return
end subroutine getso4bnd_IPCC

!###############################################################################

end module so4bnd_IPCC
