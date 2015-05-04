#include <misc.h>
#include <params.h>
module so4bnd
!----------------------------------------------------------------------- 
! 
! Purpose: 
!
! SO4 boundary module.  Deals with interpolating SO4 datasets.
! 
! Author: Brian Eaton
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!
!JR Stuck this "only" business in because Compaq compiler barfed on pcnst, pnats having dual
!JR declarations when radctl.F90 gets compiled.
!
   use pmgrid,    only: plon, plat, masterproc
   use ppgrid,    only: pcols, pver, begchunk, endchunk
   use phys_grid, only: scatter_field_to_chunk, get_ncols_p
   use shr_sys_mod, only : shr_sys_flush
   implicit none
   save
!
! Floating point data
!
   real(r8), private, allocatable, dimension(:,:,:,:) :: &
      sulfbioi    ! input sulfate bio mixing ratios (pcols,pver,begchunk:endchunk,2)
   real(r8), private, allocatable, dimension(:,:,:) :: &
      sulfbio     ! time interpolated sulfate bio mixing ratios (pcols,pver,begchunk:endchunk)
   real(r8), private, allocatable, dimension(:,:,:,:) :: &
      sulfanti    ! input sulfate ant mixing ratios (pcols,pver,begchunk:endchunk,2)
   real(r8), private, allocatable, dimension(:,:,:) :: &
      sulfant     ! time interpolated sulfate ant mixing ratios (pcols,pver,begchunk:endchunk)

   real(r8), private :: sulfscalef                  ! Sulfate scale factor (for 1870->1990 ramp) 
   real(r8), private :: cdaysulfm           ! calendar day for prv. month sulfate values read in
   real(r8), private :: cdaysulfp           ! calendar day for nxt. month sulfate values read in

   integer, private :: date_sulf(1000)              ! Date on sulfate dataset (YYYYMMDD)
   integer, private :: sec_sulf(1000)               ! seconds of date on sulfate dataset (0-86399)
!
! just check that hard-wired size is big enough
!
!
! Integer data
!
   integer, private :: nm,np      ! Array indices for prv., nxt month sulfate data
   integer, private :: np1        ! current forward time index of sulfate dataset
   integer, private :: ncid_sulf  ! sulfate dataset id
   integer, private :: sulfbio_id ! netcdf id for sulfate mmr bio variable
   integer, private :: sulfant_id ! netcdf id for sulfate mmr anth variable
   integer, private :: lonsiz     ! size of longitude dimension on sulfate dataset
   integer, private :: levsiz     ! size of level dimension on sulfate dataset
   integer, private :: latsiz     ! size of latitude dimension on sulfate dataset
   integer, private :: timsiz     ! size of time dimension on sulfate dataset
 
!
! Logical variables
!
   logical, private :: sulfcyc    ! If sulfur cycle code turned on or not

   character*80, private ::  sulfdata ! full pathname for sulfate dataset

contains

subroutine so4bndnl( xsulfdata )
!----------------------------------------------------------------------- 
! 
! Purpose: Set variables from namelist input.
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------

   character*80, intent(in):: xsulfdata ! full pathname for sulfate dataset

!-----------------------------------------------------------------------
   sulfdata = xsulfdata
   if (masterproc) &
      write(6,*)'Time-variant sulfate dataset is: ',trim(sulfdata)

   return
end subroutine so4bndnl

!###############################################################################

subroutine sulfini
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
   character(len=256) locfn      ! local filename        
!
   integer dateid                ! netcdf id for date variable
   integer secid                 ! netcdf id for seconds variable
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
   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
   real(r8) :: calday              ! current calendar day
   real(r8) caldayloc                ! calendar day (includes yr if no cycling)
   real(r8) xsulfbioi(plon,pver,plat,2)  ! input sulfate bio mixing ratios
   real(r8) xsulfanti(plon,pver,plat,2)  ! input sulfate ant mixing ratios
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
  allocate( sulfbioi(pcols,pver,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulfbioi', &
       pcols*pver*(endchunk-begchunk+1)*2 )
  allocate( sulfbio(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulfbio', &
       pcols*pver*(endchunk-begchunk+1) )
  allocate( sulfanti(pcols,pver,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulfanti', &
       pcols*pver*(endchunk-begchunk+1)*2 )
  allocate( sulfant(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'sulfini', 'sulfant', &
       pcols*pver*(endchunk-begchunk+1) )
!
! Currently assume that cycle over 12 months of data
!
  sulfcyc = .true.
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
   if (masterproc) then
!
! Obtain dataset
!
      call getfil(sulfdata, locfn)
      call wrap_open(locfn, 0, ncid_sulf)
!
! Use year information only if not cycling sulfate dataset
!
      calday = get_curr_calday()
      if ( is_perpetual() ) then
         call get_perp_date(yr, mon, day, ncsec)
      else
         call get_curr_date(yr, mon, day, ncsec)
      end if
      ncdate = yr*10000 + mon*100 + day
      if (sulfcyc) then
         caldayloc = calday
      else
         caldayloc = calday + yr*365.
      end if
!
! Obtain dimension id's
!
      call wrap_inq_dimid( ncid_sulf, 'lon', londimid)
      call wrap_inq_dimid( ncid_sulf, 'lat', latdimid)
      call wrap_inq_dimid( ncid_sulf, 'lev', levdimid)
      call wrap_inq_dimid( ncid_sulf, 'time',timeid  )
!
! Obtain size of dimensions.
! Check that horizontal and vertical dimensions are same as model's
!
      call wrap_inq_dimlen( ncid_sulf, londimid, lonsiz   )
      if (lonsiz /= plon) then
         write(6,*)'SULFINI: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf, latdimid, latsiz   )
      if (latsiz /= plat) then
         write(6,*)'SULFINI: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf, levdimid, levsiz   )
      if (levsiz /= pver) then
         write(6,*)'SULFINI: levsiz=',levsiz,' must = ',pver
         call endrun
      end if

      call wrap_inq_dimlen( ncid_sulf, timeid, timsiz   )
!
! Obtain date info id's
!
      call wrap_inq_varid( ncid_sulf, 'date'   , dateid  )
      call wrap_inq_varid( ncid_sulf, 'datesec', secid   )
! 
! Obtain sulfate mixing ratio id
!
      call wrap_inq_varid( ncid_sulf, 'sulfmmrbio' , sulfbio_id )
      call wrap_inq_varid( ncid_sulf, 'sulfmmranth', sulfant_id )

      call wrap_inq_vardimid (ncid_sulf, sulfbio_id, dimids)
      if (dimids(1) /= londimid .and. dimids(2) /= levdimid .and. dimids(3) /= latdimid) then
         write(6,*)'SULFINI: Data must be ordered lon, lev, lat, time'
         call endrun
      end if
!
! just check that hard-wired size is big enough
!
      if (timsiz > 1000) then
         write(6,*)'SO4BND: timsiz=',timsiz,' too small'
         call endrun
      end if
!
! Determine date ids
!
      call wrap_get_var_int (ncid_sulf, dateid, date_sulf)
      call wrap_get_var_int (ncid_sulf, secid, sec_sulf)
!
! If cycling data first do error checks
!
      if (sulfcyc) then
         if (timsiz.lt.12) then 
            write(6,*)'SULFINI: When cycling sulfate dataset must have 12 consecutive ', &
                      'months of data starting with Jan'
            write(6,*)'Current dataset has only ',timsiz,' months'
            call endrun
         end if
         do n = 1,12
            if (mod(date_sulf(n),10000)/100/=n) then
               write(6,*)'SULFINI: When cycling sulfate dataset must have 12 consecutive ', &
                         'months of data starting with Jan'
               write(6,*)'Month ',n,' of dataset says date= ', date_sulf(n)
               call endrun
            end if
         end do
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
! Special code for interpolation between December and January
!
      if (sulfcyc) then
         n = 12
         np1 = 1
         call bnddyi(date_sulf(n  ), sec_sulf(n  ), cdaysulfm)
         call bnddyi(date_sulf(np1), sec_sulf(np1), cdaysulfp)
         if (caldayloc.le.cdaysulfp .or. caldayloc.gt.cdaysulfm) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_sulf,sulfbio_id,strt4,cnt4,xsulfbioi(1,1,1,nm))
            call wrap_get_vara_realx (ncid_sulf,sulfant_id,strt4,cnt4,xsulfanti(1,1,1,nm))
            strt4(4) = np1
            call wrap_get_vara_realx (ncid_sulf,sulfbio_id,strt4,cnt4,xsulfbioi(1,1,1,np))
            call wrap_get_vara_realx (ncid_sulf,sulfant_id,strt4,cnt4,xsulfanti(1,1,1,np))
            goto 10
         end if
      end if
!
! Normal interpolation between consecutive time slices.
!
      do n=1,timsiz-1
         np1 = n + 1
         call bnddyi(date_sulf(n  ), sec_sulf(n  ), cdaysulfm)
         call bnddyi(date_sulf(np1), sec_sulf(np1), cdaysulfp)
         if (.not.sulfcyc) then
            yr = date_sulf(n)/10000
            cdaysulfm = cdaysulfm + yr*365.
            yr = date_sulf(np1)/10000
            cdaysulfp = cdaysulfp + yr*365.
         end if
         if (caldayloc.gt.cdaysulfm .and. caldayloc.le.cdaysulfp) then
            strt4(4) = n
            call wrap_get_vara_realx (ncid_sulf,sulfbio_id,strt4,cnt4,xsulfbioi(1,1,1,nm))
            call wrap_get_vara_realx (ncid_sulf,sulfant_id,strt4,cnt4,xsulfanti(1,1,1,nm))
            strt4(4) = np1
            call wrap_get_vara_realx (ncid_sulf,sulfbio_id,strt4,cnt4,xsulfbioi(1,1,1,np))
            call wrap_get_vara_realx (ncid_sulf,sulfant_id,strt4,cnt4,xsulfanti(1,1,1,np))
            goto 10
         end if
      end do
      write(6,*)'SULFINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
      call endrun

10    continue
      write(6,*)'SULFINI: Read sulfate data for dates ',date_sulf(n), &
           sec_sulf(n),' and ',date_sulf(np1),sec_sulf(np1)

   endif
#if (defined SPMD )
   call mpibcast( timsiz, 1, mpiint, 0, mpicom )
   call mpibcast( date_sulf, 1000, mpiint, 0, mpicom )
   call mpibcast( sec_sulf, 1000, mpiint, 0, mpicom )
   call mpibcast( cdaysulfm, 1, mpir8, 0, mpicom )
   call mpibcast( cdaysulfp, 1, mpir8, 0, mpicom )
   call mpibcast( np1, 1, mpiint, 0, mpicom )
#endif

   call scatter_field_to_chunk(1,pver,2,plon,xsulfbioi,sulfbioi)
   call scatter_field_to_chunk(1,pver,2,plon,xsulfanti,sulfanti)

   return
end subroutine sulfini

!###############################################################################

subroutine sulfint
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate sulfate mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
!-----------------------------------------------------------------------
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
   integer i, k, lchnk        ! column, level, chunk indices
   integer ncol               ! number of columns in current chunk
   integer ntmp               ! temporary
   integer  :: yr, mon, day   ! components of a date
   integer  :: ncdate         ! current date in integer format [yyyymmdd]
   integer  :: ncsec          ! current time of day [seconds]
   real(r8) :: calday         ! current calendar day
   real(r8) fact1, fact2      ! time interpolation factors
   real(r8) caldayloc         ! calendar day (includes yr if no cycling)
   real(r8) deltat            ! time (days) between interpolating sulfate data
   real(r8) xsulfbioi(plon,pver,plat)  ! input sulfate bio mixing ratios
   real(r8) xsulfanti(plon,pver,plat)  ! input sulfate ant mixing ratios

!-----------------------------------------------------------------------
!
! Use year information only if a multiyear dataset
!
   calday = get_curr_calday()
   if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
   else
      call get_curr_date(yr, mon, day, ncsec)
   end if
   ncdate = yr*10000 + mon*100 + day
   if (sulfcyc) then
      caldayloc = calday
   else
      yr = ncdate/10000
      caldayloc = calday + yr*365.
   end if
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
! timeslice for time interpolation.  Messy logic is for sulfcyc = .true. 
! interpolation between December and January (np1.eq.1).  Note that 
! np1 is never 1 when sulfcyc is .false.
!
   if (caldayloc .gt. cdaysulfp .and. .not. (np1.eq.1 .and. caldayloc.gt.cdaysulfm)) then

      if (sulfcyc) then
         np1 = mod(np1,12) + 1
      else
         np1 = np1 + 1
      end if
      if (np1.gt.timsiz) then
         write(6,*)'SULFINT: Attempt to read past end of dataset'
         call endrun
      end if
      cdaysulfm = cdaysulfp
      call bnddyi(date_sulf(np1), sec_sulf(np1), cdaysulfp)
      if (.not.sulfcyc) then
         yr = date_sulf(np1)/10000
         cdaysulfp = cdaysulfp + yr*365.
      end if
      if (np1.eq.1 .or. caldayloc.le.cdaysulfp) then
         ntmp = nm
         nm = np
         np = ntmp
         if (masterproc) then
            strt4(4) = np1
            call wrap_get_vara_realx (ncid_sulf,sulfbio_id,strt4,cnt4,xsulfbioi(1,1,1))
            call wrap_get_vara_realx (ncid_sulf,sulfant_id,strt4,cnt4,xsulfanti(1,1,1))
            write(6,*)'SULFINT: Read sulfate for date (yyyymmdd) ',date_sulf(np1),' sec ', &
                      sec_sulf(np1)
         endif
         call scatter_field_to_chunk(1,pver,1,plon,xsulfbioi,sulfbioi(1,1,begchunk,np))
         call scatter_field_to_chunk(1,pver,1,plon,xsulfanti,sulfanti(1,1,begchunk,np))
      else
         write(6,*)'SULFINT: Input sulfate for date',date_sulf(np1), &
                   ' sec ',sec_sulf(np1), 'does not exceed model date', &
                   ncdate,' sec ',ncsec,' Stopping.'
         call endrun
      end if
   end if
!
! Determine time interpolation factor.  Account for December-January 
! interpolation if cycling sulfate dataset.  Again note that np1 is never 1 
! when sulfcyc is false
!
   if (np1.eq.1) then                    ! Dec-Jan interpolation
      deltat = cdaysulfp + 365. - cdaysulfm
      if (caldayloc.gt.cdaysulfp) then   ! We're in December
         fact1 = (cdaysulfp + 365. - caldayloc)/deltat
         fact2 = (caldayloc - cdaysulfm)/deltat
      else                               ! We're in January
         fact1 = (cdaysulfp - caldayloc)/deltat
         fact2 = (caldayloc + 365. - cdaysulfm)/deltat
      end if
   else
      deltat = cdaysulfp - cdaysulfm
      fact1 = (cdaysulfp - caldayloc)/deltat
      fact2 = (caldayloc - cdaysulfm)/deltat
   end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
   if (abs(fact1+fact2-1.).gt.1.e-6 .or.  &
        fact1.gt.1.000001 .or. fact1.lt.-1.e-6 .or. &
        fact2.gt.1.000001 .or. fact2.lt.-1.e-6) then
      write(6,*)'SULFINT: Bad fact1 and/or fact2=',fact1,fact2
      call endrun
   end if
!
! Time interpolation.
!
   do lchnk=begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do k=1,pver
         do i=1,ncol
            sulfbio(i,k,lchnk) = sulfbioi(i,k,lchnk,nm)*fact1 +  &
                               sulfbioi(i,k,lchnk,np)*fact2
            sulfant(i,k,lchnk) = sulfanti(i,k,lchnk,nm)*fact1 +  &
                               sulfanti(i,k,lchnk,np)*fact2
         end do
      end do
   end do

   return
end subroutine sulfint

!###############################################################################

subroutine getso4bnd( lchnk, ncol, bio, anth )

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

   real(r8), intent(out) :: bio(pcols,pver)  ! biogenic sulfate
   real(r8), intent(out) :: anth(pcols,pver) ! anthropogenic sulfate
!
! Local variables.
!
   integer :: i, k                 ! longitude, level indices
!-----------------------------------------------------------------------
   do k = 1, pver
      do i = 1, ncol
         bio(i,k)  = sulfbio(i,k,lchnk)
         anth(i,k) = sulfant(i,k,lchnk)
      end do
   end do
   return
end subroutine getso4bnd

!###############################################################################

subroutine setso4ramp( x )

!----------------------------------------------------------------------- 
! 
! Purpose: Set so4 ramp value.
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
   real(r8),intent(in) :: x  ! sulfate scale factor computed in ramp subroutine
!-----------------------------------------------------------------------
   sulfscalef = x
   return
end subroutine setso4ramp

!###############################################################################

real*8 function so4ramp()
!----------------------------------------------------------------------- 
! 
! Purpose: Return so4 ramp value.
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   so4ramp = sulfscalef

   return
end function so4ramp

end module so4bnd
