#include <misc.h>
#include <preproc.h>

module histFileMod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varpar       !parameters
  use clm_varmap       !mapping variables
  use shr_const_mod, only: SHR_CONST_CDAY
  use fileutils, only : get_filename
  implicit none

!----------------------------------------------------------------------------
! Creating a netCDF dataset is a two phase process in which (1) dimensions
! and variables are first defined (define mode), but can not read or write
! data, and (2) variables are written (data mode), but can not create new
! dimensions or variables. The general netCDF calling sequence to do this is:
!
!            - nf_create     !create dataset: enter define mode  -
!            |    nf_def_dim !define dimensions                  |
! define mode|    nf_def_var !define variables and set id's      | histcrt
!            |    nf_put_att !assign attributes to variables     |
!            - nf_enddef     !end definitions: leave define mode -
! data mode  |    nf_put_var !provide values for variables       | histwrt
!            | nf_close      !close: save dataset                | histcls
!
! There is one call to nf_def_dim for each dimension
! There is one call to nf_def_var for each variable
! There is one call to nf_put_att for each attribute for each variable
! There is one call to nf_put_var for each variable
!
! Subroutine histcrt creates a netCDF dataset and creates dimensions/variables
! Subroutine histwrt writes data values to a netCDF dataset
! Subroutine histcls closes a netCDF dataset
!
! Every file that references a netCDF function must have the 
! include statement: include 'netcdf.inc'
!
! netCDF datasets -
!    o referenced by a dataset id that is obtained when the
!      dataset is first created or opened
!    o nf_create (fname, nf_clobber, ncid): creates a new netCDF data file
!      with name [fname], overwritting if already exists, and returning
!      a dataset id [ncid] that is used to refer to dataset in other
!      netCDF function calls
!    o nf_open (fname, nf_nowrite, ncid): opens existing netCDF file
!      [fname] in read mode only, returning id [ncid]
!    o nf_enddef (ncid): takes open netCDF dataset, referenced by [ncid],
!      out of define mode
!    o nf_close (ncid): closes an open netCDF dataset, referenced by [ncid]
!
! netCDF dimensions -
!    o has both a name and a length
!    o one dimension in the dataset can have length unlimited (e.g., time
!      dimension can be unlimited to have multiple time slices in dataset)
!
! netCDF variables -
!    o A netCDF variable has a name, type [nf_char, nf_int, nf_float, 
!      nf_double], shape (dimension), and attributes (e.g., long name, units).
!      These are specified when the variable is defined. A variable also
!      has values, which are specified in data mode. 
!    o A netCDF variable is referenced by an integer variable (1,2,3,...),
!      which is in the order in which the variables are defined
!    o Character string is treated as an array of characters
!    o A coordinate variable is a special netCDF variable that has the
!      same name as a dimensions (e.g., lat(lat), lon(lon)). It is used
!      by some appication packages to define the physical coordinate
!      corresponding to that dimension
!    o Variables are dimensioned in Fortran opposite of how listed
!      in netCDF file, e.g.:
!               Fortran             netCDF
!              --------------    --------------
!              x(lon,lat)     -> x(lat,lon)
!              x(lon,lat,lev) -> x(lev,lat,lon)
!              x(lon,lat,tim) -> x(tim,lat,lon)
!
! netCDF functions used:
!    o nf_create       (fname, nf_clobber, ncid)
!    o nf_def_dim      (ncid, dimnam, dimlen, dimid)
!    o nf_def_var      (ncid, varnam, vartyp, ndim, vdim, varid)
!    o nf_put_att_text (ncid, varid, attnam, len, text)
!    o nf_put_att_real (ncid, varid, attnam, vartyp, len, attval)
!    o nf_enddef       (ncid)
!
!    character fname  - netCDF dataset name
!    'nf_clobber'     - overwrite existing file
!    integer   ncid   - returned netCDF dataset id
!    character dimnam - dimension name 
!    character dimid  - associated dimension id 
!    character dimlen - dimension length 
!    character varnam - variable name 
!    integer   varid  - associated netCDF variable id
!              varytp - 'nf_int', 'nf_float', 'nf_double'
!    integer   ndim   - number of dimensions: 0 - scalar. 1 - vector
!    integer   vdim   - vector of ndim dimension id's corresponding to 
!                       the variables dimensions
!    character attnam - character attribute name (e.g., 'units')
!    character text   - attribute text
!    integer   len    - length of attribute text or attribute array
!    real      attval - array of len attribute values
!----------------------------------------------------------------------------
! $Id: histFileMod.F90,v 1.19.6.8.2.1 2002/10/03 20:07:37 erik Exp $
!-----------------------------------------------------------------------

! History file parameters

  real(r8), public, parameter :: spval = 1.e36     !special value for fill value

! History file structures  

  type histentry
     logical           :: active(maxflds)         !true => field is active
     character(len= 8) :: name(maxflds)           !field name
     character(len= 8) :: unit(maxflds)           !field units
     character(len= 8) :: levl(maxflds)           !field levels: single level, multi soil
     character(len= 8) :: type(maxflds)           !field time accumulation type: inst, maxi, mini, aver
     character(len=40) :: desc(maxflds)           !field description
  end type histentry

  type singl_level
     integer           :: num(maxhist)               !number of active single-level fields
     character(len= 8) :: nam(max_slevflds,maxhist)  !single-level field: name
     character(len= 8) :: uni(max_slevflds,maxhist)  !single-level field: units
     character(len= 8) :: typ(max_slevflds,maxhist)  !time accumation type: ninst, nmaxi, nmini, naver
     character(len=40) :: des(max_slevflds,maxhist)  !description of single-level fields
     integer , pointer :: count(:,:,:)               !number accumulations, single-level field
     real(r8), pointer :: value(:,:,:)               !accumulated single-lev field
  end type singl_level

  type multi_level
     integer           :: num(maxhist)               !number of active multi-level fields
     character(len= 8) :: nam(max_mlevflds,maxhist)  !multi-level field : name
     character(len= 8) :: uni(max_mlevflds,maxhist)  !multi-level field : units
     character(len= 8) :: typ(max_mlevflds,maxhist)  !time accumation type: ninst, nmaxi, nmini, naver
     character(len=40) :: des(max_mlevflds,maxhist)  !description of multi-level fields
     integer , pointer :: count(:,:,:,:)             !number accumulations, mutli-level field 
     real(r8), pointer :: value(:,:,:,:)             !accumulated multi-lev field
  end type multi_level

! History file variables

  integer :: nhist                     !actual number of history files
  integer :: ncid(maxhist)             !netCDF id from nf_open or nf_create
  logical :: ncgetid(maxhist)          !true: need to get netCDF variable id's (masterproc only)
  logical :: ehi(maxhist)              !true: current nstep is end of history interval
  integer :: ntim(maxhist)             !current number of time samples for history file
  integer :: nbeghis(maxhist)          !nbeghis=1:current nstep begins history interval

  character(len=80), public :: timcom(maxhist)           !comment: start and end of history interval 
  character(len= 8), public :: fldaux(maxalflds,maxhist) !fields for auxillary history files

! History field level types

  type(singl_level) :: slfld                      !history file
  type(multi_level) :: mlsoifld                   !history file
  character(len= 8) :: nsing = 'sing-lev'         !single-level field
  character(len= 8) :: nsoil = 'mlev_soi'         !multi-level soil field

! History field time accumulation types

  character(len= 8) :: naver = 'average'          !average field over history interval
  character(len= 8) :: nmaxi = 'maximum'          !max field value over history interval
  character(len= 8) :: nmini = 'minimum'          !min field value over history interval
  character(len= 8) :: ninst = 'instant'          !instantaneous field value
  character(len= 8) :: ncnst = 'constnt'          !instantaneous field value

! History file grid variable id's

  integer :: lonvar_id(maxhist)        !id full grid longitude coordinate variable
  integer :: latvar_id(maxhist)        !id full grid latitude  coordinate variable
  integer :: levvar_id(maxhist)        !id soil level coordinate variable
  integer :: timvar_id(maxhist)        !id timecoordinate variable
  integer :: longxy_id(maxhist)        !id 2d longitudes (longxy)
  integer :: latixy_id(maxhist)        !id 2d latitudes (latixy)
  integer :: area_id(maxhist)          !id 2d area (area)
  integer :: landfrac_id(maxhist)      !id 2d land fraction
  integer :: numlon_id(maxhist)        !id number of longitudes at each latitude
  integer :: landmask_id(maxhist)      !id 2d land/ocean mask (landmask)
#if (defined OFFLINE)
  integer :: edgen_id(maxhist)         !id northern edge of grid (lsmedge(1))
  integer :: edgee_id(maxhist)         !id eastern  edge of grid (lsmedge(2))
  integer :: edges_id(maxhist)         !id southern edge of grid (lsmedge(3))
  integer :: edgew_id(maxhist)         !id western  edge of grid (lsmedge(4))
#endif

! History file time variant variable id's

  integer :: slfld_id(max_slevflds,maxhist)    !id single-level fields (slfld%value)
  integer :: mlsoifld_id(max_mlevflds,maxhist) !id multi-level fields (mlsoifld%value)
  integer :: mcdate_id(maxhist)                !id current date, yyyymmdd format (mcdate)
  integer :: mcsec_id(maxhist)                 !id current seconds in day (mcsec)
  integer :: mdcur_id(maxhist)                 !id current day (from base day) (mdcur)
  integer :: mscur_id(maxhist)                 !id current seconds of current day (mdcur)
  integer :: nstep_id(maxhist)                 !id current nstep 
  integer :: timcom_id(maxhist)                !id time comment (timcom)

  SAVE

!=======================================================================
CONTAINS
!=======================================================================

  subroutine histini ()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize variables for history files
!
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

    use clm_varctl
    use spmdMod, only : masterproc

! ------------------------ local variables ------------------------
    integer :: i                   !loop index
! -----------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Initializing variables for history files .....'
       write(6,'(72a1)') ("-",i=1,60)
    endif

! -----------------------------------------------------------------
! Initialize active history fields
! -----------------------------------------------------------------

    call histlst 

! -----------------------------------------------------------------
! Initialize variables for initial or branch runs
! -----------------------------------------------------------------

    if (nsrest==0 .or. nsrest==3) then

! nbeghis = 1 indicates the current time step is start of a history
! interval. This is part of the restart file if continuation run

       nbeghis(:) = 1

! Set accumulation counters to zero: only if current time step
! start of history interval. Otherwise read in from restart file

       slfld%count(:,:,:) = 0
       mlsoifld%count(:,:,:,:) = 0

!  Initialize local file name for history files

       locfnh(:) = ' '

! Set current number of time samples in history file and current 
! history file counter

       ntim(:) = 0

! No need to obtain time dependent netCDF variable id's from history
! file because a new history file will be created

       ncgetid(:) = .false.

    end if

    if (masterproc) then
       write(6,'(72a1)') ("-",i=1,60)
       write(6,*) 'Successfully initialized history files'
       write(6,*)
    endif

    return
  end subroutine histini

!=======================================================================

  subroutine histlst 

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize active field list for history files
!
! Method: 
! This subroutine sets for both primary and auxillary history files:
!    o number      of active single-level and multi-level fields
!    o names       of active single-level and multi-level fields
!    o units       of active single-level and multi-level fields
!    o type        of active single-level and multi-level fields
!    o description of active single-level and multi-level fields
!
! The field types, which are set for each active field, are:
!    o average over history interval
!    o maximum in history interval
!    o minimum in history interval
!    o instantaneous when history file written
!
! Default inactive fields can be made active by setting the [hist_fldadd] 
! variable to the appropriate field name via the namelist input
!
! Field type can be overridden by setting the [hist_chntyp] variable to the
! appropriate field name and new field type via the namelist input
!
! Fields for auxillary files are read from namelist and must be
! a subset of the primary history file fields
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

    use clm_varctl
    use spmdMod, only : masterproc	

! ------------------------ local variables ------------------------
    integer :: i,j,k,n                      !indices
    integer :: nflds = 0                    !number of declared fields (active+inactive)
    integer :: nacti = 0                    !number of active fields
    integer :: ind(maxflds)                 !index to active fields
    type(histentry) histfld                 !primary field names
    type(histentry) tempfld                 !temporary field name
! -----------------------------------------------------------------

! set default fields for primary history files:

! snow properties (will be vertically averaged over the snow profile)

    call histfldini(nflds, 'SNOWDP  ', 'm       ', nsing, naver,  &
         'snow height', .true., histfld)

    call histfldini(nflds, 'SNOWAGE ', 'unitless', nsing, naver,  &
         'snow age', .true., histfld)

! soil time invariant properties (note that lake levels and
! depths do not appear in the fields DZSOI and ZSOI)

    call histfldini(nflds, 'WATSAT  ', 'mm3/mm3 ', nsoil, ncnst,  &
         'saturated soil water content (porosity)', .true., histfld)

    call histfldini(nflds, 'SUCSAT  ', 'mm      ', nsoil, ncnst,  &
         'saturated soil matric potential', .true., histfld) 

    call histfldini(nflds, 'BSW     ', 'unitless', nsoil, ncnst,  &
         'slope of soil water retention curve', .true., histfld) 

    call histfldini(nflds, 'ZSOI    ', 'm       ', nsoil, ncnst,  &
         'soil depth', .true., histfld)

    call histfldini(nflds, 'DZSOI   ', 'm       ', nsoil, ncnst,  &
         'soil thickness', .true., histfld) 

! temperatures

    call histfldini(nflds, 'TSA     ', 'K       ', nsing, naver,  &
         '2 m air temperature', .true., histfld)

    call histfldini(nflds, 'TV      ', 'K       ', nsing, naver,  &
         'vegetation temperature', .true., histfld)

    call histfldini(nflds, 'TG      ', 'K       ', nsing, naver,  &
         'ground temperature', .true., histfld)

    call histfldini(nflds, 'TSOI    ', 'K       ', nsoil, naver,  &  
         'soil temperature', .true., histfld)

    call histfldini(nflds, 'TLAKE   ', 'K       ', nsoil, naver,  &  
         'lake temperature', .true., histfld)

    call histfldini(nflds, 'TSNOW   ', 'K       ', nsing, naver,  &
        'snow temperature', .true., histfld)

! surface radiation                                          

    call histfldini(nflds, 'FSA     ', 'watt/m^2', nsing, naver,  &
         'absorbed solar radiation', .true., histfld)

    call histfldini(nflds, 'FSR     ', 'watt/m^2', nsing, naver,  &
         'reflected solar radiation', .true., histfld)

    call histfldini(nflds, 'NDVI    ', 'unitless', nsing, naver,  &
         'surface ndvi',     .true., histfld)

    call histfldini(nflds, 'FIRA    ', 'watt/m^2', nsing, naver,  &
         'net infrared (longwave) radiation', .true., histfld)

    call histfldini(nflds, 'FIRE    ', 'watt/m^2', nsing, naver,  &
         'emitted infrared (longwave) radiation',.true., histfld)

! surface energy fluxes                                      

    call histfldini(nflds, 'FCTR    ', 'watt/m^2', nsing, naver,  &
         'canopy transpiration',  .true., histfld)

    call histfldini(nflds, 'FCEV    ', 'watt/m^2', nsing, naver,  &
         'canopy evaporation',  .true., histfld)

    call histfldini(nflds, 'FGEV    ', 'watt/m^2', nsing, naver,  &
         'ground evaporation',  .true., histfld)

    call histfldini(nflds, 'FSH     ', 'watt/m^2', nsing, naver,  &
         'sensible heat',  .true., histfld)

    call histfldini(nflds, 'FGR     ', 'watt/m^2', nsing, naver,  &
         'heat flux into soil',  .true., histfld)

    call histfldini(nflds, 'FSM     ', 'watt/m^2', nsing, naver,  &
         'snow melt heat flux',  .true., histfld)

    call histfldini(nflds, 'TAUX    ', 'kg/m/s^2', nsing, naver,  &
         'zonal surface stress',      .true., histfld)

    call histfldini(nflds, 'TAUY    ', 'kg/m/s^2', nsing, naver,  &
         'meridional surface stress',      .true., histfld)

! vegetation phenology                                       

    call histfldini(nflds, 'ELAI    ', 'm^2/m^2 ', nsing, naver,  &
         'exposed one-sided leaf area index',  .true., histfld)

    call histfldini(nflds, 'ESAI    ', 'm^2/m^2 ', nsing, naver,  &
         'exposed one-sided stem area index',  .true., histfld)

! canopy physiology                                          

    call histfldini(nflds, 'RSSUN   ', 's/m     ', nsing, nmini,  &
         'sunlit leaf stomatal resistance',  .true., histfld)

    call histfldini(nflds, 'RSSHA   ', 's/m     ', nsing, nmini,  &
         'shaded leaf stomatal resistance',  .true., histfld)

    call histfldini(nflds, 'BTRAN   ', 'unitless', nsing, naver,  &
         'transpiration beta factor', .true., histfld)       

    call histfldini(nflds, 'FPSN    ', 'umol/m2s', nsing, naver,  &
         'photosynthesis',  .true., histfld)

! hydrology

    call histfldini(nflds, 'H2OSOI  ', 'mm3/mm3 ', nsoil, naver,  &
         'volumetric soil water', .true., histfld)           

    call histfldini(nflds, 'H2OSNO  ', 'mm      ', nsing, naver,  &
         'snow depth (liquid water)', .true., histfld)       

    call histfldini(nflds, 'H2OCAN  ', 'mm      ', nsing, naver,  &
         'intercepted water', .true., histfld)               

    call histfldini(nflds, 'SOILLIQ', 'kg/m2    ', nsoil, naver,  &
        'soil liquid water', .true., histfld)

    call histfldini(nflds, 'SOILICE', 'kg/m2    ', nsoil, naver,  &
        'soil ice', .true., histfld)

    call histfldini(nflds, 'SNOWLIQ', 'kg/m2    ', nsing, naver,  &
        'snow liquid water', .true., histfld)

    call histfldini(nflds, 'SNOWICE', 'kg/m2    ', nsing, naver,  &
        'snow ice', .true., histfld)

    call histfldini(nflds, 'QINFL   ', 'mm/S    ', nsing, naver,  &
         'infiltration', .true., histfld)                    

    call histfldini(nflds, 'QOVER   ', 'mm/S    ', nsing, naver,  &
         'surface runoff', .true., histfld)                  

    call histfldini(nflds, 'QRGWL   ', 'mm/S    ', nsing, naver,  &
         'surf runoff at glaciers, wetlands, lakes', .true., histfld)

    call histfldini(nflds, 'QDRAI   ', 'mm/S    ', nsing, naver,  &
         'sub-surface drainage', .true., histfld)

    call histfldini(nflds, 'QINTR   ', 'mm/S    ', nsing, naver,  &
         'interception',      .true., histfld)        

    call histfldini(nflds, 'QDRIP   ', 'mm/S    ', nsing, naver,  &
         'throughfall',      .true., histfld)         

    call histfldini(nflds, 'QMELT   ', 'mm/S    ', nsing, naver,  &
         'snow melt',      .true., histfld)            

    call histfldini(nflds, 'QSOIL   ', 'mm/S    ', nsing, naver,  &
         'ground evaporation',      .true., histfld)   

    call histfldini(nflds, 'QVEGE   ', 'mm/S    ', nsing, naver,  &
         'canopy evaporation',      .true., histfld)   

    call histfldini(nflds, 'QVEGT   ', 'mm/S    ', nsing, naver,  &
         'canopy transpiration',      .true., histfld) 

#ifdef RTM
    call histfldini(nflds, 'QCHOCNR ', 'm3/s    ', nsing, naver,  &
         'RTM river discharge into ocean', .true., histfld)

    call histfldini(nflds, 'QCHANR  ', 'M3/S    ', nsing, naver,  &
         'RTM river flow (maximum subgrid flow)', .true., histfld)     
#endif

! water and energy balance checks

    call histfldini(nflds, 'ERRSOI  ', 'watt/m^2', nsing, naver,  &
         'soil/lake energy conservation error', .true., histfld)

    call histfldini(nflds, 'ERRSEB  ', 'watt/m^2', nsing, naver,  &
         'surface energy conservation error',  .true., histfld)

    call histfldini(nflds, 'ERRSEBMX', 'watt/m^2', nsing, nmaxi,  &
         'maximum of surface energy conservation error',  .true., histfld)

    call histfldini(nflds, 'ERRSOL  ', 'watt/m^2', nsing, naver,  &
         'solar radiation conservation error',  .true., histfld)

    call histfldini(nflds, 'ERRH2O  ', 'mm      ', nsing, naver,  &
         'total water conservation error',  .true., histfld)

! atmospheric forcing                                

    call histfldini(nflds, 'RAIN    ', 'mm/S    ', nsing, naver,  &
         'rain',      .true., histfld)

    call histfldini(nflds, 'SNOW    ', 'mm/S    ', nsing, naver,  &
         'snow',      .true., histfld)

    call histfldini(nflds, 'TBOT    ', 'K       ', nsing, naver,  &
         'atmospheric air temperature',  .true., histfld)

    call histfldini(nflds, 'WIND    ', 'm/s     ', nsing, naver,  &
         'atmospheric wind velocity magnitude',  .true., histfld)

    call histfldini(nflds, 'THBOT   ', 'K       ', nsing, naver,  &
         'atmospheric air potential temperature',  .true., histfld)

    call histfldini(nflds, 'QBOT    ', 'kg/kg   ', nsing, naver,  &
         'atmospheric specific humidity',  .true., histfld)

    call histfldini(nflds, 'ZBOT    ', 'M       ', nsing, naver,  &
         'atmospheric reference height',  .true., histfld)

    call histfldini(nflds, 'FLDS    ', 'watt/m^2', nsing, naver,  &
         'incident longwave radiation',  .true., histfld)

    call histfldini(nflds, 'FSDS    ', 'watt/m^2', nsing, naver,  &
         'incident solar radiation',  .true., histfld)

#if (defined OFFLINE)
! accumulation variables                             

    call histfldini(nflds, 'TDA     ', 'K       ', nsing, ninst,  &
         'daily average 2-m temperature', .false., histfld)

    call histfldini(nflds, 'T15     ', 'K       ', nsing, ninst,  &
         '15-day running mean of 2-m temperature', .false., histfld)

    call histfldini(nflds, 'AGDD0   ', 'K       ', nsing, ninst,  &
         'growing degree-days base 0C', .false., histfld)

    call histfldini(nflds, 'AGDD5   ', 'K       ', nsing, ninst,  &
         'growing degree-days base 5C', .false., histfld)

#endif

! use namelist variable [hist_chntyp] to override default [histfld%type] variable 

    do i = 1, nflds
       do j = 1, maxalflds
          if (histfld%name(i) == hist_chntyp(1,j)) then
             if (hist_chntyp(2,j)==naver .or. hist_chntyp(2,j)==nmaxi .or. &
                 hist_chntyp(2,j)==nmini .or. hist_chntyp(2,j)==ninst) then
                histfld%type(i) = hist_chntyp(2,j)
             else
                write(6,*) 'HISTLST error: attempting to change', &
                     ' field type for field = ',histfld%name(i) 
                write(6,*) 'to inaccurate type = ',hist_chntyp(2,j), &
                     '. valid types are: ',naver,nmaxi,nmini,ninst
                call endrun
             end if
          end if
       end do
    end do

! use namelist variable [hist_fldadd] to add to default active fields

    do i = 1, nflds
       do j = 1, maxalflds
          if (histfld%name(i) == hist_fldadd(j)) then
             histfld%active(i) = .true.
          endif
       end do
    end do

! number of active fields (nacti <= nflds) and pointer (1 to nflds)

    nacti = 0
    do i = 1, nflds
       if (histfld%active(i)) then
          nacti = nacti + 1
          ind(nacti) = i
       end if
    end do

! re-order fields from 1 -> nflds to 1 -> nacti

    do i = 1, nacti
       tempfld%name(i)   = histfld%name(ind(i))
       tempfld%unit(i)   = histfld%unit(ind(i))
       tempfld%levl(i)   = histfld%levl(ind(i))
       tempfld%type(i)   = histfld%type(ind(i))
       tempfld%desc(i)   = histfld%desc(ind(i))
    end do
    do i = 1, nacti
       histfld%name(i)   = tempfld%name(i) 
       histfld%unit(i)   = tempfld%unit(i) 
       histfld%levl(i)   = tempfld%levl(i) 
       histfld%type(i)   = tempfld%type(i) 
       histfld%desc(i)   = tempfld%desc(i)   
    end do

! separate single-level and multi-level fields for primary history files

    slfld%num(1) = 0
    mlsoifld%num(1) = 0
    do i = 1, nacti
       if (histfld%levl(i) == nsing) then
          slfld%num(1) = slfld%num(1) + 1
          n = slfld%num(1)
          slfld%nam(n,1) = histfld%name(i)
          slfld%uni(n,1) = histfld%unit(i)
          slfld%typ(n,1) = histfld%type(i)
          slfld%des(n,1) = histfld%desc(i)
       else if (histfld%levl(i) == nsoil) then
          mlsoifld%num(1) = mlsoifld%num(1) + 1
          n = mlsoifld%num(1)
          mlsoifld%nam(n,1) = histfld%name(i)
          mlsoifld%uni(n,1) = histfld%unit(i)
          mlsoifld%typ(n,1) = histfld%type(i)
          mlsoifld%des(n,1) = histfld%desc(i)
       end if
    end do
    if (slfld%num(1) > max_slevflds) then
       write(6,*) 'HISTLST error: number single-level fields', &
            ' for primary files > parameter max_slevflds'
       call endrun
    endif
    if (mlsoifld%num(1) > max_mlevflds) then
       write(6,*) 'HISTLST error: number multi-level soil fields', &
            ' for primary files > parameter max_mlevflds'
       call endrun
    endif

! now make list of single-level and multi-level fields for auxillary files

    slfld%num(2:nhist) = 0
    mlsoifld%num(2:nhist) = 0
    do k = 2, nhist
       do j = 1, nacti
          do i = 1, slfld%num(1)
             if (fldaux(j,k-1) == slfld%nam(i,1)) then
                slfld%num(k) = slfld%num(k) + 1
                n = slfld%num(k)
                slfld%nam(n,k) = slfld%nam(i,1)
                slfld%uni(n,k) = slfld%uni(i,1)
                slfld%typ(n,k) = slfld%typ(i,1)
                slfld%des(n,k) = slfld%des(i,1)
             end if
          end do
          do i = 1, mlsoifld%num(1)
             if (fldaux(j,k-1) == mlsoifld%nam(i,1)) then
                mlsoifld%num(k) = mlsoifld%num(k) + 1
                n = mlsoifld%num(k)
                mlsoifld%nam(n,k) = mlsoifld%nam(i,1)
                mlsoifld%uni(n,k) = mlsoifld%uni(i,1)
                mlsoifld%typ(n,k) = mlsoifld%typ(i,1)
                mlsoifld%des(n,k) = mlsoifld%des(i,1)
             end if
          end do
       end do
       if (slfld%num(k) > max_slevflds) then
          write(6,*) 'HISTLST error: number single-level fields', &
               ' for auxillary files > parameter max_slevflds '
          call endrun
       endif
       if (mlsoifld%num(k) > max_mlevflds) then
          write(6,*) 'HISTLST error: number multi-level fields', &
               ' for auxillary files > parameter max_mlevflds '
          call endrun
       endif
    end do

    j = 0
    do i = 1, nhist
       if ((slfld%num(i)+mlsoifld%num(i)) > 0) j = j + 1
    end do
    if (j /= nhist) then
       write(6,*) 'HISTLST error: number of history files = ', nhist
       write(6,*) 'but number of files based on active fields = ',j
       call endrun
    end if

! echo active fields 

    if (masterproc) then
       do j = 1, nhist
          if (slfld%num(j) > 0) then
             write(6,*)
             write(6,*) 'History file ',j,': Active single-level fields'
             write(6,1002) 
             write(6,'(72a1)') ("_",i=1,71)
             do i = 1, slfld%num(j)
                write(6,1003)i,slfld%nam(i,j),&
                     slfld%uni(i,j),slfld%typ(i,j),slfld%des(i,j)
             end do
             write(6,'(72a1)') ("_",i=1,71)
             write(6,*)
          end if
          if (mlsoifld%num(j) > 0) then
             write(6,*) 'History file ',j,': Active multi-level soil fields'
             write(6,1002) 
             write(6,'(72a1)') ("_",i=1,71)
             do i = 1, mlsoifld%num(j)
                write(6,1003)i,mlsoifld%nam(i,j),&
                     mlsoifld%uni(i,j),mlsoifld%typ(i,j),mlsoifld%des(i,j)
             end do
             write(6,'(72a1)') ("_",i=1,71)
             write(6,*)
          end if
       end do
    endif
1002 format(' No',' Name    ',' Units   ',' Type    ',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),(1x,a40))

    return
  end subroutine histlst

!=======================================================================

  subroutine histfldini (nflds, name, unit, levl, type, &
                         desc, active, histfld)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set up history file field (active or inactive)
!
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

! ------------------------ arguments ---------------------------------
    integer, intent(inout)       :: nflds   !number of fields       
    character(len=*), intent(in) :: name    !field name
    character(len=*), intent(in) :: unit    !field units 
    character(len=*), intent(in) :: levl    !field level type 
    character(len=*), intent(in) :: type    !field time averaging type
    character(len=*), intent(in) :: desc    !field description
    logical         , intent(in) :: active  !true=> field is active
    type(histentry) , intent(out):: histfld !history field entry
! --------------------------------------------------------------------

    nflds = nflds + 1
    histfld%name(nflds)   = name
    histfld%unit(nflds)   = unit
    histfld%levl(nflds)   = levl
    histfld%type(nflds)   = type
    histfld%desc(nflds)   = desc
    histfld%active(nflds) = active

    return
  end subroutine histfldini

!=======================================================================

  subroutine histcrt (nf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! create netCDF history file
!
! Method: 
! This subroutine opens a new netCDF data file. Global attributes
! and variables are defined in define mode. Upon exiting this
! routine, define mode is exited and the file is ready to write.
!
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use clm_varsur  , only : fullgrid, offline_rdgrid 	
    use time_manager, only : get_ref_date 
    use clm_varctl        
    include 'netcdf.inc'

! ------------------------ arguments ---------------------------------
    integer, intent(in) :: nf  !current history file: 1=primary. 2,3,4=auxillary
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
    integer i           !field do loop index
    integer status      !netCDF error status
    integer tim_id      !netCDF id for time dimension
    integer lon_id      !netCDF id for longitude dimension
    integer lat_id      !netCDF id for latitude dimension
    integer levsoi_id   !netCDF id for soil layer dimension
    integer patch_id    !netCDF id for total number subgrid patches
    integer strlen_id   !netCDF id for character string variables
    integer dim1_id(1)  !netCDF dimension id for 1-d variables
    integer dim2_id(2)  !netCDF dimension id for 2-d variables
    integer dim3_id(3)  !netCDF dimension id for 3-d variables
    integer dim4_id(4)  !netCDF dimension id for 4-d variables
    integer omode       !netCDF dummy variable
    character(len=256) name     !name of attribute
    character(len=256) unit     !units of attribute
    character(len=256) mode     !field mode (aver, inst, max, min, etc)
    character(len=256) str      !global attribute string 
    character(len=  8) curdate  !current date
    character(len=  8) curtime  !current time 
    character(len= 10) basedate !base date (yyyymmdd)
    character(len=  8) basesec  !base seconds
    integer yr,mon,day,nbsec    !year,month,day,seconds components of a date
    integer hours,minutes,secs  !hours,minutes,seconds of hh:mm:ss
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Create new netCDF file. File will be in define mode
! --------------------------------------------------------------------

    call wrap_create (trim(locfnh(nf)), nf_clobber, ncid(nf))

! set fill mode to "no fill" to optimize performance

    status = nf_set_fill (ncid(nf), nf_nofill, omode)
    if (status /= nf_noerr) then
       write(6,*) ' netCDF error = ',nf_strerror(status)
       call endrun
    end if

! --------------------------------------------------------------------
! Create global attributes. Attributes are used to store information
! about the data set. Global attributes are information about the
! data set as a whole, as opposed to a single variable
! --------------------------------------------------------------------

    str = 'CF1.0'
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'conventions', trim(str))
    
    call datetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call wrap_put_att_text(ncid(nf), NF_GLOBAL,'history', trim(str))

!!  call getenv ('LOGNAME', str)      !!(2003.10.22)
    call getenv ('HOME'   , str)
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'logname',trim(str))
    
    call getenv ('HOST', str)
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'host', trim(str))
    
    str = 'Community Land Model: CLM2'
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'source', trim(str))
    
    str = '$Name: GAMIL1.0 $'
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'version', trim(str))
    
    str = '$Id: histFileMod.F90,v 1.19.6.8.2.1 2002/10/03 20:07:37 erik Exp $'
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'revision_id', trim(str))
    
    str = ctitle 
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'case_title', trim(str))

    str = caseid
    call wrap_put_att_text (ncid(nf), NF_GLOBAL, 'case_id', trim(str))
    
    if (fsurdat == ' ') then
       str = 'created at run time'
    else
       str = get_filename(fsurdat)
    endif
    call wrap_put_att_text(ncid(nf), NF_GLOBAL, 'Surface_dataset', trim(str))

    if (finidat == ' ') then
       str = 'arbitrary initialization'
    else
       str = get_filename(finidat)
    endif
    call wrap_put_att_text(ncid(nf), NF_GLOBAL, 'Initial_conditions_dataset', trim(str))

    str = get_filename(fpftcon)
    call wrap_put_att_text(ncid(nf), NF_GLOBAL, 'PFT_physiological_constants_dataset', trim(str))

    if (frivinp_rtm /= ' ') then
       str = get_filename(frivinp_rtm)
       call wrap_put_att_text(ncid(nf), NF_GLOBAL, 'RTM_input_datset', trim(str))
    endif

! --------------------------------------------------------------------
! Define dimensions. Array dimensions are referenced by an
! associated dimenision id: e.g., lon_id -> lon.
! o Time is an unlimited dimension.
! o Character string is treated as an array of characters. 
! --------------------------------------------------------------------

    call wrap_def_dim (ncid(nf), 'lon'   , lsmlon , lon_id)
    call wrap_def_dim (ncid(nf), 'lat'   , lsmlat , lat_id)
    call wrap_def_dim (ncid(nf), 'levsoi', nlevsoi, levsoi_id)
    if (.not. hist_dov2xy(nf)) then
       call wrap_def_dim (ncid(nf), 'patch', numpatch, patch_id)
    end if
    call wrap_def_dim (ncid(nf), 'time'  , nf_unlimited, tim_id)
    call wrap_def_dim (ncid(nf), 'string_length', 80, strlen_id)

! --------------------------------------------------------------------
! Define time-independent grid variables 
! --------------------------------------------------------------------

    mode = 'time-invariant'

! coordinate variables (including time)

    if (fullgrid) then
       dim1_id(1) = lon_id
       name = 'coordinate longitude'
       unit = 'degrees_east'
       call wrap_def_var (ncid(nf), 'lon' , ncprec, 1, dim1_id, lonvar_id(nf))
       call wrap_put_att_text (ncid(nf), lonvar_id(nf), 'long_name',name)
       call wrap_put_att_text (ncid(nf), lonvar_id(nf), 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), lonvar_id(nf), 'mode'     ,mode)
       
       dim1_id(1) = lat_id
       name = 'coordinate latitude'
       unit = 'degrees_north'
       call wrap_def_var (ncid(nf), 'lat' , ncprec, 1, dim1_id, latvar_id(nf))
       call wrap_put_att_text (ncid(nf), latvar_id(nf), 'long_name',name)
       call wrap_put_att_text (ncid(nf), latvar_id(nf), 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), latvar_id(nf), 'mode'     ,mode)
    endif

    dim1_id(1) = levsoi_id
    name = 'coordinate soil levels'
    unit = 'm'
    call wrap_def_var (ncid(nf), 'levsoi' , ncprec, 1, dim1_id, levvar_id(nf))
    call wrap_put_att_text (ncid(nf), levvar_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), levvar_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), levvar_id(nf), 'mode'     ,mode)

    dim1_id(1) = tim_id
    name = 'time'
    call get_ref_date(yr, mon, day, nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,80) yr,mon,day
80  format(i4.4,'-',i2.2,'-',i2.2)
    write(basesec ,90) hours, minutes, secs
90  format(i2.2,':',i2.2,':',i2.2)
    unit = 'days since ' // basedate // " " // basesec
    call wrap_def_var (ncid(nf), 'time', ncprec, 1, dim1_id, timvar_id(nf))
    call wrap_put_att_text (ncid(nf), timvar_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), timvar_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), timvar_id(nf), 'calendar' ,'noleap')

#if (defined OFFLINE)

! surface grid edges

    if (.not. offline_rdgrid) then
       name = 'northern edge of surface grid'
       unit = 'degrees_north'
       call wrap_def_var (ncid(nf) , 'edgen', ncprec, 0, 0, edgen_id(nf))
       call wrap_put_att_text (ncid(nf), edgen_id(nf), 'long_name',name)
       call wrap_put_att_text (ncid(nf), edgen_id(nf), 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), edgen_id(nf), 'mode'     ,mode)
       
       name = 'eastern edge of surface grid'
       unit = 'degrees_east'
       call wrap_def_var (ncid(nf), 'edgee', ncprec,0, 0, edgee_id(nf))
       call wrap_put_att_text (ncid(nf), edgee_id(nf), 'long_name',name)
       call wrap_put_att_text (ncid(nf), edgee_id(nf), 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), edgee_id(nf), 'mode'     ,mode)
       
       name = 'southern edge of surface grid'
       unit = 'degrees_north'
       call wrap_def_var (ncid(nf), 'edges', ncprec, 0, 0, edges_id(nf))
       call wrap_put_att_text (ncid(nf), edges_id(nf), 'long_name',name)
       call wrap_put_att_text (ncid(nf), edges_id(nf), 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), edges_id(nf), 'mode'     ,mode)
       
       name = 'western edge of surface grid'
       unit = 'degrees_east'
       call wrap_def_var (ncid(nf), 'edgew', ncprec, 0, 0, edgew_id(nf))
       call wrap_put_att_text (ncid(nf), edgew_id(nf) , 'long_name',name)
       call wrap_put_att_text (ncid(nf), edgew_id(nf) , 'units'    ,unit)
       call wrap_put_att_text (ncid(nf), edgew_id(nf) , 'mode'     ,mode)
    endif
       
#endif

! longitude, latitude, surface type: real (lsmlon x lsmlat)

    dim2_id(1) = lon_id
    dim2_id(2) = lat_id

    if (fullgrid) then
       name = 'longitude'
       unit = 'degrees_east'
       call wrap_def_var (ncid(nf), 'longxy' , ncprec, 2, dim2_id, longxy_id(nf))
    else
       name = 'rlongitude'
       unit = 'degrees_east'
       call wrap_def_var (ncid(nf), 'rlongxy', ncprec, 2, dim2_id, longxy_id(nf))
    endif
    call wrap_put_att_text (ncid(nf), longxy_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), longxy_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), longxy_id(nf), 'mode'     ,mode)

    name = 'latitude'
    unit = 'degrees_north'
    call wrap_def_var (ncid(nf), 'latixy', ncprec, 2, dim2_id, latixy_id(nf))
    call wrap_put_att_text (ncid(nf), latixy_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), latixy_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), latixy_id(nf), 'mode'     ,mode)

    name = 'grid cell areas'
    unit = 'km^2'
    call wrap_def_var (ncid(nf), 'area', ncprec, 2, dim2_id, area_id(nf))
    call wrap_put_att_text (ncid(nf), area_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), area_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), area_id(nf), 'mode'     ,mode)

    name = 'land fraction'
    call wrap_def_var (ncid(nf), 'landfrac', ncprec, 2, dim2_id, landfrac_id(nf))
    call wrap_put_att_text (ncid(nf), landfrac_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), landfrac_id(nf), 'mode'     ,mode)

! number of longitudes per latitude (reduced grid only)

    dim1_id(1) = lat_id
    name = 'number of longitudes at each latitude'
    call wrap_def_var (ncid(nf), 'numlon', nf_int, 1, dim1_id, numlon_id(nf))
    call wrap_put_att_text (ncid(nf), numlon_id(nf), 'long_name', name)

! Surface type

    name = 'land/ocean mask (0.=ocean and 1.=land)'
    call wrap_def_var (ncid(nf), 'landmask', nf_int, 2, dim2_id, landmask_id(nf))
    call wrap_put_att_text (ncid(nf),landmask_id(nf),'long_name',name)
    call wrap_put_att_text (ncid(nf),landmask_id(nf),'mode'     ,mode)

! --------------------------------------------------------------------
! Define time-dependent variables: time information
! --------------------------------------------------------------------

    mode = trim(ninst)

! current date, day and time step

    dim1_id(1) = tim_id

    name = 'current date (YYYYMMDD)'
    call wrap_def_var (ncid(nf) , 'mcdate', nf_int, 1, dim1_id  , mcdate_id(nf))
    call wrap_put_att_text (ncid(nf), mcdate_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), mcdate_id(nf), 'mode'     ,mode)

    name = 'current seconds of current date'
    unit = 's'
    call wrap_def_var (ncid(nf) , 'mcsec' , nf_int, 1, dim1_id , mcsec_id(nf))
    call wrap_put_att_text (ncid(nf), mcsec_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), mcsec_id(nf), 'units'    ,unit)
    call wrap_put_att_text (ncid(nf), mcsec_id(nf), 'mode'     ,mode)

    name = 'current day (from base day)'
    call wrap_def_var (ncid(nf) , 'mdcur' , nf_int, 1, dim1_id , mdcur_id(nf))
    call wrap_put_att_text (ncid(nf), mdcur_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), mdcur_id(nf), 'mode'     ,mode)

    name = 'current seconds of current day'
    call wrap_def_var (ncid(nf) , 'mscur' , nf_int, 1, dim1_id , mscur_id(nf))
    call wrap_put_att_text (ncid(nf), mscur_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), mscur_id(nf), 'mode'     ,mode)

    name = 'time step'
    call wrap_def_var (ncid(nf) , 'nstep' , nf_int, 1, dim1_id , nstep_id(nf))
    call wrap_put_att_text (ncid(nf), nstep_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), nstep_id(nf), 'mode'     ,mode)

! character time comment: character (80 x time)

    dim2_id(1) = strlen_id
    dim2_id(2) = tim_id

    name = 'history interval for time slice'
    call wrap_def_var (ncid(nf) , 'time_comment', nf_char, 2, dim2_id, timcom_id(nf))
    call wrap_put_att_text (ncid(nf), timcom_id(nf), 'long_name',name)
    call wrap_put_att_text (ncid(nf), timcom_id(nf), 'mode'     ,mode)

! --------------------------------------------------------------------
! Define time-dependent variables: active history file fields.
! Array dimensions depend on whether it is
!
! single-level
! o 1-d vector   (hist_dov2xy = false): numpatch x time
! o grid average (hist_dov2xy = true ): lsmlon x lsmlat x time
! 
! multi-level soil (static levels)
! o 1-d vector   (hist_dov2xy = false): numpatch x nlevsoi x time
! o grid average (hist_dov2xy = true ): lsmlon x lsmlat x nlevsoi x time
! --------------------------------------------------------------------

! single level fields

    do i = 1, slfld%num(nf)
       if (hist_dov2xy(nf)) then
          dim3_id(1) = lon_id         
          dim3_id(2) = lat_id         
          dim3_id(3) = tim_id         
          call wrap_def_var (ncid(nf), slfld%nam(i,nf),  ncprec, 3, dim3_id, slfld_id(i,nf))
       else                                                                
          dim2_id(1) = patch_id                                            
          dim2_id(2) = tim_id                                              
          call wrap_def_var (ncid(nf), slfld%nam(i,nf),  ncprec, 2, dim2_id, slfld_id(i,nf))
       end if
    end do

! multi-level soil fields 

    do i = 1, mlsoifld%num(nf)
       if (hist_dov2xy(nf)) then
          dim4_id(1) = lon_id      
          dim4_id(2) = lat_id   
          dim4_id(3) = levsoi_id   
          dim4_id(4) = tim_id      
          call wrap_def_var (ncid(nf), mlsoifld%nam(i,nf), ncprec, 4, dim4_id, mlsoifld_id(i,nf))
       else                                                               
          dim3_id(1) = patch_id                                           
          dim3_id(2) = levsoi_id                                             
          dim3_id(3) = tim_id                                             
          call wrap_def_var (ncid(nf), mlsoifld%nam(i,nf), ncprec, 3, dim3_id, mlsoifld_id(i,nf))
       endif
    end do

! define attributes for each field: long name, units, 
! mode (inst, aver, etc), and fill value (spval)

    do i = 1, slfld%num(nf)
       call wrap_put_att_text (ncid(nf), slfld_id(i,nf), 'long_name' , slfld%des(i,nf))
       call wrap_put_att_text (ncid(nf), slfld_id(i,nf), 'units'     , slfld%uni(i,nf))
       call wrap_put_att_text (ncid(nf), slfld_id(i,nf), 'mode'      , slfld%typ(i,nf))
       call wrap_put_att_realx(ncid(nf), slfld_id(i,nf), '_FillValue', ncprec,1 ,spval)
       call wrap_put_att_realx(ncid(nf), slfld_id(i,nf), 'missing_value', ncprec,1 ,spval)
    end do

    do i = 1, mlsoifld%num(nf)
       call wrap_put_att_text (ncid(nf), mlsoifld_id(i,nf), 'long_name' , mlsoifld%des(i,nf))
       call wrap_put_att_text (ncid(nf), mlsoifld_id(i,nf), 'units'     , mlsoifld%uni(i,nf))
       call wrap_put_att_text (ncid(nf), mlsoifld_id(i,nf), 'mode'      , mlsoifld%typ(i,nf))
       call wrap_put_att_realx(ncid(nf), mlsoifld_id(i,nf), '_FillValue', ncprec,1 ,spval)
       call wrap_put_att_realx(ncid(nf), mlsoifld_id(i,nf), 'missing_value', ncprec,1 ,spval)
    end do

! --------------------------------------------------------------------
! Finish creating netCDF file (end define mode)
! --------------------------------------------------------------------

    status = nf_enddef(ncid(nf))

    return
  end subroutine histcrt

!=======================================================================

  subroutine histwrt (nf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! write to netCDF history file
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8	
    use clm_varder
    use clm_varsur        !surface data   
    use clm_varctl        !run control variables 
#if (defined SPMD)
    use spmdMod     , only : masterproc, npes, compute_mpigs_patch
    use mpishorthand, only : mpir8, mpilog, mpiint, mpicom 
#else
    use spmdMod     , only : masterproc
#endif
    use time_manager, only : get_nstep, get_curr_date, get_curr_time
    use shr_sys_mod, only : shr_sys_flush
    implicit none

! ------------------------ includes ----------------------------------
    include 'netcdf.inc'
! --------------------------------------------------------------------

! ------------------------ arguments ---------------------------------
    integer, intent(in) :: nf   !current history file: 
                                !1 = primary  2,3 = auxillary
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
    integer  :: i,j,k,l,m,n                    !do loop indices
    integer  :: beg1d(1)                       !netCDF 1-d start index 
    integer  :: len1d(1)                       !netCDF 1-d count index 
    integer  :: beg2d(2)                       !netCDF 2-d start index 
    integer  :: len2d(2)                       !netCDF 2-d count index 
    integer  :: beg3d(3)                       !netCDF 3-d start index 
    integer  :: len3d(3)                       !netCDF 3-d count index 
    integer  :: beg4d(4)                       !netCDF 4-d start index 
    integer  :: len4d(4)                       !netCDF 4-d count index 
    real(r8) :: slfxy(lsmlon,lsmlat)           !grid-average single-level field
    real(r8) :: mlsoifxy(lsmlon,lsmlat,nlevsoi)!grid-average multi-level soil field
    real(r8) :: lonvar(lsmlon)                 !only used for full grid 
    real(r8) :: latvar(lsmlat)                 !only used for full grid
    real(r8) :: time                           !current time
    integer  :: mdcur, mscur                   !outputs from get_curr_time
    integer  :: yr,mon,day,mcsec               !outputs from get_curr_date
    integer  :: mcdate                         !current date 
    integer  :: nstep                          !time step
#if (defined SPMD)
    integer :: numrecvv(0:npes-1)              !vector of items to be received  
    integer :: displsv(0:npes-1)               !displacement vector
    integer :: numsend                         !number of items to be sent
    integer :: ier                             !MPI error status
    real(r8), pointer :: buf1d(:)              !temporary for MPI gatherv
    real(r8), pointer :: buf2d(:,:)            !temporary for MPI gatherv
#endif
    real(r8), pointer :: gather1d(:)           !temporary for MPI gatherv
    real(r8), pointer :: output2d(:,:)         !temporary for MPI gatherv
    real(r8), pointer :: gather2d(:,:)         !temporary for MPI gatherv
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! Write out time-invariant variables. Do once, at first write to file.
! --------------------------------------------------------------------

    if (ntim(nf) == 1) then

#if (defined OFFLINE)
       if (masterproc) then
          if (.not. offline_rdgrid) then
             call wrap_put_var_realx (ncid(nf), edgen_id(nf), lsmedge(1))
             call wrap_put_var_realx (ncid(nf), edgee_id(nf), lsmedge(2))
             call wrap_put_var_realx (ncid(nf), edges_id(nf), lsmedge(3))
             call wrap_put_var_realx (ncid(nf), edgew_id(nf), lsmedge(4))
          endif
       endif
#endif

! Surface grid (coordinate variables, latitude, longitude, surface type). 

       if (masterproc) then
          if (fullgrid) then
             lonvar(1:lsmlon) = longxy(1:lsmlon,1)
             call wrap_put_var_realx (ncid(nf), lonvar_id(nf), lonvar)
             latvar(1:lsmlat) = latixy(1,1:lsmlat)
             call wrap_put_var_realx (ncid(nf), latvar_id(nf), latvar)
          endif
          call wrap_put_var_realx (ncid(nf), levvar_id(nf)  , zsoi)
          call wrap_put_var_realx (ncid(nf), longxy_id(nf)  , longxy)
          call wrap_put_var_realx (ncid(nf), latixy_id(nf)  , latixy)
          call wrap_put_var_realx (ncid(nf), area_id(nf)    , area) 
          call wrap_put_var_realx (ncid(nf), landfrac_id(nf), landfrac) 
          call wrap_put_var_int   (ncid(nf), landmask_id(nf), landmask)
          call wrap_put_var_int   (ncid(nf), numlon_id(nf)  , numlon)
       endif

    end if !end of write of time constant variables

! --------------------------------------------------------------------
! Get variable id's for time-varying variables if restart and
! current history file is not full. Needs to be done so that 
! non-full history files can be filled before a new file is created
! --------------------------------------------------------------------

    if (masterproc) then
       if (ncgetid(nf)) then
          call wrap_inq_varid (ncid(nf), 'mcdate', mcdate_id(nf))
          call wrap_inq_varid (ncid(nf), 'mcsec' , mcsec_id(nf))
          call wrap_inq_varid (ncid(nf), 'mdcur' , mdcur_id(nf))
          call wrap_inq_varid (ncid(nf), 'mscur' , mscur_id(nf))
          call wrap_inq_varid (ncid(nf), 'nstep' , nstep_id(nf))
          call wrap_inq_varid (ncid(nf), 'time'  , timvar_id(nf))
          call wrap_inq_varid (ncid(nf), 'time_comment', timcom_id(nf))
          do i = 1, slfld%num(nf)
             call wrap_inq_varid (ncid(nf), slfld%nam(i,nf), slfld_id(i,nf))
          end do
          do i = 1, mlsoifld%num(nf)
             call wrap_inq_varid (ncid(nf), mlsoifld%nam(i,nf), mlsoifld_id(i,nf))
          end do
          ncgetid(nf) = .false.
       end if
    endif

! --------------------------------------------------------------------
! Write time-varying variables
! --------------------------------------------------------------------

! current date, seconds, day and nstep

    if (masterproc) then
       beg1d(1) = ntim(nf) ; len1d(1) = 1
       call get_curr_date(yr, mon, day, mcsec)
       mcdate = yr*10000 + mon*100 + day
       call get_curr_time(mdcur,mscur)  
       time = mdcur + mscur/SHR_CONST_CDAY
       nstep = get_nstep()
       call wrap_put_vara_int (ncid(nf), mcdate_id(nf), beg1d, len1d, mcdate)
       call wrap_put_vara_int (ncid(nf), mcsec_id(nf) , beg1d, len1d, mcsec)
       call wrap_put_vara_int (ncid(nf), mdcur_id(nf) , beg1d, len1d, mdcur)
       call wrap_put_vara_int (ncid(nf), mscur_id(nf) , beg1d, len1d, mscur)
       call wrap_put_vara_int (ncid(nf), nstep_id(nf) , beg1d, len1d, nstep)
       call wrap_put_vara_realx (ncid(nf), timvar_id(nf), beg1d, len1d, time)
    endif

! time comment for history interval

    if (masterproc) then
       beg2d(1) = 1        ; len2d(1) = len_trim(timcom(nf))
       beg2d(2) = ntim(nf) ; len2d(2) = 1
       call wrap_put_vara_text (ncid(nf), timcom_id(nf), beg2d, len2d, timcom(nf))
    endif

! set beginning and ending output indices

    if (slfld%num(nf) > 0) then
       if (hist_dov2xy(nf)) then
          beg3d(1) = 1        ; len3d(1) = lsmlon
          beg3d(2) = 1        ; len3d(2) = lsmlat
          beg3d(3) = ntim(nf) ; len3d(3) = 1
       else
          beg2d(1) = 1        ; len2d(1) = numpatch
          beg2d(2) = ntim(nf) ; len2d(2) = 1
       endif
    endif

    if (mlsoifld%num(nf) > 0) then
       if (hist_dov2xy(nf)) then
          beg4d(1) = 1       ; len4d(1) = lsmlon
          beg4d(2) = 1       ; len4d(2) = lsmlat
          beg4d(3) = 1       ; len4d(3) = nlevsoi
          beg4d(4) = ntim(nf); len4d(4) = 1
       else
          beg3d(1) = 1        ; len3d(1) = numpatch
          beg3d(2) = 1        ; len3d(2) = nlevsoi
          beg3d(3) = ntim(nf) ; len3d(3) = 1
       endif
    endif
       
! active single-level fields (either grid averages or 1-d vectors)

    if (slfld%num(nf) > 0) then
#if (defined SPMD)
       allocate (buf1d(begpatch:endpatch))
       allocate (gather1d(numpatch))
       call compute_mpigs_patch(1, numsend, numrecvv, displsv)
#endif
       do n = 1, slfld%num(nf)
#if (defined SPMD)
          do k = begpatch, endpatch
             buf1d(k) = slfld%value(k,n,nf)
          end do
          call mpi_gatherv (buf1d(begpatch), numsend , mpir8, &
               gather1d, numrecvv, displsv, mpir8, 0, mpicom, ier)
#else
          gather1d => slfld%value(:,n,nf)
#endif          
          if (masterproc) then
             if (hist_dov2xy(nf)) then
                call v2xy (gather1d, spval, slfxy)
                call wrap_put_vara_realx (ncid(nf), slfld_id(n,nf), beg3d, len3d, slfxy)
             else
                call wrap_put_vara_realx (ncid(nf), slfld_id(n,nf), beg2d, len2d, gather1d)
             endif
          endif
       end do
#if (defined SPMD)
       deallocate (buf1d)
       deallocate (gather1d)
#endif
    endif

! active multi-level soil fields (either grid averages or 1-d vectors)

    if (mlsoifld%num(nf) > 0) then
#if (defined SPMD)
       allocate (buf2d(nlevsoi,begpatch:endpatch))
       allocate (gather2d(nlevsoi,numpatch))
       allocate (output2d(numpatch,nlevsoi))
       call compute_mpigs_patch(nlevsoi, numsend, numrecvv, displsv)
#endif
       do n = 1, mlsoifld%num(nf)
#if (defined SPMD)
          do l = 1, nlevsoi
             do k = begpatch, endpatch
                buf2d(l,k) = mlsoifld%value(k,l,n,nf)
             end do
          end do
          call mpi_gatherv (buf2d(1,begpatch), numsend , mpir8, &
               gather2d, numrecvv, displsv, mpir8, 0, mpicom, ier)
          if (masterproc) then
             do l = 1, nlevsoi
                do k = 1,numpatch
                   output2d(k,l) = gather2d(l,k)
                end do
             end do
          endif
#else
          output2d => mlsoifld%value(:,:,n,nf)
#endif
          if (masterproc) then
             if (hist_dov2xy(nf)) then
                do l = 1, nlevsoi
                   call v2xy (output2d(1,l), spval, mlsoifxy(1,1,l))
                end do
                call wrap_put_vara_realx (ncid(nf), mlsoifld_id(n,nf),beg4d, len4d, mlsoifxy)
             else
                call wrap_put_vara_realx (ncid(nf), mlsoifld_id(n,nf),beg3d, len3d, output2d)
             endif
          endif
       end do
#if (defined SPMD)
       deallocate (buf2d)
       deallocate (gather2d)
       deallocate (output2d)
#endif
    endif

    return 
  end subroutine histwrt

!=======================================================================

  subroutine histcls (nf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! close netCDF file 
!
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------

    include 'netcdf.inc'

! ------------------------ arguments ---------------------------------
    integer, intent(in) :: nf           !history file number
! --------------------------------------------------------------------

    call wrap_close(ncid(nf))

    return
  end subroutine histcls

!=======================================================================

  subroutine histslf (name, fld)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! accumulate single-level field over history time interval
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varmap, only : begpatch, endpatch

! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: name            !field name
    real(r8), intent(in) :: fld(begpatch:endpatch)  !field values for current time step
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i,n,m,k     !loop indices
    character(len= 8) :: type
! -----------------------------------------------------------------

    do m = 1, nhist

       ! find field index. return if "name" is not on active list
       n = 0
       do i = 1, slfld%num(m)
          if (name == slfld%nam(i,m)) n = i
       end do
       if (n == 0) go to 1000
     
       ! determine field attributes
       type  =  slfld%typ(n,m)
       
!$OMP PARALLEL DO PRIVATE (K)
       do k = begpatch,endpatch
        
          ! accumulate field 
          if (fld(k) /= spval) then
             if (type == naver) then           !time average field
                if (slfld%count(k,n,m) == 0) slfld%value(k,n,m) = 0.  
                slfld%value(k,n,m) = slfld%value(k,n,m) + fld(k)
                slfld%count(k,n,m) = slfld%count(k,n,m) + 1
             else if (type == ncnst) then      !constant field value
                if (slfld%count(k,n,m) == 0) then
                   slfld%value(k,n,m) = fld(k)
                   slfld%count(k,n,m) = 1
                endif
             else if (type == ninst) then      !instantaneous field value
                slfld%value(k,n,m) = fld(k)
                slfld%count(k,n,m) = 1
             else if (type == nmaxi) then      !maximum field value
                if (slfld%count(k,n,m) == 0) slfld%value(k,n,m) = -1.e50
                slfld%value(k,n,m) = max( slfld%value(k,n,m), fld(k) )
                slfld%count(k,n,m) = 1
             else if (type == nmini) then      !minimum field value
                if (slfld%count(k,n,m) == 0) slfld%value(k,n,m) = +1.e50
                slfld%value(k,n,m) = min( slfld%value(k,n,m), fld(k) )
                slfld%count(k,n,m) = 1
             end if
          else
             if (slfld%count(k,n,m)== 0) slfld%value(k,n,m) = fld(k)
          endif
          
          ! end of history interval: normalize accumulated values 
          if (ehi(m)) then
             if (type == naver .and. slfld%count(k,n,m)/=0) then  
                slfld%value(k,n,m) = slfld%value(k,n,m) / float(slfld%count(k,n,m))
             end if
          endif
          
       end do
!$OMP END PARALLEL DO
     
1000   continue
    end do

    return
  end subroutine histslf

!=======================================================================

  subroutine histmlf (name, fld, nlev)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! accumulate multi-level field over history time interval
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varmap , only : begpatch, endpatch
    use clm_varpar , only : nlevsoi
  
! ------------------------ arguments ------------------------------
    character(len=*), intent(in) :: name                !field name
    integer , intent(in) :: nlev                        !number of levels
    real(r8), intent(in) :: fld(begpatch:endpatch,nlev) !field values for current time step
! -----------------------------------------------------------------

! ------------------------ local variables ------------------------
    integer i,j,n,m,k         !do loop indices
    character(len= 8) :: type
! -----------------------------------------------------------------

! loop over history tapes
    
    do m = 1, nhist
       
       ! find field index. return if "name" is not on active list
       n = 0
       do i = 1, mlsoifld%num(m)
          if (name == mlsoifld%nam(i,m)) n = i
       end do
       if (n == 0) go to 1000
       
       ! initialize field attributes
       type  =  mlsoifld%typ(n,m)
       
!$OMP PARALLEL DO PRIVATE (J,K)
       do k = begpatch,endpatch
          do j = 1, nlev
           
             ! accumulate field 
             if (fld(k,j) /= spval) then
                if (type == naver) then           !time average field
                   if (mlsoifld%count(k,j,n,m) == 0) mlsoifld%value(k,j,n,m) = 0.  
                   mlsoifld%value(k,j,n,m) = mlsoifld%value(k,j,n,m) + fld(k,j)
                   mlsoifld%count(k,j,n,m) = mlsoifld%count(k,j,n,m) + 1
                else if (type == ncnst) then      !constant field value
                   if (mlsoifld%count(k,j,n,m) == 0) then
                      mlsoifld%value(k,j,n,m) = fld(k,j)
                      mlsoifld%count(k,j,n,m) = 1
                   endif
                else if (type == ninst) then      !instantaneous field value
                   mlsoifld%value(k,j,n,m) = fld(k,j)
                   mlsoifld%count(k,j,n,m) = 1
                else if (type == nmaxi) then      !maximum field value
                   if (mlsoifld%count(k,j,n,m) == 0) mlsoifld%value(k,j,n,m) = -spval
                   mlsoifld%value(k,j,n,m) = max(mlsoifld%value(k,j,n,m),fld(k,j))
                   mlsoifld%count(k,j,n,m) = 1
                else if (type == nmini) then      !minimum field value
                   if (mlsoifld%count(k,j,n,m) == 0) mlsoifld%value(k,j,n,m) = +spval
                   mlsoifld%value(k,j,n,m) = min(mlsoifld%value(k,j,n,m),fld(k,j))
                   mlsoifld%count(k,j,n,m) = 1
                end if
             else
                if (mlsoifld%count(k,j,n,m)== 0) mlsoifld%value(k,j,n,m) = fld(k,j)
             endif
             
             ! end of history interval, normalize accumulated values 
             if (ehi(m)) then
                if (type==naver .and. mlsoifld%count(k,j,n,m)/=0) then
                   mlsoifld%value(k,j,n,m) = mlsoifld%value(k,j,n,m) / float(mlsoifld%count(k,j,n,m))
                endif
             endif
             
          end do
       end do
!$OMP END PARALLEL DO
     
1000   continue
    end do
  
    return
  end subroutine histmlf

!=======================================================================

  subroutine histzero(nfile)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! zero out history counters
!
! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------

    integer, intent(in) :: nfile ! history file index

! Reset counters to zero

    slfld%count(:,:,nfile) = 0
    mlsoifld%count(:,:,:,nfile) = 0

  end subroutine histzero

!=======================================================================

  subroutine histvar_ini

    use infnan

! Allocate dynamic history variables

    allocate (slfld%value(begpatch:endpatch,max_slevflds,nhist))
    allocate (slfld%count(begpatch:endpatch,max_slevflds,nhist))

    allocate (mlsoifld%value(begpatch:endpatch,nlevsoi,max_mlevflds,nhist))
    allocate (mlsoifld%count(begpatch:endpatch,nlevsoi,max_mlevflds,nhist))

! Set variables to infinity 

    slfld%value(:,:,:) = spval
    slfld%count(:,:,:) = 0

    mlsoifld%value(:,:,:,:) = spval
    mlsoifld%count(:,:,:,:) = 0

  end subroutine histvar_ini

!=======================================================================

end module histfileMod
