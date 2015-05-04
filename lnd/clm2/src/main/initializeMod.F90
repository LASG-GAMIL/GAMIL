#include <misc.h>
#include <preproc.h>

module initializeMod

    use spmdMod, only : masterproc

!=======================================================================
CONTAINS
!=======================================================================

#if (defined OFFLINE) || (defined COUP_CSM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp    )
#elif (defined COUP_CAM)
  subroutine initialize(eccen       , obliqr      , lambm0    , mvelpp    , cam_caseid , &
                        cam_ctitle  , cam_nsrest  , cam_nstep , cam_irad  , cam_crtinic, &
                        cam_nhtfrq  , cam_mfilt   , cam_longxy, cam_latixy, cam_numlon , &
                        cam_landmask, cam_landfrac, cam_irt   )
#endif

!-----------------------------------------------------------------------
!
! Purpose:
! Land model initialization.
!
! Method:
! Initialization routine for land surface model. This subroutine:
!   o Initializes run control variables via the [clmexp] namelist
!   o Reads in and/or creates surface data on [lsmlon] x [lsmlat] grid
!   o Defines the multiple plant types and fraction areas for each surface type
!   o Builds the appropriate subgrid <-> grid mapping indices and weights
!   o Assigns subgrid patches the appropriate time-constant data
!   o Initializes history file variables
!   o Initializes RTM river routing model (if cpp varialbe RTM is defined)
!   o Reads in restart data (for a restart or branch run)
!   o Reads in initial data and initializes the time variant variables
!     (for an initial run)
!
! Author: Gordon Bonan
!
!-----------------------------------------------------------------------
! $Id: initializeMod.F90,v 1.7.2.8.6.1 2002/10/03 20:07:41 erik Exp $
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varder
    use clm_varpar     , only : lsmlon, lsmlat, maxpatch
    use clm_varmap     , only : numland, numpatch, begpatch, endpatch
    use clm_varsur     , only : varsur_alloc, varsur_dealloc
    use clm_varctl     , only : fsurdat, finidat, nsrest, irad, &
                                mksrf_offline_fgrid, mksrf_offline_fnavyoro
    use controlMod     , only : control_init, control_print
    use fileutils      , only : lsmiou
    use mksrfdatMod    , only : mksrfdat
    use surfFileMod    , only : surfrd
    use pftcFileMod    , only : pftconrd
    use mvegFileMod    , only : interpMonthlyVeg
    use histFileMod    , only : histini
    use restFileMod    , only : restrd
#if (defined COUP_CAM)
    use time_manager   , only:  get_curr_date, get_nstep
#else
    use time_manager   , only:  get_curr_date, timemgr_init, get_nstep, &
	                        advance_timestep
#endif
#if (defined OFFLINE)
    use atmdrvMod      , only : atm_getgrid
#endif
#if (defined RTM)
    use RtmMod         , only : Rtmgridini, Rtmlandini
#endif
#ifdef COUP_CSM
    use clm_csmMod
#endif
    implicit none

! ------------------------ arguments ----------------------------------
    real(r8), intent(inout) :: eccen    !Earth's orbital eccentricity
    real(r8), intent(inout) :: obliqr   !Earth's obliquity in radians
    real(r8), intent(inout) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(inout) :: mvelpp   !Earth's moving vernal equinox longitude of perihelion + pi (radians)
#if (defined COUP_CAM)
    character(len=*),  intent(in) :: cam_caseid   !cam caseid
    character(len=*),  intent(in) :: cam_ctitle   !cam title
    character(len=*),  intent(in) :: cam_crtinic  !cam initial dataset generation frequency
    integer ,  intent(in) :: cam_irad             !cam radiation frequency
    integer ,  intent(in) :: cam_nsrest           !cam 0=initial run, > 0=continuation run
    integer ,  intent(in) :: cam_nstep            !cam current time index
    integer ,  intent(in) :: cam_nhtfrq           !cam history write freq for tape 1
    integer ,  intent(in) :: cam_mfilt            !cam number of files per tape for tape 1
    integer ,  intent(in) :: cam_irt              !cam mss retention time
    integer ,  intent(in) :: cam_numlon(:)        !cam number of longitudes
    real(r8),  intent(in) :: cam_longxy(:,:)      !cam lon values
    real(r8),  intent(in) :: cam_latixy(:,:)      !cam lat values
    real(r8),  intent(in) :: cam_landfrac(:,:)    !cam fractional land
    integer ,  intent(in) :: cam_landmask(:,:)    !cam land mask
#endif
! ----------------------------------------------------------------------

! ------------------------ local variables -----------------------------
    logical  :: readini               !true if read in initial data set
    integer  :: i,j,k                 !indices
    integer  :: yr                    !current year (0, ...)
    integer  :: mon                   !current month (1 -> 12)
    integer  :: day                   !current day (1 -> 31)
    integer  :: ncsec                 !current time of day [seconds]
    integer  :: vegxy(lsmlon,lsmlat,maxpatch) !vegetation type
    real(r8) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid weights
#ifdef COUP_CSM
    integer  :: cam_numlon(lsmlat)           !cam number of longitudes
    real(r8) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8) :: cam_latixy(lsmlon,lsmlat)    !cam lat values
    real(r8) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer  :: cam_landmask(lsmlon,lsmlat)  !cam land mask
#endif
    integer  :: nstep                        !current time step
! ----------------------------------------------------------------------

! Echo initialization to standard output

    call header()

    if (masterproc) then
       write (6,*) 'Attempting to initialize the land model .....'
       write (6,*)
    endif

! ----------------------------------------------------------------------
! Initialize run control variables, time manager, timestep
! ----------------------------------------------------------------------

#if (defined COUP_CAM)
    call control_init (cam_caseid , cam_ctitle, cam_irad , cam_nsrest, &
                       cam_crtinic, cam_nhtfrq, cam_mfilt, cam_irt )
#else
    call control_init ()
#endif
    if (masterproc) call control_print()

#if (defined OFFLINE) || (defined COUP_CSM)

! Initialize time manager for initial run

    if (nsrest == 0) call timemgr_init()

#endif

! Allocate surface grid dynamic memory (only needed for initialization)

    call varsur_alloc()

#if (defined OFFLINE)

! start at nstep = 1 for an initial offline run

  if (nsrest == 0) call advance_timestep()

#endif

! ----------------------------------------------------------------------
! Initialize Fortran unit numbers 1 to 99 to inactive: except standard
! input (5) and standard output (6)
! ----------------------------------------------------------------------

    lsmiou(:) = .false.
    lsmiou(5) = .true.
    lsmiou(6) = .true.

    if (masterproc) then
       write (6,*) 'Preset Fortran unit numbers:'
       write (6,*) '   unit  5 = standard input'
       write (6,*) '   unit  6 = standard output'
    endif

#if (defined RTM)
! ----------------------------------------------------------------------
! Initialize RTM river routing grid and mask
! ----------------------------------------------------------------------

     call Rtmgridini ()

#endif

#ifdef COUP_CSM
! ----------------------------------------------------------------------
! Initial flux coupler communication
! ----------------------------------------------------------------------

! Receive orbital parameters from flux coupler

    call csm_recvorb (eccen, obliqr, lambm0, mvelpp)

! Send control data to flux coupler

    call csm_sendcontrol (irad)

! Get grid and land mask back from flux coupler

    call csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

#endif

! ----------------------------------------------------------------------
! If no surface dataset name is specified then make surface dataset
! from original data sources. Always read surface boundary data in.
! This insures that bit for bit results are obtained for a run where a
! surface dataset file is generated and a run where a surface dataset
! is specified and read in. Set up vegetation type [veg] and weight [wt]
! arrays for [maxpatch] subgrid patches.
! ----------------------------------------------------------------------

#if (defined OFFLINE)

    if (fsurdat == ' ') then
       call mksrfdat ()
    endif
    call surfrd (vegxy, wtxy)

#else

    if (fsurdat == ' ') then
       call mksrfdat (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
    endif
    call surfrd (vegxy, wtxy, &
         cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)

#endif

! ----------------------------------------------------------------------
! Read list of PFTs and their corresponding parameter values
! ----------------------------------------------------------------------

    call pftconrd ()

! ----------------------------------------------------------------------
! Use [wtxy] to build mapping indices and weights:
! [lsmlon] x [lsmlat] grid <-> [numland] vector of land points <->
! [numpatch] vector of subgrid points
! ----------------------------------------------------------------------

    call clm_map (wtxy)

! ----------------------------------------------------------------------
! Initialize time invariant variables as subgrid vectors of length [numpatch]
! ----------------------------------------------------------------------

    if (masterproc) then
       write (6,*) ('Attempting to initialize time invariant variables')
    endif

    call iniTimeConst (vegxy)

    if (masterproc) then
       write (6,*) ('Successfully initialized time invariant variables')
       write (6,*)
    endif

! ----------------------------------------------------------------------
! Initialize river routing model(s)
! ----------------------------------------------------------------------

#if (defined RTM)

    call Rtmlandini ()

#endif


#ifdef COUP_CSM
! Send valid ocean runoff points to coupler


    call csm_sendrunoff ()

#endif

! ----------------------------------------------------------------------
! Read restart files if continuation run
! ----------------------------------------------------------------------

    if (nsrest > 0) call restrd ()

#if (defined OFFLINE)
! ----------------------------------------------------------------------
! Read atmospheric forcing dataset one time to obtain the longitudes
! and latitudes of the atmospheric dataset, as well as the edges. When
! coupled to atm model, these are input variables. If no
! atmospheric data files are provided, model uses dummy atmospheric
! forcing and sets atmospheric grid to land grid.
! ----------------------------------------------------------------------

    if (masterproc) then
       write (6,*) 'Attempting to set up atmospheric grid .....'
    endif

    call atm_getgrid ()

    if (masterproc) then
       write (6,*) 'Succesfully set up atmospheric grid .....'
       write (6,*)
    endif
#endif

! ----------------------------------------------------------------------
! Get current date
! ----------------------------------------------------------------------

    call get_curr_date(yr, mon, day, ncsec)

! ----------------------------------------------------------------------
! Initialize variables for history files
! ----------------------------------------------------------------------

    call histini ()

! ----------------------------------------------------------------------
! Read monthly vegetation data for interpolation to daily values
! Note: routine mapinit needs to be called first since interpMonthlyVeg
! needs the array landvec%wtxy
! ----------------------------------------------------------------------

    call interpMonthlyVeg (fsurdat, mon, day)

! ----------------------------------------------------------------------
! If initial run: initialize time-varying data
! If continuation run: end of initialization because time varying
! read in from restart file
! ----------------------------------------------------------------------

    if (nsrest == 0) then

       if (masterproc) then
          write (6,*) ('Attempting to initialize time variant variables .....')
       endif

       if (finidat == ' ') then
          readini = .false.
       else
          readini = .true.
       end if

       call iniTimeVar (readini, eccen, obliqr, lambm0 , mvelpp)

       if (masterproc) then
          write (6,*) ('Successfully initialized time variant variables')
          write (6,*)
       endif

    endif

#ifdef COUP_CSM
! ----------------------------------------------------------------------
! Send first land model data to flux coupler.
! ----------------------------------------------------------------------

  call csm_sendalb ()

#endif

! Deallocate surface grid dynamic memory

    call varsur_dealloc ()

! ----------------------------------------------------------------------
! End initialization
! ----------------------------------------------------------------------

    if (masterproc) then
       write (6,*) ('Successfully initialized the land model')
       if (nsrest == 0) then
          write (6,*) 'begin initial run at: '
       else
          write (6,*) 'begin continuation run at:'
       end if
       write (6,*) '   nstep= ',get_nstep(), &
            ' year= ',yr,' month= ',mon,' day= ',day,' seconds= ',ncsec
       write (6,*)
       write (6,'(72a1)') ("*",i=1,60)
       write (6,*)
    endif

    return
  end subroutine initialize

!=======================================================================

  subroutine header

!-----------------------------------------------------------------------
!
! Purpose:
! echo and save model version number
!
! Method:
!
! Author: Gordon Bonan
!
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl, only : version
    implicit none

    version = 'CLM MODEL version 2.0'

    if ( masterproc )then
      write (6,*) trim(version)
      write (6,*)
    end if

    return
  end subroutine header

!=======================================================================

end module initializeMod
