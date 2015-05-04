#include <misc.h>
#include <params.h>

! **************************************************************************** !
! comsrf module                                                                !
!                                                                              !
! Description:                                                                 !
!                                                                              !
!   Handle surface fluxes for the subcomponents of cam/csm                     !
!                                                                              !
! **************************************************************************** !

module comsrf

    use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
    use constituents, only: pcnst, pnats
    use ppgrid,       only: pcols, begchunk, endchunk, pvermx, pverp
    use phys_grid,    only: get_ncols_p
    use infnan,       only: inf, uninit_r8
    use pmgrid,       only: plon, plat

    implicit none

    public initialize_comsrf ! set the surface temperature and sea-ice fraction
    public srfflx_parm_reset
    public srfflx_state_reset
    public update_srf_fractions
    public update_srf_fluxes
    public verify_fractions

    ! DONG Li: Verify the public access permission
    !public   ! By default all data is public to this module

    public surface_state
    public srfflx_state
    public srfflx_parm

    integer, parameter :: plevmx = 4            ! number of subsurface levels

    character(8) tsnam(plevmx)                  ! names of sub-surface temperature fields

    real(r8), allocatable :: landm(:,:)         ! land/ocean/sea ice flag
    real(r8), allocatable :: sgh(:,:)           ! land/ocean/sea ice flag
    real(r8), allocatable :: sicthk(:,:)        ! cam sea-ice thickness (m)
    real(r8), allocatable :: snowhice(:,:)      ! snow depth (liquid water) ovr ice
    real(r8), allocatable :: snowhland(:,:)     !snow depth (liquid water) ovr lnd
    real(r8), allocatable :: fv(:,:)            ! needed for dry dep velocities (over land) ! For FGOALS2.0
    real(r8), allocatable :: ram1(:,:)          ! needed for dry dep velocities (over land) !
    real(r8), allocatable :: fsns(:,:)          ! surface absorbed solar flux
    real(r8), allocatable :: fsnt(:,:)          ! Net column abs solar flux at model top
    real(r8), allocatable :: flns(:,:)          ! Srf longwave cooling (up-down) flux
    real(r8), allocatable :: flnt(:,:)          ! Net outgoing lw flux at model top
    real(r8), allocatable :: srfrpdel(:,:)      ! 1./(pint(k+1)-pint(k))
    real(r8), allocatable :: psm1(:,:)          ! surface pressure
    real(r8), allocatable :: absorb(:,:)        ! cam surf absorbed solar flux (W/m2)
    real(r8), allocatable :: prcsnw(:,:)        ! cam tot snow precip
    real(r8), allocatable :: landfrac(:,:)      ! fraction of sfc area covered by land
    real(r8), allocatable :: landfrac_field(:,:)! fraction of sfc area covered by land (needed globally) ! For FGOALS2.0
    real(r8), allocatable :: ocnfrac(:,:)       ! frac of sfc area covered by open ocean
    real(r8), allocatable :: icefrac(:,:)       ! fraction of sfc area covered by seaice
    real(r8), allocatable :: previcefrac(:,:)   ! fraction of sfc area covered by seaice
    real(r8), allocatable :: trefmxav(:,:)      ! diagnostic: tref max over the day
    real(r8), allocatable :: trefmnav(:,:)      ! diagnostic: tref min over the day

    real(r8), allocatable:: asdirice(:,:)
    real(r8), allocatable:: aldirice(:,:)
    real(r8), allocatable:: asdifice(:,:)
    real(r8), allocatable:: aldifice(:,:)

    real(r8), allocatable:: asdirlnd(:,:)
    real(r8), allocatable:: aldirlnd(:,:)
    real(r8), allocatable:: asdiflnd(:,:)
    real(r8), allocatable:: aldiflnd(:,:)

    real(r8), allocatable:: asdirocn(:,:)
    real(r8), allocatable:: aldirocn(:,:)
    real(r8), allocatable:: asdifocn(:,:)
    real(r8), allocatable:: aldifocn(:,:)

    real(r8), allocatable:: lwuplnd(:,:)
    real(r8), allocatable:: lwupocn(:,:)
    real(r8), allocatable:: lwupice(:,:)

    real(r8), allocatable:: tsice(:,:)
    real(r8), allocatable:: tsice_rad(:,:)  ! equivalent upward long wave temperature
    real(r8), allocatable:: fld_kvh(:,:,:)  ! added by LIU Li

    !JR Renamed sst to tsocn to avoid potential conflict with other CAM variable of the same name
    real(r8), allocatable :: tsocn(:,:)
    real(r8), allocatable :: Focn(:,:)
    real(r8), allocatable :: frzmlt(:,:)
    real(r8), allocatable :: aice(:,:)      ! CSIM ice fraction

    type surface_state
        real(r8) tbot(pcols)         ! bot level temperature
        real(r8) zbot(pcols)         ! bot level height above surface
        real(r8) ubot(pcols)         ! bot level u wind
        real(r8) vbot(pcols)         ! bot level v wind
        real(r8) qbot(pcols)         ! bot level specific humidity
        real(r8) pbot(pcols)         ! bot level pressure
        real(r8) flwds(pcols)        !
        real(r8) precsc(pcols)       !
        real(r8) precsl(pcols)       !
        real(r8) precc(pcols)        !
        real(r8) precl(pcols)        !
        real(r8) soll(pcols)         !
        real(r8) sols(pcols)         !
        real(r8) solld(pcols)        !
        real(r8) solsd(pcols)        !
        real(r8) srfrad(pcols)       !
        real(r8) thbot(pcols)        !
        real(r8) tssub(pcols,plevmx) ! cam surface/subsurface temperatures
    end type surface_state

    type srfflx_state
        integer lchnk
        integer ncol

        real(r8) asdir(pcols)            ! albedo: shortwave, direct
        real(r8) asdif(pcols)            ! albedo: shortwave, diffuse
        real(r8) aldir(pcols)            ! albedo: longwave, direct
        real(r8) aldif(pcols)            ! albedo: longwave, diffuse
        real(r8) lwup(pcols)             ! longwave up radiative flux
        real(r8) lhf(pcols)              ! latent heat flux
        real(r8) shf(pcols)              ! sensible heat flux
        real(r8) wsx(pcols)              ! surface u-stress (N)
        real(r8) wsy(pcols)              ! surface v-stress (N)
        real(r8) tref(pcols)             ! ref height surface air temp
        real(r8) qref(pcols)             ! ref height surface specific humidity ! For FGOALS2.0
        real(r8) rhref(pcols)            ! ref height surface relative humidity ! For FGOALS2.0
        real(r8) ts(pcols)               ! sfc temp (merged w/ocean if coupled)
        real(r8) sst(pcols)              ! surface sea temperature              ! For FGOALS2.0
        real(r8) cflx(pcols,pcnst+pnats) ! constituent flux (evap)
    end type srfflx_state

    !---------------------------------------------------------------------------
    ! This is for surface fluxes returned from individual parameterizations
    !---------------------------------------------------------------------------

    type srfflx_parm
        character(24) name               ! name of parameterization which produced tendencies
        real(r8) asdir(pcols)            ! albedo: shortwave, direct
        real(r8) asdif(pcols)            ! albedo: shortwave, diffuse
        real(r8) aldir(pcols)            ! albedo: longwave, direct
        real(r8) aldif(pcols)            ! albedo: longwave, diffuse
        real(r8) lwup(pcols)             ! longwave up radiative flux
        real(r8) lhf(pcols)              ! latent heat flux
        real(r8) shf(pcols)              ! sensible heat flux
        real(r8) wsx(pcols)              ! surface u-stress (N)
        real(r8) wsy(pcols)              ! surface v-stress (N)
        real(r8) tref(pcols)             ! ref height surface air temp
        real(r8) ts(pcols)               ! sfc temp (merged w/ocean if coupled)
        real(r8) cflx(pcols,pcnst+pnats) ! pcnst+pnats constituent flux (evap)
    end type srfflx_parm

    type(surface_state), allocatable :: surface_state2d(:)
    type(srfflx_state), allocatable :: srfflx_state2d(:)
    type(srfflx_parm), allocatable :: srfflx_parm2d(:)
    type(srfflx_parm), allocatable :: srfflx_parm2d_ocn(:)

contains

!-----------------------------------------------------------------------
!
! BOP
!
! !IROUTINE: initialize_comsrf
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying ice boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method:
!
! Author:
!
!-----------------------------------------------------------------------
!
! !INTERFACE
!

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Initialize surface data
    !
    ! Method:
    !
    ! Author: Mariana Vertenstein
    !
    !-----------------------------------------------------------------------

    subroutine initialize_comsrf
        integer k, c      ! level, constituent indices

        allocate(landm   (pcols,begchunk:endchunk))
        allocate(sgh     (pcols,begchunk:endchunk))
        allocate(sicthk  (pcols,begchunk:endchunk))
        allocate(snowhice(pcols,begchunk:endchunk))
        allocate(snowhland(pcols,begchunk:endchunk))
        allocate(fv(pcols,begchunk:endchunk))
        allocate(ram1(pcols,begchunk:endchunk))
        allocate(fsns    (pcols,begchunk:endchunk))
        allocate(fsnt    (pcols,begchunk:endchunk))
        allocate(flns    (pcols,begchunk:endchunk))
        allocate(flnt    (pcols,begchunk:endchunk))
        allocate(lwupocn (pcols,begchunk:endchunk))
        allocate(lwuplnd (pcols,begchunk:endchunk))
        allocate(lwupice (pcols,begchunk:endchunk))
        allocate(srfrpdel(pcols,begchunk:endchunk))
        allocate(psm1    (pcols,begchunk:endchunk))
        allocate(absorb  (pcols,begchunk:endchunk))
        allocate(prcsnw  (pcols,begchunk:endchunk))
        allocate(landfrac(pcols,begchunk:endchunk))
        allocate(ocnfrac (pcols,begchunk:endchunk))
        allocate(icefrac (pcols,begchunk:endchunk))
        allocate(previcefrac(pcols,begchunk:endchunk))
        allocate(trefmxav(pcols,begchunk:endchunk))
        allocate(trefmnav(pcols,begchunk:endchunk))

        allocate(tsice_rad(pcols,begchunk:endchunk))
        allocate(asdirice(pcols,begchunk:endchunk))
        allocate(aldirice(pcols,begchunk:endchunk))
        allocate(asdifice(pcols,begchunk:endchunk))
        allocate(aldifice(pcols,begchunk:endchunk))

        allocate(asdirlnd(pcols,begchunk:endchunk))
        allocate(aldirlnd(pcols,begchunk:endchunk))
        allocate(asdiflnd(pcols,begchunk:endchunk))
        allocate(aldiflnd(pcols,begchunk:endchunk))

        allocate(asdirocn(pcols,begchunk:endchunk))
        allocate(aldirocn(pcols,begchunk:endchunk))
        allocate(asdifocn(pcols,begchunk:endchunk))
        allocate(aldifocn(pcols,begchunk:endchunk))

        allocate(tsice(pcols,begchunk:endchunk))
        allocate(tsocn(pcols,begchunk:endchunk))
        allocate(Focn(pcols,begchunk:endchunk))
        allocate(frzmlt(pcols,begchunk:endchunk))
        allocate(aice(pcols,begchunk:endchunk))

        allocate(srfflx_state2d(begchunk:endchunk))
        allocate(srfflx_parm2d(begchunk:endchunk))
        allocate(surface_state2d(begchunk:endchunk))
        allocate(srfflx_parm2d_ocn(begchunk:endchunk))
        allocate(fld_kvh(pcols,pverp,begchunk:endchunk)) ! added by LIU Li
        !
        ! Initialize to NaN or Inf
        !
        landm (:,:) = inf
        sgh   (:,:) = inf
        sicthk(:,:) = inf
        !
        ! snowhice and snowhland are only defined for a single surface type
        ! elements of the array outside valid surface points must be set to
        ! zero if these fields are to be written to netcdf history files.
        !
        snowhice (:,:) = 0.
        snowhland(:,:) = 0.
        fsns     (:,:) = inf
        fsnt     (:,:) = inf
        flns     (:,:) = inf
        flnt     (:,:) = inf
        lwupocn  (:,:) = 0.
        lwuplnd  (:,:) = 0.
        lwupice  (:,:) = 0.
        srfrpdel (:,:) = inf
        psm1     (:,:) = inf
        absorb   (:,:) = inf
        prcsnw   (:,:) = inf
        landfrac (:,:) = inf
        ocnfrac  (:,:) = inf
        icefrac  (:,:) = inf
        previcefrac(:,:) = uninit_r8
        trefmxav (:,:) = -1.0e36
        trefmnav (:,:) =  1.0e36

        asdirice (:,:) = 0.
        aldirice (:,:) = 0.
        asdifice (:,:) = 0.
        aldifice (:,:) = 0.

        asdirlnd (:,:) = 0.
        aldirlnd (:,:) = 0.
        asdiflnd (:,:) = 0.
        aldiflnd (:,:) = 0.

        asdirocn (:,:) = 0.
        aldirocn (:,:) = 0.
        asdifocn (:,:) = 0.
        aldifocn (:,:) = 0.

        tsice (:,:) = inf
        tsice_rad (:,:) = inf

        tsocn (:,:) = inf
        Focn (:,:) = inf
        frzmlt (:,:) = inf
        aice (:,:) = inf
        !
        ! Sub-surface temperatures
        !
        if (plevmx > 9) then
            call endrun('comsrf: initialize_comsrf: Cannot handle more than 9 subsurface levels')
        end if
        do k = 1, plevmx
            write(unit=tsnam(k), fmt='(''TS'',i1,''     '')') k
        end do

        do c = begchunk, endchunk
            srfflx_state2d(c)%lchnk = c
            srfflx_state2d(c)%ncol  = get_ncols_p(c)
            srfflx_state2d(c)%asdir  (:) = 0.
            srfflx_state2d(c)%asdif  (:) = 0.
            srfflx_state2d(c)%aldir  (:) = 0.
            srfflx_state2d(c)%aldif  (:) = 0.
            srfflx_state2d(c)%lwup   (:) = 0.
            srfflx_state2d(c)%lhf    (:) = 0.
            srfflx_state2d(c)%shf    (:) = 0.
            srfflx_state2d(c)%cflx   (:,:) = 0.
            srfflx_state2d(c)%wsx    (:) = 0.
            srfflx_state2d(c)%wsy    (:) = 0.
            srfflx_state2d(c)%tref   (:) = 0.
            srfflx_state2d(c)%ts     (:) = 0.
            srfflx_state2d(c)%sst    (:) = 0.

            surface_state2d(c)%tbot  (:) = 0.
            surface_state2d(c)%zbot  (:) = 0.
            surface_state2d(c)%ubot  (:) = 0.
            surface_state2d(c)%vbot  (:) = 0.
            surface_state2d(c)%qbot  (:) = 0.
            surface_state2d(c)%pbot  (:) = 0.
            surface_state2d(c)%thbot (:) = 0.
            surface_state2d(c)%tssub (:,:) = 0.
        end do

        call srfflx_parm_reset(srfflx_parm2d)
        call srfflx_parm_reset(srfflx_parm2d_ocn)

        allocate (landfrac_field(plon,plat))
        landfrac_field(:,:) = inf

        return
    end subroutine initialize_comsrf

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Zero fluxes that are update by land,ocn,ice sub processes
    !
    ! Method:
    !
    ! Author: John Truesdale
    !
    !-----------------------------------------------------------------------

    subroutine srfflx_parm_reset(parm)
        ! DONG Li: "parm" should be called changes, not tendencies?
        type(srfflx_parm), intent(inout) :: parm(begchunk:endchunk) ! parameterization tendencies
        integer c        ! chunk index

        call t_startf('srfflx_rst_st')
        do c = begchunk, endchunk
            parm(c)%asdir(:) = 0.
            parm(c)%asdif(:) = 0.
            parm(c)%aldir(:) = 0.
            parm(c)%aldif(:) = 0.
            parm(c)%lwup (:) = 0.
            parm(c)%lhf  (:) = 0.
            parm(c)%shf  (:) = 0.
            parm(c)%cflx (:,:) = 0.
            parm(c)%wsx  (:) = 0.
            parm(c)%wsy  (:) = 0.
            parm(c)%tref (:) = 0.
            parm(c)%ts   (:) = 0.
        end do
        call t_stopf ('srfflx_rst_st')

        return
    end subroutine srfflx_parm_reset

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Zero fluxes that are update by land,ocn,ice sub processes
    !
    ! Method:
    !
    ! Author: John Truesdale
    !
    !-----------------------------------------------------------------------

    subroutine srfflx_state_reset(state)
        type(srfflx_state), intent(inout) :: state(begchunk:endchunk)   ! srfflx state variables
        integer c

        call t_startf ('srfflx_rst_st')
        do c = begchunk, endchunk
            state(c)%asdir(:) = 0.
            state(c)%asdif(:) = 0.
            state(c)%aldir(:) = 0.
            state(c)%aldif(:) = 0.
            state(c)%lwup(:) = 0.
            state(c)%lhf(:) = 0.
            state(c)%shf(:) = 0.
            state(c)%cflx(:,:) = 0.
            state(c)%wsx(:) = 0.
            state(c)%wsy(:) = 0.
            state(c)%tref(:) = 0.
            state(c)%qref(:) = 0. ! added by DONG Li
            state(c)%ts(:) = 0.
            state(c)%sst(:) = 0. ! added by DONG Li
        end do
        call t_stopf ('srfflx_rst_st')

        return
    end subroutine srfflx_state_reset

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! update surface fluxes
    !
    ! Method:
    !
    ! Author: John Truesdale
    !
    !-----------------------------------------------------------------------

    subroutine update_srf_fluxes(state, parm, frac)

        use physconst, only: stebol
        use phys_grid, only: get_ncols_p

        type(srfflx_parm), intent(inout)  :: parm(begchunk:endchunk)
        type(srfflx_state), intent(inout) :: state(begchunk:endchunk)
        real(r8), intent(in) :: frac(pcols,begchunk:endchunk)

        integer i, c, m     ! longitude, chunk, constituent indices
        integer ncol        ! number of longitudes this chunk

        call t_startf ('update_srf_st')
        do c = begchunk, endchunk
            ncol = get_ncols_p(c)
            do i = 1, ncol
                if (frac(i,c) > 0.) then
                    state(c)%asdir(i) = state(c)%asdir(i)+parm(c)%asdir(i)*frac(i,c)
                    state(c)%asdif(i) = state(c)%asdif(i)+parm(c)%asdif(i)*frac(i,c)
                    state(c)%aldir(i) = state(c)%aldir(i)+parm(c)%aldir(i)*frac(i,c)
                    state(c)%aldif(i) = state(c)%aldif(i)+parm(c)%aldif(i)*frac(i,c)
                    state(c)%lwup(i)  = state(c)%lwup(i) +parm(c)%lwup(i) *frac(i,c)
                    state(c)%lhf(i)   = state(c)%lhf(i)  +parm(c)%lhf(i)  *frac(i,c)
                    state(c)%shf(i)   = state(c)%shf(i)  +parm(c)%shf(i)  *frac(i,c)
                    state(c)%wsx(i)   = state(c)%wsx(i)  +parm(c)%wsx(i)  *frac(i,c)
                    state(c)%wsy(i)   = state(c)%wsy(i)  +parm(c)%wsy(i)  *frac(i,c)
                    state(c)%tref(i)  = state(c)%tref(i) +parm(c)%tref(i) *frac(i,c)
                    !
                    ! if we are calculating ts for a non-fractional grid box (ie all land or
                    ! all ocean or all ice then use the ts given by the parameterization)
                    ! otherwise calculate ts based on the grid averaged lwup
                    !
                    ! jt pull this code out after bit for bit testing
                    !
                    if (frac(i,c) == 1.) then
                        state(c)%ts(i) = state(c)%ts(i)+parm(c)%ts(i)*frac(i,c)
                    else
                        state(c)%ts(i) = sqrt(sqrt(state(c)%lwup(i)/stebol))
                    end if

                    do m = 1, pcnst+pnats
                        state(c)%cflx(i,m) = state(c)%cflx(i,m)+parm(c)%cflx(i,m)*frac(i,c)
                    end do
                end if
            end do
        end do
        !
        ! zero srfflx parameterization
        !
        call t_stopf ('update_srf_st')
        call srfflx_parm_reset(parm)

        return
    end subroutine update_srf_fluxes

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! DONG Li: This comment is wrong.
    ! update surface fluxes
    !
    ! Method:
    !
    ! Make land, ocean and ice concentrations consistant - after advnce
    ! is called the ice fraction will change due to interpolation in
    ! iceint so we must reapportion the ocean fraction accordingly.
    !
    ! 1) land fraction is taken as given
    ! 2) ice percentage is allowed to occupy what is not land
    !    ice percentage may be modified so that the sum of land and ice
    !        fractions is <= 1.
    ! 3) ocean fraction is determined last.  Ocean occupies whatever is
    !    not land or ice
    !
    ! Author: John Truesdale
    !
    !-----------------------------------------------------------------------

    subroutine update_srf_fractions

        use phys_grid, only: get_ncols_p

        integer j, i, ncol

        call t_startf ('update_srf_st')
        do j = begchunk, endchunk
            ncol = get_ncols_p(j)
            do i = 1, ncol
            ! fix inconsistancies between ice and land fractions, use land as truth
                if ((icefrac(i,j) > 1.0 .or. icefrac(i,j) < 0.) .or. &
                    (landfrac(i,j) > 1.0.or.landfrac(i,j) < 0.)) then
                    write(6, "('Error: comsrf: update_srf_fluxes: ')", advance="no")
                    write(6, "('In orography fractions, encountered a fraction greater than 1 or less than 0.')")
                    call endrun
                end if
                if (icefrac(i,j)+landfrac(i,j) > 1.0_r8) then
                    icefrac(i,j) = 1.0_r8-landfrac(i,j)
                end if
                ocnfrac(i,j) = 1.0_r8-landfrac(i,j)-icefrac(i,j)
            end do
        end do
        call t_stopf ('update_srf_st')

        return
    end subroutine update_srf_fractions

    !-----------------------------------------------------------------------
    !
    ! Purpose: Ensure that input array frac is bounded between 0 and 1 to
    !          within roundoff.  Only SOM currently requires the epsilon ()
    !          application.
    !
    ! Method:
    !
    ! Author:
    !
    !-----------------------------------------------------------------------

    subroutine bounding(frac, c, ncol, string)
        real(r8), intent(in) :: frac(pcols)      ! surface fraction (land, ocn, or ice)
        integer, intent(in) :: c                 ! chunk index
        integer, intent(in) :: ncol              ! number of columns
        character(*), intent(in) :: string       ! string for error print

        integer i, isave   ! loop index
        integer csave      ! chunk index of bad point
        logical bad        ! flag indicates some bad points found
        real(r8) fracsave  ! bad fraction value

        bad = .false.
        do i = 1, ncol
            if (frac(i) < -10.0*epsilon (1.0_r8) .or. frac(i) > 1.0+10.0*epsilon(1.0_r8)) then
                bad = .true.
                fracsave = frac(i)
                isave = i
                csave = c
            end if
        end do

        if (bad) then
            write(6, "('Error: comsrf: bounding: ')", advance="no")
            write(6, "('bad fraction ', F, 'of ', A, ' at column ', I5, ', chunk ', I5)") string, isave, csave
            call endrun
        end if

        return
    end subroutine bounding

    !-----------------------------------------------------------------------
    !
    ! Purpose: Ensure that surface fractions are valid
    !
    ! Method:
    !
    ! Author:
    !
    !-----------------------------------------------------------------------

    subroutine verify_fractions(c, ncol)

        use phys_grid, only: get_ncols_p

        integer, intent(in) :: c     ! chunk index
        integer, intent(in) :: ncol  ! number of columns

        integer i                      ! loop index
        logical bad                    ! flag
        real(r8) icesv, ocnsv, lndsv   ! saved values
        real(r8) sfcfrac               ! icefrac + ocnfrac + landfrac
        real(r8) delta                 ! 1. - sum of sfc fractions

        !call t_startf ('verify_fractions')

        call bounding(icefrac(1,c), c, ncol, 'ICE')
        call bounding(ocnfrac(1,c), c, ncol, 'OCN')
        call bounding(landfrac(1,c), c, ncol, 'LND')

        bad = .false.
        do i = 1, ncol
            sfcfrac = icefrac(i,c)+ocnfrac(i,c)+landfrac(i,c)
            delta   = 1.0-sfcfrac
            if (abs(delta) > 10.0*epsilon(1.0_r8)) then
                bad = .true.
                icesv = icefrac(i,c)
                ocnsv = ocnfrac(i,c)
                lndsv = landfrac(i,c)
            end if
        end do

        if (bad) then
            write(6,*)'VERIFY_FRACTIONS: total surface fraction: ', icesv, ocnsv, lndsv
            call endrun ()
        end if

        !call t_stopf ('verify_fractions')

        return
    end subroutine verify_fractions

end module comsrf
