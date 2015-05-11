!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use ppgrid,         only: pcols, pver, pverp
    use constituents,   only: pcnst, pnats, qmin, cnst_name
    use tracers,        only: ixcldw
    use phys_grid,      only: get_ncols_p, get_rlon_all_p, get_rlat_all_p, get_gcol_all_p

    implicit none

    private

    public physics_state
    public physics_tend
    public physics_ptend

    public physics_update
    public physics_ptend_reset
    public physics_ptend_init
    public physics_ptend_sum      ! added by SHI Xiangjun
    public physics_tend_init      !
    public physics_state_copy     !
    public physics_state_set_grid


    ! DONG Li: What is the difference between "active" columns (ncol) and general columns (pcols)?
    type physics_state
        integer lchnk             ! chunk index
        integer ncol              ! number of active columns
        integer ulatcnt           ! number of unique lats in chunk                ! added by SHI Xiangjun
        integer uloncnt           ! number of unique lons in chunk                !
        ! 1D variables
        real(r8) lat(pcols)
        real(r8) lon(pcols)
        real(r8) ulat(pcols)      ! unique latitudes                              ! added by SHI Xiangjun
        real(r8) ulon(pcols)      ! unique longitudes                             !
        integer latmapback(pcols) ! map from column to unique lat for that column !
        integer lonmapback(pcols) ! map from column to unique lon for that column !
        integer cid(pcols)        ! unique column id                              !
        ! 2D variables
        real(r8) ps(pcols)
        real(r8) phis(pcols)
        ! 3D variables on the model layers
        real(r8) t(pcols,pver)
        real(r8) t1(pcols,pver) ! added by LI Lijuan
        real(r8) t2(pcols,pver) !
        real(r8) u(pcols,pver)
        real(r8) v(pcols,pver)
        real(r8) s(pcols,pver)  ! dry static energy
        real(r8) omega(pcols,pver)
        real(r8) pmid(pcols,pver)
        real(r8) pdel(pcols,pver)
        real(r8) rpdel(pcols,pver)
        real(r8) lnpmid(pcols,pver) ! ln(pmid)
        real(r8) exner(pcols,pver)  ! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
        real(r8) zm(pcols,pver)
        ! 3D variables on the interface layers
        real(r8) pint(pcols,pverp)
        real(r8) lnpint(pcols,pverp)
        real(r8) zi(pcols,pverp)
        ! constituent variables
        real(r8) q(pcols,pver,pcnst+pnats)
        real(r8) q1(pcols,pver) ! added by LI Lijuan
        real(r8) q2(pcols,pver) !
    end type physics_state

    type physics_tend
        real(r8) dtdt(pcols,pver)
        real(r8) dudt(pcols,pver)
        real(r8) dvdt(pcols,pver)
        real(r8) flx_net(pcols)
    end type physics_tend


    ! This is for tendencies returned from individual parameterizations
    type physics_ptend
        character(24) :: name       ! name of parameterization which produced tendencies.

        logical ls                  ! true if dsdt is returned
        logical lu                  ! true if dudt is returned
        logical lv                  ! true if dvdt is returned
        logical lq(pcnst+pnats)     ! true if dqdt() is returned

        integer top_level           ! top level index for which nonzero tendencies have been set
        integer bot_level           ! bottom level index for which nonzero tendencies have been set

        real(r8) s(pcols,pver)              ! heating rate (J/kg/s)
        real(r8) u(pcols,pver)              ! u momentum tendency (m/s/s)
        real(r8) v(pcols,pver)              ! v momentum tendency (m/s/s)
        real(r8) q(pcols,pver,pcnst+pnats)  ! consituent tendencies (kg/kg/s)
    end type physics_ptend

contains
!===============================================================================

!===============================================================================
  subroutine physics_update(state, tend, ptend, dt)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use geopotential, only: geopotential_dse
    use physconst,    only: cpair, gravit, rair, zvir
!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables
    type(physics_tend ), intent(inout)  :: tend    ! Physics tendencies

    real(r8), intent(in) :: dt                     ! time step
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3
!-----------------------------------------------------------------------
    ncol = state%ncol

! Update u,v fields
    if(ptend%lu) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%u  (i,k) = state%u  (i,k) + ptend%u(i,k) * dt
             tend%dudt(i,k) = tend%dudt(i,k) + ptend%u(i,k)
          end do
       end do
    end if

    if(ptend%lv) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%v  (i,k) = state%v  (i,k) + ptend%v(i,k) * dt
             tend%dvdt(i,k) = tend%dvdt(i,k) + ptend%v(i,k)
          end do
       end do
    end if

! Update dry static energy
    if(ptend%ls) then
       do k = ptend%top_level, ptend%bot_level
          do i = 1, ncol
             state%s(i,k)   = state%s(i,k)   + ptend%s(i,k) * dt
             tend%dtdt(i,k) = tend%dtdt(i,k) + ptend%s(i,k)/cpair
          end do
       end do
    end if

! Update constituents, all schemes use time split q: no tendency kept
    do m = 1, pcnst+pnats
       if(ptend%lq(m)) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1,ncol
                state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
             end do
          end do
! special test for cloud water (zero +/- 1.e-12 = 0)
          if (ptend%name == 'pcond' .and. m == ixcldw) then
             do k = 1,pver
                do i = 1,ncol
                   if (abs(state%q(i,k,m)) < 1.e-12) state%q(i,k,m) = 0.
                end do
             end do
          end if
! now test for mixing ratios which are too small
          name = trim(ptend%name) // '/' // trim(cnst_name(m))
          call qneg3(trim(name), state%lchnk, ncol, pcols, pver, 1, qmin(m), state%q(1,1,m))
       end if
    end do

! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(                                                                    &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(1,1,1), rair      , gravit      , cpair       , zvir         , &
            state%t     , state%zi    , state%zm    , ncol         )
    end if

! Reset all parameterization tendency flags to false
    call physics_ptend_reset(ptend)

  end subroutine physics_update

!===============================================================================
  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) ptend%s = 0.
    if(ptend%lu) ptend%u = 0.
    if(ptend%lv) ptend%v = 0.
    do m = 1, pcnst+pnats
       if(ptend%lq(m)) ptend%q(:,:,m) = 0.
    end do

    ptend%name  = "none"
    ptend%lq(:) = .FALSE.
    ptend%ls    = .FALSE.
    ptend%lu    = .FALSE.
    ptend%lv    = .FALSE.

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend)
!-----------------------------------------------------------------------
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    ptend%name  = "none"
    ptend%lq(:) = .true.
    ptend%ls    = .true.
    ptend%lu    = .true.
    ptend%lv    = .true.

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!
  subroutine physics_ptend_sum(ptend, ptend_sum, state)         !!sxj
    type(physics_ptend), intent(in)     :: ptend   ! New parameterization tendencies
    type(physics_ptend), intent(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    type(physics_state), intent(in)     :: state   ! New parameterization tendencies
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ncol                                ! number of columns
    ncol = state%ncol
! Update u,v fields
    if(ptend%lu) then
       ptend_sum%lu = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%u(i,k) = ptend_sum%u(i,k) + ptend%u(i,k)
          end do
       end do
    end if
    if(ptend%lv) then
       ptend_sum%lv = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%v(i,k) = ptend_sum%v(i,k) + ptend%v(i,k)
          end do
       end do
    end if
    if(ptend%ls) then
       ptend_sum%ls = .true.
       do i = 1, ncol
          do k = ptend%top_level, ptend%bot_level
             ptend_sum%s(i,k) = ptend_sum%s(i,k) + ptend%s(i,k)
          end do
       end do
    end if
! Update constituents
    do m = 1, pcnst+pnats
       if(ptend%lq(m)) then
          ptend_sum%lq(m) = .true.
          do i = 1,ncol
             do k = ptend%top_level, ptend%bot_level
                ptend_sum%q(i,k,m) = ptend_sum%q(i,k,m) + ptend%q(i,k,m)
             end do
          end do
       end if
    end do
  end subroutine physics_ptend_sum

!
  subroutine physics_tend_init(tend)   !!sxj
    implicit none
    type(physics_tend), intent(inout) :: tend
    tend%dtdt  = 0._r8
    tend%dudt = 0._r8
    tend%dvdt  = 0._r8
    tend%flx_net   = 0._r8
  end subroutine physics_tend_init

    subroutine physics_state_copy(state_in, state_out)
        type(physics_state), intent(in) :: state_in
        type(physics_state), intent(out) :: state_out

        integer i, k, m, ncol

        ncol = state_in%ncol
        state_out%lchnk = state_in%lchnk
        state_out%ncol  = state_in%ncol
        do i = 1, ncol
            state_out%ps(i)   = state_in%ps(i)
            state_out%phis(i) = state_in%phis(i)
        end do
        do k = 1, pver
            do i = 1, ncol
                state_out%t(i,k)      = state_in%t(i,k)
                state_out%u(i,k)      = state_in%u(i,k)
                state_out%v(i,k)      = state_in%v(i,k)
                state_out%s(i,k)      = state_in%s(i,k)
                state_out%omega(i,k)  = state_in%omega(i,k)
                state_out%pmid(i,k)   = state_in%pmid(i,k)
                state_out%pdel(i,k)   = state_in%pdel(i,k)
                state_out%rpdel(i,k)  = state_in%rpdel(i,k)
                state_out%lnpmid(i,k) = state_in%lnpmid(i,k)
                state_out%exner(i,k)  = state_in%exner(i,k)
                state_out%zm(i,k)     = state_in%zm(i,k)
            end do
        end do
        do k = 1, pverp
            do i = 1, ncol
                state_out%pint(i,k)      = state_in%pint(i,k)
                state_out%lnpint(i,k)    = state_in%lnpint(i,k)
                state_out%zi(i,k)        = state_in% zi(i,k)
            end do
        end do
        do m = 1, pcnst+pnats
            do k = 1, pver
                do i = 1, ncol
                    state_out%q(i,k,m) = state_in%q(i,k,m)
                end do
            end do
        end do

    end subroutine physics_state_copy

!!! sxj--2009-0309
  subroutine physics_state_set_grid(lchnk, phys_state)
    !use abortutils, only : endrun
!-----------------------------------------------------------------------
! Set the grid components of the physics_state object
!-----------------------------------------------------------------------

    integer,             intent(in)    :: lchnk
    type(physics_state), intent(inout) :: phys_state

    ! local variables
    integer  :: i, ncol
    real(r8) :: rlon(pcols)
    real(r8) :: rlat(pcols)
    !-----------------------------------------------------------------------

    ncol = get_ncols_p(lchnk)

    if(ncol<=0) then
       print *, lchnk, ncol
       call endrun('physics_state_set_grid')
    end if

    call get_rlon_all_p(lchnk, ncol, rlon)
    call get_rlat_all_p(lchnk, ncol, rlat)
    phys_state%ncol  = ncol
    phys_state%lchnk = lchnk
    do i=1,ncol
       phys_state%lat(i) = rlat(i)
       phys_state%lon(i) = rlon(i)
    end do
    call init_geo_unique(phys_state,ncol)

  end subroutine physics_state_set_grid


  subroutine init_geo_unique(phys_state,ncol)
    integer,             intent(in)    :: ncol
    type(physics_state), intent(inout) :: phys_state
    logical :: match
    integer :: i, j, ulatcnt, uloncnt

    phys_state%ulat=-999.
    phys_state%ulon=-999.
    phys_state%latmapback=0
    phys_state%lonmapback=0
    match=.false.
    ulatcnt=0
    uloncnt=0
    match=.false.
    do i=1,ncol
       do j=1,ulatcnt
          if(phys_state%lat(i) .eq. phys_state%ulat(j)) then
             match=.true.
             phys_state%latmapback(i)=j
          end if
       end do
       if(.not. match) then
          ulatcnt=ulatcnt+1
          phys_state%ulat(ulatcnt)=phys_state%lat(i)
          phys_state%latmapback(i)=ulatcnt
       end if

       match=.false.
       do j=1,uloncnt
          if(phys_state%lon(i) .eq. phys_state%ulon(j)) then
             match=.true.
             phys_state%lonmapback(i)=j
          end if
       end do
       if(.not. match) then
          uloncnt=uloncnt+1
          phys_state%ulon(uloncnt)=phys_state%lon(i)
          phys_state%lonmapback(i)=uloncnt
       end if
       match=.false.

    end do
    phys_state%uloncnt=uloncnt
    phys_state%ulatcnt=ulatcnt

    call get_gcol_all_p(phys_state%lchnk,pcols,phys_state%cid)

  end subroutine init_geo_unique




end module physics_types
