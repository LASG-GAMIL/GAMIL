#include <misc.h>
#include <params.h>

module vertical_diffusion

    !---------------------------------------------------------------------------------
    ! Module to compute vertical diffusion of momentum,  moisture, trace constituents
    ! and temperature. Uses turbulence module to get eddy diffusivities, including boundary
    ! layer diffusivities and nonlocal transport terms. Computes and applies molecular
    ! diffusion, if model top is high enough (above ~90 km).
    !
    ! calling sequence:
    !
    !    vd_inti        initializes vertical diffustion constants
    !      trbinti      initializes turblence/pbl constants
    !     .
    !     .
    !    vd_intr            interface for vertical diffusion and pbl scheme
    !      vdiff            performs vert diff and pbl
    !        trbintr        compute turbulent diffusivities and boundary layer cg terms
    !        vd_add_cg      add boudary layer cg term to profiles
    !        vd_add_ent     add boundary layer entrainment term to profiles
    !        vd_lu_decomp   perform lu decomposition of vertical diffusion matrices
    !        vd_lu_qmolec   update decomposition for molecular diffusivity of a constituent
    !        vd_lu_solve    solve the euler backward problem for vertical diffusion 
    !
    !---------------------------Code history--------------------------------
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          P. Rasch, B. Boville, August 1992
    ! Reviewed:          P. Rasch, April 1996
    ! Reviewed:          B. Boville, April 1996
    ! rewritten:         B. Boville, May 2000
    ! more changes       B. Boville, May-Aug 2001
    !---------------------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,       only: pcols, pver, pverp
    use constituents, only: pcnst, pnats, cnst_name, cnst_add, qmincg
    use tracers,      only: ixcldw
    use pmgrid,       only: masterproc

    implicit none

    private
    save
    !
    ! Public interfaces
    !
    public vd_inti     ! Initialization
    public vd_register ! register pbuf  !!sxj
    public vd_intr     ! full routine

#include<RK_or_MG.h>
    !
    ! Private data
    !
    real(r8), private :: cpair      ! Specific heat of dry air
    real(r8), private :: gravit     ! Acceleration due to gravity
    real(r8), private :: rair       ! Gas const for dry air
    real(r8), private :: zvir       ! rh2o/rair - 1

    real(r8), parameter :: km_fac = 3.55E-7 ! molecular viscosity constant
    real(r8), parameter :: pr_num = 1.       ! Prandtl number
    real(r8), parameter :: d0     = 1.52E20  ! diffusion factor m-1 s-1 molec sqrt(kg/kmol/K)
    ! Aerononmy, Part B, Banks and Kockarts (1973), p39
    ! note text cites 1.52E18 cm-1 ...
    real(r8), parameter :: n_avog = 6.022E26 ! Avogadro's number (molec/kmol)

    real(r8), parameter :: kq_fac = 2.52E-13 ! molecular constituent diffusivity constant
    real(r8), parameter :: mw_dry = 29.      ! molecular weight of dry air ****REMOVE THIS***
    real(r8) :: mw_fac(pcnst+pnats)          ! sqrt(1/M_q + 1/M_d) in constituent diffusivity

    integer, private :: ntop        ! Top level to which vertical diffusion is applied (=1).
    integer, private :: nbot        ! Bottom level to which vertical diffusion is applied (=pver).
    integer, private :: ntop_eddy   ! Top level to which eddy vertical diffusion is applied.
    integer, private :: nbot_eddy   ! Bottom level to which eddy vertical diffusion is applied (=pver).
    integer, private :: ntop_molec  ! Top level to which molecular vertical diffusion is applied (=1).
    integer, private :: nbot_molec  ! Bottom level to which molecular vertical diffusion is applied.
    integer, private :: kvh_idx      ! sxj

    character(len=8), private :: vdiffnam(pcnst+pnats) ! names of v-diff tendencies

contains

    subroutine vd_register
        !-----------------------------------------------------------------------
        ! Register physics buffer fields and constituents
        !-----------------------------------------------------------------------
        use phys_buffer, only: pbuf_times, pbuf_add

        implicit none

        call pbuf_add('kvh', 'global', 1, pverp, pbuf_times, kvh_idx)  

    end subroutine vd_register

    subroutine vd_inti(cpairx, cpwvx, gravx, rairx, zvirx, hypm, vkx)
        !-----------------------------------------------------------------------
        ! Initialization of time independent fields for vertical diffusion.
        ! Calls initialization routine for boundary layer scheme.
        !-----------------------------------------------------------------------
        use history,    only: addfld, add_default, phys_decomp
        use turbulence, only: trbinti
        !------------------------------Arguments--------------------------------
        real(r8), intent(in) :: cpairx                 ! specific heat of dry air
        real(r8), intent(in) :: cpwvx                  ! spec. heat of water vapor at const. pressure
        real(r8), intent(in) :: gravx                  ! acceleration due to gravity
        real(r8), intent(in) :: rairx                  ! gas constant for dry air
        real(r8), intent(in) :: zvirx                  ! rh2o/rair - 1
        real(r8), intent(in) :: hypm(pver)             ! reference pressures at midpoints
        real(r8), intent(in) :: vkx                    ! Von Karman's constant

        !---------------------------Local workspace-----------------------------
        integer :: k                                   ! vertical loop index

        !-----------------------------------------------------------------------
        ! Set physical constants for vertical diffusion and pbl:
        cpair  = cpairx
        gravit = gravx
        rair   = rairx
        zvir   = zvirx

        ! Range of levels to which v-diff is applied.
        ntop_eddy  = 1       ! no reason not to make this 1, if >1, must be <= nbot_molec
        nbot_eddy  = pver    ! should always be pver
        ntop_molec = 1       ! should always be 1
        nbot_molec = 0       ! should be set below about 70 km

        ! Molecular diffusion turned on above ~60 km (50 Pa) if model top is above ~90 km (.1 Pa).
        ! Note that computing molecular diffusivities is a trivial expense, but constituent
        ! diffusivities depend on their molecular weights. Decomposing the diffusion matric
        ! for each constituent is a needless expense unless the diffusivity is significant.

        if (hypm(1) .lt. 0.1) then
            write(6,*) 'hypm(1)=',hypm(1)
            do k = 1, pver
                if (hypm(k) .lt. 50.) nbot_molec  = k
            end do
            if (masterproc) then
                write (6,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =',ntop_molec, 'NBOT_MOLEC =',nbot_molec
                write (6,fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =',ntop_eddy,  'NBOT_EDDY  =',nbot_eddy
            endif
        end if

        ! The vertical diffusion solver must operate over the full range of molecular and eddy diffusion
        ntop = 1
        nbot = pver

        ! Molecular weight factor in constitutent diffusivity
        ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****
        do k = 1, pcnst+pnats
            mw_fac(k) = d0 * mw_dry * sqrt(1./mw_dry + 1./mw_dry) / n_avog
        end do

        ! Initialize turbulence variables
        call trbinti(gravit, cpair, rair, zvir, ntop_eddy, nbot_eddy,   hypm, vkx)
        ! Set names of diffused variable tendencies and declare them as history variables
        do k = 1, pcnst+pnats
            vdiffnam(k) = 'VD'//cnst_name(k)
            if(k==1)vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
            call addfld (vdiffnam(k),'kg/kg/s ',pver, 'A','Vertical diffusion of '//cnst_name(k),phys_decomp)
        end do

        return
    end subroutine vd_inti

    subroutine vd_intr(                                    &
        ztodt    ,state    ,                               &
        taux     ,tauy     ,shflx    ,cflx     ,pblh     , &
        tpert    ,qpert    ,we       ,ehf      ,eqf      , &
        euf      ,evf      ,esf      ,ustar    ,obklen   ,ptend    , &
        khfs_tot ,kqfs_tot ,cldn) 
        !-----------------------------------------------------------------------
        ! interface routine for vertical diffusion and pbl scheme
        !-----------------------------------------------------------------------
        use physics_types,  only: physics_state, physics_ptend
        use history,        only: outfld
        use phys_buffer,    only: pbuf
        !------------------------------Arguments--------------------------------
        real(r8), intent(in) :: taux(pcols)            ! x surface stress (N/m2)
        real(r8), intent(in) :: tauy(pcols)            ! y surface stress (N/m2)
        real(r8), intent(in) :: shflx(pcols)           ! surface sensible heat flux (w/m2)
        real(r8), intent(in) :: cflx(pcols,pcnst+pnats)! surface constituent flux (kg/m2/s)
        real(r8), intent(in) :: ztodt                  ! 2 delta-t
        real(r8), intent(in) :: cldn(pcols,pver)       ! cloud fraction [fraction]
        type(physics_state), intent(in)  :: state      ! Physics state variables
        type(physics_ptend), intent(inout)  :: ptend   ! indivdual parameterization tendencies

        real(r8), intent(out) :: pblh(pcols)           ! planetary boundary layer height
        real(r8), intent(out) :: tpert(pcols)          ! convective temperature excess
        real(r8), intent(out) :: qpert(pcols,pcnst+pnats) ! convective humidity and constituent excess
        real(r8), intent(out) :: we(pcols)                ! entrainment velocity
        real(r8), intent(out) :: ehf(pcols)               ! entrainment heat flux
        real(r8), intent(out) :: esf(pcols)               ! entrainment dse flux
        real(r8), intent(out) :: euf(pcols)               ! entrainment U flux
        real(r8), intent(out) :: evf(pcols)               ! entrainment V flux
        real(r8), intent(out) :: eqf(pcols,pcnst+pnats)   ! entrainment constituent flux
        real(r8), intent(out) :: ustar(pcols)          ! surface friction velocity
        real(r8), intent(out) :: obklen(pcols)         ! Obukhov length
        real(r8), intent(out) :: khfs_tot(pcols,pverp)
        real(r8), intent(out) :: kqfs_tot(pcols,pverp)
        !
        !---------------------------Local storage-------------------------------
        integer :: lchnk                               ! chunk identifier
        integer :: ncol                                ! number of atmospheric columns
        integer :: i,k,m                               ! longitude,level,constituent indices

        real(r8) :: dtk(pcols,pver)                    ! T tendency from KE dissipation
        real(r8) :: kvh(pcols,pverp)                   ! diffusion coefficient for heat
        real(r8) :: kvm(pcols,pverp)                   ! diffusion coefficient for momentum
        real(r8) :: cgs(pcols,pverp)                   ! counter-gradient star (cg/flux)
        real(r8) :: rztodt                             ! 1./ztodt
        
        real(r8) :: rhoi(pcols,pverp)
        real(r8) :: khfs_1(pcols)
        real(r8) :: kqfs_1(pcols)
        
        real(r8) :: coun_h(pcols,pverp)
        real(r8) :: coun_q(pcols,pverp)
        real(r8) :: entr_h(pcols,pverp)
        real(r8) :: entr_q(pcols,pverp)
        real(r8) :: entr_u(pcols,pverp)
        real(r8) :: entr_v(pcols,pverp)
        !-----------------------------------------------------------------------
        ! local constants
        rztodt = 1./ztodt

        lchnk = state%lchnk
        ncol  = state%ncol

        ! Operate on copies of the input states, convert to tendencies at end.
        ! The input values of u,v are required for computing the kinetic energy dissipation.

        ptend%q(:ncol,:,:) = state%q(:ncol,:,:)
        ptend%s(:ncol,:) = state%s(:ncol,:)
        ptend%u(:ncol,:) = state%u(:ncol,:)
        ptend%v(:ncol,:) = state%v(:ncol,:)

        ! All variables are modified by vertical diffusion

        ptend%name  = "vertical diffusion"
        ptend%lq(:) = .TRUE.
        ptend%ls    = .TRUE.
        ptend%lu    = .TRUE.
        ptend%lv    = .TRUE.

        ! Call the real vertical diffusion code.
        call vdiff(                                                           &
            lchnk       ,ncol        ,                                        &
            ptend%u     ,ptend%v     ,state%t     ,ptend%q     ,ptend%s     , &
            state%pmid  ,state%pint  ,state%lnpint,state%lnpmid,state%pdel  , &
            state%rpdel ,state%zm    ,state%zi    ,ztodt       ,taux        , &
            tauy        ,shflx       ,cflx        ,pblh        ,ustar       , &
            kvh         ,kvm         ,tpert       ,qpert(1,1)  ,we          , &
            ehf         ,eqf         ,euf         ,evf         ,esf ,cgs         , &
            obklen      ,state%exner , &
            dtk         ,coun_h,coun_q,entr_h,entr_q,entr_u,entr_v,cldn)

        ! Convert the new profiles into vertical diffusion tendencies.
        ! Convert KE dissipative heat change into "temperature" tendency.
        
        if (RK_or_MG == 'MG') then
            pbuf(kvh_idx)%fld_ptr(1,1:ncol,1:pverp,lchnk,1) = kvh(:ncol,:pverp)
        end if

        do k = 1, pver
            do i = 1, ncol
                ptend%s(i,k) = (ptend%s(i,k) - state%s(i,k)) * rztodt
                ptend%u(i,k) = (ptend%u(i,k) - state%u(i,k)) * rztodt
                ptend%v(i,k) = (ptend%v(i,k) - state%v(i,k)) * rztodt
            end do
            do m = 1, pcnst+pnats
                do i = 1, ncol
                    ptend%q(i,k,m) = (ptend%q(i,k,m) - state%q(i,k,m))*rztodt
                end do
            end do
        end do
#ifdef PERGRO
        ! For pergro case, do not diffuse cloud water: replace with input values
        ptend%lq(ixcldw) = .FALSE.
        ptend%q(:ncol,:,ixcldw) = 0.
#endif
        do k = 2, pver
           do i = 1, ncol
              rhoi(i,k) = state%pint(i,k) * 2. / (rair*(state%t(i,k) + state%t(i,k-1)))
           end do
        end do
        do i = 1, ncol
           rhoi(i,pverp) = state%pmid(i,pver) / (rair*state%t(i,pver))
        end do

        do i = 1, ncol
           khfs_1(i) = shflx(i)/(rhoi(i,pverp)*cpair)
           kqfs_1(i) = cflx(i,1)/rhoi(i,pverp)
        end do
        
        khfs_tot = 0.
        kqfs_tot = 0.
        do i = 1, ncol
           khfs_tot(i,pverp) = khfs_1(i)
           kqfs_tot(i,pverp) = kqfs_1(i)
        end do   
        do k = 2, pver
           do i = 1,ncol
              khfs_tot(i,k) = rhoi(i,k+1)*khfs_tot(i,k+1)/rhoi(i,k) - (ptend%s(i,k)/(state%rpdel(i,k)*cpair*gravit*rhoi(i,k)))
              kqfs_tot(i,k) = rhoi(i,k+1)*kqfs_tot(i,k+1)/rhoi(i,k) - (ptend%q(i,k,1)/(state%rpdel(i,k)*gravit*rhoi(i,k)))
           end do
        end do
        ! Save the vertical diffusion variables on the history file
        dtk(:ncol,:) = dtk(:ncol,:)/cpair                ! normalize heating for history
        call outfld('DTVKE   ', dtk,            pcols, lchnk)
        dtk(:ncol,:) = ptend%s(:ncol,:)/cpair            ! normalize heating for history using dtk
        call outfld('DTV     ', dtk,            pcols, lchnk)
        call outfld('DQV     ', ptend%q(:,:,1), pcols, lchnk)
        call outfld('DUV     ', ptend%u,        pcols, lchnk)
        call outfld('DVV     ', ptend%v,        pcols, lchnk)
        call outfld('khfs_tot', khfs_tot,       pcols, lchnk)
        call outfld('kqfs_tot', kqfs_tot,       pcols, lchnk)
        do m = 1, pcnst+pnats
            call outfld(vdiffnam(m),ptend%q(1,1,m),pcols,lchnk)
        end do

        return
    end subroutine vd_intr

    subroutine vdiff (lchnk      ,ncol     ,                         &
        u          ,v          ,t          ,q          ,dse        , &
        pmid       ,pint       ,piln       ,pmln       ,pdel       , &
        rpdel      ,zm         ,zi         ,ztodt      ,taux       , &
        tauy       ,shflx      ,cflx       ,pblh       ,ustar      , &
        kvh        ,kvm        ,tpert      ,qpert      ,we         , &
        ehf        ,eqf        ,euf        ,evf        ,esf   ,cgs        , &
        obklen     ,exner      , &
        dtk        ,coun_h,coun_q,entr_h,entr_q,entr_u,entr_v,cldn)
        !-----------------------------------------------------------------------
        ! Driver routine to compute vertical diffusion of momentum, moisture, trace 
        ! constituents and dry static energy. The new temperature is computed from
        ! the diffused dry static energy.

        ! Turbulent diffusivities and boundary layer nonlocal transport terms are 
        ! obtained from the turbulence module.
        !-----------------------------------------------------------------------
        use turbulence,   only: trbintr
        use history,      only: outfld
        use phys_grid   !ljli2009
        use commap,       only: clat
        !------------------------------Arguments--------------------------------
        integer, intent(in) :: lchnk                   ! chunk identifier
        integer, intent(in) :: ncol                    ! number of atmospheric columns

        real(r8), intent(in) :: exner(pcols,pver)      ! (ps/pmid)**cappa
        real(r8), intent(in) :: pmid(pcols,pver)       ! midpoint pressures
        real(r8), intent(in) :: pint(pcols,pverp)      ! interface pressures
        real(r8), intent(in) :: piln (pcols,pverp)     ! Log interface pressures
        real(r8), intent(in) :: pmln (pcols,pver)      ! Log midpoint pressures
        real(r8), intent(in) :: pdel(pcols,pver)       ! thickness between interfaces
        real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel
        real(r8), intent(in) :: t(pcols,pver)          ! temperature input
        real(r8), intent(in) :: ztodt                  ! 2 delta-t
        real(r8), intent(in) :: taux(pcols)            ! x surface stress (N/m2)
        real(r8), intent(in) :: tauy(pcols)            ! y surface stress (N/m2)
        real(r8), intent(in) :: shflx(pcols)           ! surface sensible heat flux (W/m2)
        real(r8), intent(in) :: cflx(pcols,pcnst+pnats)! surface constituent flux (kg/m2/s)
        real(r8), intent(in) :: cldn(pcols,pver)       ! cloud fraction [fraction]

        real(r8), intent(inout) :: u(pcols,pver)       ! u wind
        real(r8), intent(inout) :: v(pcols,pver)       ! v wind
        real(r8), intent(inout) :: q(pcols,pver,pcnst+pnats)! moisture and trace constituents
        real(r8), intent(inout) :: dse(pcols,pver)     ! dry static energy [J/kg]
        real(r8), intent(in) :: zm(pcols,pver)         ! midpoint geoptl height above sfc
        real(r8), intent(in) :: zi(pcols,pverp)        ! interface geoptl height above sfc

        real(r8), intent(out) :: pblh(pcols)           ! planetary boundary layer height
        real(r8), intent(out) :: ustar(pcols)          ! surface friction velocity
        real(r8), intent(out) :: kvh(pcols,pverp)      ! diffusivity for heat
        real(r8), intent(out) :: kvm(pcols,pverp)      ! viscosity (diffusivity for momentum)
        real(r8), intent(out) :: tpert(pcols)          ! convective temperature excess
        real(r8), intent(out) :: qpert(pcols)          ! convective humidity excess
        real(r8), intent(out) :: we(pcols)             ! entrainment velocity
        real(r8), intent(out) :: ehf(pcols)            ! entrainment heat flux
        real(r8), intent(out) :: esf(pcols)            ! entrainment dse flux
        real(r8), intent(out) :: euf(pcols)            ! entrainment U flux
        real(r8), intent(out) :: evf(pcols)            ! entrainment V flux
        real(r8), intent(out) :: eqf(pcols,pcnst+pnats)! entrainment constituent flux
        real(r8), intent(out) :: cgs(pcols,pverp)      ! counter-grad star (cg/flux)
        real(r8), intent(out) :: obklen(pcols)         ! Obukhov length
        real(r8), intent(out) :: dtk(pcols,pver)       ! T tendency from KE dissipation
        real(r8), intent(out) :: coun_h(pcols,pverp)
        real(r8), intent(out) :: coun_q(pcols,pverp)
        real(r8), intent(out) :: entr_h(pcols,pverp)
        real(r8), intent(out) :: entr_q(pcols,pverp)
        real(r8), intent(out) :: entr_u(pcols,pverp)
        real(r8), intent(out) :: entr_v(pcols,pverp)
        !
        !---------------------------Local workspace-----------------------------
        real(r8) :: tmpm(pcols,pver)                   ! potential temperature, ze term in tri-diag sol'n

        real(r8) :: ca(pcols,pver)                     ! -upper diag of tri-diag matrix
        real(r8) :: cc(pcols,pver)                     ! -lower diag of tri-diag matrix
        real(r8) :: dnom(pcols,pver)                   ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))

        real(r8) :: kqfs(pcols,pcnst+pnats)            ! kinematic surf constituent flux (kg/m2/s)
        real(r8) :: kvq(pcols,pverp)                   ! diffusivity for constituents
        real(r8) :: kq_scal(pcols,pverp)               ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
        real(r8) :: mkvisc                             ! molecular kinematic viscosity c*tint**(2/3)/rho
        real(r8) :: tint                               ! interface temperature
        real(r8) :: tmp1(pcols)                        ! temporary storage
        real(r8) :: tmpi1(pcols,pverp)                 ! heat cg term [J/kg/m], or interface KE dissipation
        real(r8) :: tmpi2(pcols,pverp)                 ! rho or dt*(g*rho)**2/dp at interfaces

        real(r8) :: dinp_u(pcols,pverp)                ! vertical difference at interfaces, input u
        real(r8) :: dinp_v(pcols,pverp)                ! vertical difference at interfaces, input v
        real(r8) :: dout_u                             ! vertical difference at interfaces, output u
        real(r8) :: dout_v                             ! vertical difference at interfaces, output v
        real(r8) :: latitude
        integer :: i,k,m                               ! longitude,level,constituent indices
        integer :: lats(pcols)
        integer :: ktopbl(pcols)                       ! index of first midpoint inside pbl
        integer :: ktopblmn                            ! min value of ktopbl

        ! Compute the vertical differences of the input u,v for KE dissipation, interface k-
        ! Note, velocity=0 at surface, so then difference at the bottom interface is -u,v(pver)
        do i = 1, ncol
            dinp_u(i,1) = 0.
            dinp_v(i,1) = 0.
            dinp_u(i,pverp) = -u(i,pver)
            dinp_v(i,pverp) = -v(i,pver)
        end do
        do k = 2, pver
            do i = 1, ncol
                dinp_u(i,k) = u(i,k) - u(i,k-1)
                dinp_v(i,k) = v(i,k) - v(i,k-1)
            end do
        end do

        ! Compute potential temperature
        tmpm(:ncol,:pver) = t(:ncol,:pver) * exner(:ncol,:pver)

        ! Get diffusivities and c-g terms from turbulence module
        call trbintr(lchnk   ,ncol    ,                   &
            tmpm    ,t       ,q       ,zm      ,zi      , &
            pmid    ,u       ,v       ,taux    ,tauy    , &
            shflx   ,cflx    ,obklen  ,ustar   ,pblh    , &
            kvm     ,kvh     ,tmpi1   ,cgs     ,kqfs    , &
            tpert   ,qpert   ,we      ,ehf     ,eqf     , &
            euf     ,evf     ,ktopbl  ,ktopblmn,cldn)

        ! Compute rho at interfaces p/RT,  Ti = (Tm_k + Tm_k-1)/2,  interface k-
        do k = 2, pver
            do i = 1, ncol
                tmpi2(i,k) = pint(i,k) * 2. / (rair*(t(i,k) + t(i,k-1)))
            end do
        end do
        do i = 1, ncol
            tmpi2(i,pverp) = pint(i,pverp) / (rair*t(i,pver))
        end do

        ! Copy turbulent diffusivity for constituents from heat diffusivity

        kvq(:ncol,:) = kvh(:ncol,:)

        ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
        ! For the moment, keep molecular diffusivities all the same 

        kq_scal(:ncol,ntop_molec) = 0.
        do k = ntop_molec+1, nbot_molec
            do i = 1, ncol
                tint     = 0.5*(t(i,k) + t(i,k-1))
                mkvisc   = km_fac * tint**(2./3.) / tmpi2(i,k)
                kvm(i,k) = kvm(i,k) + mkvisc
                kvh(i,k) = kvh(i,k) + mkvisc*pr_num
                kq_scal(i,k) = sqrt(tint) / tmpi2(i,k)
            end do
        end do
        !ljli2009
         call get_lat_all_p(lchnk,ncol,lats)
         do k=1,pver
            do i=1,ncol
                latitude = clat(lats(i))
                kvm(i,k) = max(5.0*cos(latitude),kvm(i,k))
                kvh(i,k) = max(1.0*cos(latitude),kvh(i,k))
            end do
         end do
        !ljli2009

        ! Add the nonlocal transport terms to dry static energy, specific humidity and 
        ! other constituents in the boundary layer.

        call vd_add_cg(                                                   &
            lchnk      ,ncol       ,                                     &
            ktopblmn   ,q   ,dse        ,tmpi2      ,kvh        , &
            tmpi1      ,cgs       ,kqfs       ,rpdel      ,ztodt,coun_h,coun_q )
        
        ! Add the entrainment flux terms to dry static energy, specific humidity and other
        ! constituents in the boundary layer.

        call vd_add_ent(                                                 &
            lchnk     ,ncol         ,we         ,ktopbl                , &
            zm        ,ktopblmn     ,q          ,dse       ,u  ,v      , &
            tmpi2     ,ehf          ,eqf        , &
            euf       ,evf        ,pblh       ,rpdel       ,ztodt,esf,entr_h,entr_q,entr_u,entr_v )

        ! Add the (explicit) surface fluxes to the lowest layer

        do i = 1, ncol
            tmp1(i)  = ztodt*gravit*rpdel(i,pver)
            u(i,pver) = u(i,pver) + tmp1(i) * taux(i)
            v(i,pver) = v(i,pver) + tmp1(i) * tauy(i)
            dse(i,pver) = dse(i,pver) + tmp1(i) * shflx(i)
        end do
        do m = 1, pcnst+pnats
            do i = 1, ncol
                q(i,pver,m) = q(i,pver,m) + tmp1(i) * cflx(i,m)
            end do
        end do

        ! Calculate dt * (g*rho)^2/dp at interior interfaces,  interface k-
        do k = 2, pver
            do i = 1, ncol
                tmpi2(i,k) = ztodt * (gravit*tmpi2(i,k))**2 / (pmid(i,k) - pmid(i,k-1))
            end do
        end do

        ! Diffuse momentum
        call vd_lu_decomp(                                        &
            lchnk    ,ncol     ,                                 &
            kvm      ,tmpi2    ,rpdel    ,ztodt    ,ca         , &
            cc       ,dnom     ,tmpm     ,ntop     ,nbot       )
        call vd_lu_solve(                                                &
            lchnk    ,ncol     ,                                        &
            u        ,ca       ,tmpm     ,dnom     ,ntop     ,nbot     )
        call vd_lu_solve(                                                &
            lchnk    ,ncol     ,                                        &
            v        ,ca       ,tmpm     ,dnom     ,ntop     ,nbot     )

        ! Calculate kinetic energy dissipation
        ! 1. compute dissipation term at interfaces
        k = pverp
        do i = 1, ncol
            tmpi1(i,1) = 0.
            tmpi1(i,k) = 0.5 * ztodt * gravit * &
                ( (-u(i,k-1) + dinp_u(i,k))*taux(i) + (-v(i,k-1)+dinp_v(i,k))*tauy(i) )
        end do
        do k = 2, pver
            do i = 1, ncol
                dout_u = u(i,k) - u(i,k-1)
                dout_v = v(i,k) - v(i,k-1)
                tmpi1(i,k) = 0.25 * tmpi2(i,k) * kvm(i,k) *   &
                    (dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k))
            end do
        end do
        ! 2. Compute dissipation term at midpoints, add to dry static energy
        do k = 1, pver
            do i = 1, ncol
                dtk(i,k) = (tmpi1(i,k+1) + tmpi1(i,k)) * rpdel(i,k)
                dse(i,k) = dse(i,k) + dtk(i,k)
            end do
        end do

        ! Diffuse dry static energy
        call vd_lu_decomp(                                      &
            lchnk    ,ncol     ,                               &
            kvh      ,tmpi2    ,rpdel    ,ztodt    ,ca       , &
            cc       ,dnom     ,tmpm     ,ntop     ,nbot     )
        call vd_lu_solve(                                               &
            lchnk    ,ncol     ,                                       &
            dse      ,ca       ,tmpm     ,dnom     ,ntop     ,nbot    )

        ! Diffuse constituents.
        ! New decomposition required only if kvq != kvh (gw or molec diffusion)
        call vd_lu_decomp(                                                &
            lchnk      ,ncol       ,                                     &
            kvq        ,tmpi2      ,rpdel      ,ztodt      ,ca         , &
            cc         ,dnom       ,tmpm       ,ntop_eddy  ,nbot_eddy  )
        do m = 1, pcnst+pnats
            ! update decomposition in molecular diffusion range
            if (ntop_molec < nbot_molec) then
                call vd_lu_qmolec(                                                &
                    ncol       ,                                                 &
                    kvq        ,kq_scal    ,mw_fac(m)  ,tmpi2      ,rpdel      , &
                    ztodt      ,ca         ,cc         ,dnom       ,tmpm       , &
                    ntop_molec ,nbot_molec )
            end if
            call vd_lu_solve(                                              &
                lchnk   ,ncol     ,                                       &
                q(1,1,m),ca       ,tmpm     ,dnom     ,ntop    ,nbot     )
        end do
        call outfld ('coun_h',coun_h,pcols,lchnk)
        call outfld ('coun_q',coun_q,pcols,lchnk)
        call outfld ('entr_h',entr_h,pcols,lchnk)
        call outfld ('entr_q',entr_q,pcols,lchnk)
        call outfld ('entr_u',entr_u,pcols,lchnk)
        call outfld ('entr_v',entr_v,pcols,lchnk)
        call outfld ('ESF   ',esf/cpair,pcols,lchnk)
        return
    end subroutine vdiff

    subroutine vd_add_cg (                                &
        lchnk    ,ncol     ,                             &
        ktopblmn ,q        ,dse      ,rhoi    ,kvh     , &
        cgh      ,cgs      ,kqfs     ,rpdel   ,ztodt,coun_h,coun_q   )
        !-----------------------------------------------------------------------
        ! Add the "counter-gradient" term to the dry static energy and tracer fields
        ! in the boundary layer.
        ! Note, ktopblmn gives the minimum vertical index of the first midpoint
        ! within the boundary layer.

        ! This is really a boundary layer routine, but, since the pbl also supplies the 
        ! diffusivities in the boundary layer, there is not a clear divsion between
        ! pbl and eddy vertical diffusion code.
        !------------------------------Arguments--------------------------------

        
        integer, intent(in) :: lchnk                   ! chunk identifier
        integer, intent(in) :: ncol                    ! number of atmospheric columns
        integer, intent(in) :: ktopblmn                ! min value of ktopbl (highest bl top)

        real(r8), intent(in) :: rhoi(pcols,pverp)      ! density (p/RT) at interfaces
        real(r8), intent(in) :: kvh(pcols,pverp)       ! coefficient for heat and tracers
        real(r8), intent(in) :: cgh(pcols,pverp)       ! countergradient term for heat [J/kg/m]
        real(r8), intent(in) :: cgs(pcols,pverp)       ! unscaled cg term for constituents
        real(r8), intent(in) :: kqfs(pcols,pcnst+pnats)! kinematic surf constituent flux (kg/m2/s)
        real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel  (thickness bet interfaces)
        real(r8), intent(in) :: ztodt                  ! 2 delta-t

        real(r8), intent(inout) :: q(pcols,pver,pcnst+pnats)  ! moisture and trace constituent input
        real(r8), intent(inout) :: dse(pcols,pver)     ! dry static energy

        real(r8),intent(out) :: coun_h(pcols,pverp)
        real(r8),intent(out) :: coun_q(pcols,pverp)
        !---------------------------Local workspace-----------------------------

        real(r8) :: qtm(pcols,pver)                    ! temporary trace constituent input
        logical lqtst(pcols)                            ! adjust vertical profiles
        integer i                                      ! longitude index
        integer k                                      ! vertical index
        integer m                                      ! constituent index

        !-----------------------------------------------------------------------

        ! Add counter-gradient to input static energy profiles

        do k = ktopblmn-1, pver
            do i = 1, ncol
                dse(i,k) = dse(i,k)+ztodt*rpdel(i,k)*gravit &
                    *(rhoi(i,k+1)*kvh(i,k+1)*cgh(i,k+1)     &
                    - rhoi(i,k  )*kvh(i,k  )*cgh(i,k))
            end do
        end do

        ! Add counter-gradient to input mixing ratio profiles.
        ! Check for neg q's in each constituent and put the original vertical
        ! profile back if a neg value is found. A neg value implies that the
        ! quasi-equilibrium conditions assumed for the countergradient term are
        ! strongly violated.

        do m = 1, pcnst+pnats
            qtm(:ncol,ktopblmn-1:pver) = q(:ncol,ktopblmn-1:pver,m)
            do k = ktopblmn-1, pver
                do i = 1, ncol
                    q(i,k,m) = q(i,k,m) + ztodt*rpdel(i,k)*gravit*kqfs(i,m)     &
                        *(rhoi(i,k+1)*kvh(i,k+1)*cgs(i,k+1)                    &
                        - rhoi(i,k  )*kvh(i,k  )*cgs(i,k  ))
                end do
            end do
            lqtst(:ncol) = all(q(:ncol,ktopblmn-1:pver,m) >= qmincg(m), 2)
            do k = ktopblmn-1, pver
                q(:ncol,k,m) = merge (q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol))
            end do
        end do
        do k = 1, pverp
           do i = 1,ncol
              coun_h(i,k) = kvh(i,k)*cgh(i,k)/cpair
              coun_q(i,k) = kvh(i,k)*cgs(i,k)*kqfs(i,1)
           end do
        end do
        return
    end subroutine vd_add_cg

    subroutine vd_add_ent(                                               &
            lchnk     ,ncol         ,we         ,ktopbl                , &
            z         ,ktopblmn     ,q          ,dse       ,u    ,v    , &
            rhoi      ,ehf          ,eqf      , &
            euf       ,evf        ,pblh       ,rpdel       ,ztodt, esf, entr_h,entr_q,entr_u,entr_v )
        !-----------------------------------------------------------------------
        ! Add the "Entrainment Flux" term to the dry static energy and tracer fields
        ! in the boundary layer.
        ! Note, ktopbl gives the vertical index of the first midpoint
        ! within the boundary layer.
        !------------------------------Arguments--------------------------------
        use history, only: outfld 
        
        integer, intent(in) :: lchnk                   ! chunk identifier
        integer, intent(in) :: ncol                    ! number of atmospheric columns
        real(r8),intent(in) :: we(pcols)
        integer, intent(in) :: ktopbl(pcols)           ! index of first midpoint inside pbl
        integer, intent(in) :: ktopblmn           ! the highest midpoint in bl 

        real(r8), intent(in) :: z(pcols,pver)          !  mid-point geoptl height above sfc
        real(r8), intent(in) :: rhoi(pcols,pverp)      ! density (p/RT) at interfaces
        real(r8), intent(in) :: ehf(pcols)             ! entrainment term for heat [J/kg/m]
        real(r8), intent(in) :: euf(pcols)             ! entrainment term for U
        real(r8), intent(in) :: evf(pcols)             ! entrainment term for V
        real(r8), intent(in) :: eqf(pcols,pcnst+pnats) ! entrainment term for constituents
        real(r8), intent(in) :: pblh(pcols)            ! PBL height
        real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel  (thickness bet interfaces)
        real(r8), intent(in) :: ztodt                  ! 2 delta-t

        real(r8), intent(inout) :: q(pcols,pver,pcnst+pnats)  ! moisture and trace constituent input
        real(r8), intent(inout) :: dse(pcols,pver)     ! dry static energy
        real(r8), intent(inout) :: u(pcols,pver)     ! Velocity U
        real(r8), intent(inout) :: v(pcols,pver)     ! Velocity V

        real(r8), intent(out) :: esf(pcols)          ! entrainment term for dse[J/kg/m]
        real(r8),intent(out) :: entr_h(pcols,pverp)
        real(r8),intent(out) :: entr_q(pcols,pverp)
        real(r8),intent(out) :: entr_u(pcols,pverp)
        real(r8),intent(out) :: entr_v(pcols,pverp)
        !---------------------------Local workspace-----------------------------

        real(r8) :: qtm(pcols,pver)                    ! temporary trace constituent input
        real(r8) :: ehfz(pcols,pverp)
        real(r8) :: eufz(pcols,pverp)
        real(r8) :: evfz(pcols,pverp)
        real(r8) :: eqfz(pcols,pverp,pcnst+pnats)
        real(r8) :: dtent(pcols,pver)                    ! Temperature tendency for entrainment
        real(r8) :: duent(pcols,pver)                    ! U tendency for entrainment
        real(r8) :: dvent(pcols,pver)                    ! V tendency for entrainment
        real(r8) :: dqent(pcols,pver,pcnst+pnats)        ! moisture tendency for entrainment
        real(r8) :: dqent1(pcols,pver,pcnst+pnats)       ! moisture tendency for entrainment(after adjustment)
        real(r8) :: zp(pcols)

        logical lqtst(pcols)                            ! adjust vertical profiles(across pbl)
        real(r8) zmzp
        real(r8) minlevel
        real(r8) delta_s
        integer i                                      ! longitude index
        integer k                                      ! vertical index
        integer m                                      ! constituent index
        
        !-----------------------------------------------------------------------
        do i = 1, ncol
           minlevel = ktopbl(i)
           delta_s = dse(i, minlevel-1) - dse(i, minlevel)
           esf(i) = we(i)*delta_s
        end do
        
        do i = 1, ncol
           do k = 1, pverp
              ehfz(i,k) = 0.0
              eufz(i,k) = 0.0
              evfz(i,k) = 0.0
              do m = 1, pcnst+pnats
                 eqfz(i,k,m) = 0.0
              end do
           end do
        end do
        
        do i = 1, ncol
           do k = ktopblmn, pver
              if(z(i,k)<=pblh(i)) then
                 zp(i) = z(i, k-1)
                 if(zp(i) > pblh(i)) zp(i) = pblh(i)
                 zmzp = (z(i,k)+zp(i))/2
                 ehfz(i,k) = esf(i)*(zmzp/pblh(i))**3
                 eufz(i,k) = euf(i)*(zmzp/pblh(i))**3
                 evfz(i,k) = evf(i)*(zmzp/pblh(i))**3
                ! do m = 1, pcnst+pnats
                !    eqfz(i,k,m) = eqf(i,m)*(zmzp/pblh(i))**3
                ! end do
              end if
           end do
        end do
        ! Add entrainment heat flux  to input static energy profiles
        do k = ktopblmn-1, pver
            do i = 1, ncol
                dtent(i,k) = rpdel(i,k)*gravit &
                             *(rhoi(i,k+1)*ehfz(i,k+1)     &
                               - rhoi(i,k)*ehfz(i,k))/cpair
                
                dse(i,k) = dse(i,k)+ztodt*rpdel(i,k)*gravit &
                    *(rhoi(i,k+1)*ehfz(i,k+1)     &
                    - rhoi(i,k  )*ehfz(i,k))
            end do
        end do

        ! Add entrainment velocity flux  to input velocity profiles
        do k = ktopblmn-1, pver
            do i = 1, ncol
                duent(i,k) = rpdel(i,k)*gravit &
                             *(rhoi(i,k+1)*eufz(i,k+1)     &
                               - rhoi(i,k)*eufz(i,k))

                u(i,k) = u(i,k)+ztodt*rpdel(i,k)*gravit &
                    *(rhoi(i,k+1)*eufz(i,k+1)     &
                    - rhoi(i,k  )*eufz(i,k))
            end do
        end do
        do k = ktopblmn-1, pver
            do i = 1, ncol
                dvent(i,k) = rpdel(i,k)*gravit &
                             *(rhoi(i,k+1)*evfz(i,k+1)     &
                               - rhoi(i,k  )*evfz(i,k))

                v(i,k) = v(i,k)+ztodt*rpdel(i,k)*gravit &
                    *(rhoi(i,k+1)*evfz(i,k+1)     &
                    - rhoi(i,k  )*evfz(i,k))
            end do
        end do
        
        ! Add entrainment constituent to input constituent profiles.
        ! Check for neg q's in each constituent and put the original vertical
        ! profile back if a neg value is found. 

        do m = 1, pcnst+pnats
            qtm(:ncol,ktopblmn-1:pver) = q(:ncol,ktopblmn-1:pver,m)
            do k = ktopblmn-1, pver
                do i = 1, ncol
                    dqent(i,k,m) = rpdel(i,k)*gravit*     &
                                  (rhoi(i,k+1)*eqfz(i,k+1,m)                    &
                                   - rhoi(i,k)*eqfz(i,k,m))

                    q(i,k,m) = q(i,k,m) + ztodt*rpdel(i,k)*gravit*     &
                         (rhoi(i,k+1)*eqfz(i,k+1,m)                    &
                        - rhoi(i,k  )*eqfz(i,k,m))
                end do
            end do
            lqtst(:ncol) = all(q(:ncol,ktopblmn-1:pver,m) >= qmincg(m), 2)
            do k = ktopblmn-1, pver
                q(:ncol,k,m) = merge (q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol))
            end do
            dqent1(:,:,m) = (q(:,:,m)-qtm)/ztodt
        end do
        call outfld('DTENT     ',dtent    ,pcols,lchnk)
        call outfld('DUENT     ',duent    ,pcols,lchnk)
        call outfld('DVENT     ',dvent    ,pcols,lchnk)
        call outfld('DQENT     ',dqent(:,:,1)    ,pcols,lchnk)
        call outfld('DQENT1    ',dqent1(:,:,1)   ,pcols,lchnk)
        do k = 1, pverp
           do i = 1,ncol
              entr_h(i,k) = ehfz(i,k)
              entr_q(i,k) = eqfz(i,k,1)
              entr_u(i,k) = eufz(i,k)
              entr_v(i,k) = evfz(i,k)
           end do
        end do
        return
    end subroutine vd_add_ent
    
    !==============================================================================
    subroutine vd_lu_decomp(                                          &
        lchnk      ,ncol       ,                                     &
        kv         ,tmpi       ,rpdel      ,ztodt      ,ca         , &
        cc         ,dnom       ,ze         ,ntop       ,             &
        nbot       )
        !-----------------------------------------------------------------------
        ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
        ! tridiagonal diffusion matrix. 
        ! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver.
        ! Also determine ze factor and denominator for ze and zf (see solver).
        !------------------------------Arguments--------------------------------
        integer, intent(in) :: lchnk                   ! chunk identifier
        integer, intent(in) :: ncol                    ! number of atmospheric columns
        integer, intent(in) :: ntop                    ! top level to operate on
        integer, intent(in) :: nbot                    ! bottom level to operate on

        real(r8), intent(in) :: kv(pcols,pverp)        ! vertical diffusion coefficients
        real(r8), intent(in) :: tmpi(pcols,pverp)      ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
        real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel  (thickness bet interfaces)
        real(r8), intent(in) :: ztodt                  ! 2 delta-t

        real(r8), intent(out) :: ca(pcols,pver)        ! -upper diagonal
        real(r8), intent(out) :: cc(pcols,pver)        ! -lower diagonal
        real(r8), intent(out) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
        ! 1./(b(k) - c(k)*e(k-1))
        real(r8), intent(out) :: ze(pcols,pver)        ! term in tri-diag. matrix system

        !---------------------------Local workspace-----------------------------
        integer :: i                                   ! longitude index
        integer :: k                                   ! vertical index

        !-----------------------------------------------------------------------
        ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
        ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
        ! a combination of ca and cc; they are not required by the solver.

        do k = nbot-1, ntop, -1
            do i = 1, ncol
                ca(i,k  ) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k  )
                cc(i,k+1) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k+1)
            end do
        end do

        ! The bottom element of the upper diagonal (ca) is zero (element not used).
        ! The subdiagonal (cc) is not needed in the solver.

        do i=1,ncol
            ca(i,nbot) = 0.
        end do

        ! Calculate e(k).  This term is 
        ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

        do i = 1, ncol
            dnom(i,nbot) = 1./(1. + cc(i,nbot))
            ze(i,nbot)   = cc(i,nbot)*dnom(i,nbot)
        end do
        do k = nbot-1, ntop+1, -1
            do i=1,ncol
                dnom(i,k) = 1./ (1. + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
                ze(i,k) = cc(i,k)*dnom(i,k)
            end do
        end do
        do i=1,ncol
            dnom(i,ntop) = 1./ (1. + ca(i,ntop) - ca(i,ntop)*ze(i,ntop+1))
        end do

        return
    end subroutine vd_lu_decomp

    !==============================================================================
    subroutine vd_lu_qmolec(                                          &
        ncol       ,                                                 &
        kv         ,kq_scal    ,mw_facm    ,tmpi       ,rpdel      , &
        ztodt      ,ca         ,cc         ,dnom       ,ze         , &
        ntop       ,nbot       )
        !-----------------------------------------------------------------------
        ! Add the molecular diffusivity to the turbulent and gw diffusivity for a 
        ! consitutent.
        ! Update the superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
        ! tridiagonal diffusion matrix, also ze and denominator.
        !------------------------------Arguments--------------------------------
        integer, intent(in)     :: ncol                 ! number of atmospheric columns
        integer, intent(in)     :: ntop                 ! top level to operate on
        integer, intent(in)     :: nbot                 ! bottom level to operate on

        real(r8), intent(in)    :: kv(pcols,pverp)      ! vertical diffusion coefficients
        real(r8), intent(in)    :: kq_scal(pcols,pverp) ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
        real(r8), intent(in)    :: mw_facm              ! sqrt(1/M_q + 1/M_d) for this constituent
        real(r8), intent(in)    :: tmpi(pcols,pverp)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
        real(r8), intent(in)    :: rpdel(pcols,pver)    ! 1./pdel  (thickness bet interfaces)
        real(r8), intent(in)    :: ztodt                ! 2 delta-t

        real(r8), intent(inout) :: ca(pcols,pver)       ! -upper diagonal
        real(r8), intent(inout) :: cc(pcols,pver)       ! -lower diagonal
        real(r8), intent(inout) :: dnom(pcols,pver)     ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
        ! 1./(b(k) - c(k)*e(k-1))
        real(r8), intent(inout) :: ze(pcols,pver)       ! term in tri-diag. matrix system

        !---------------------------Local workspace-----------------------------
        integer :: i                                    ! longitude index
        integer :: k                                    ! vertical index

        real(r8) :: kmq                                 ! molecular diffusivity for constituent

        !-----------------------------------------------------------------------
        ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
        ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
        ! a combination of ca and cc; they are not required by the solver.

        do k = nbot-1, ntop, -1
            do i = 1, ncol
                kmq = kq_scal(i,k+1) * mw_facm
                ca(i,k  ) = (kv(i,k+1)+kmq) * tmpi(i,k+1) * rpdel(i,k  )
                cc(i,k+1) = (kv(i,k+1)+kmq) * tmpi(i,k+1) * rpdel(i,k+1)
            end do
        end do

        ! Calculate e(k).  This term is 
        ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

        do k = nbot, ntop+1, -1
            do i=1,ncol
                dnom(i,k) = 1./ (1. + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
                ze(i,k) = cc(i,k)*dnom(i,k)
            end do
        end do
        do i=1,ncol
            dnom(i,ntop) = 1./ (1. + ca(i,ntop) - ca(i,ntop)*ze(i,ntop+1))
        end do

        return
    end subroutine vd_lu_qmolec

    !==============================================================================
    subroutine vd_lu_solve(                                           &
        lchnk      ,ncol       ,                                     &
        q          ,ca         ,ze         ,dnom       ,             &
        ntop       ,nbot       )
        !-----------------------------------------------------------------------
        ! Solve the implicit vertical diffusion equation with zero flux boundary conditions.
        ! Actual surface fluxes are explicit (rather than implicit) and are applied separately. 
        ! Procedure for solution of the implicit equation follows Richtmyer and 
        ! Morton (1967,pp 198-200).

        ! The equation solved is
        !
        ! -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k), 
        !
        ! where d(k) is the input profile and q(k) is the output profile
        !
        ! The solution has the form
        !
        ! q(k)  = ze(k)*q(k-1) + zf(k)
        !
        ! ze(k) = cc(k) * dnom(k)
        !
        ! zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)
        !
        ! dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]

        ! Note that the same routine is used for temperature, momentum and tracers,
        ! and that input variables are replaced.
        !------------------------------Arguments--------------------------------
        integer, intent(in) :: lchnk                   ! chunk identifier
        integer, intent(in) :: ncol                    ! number of atmospheric columns
        integer, intent(in) :: ntop                    ! top level to operate on
        integer, intent(in) :: nbot                    ! bottom level to operate on

        real(r8), intent(in) :: ca(pcols,pver)         ! -upper diag coeff.of tri-diag matrix
        real(r8), intent(in) :: ze(pcols,pver)         ! term in tri-diag solution
        real(r8), intent(in) :: dnom(pcols,pver)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))

        real(r8), intent(inout) :: q(pcols,pver)       ! constituent field

        !---------------------------Local workspace-----------------------------
        real(r8) :: zf(pcols,pver)                     ! term in tri-diag solution

        integer i,k                                    ! longitude, vertical indices

        !-----------------------------------------------------------------------

        ! Calculate zf(k).  Terms zf(k) and ze(k) are required in solution of 
        ! tridiagonal matrix defined by implicit diffusion eqn.
        ! Note that only levels ntop through nbot need be solved for.

        do i = 1, ncol
            zf(i,nbot) = q(i,nbot)*dnom(i,nbot)
        end do
        do k = nbot-1, ntop, -1
            do i=1,ncol
                zf(i,k) = (q(i,k) + ca(i,k)*zf(i,k+1))*dnom(i,k)
            end do
        end do

        ! Perform back substitution

        do i=1,ncol
            q(i,ntop) = zf(i,ntop)
        end do
        do k=ntop+1, nbot, +1
            do i=1,ncol
                q(i,k) = zf(i,k) + ze(i,k)*q(i,k-1)
            end do
        end do

        return
    end subroutine vd_lu_solve

end module vertical_diffusion
