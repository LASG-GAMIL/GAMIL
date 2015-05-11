#include<misc.h>   ! SPMD

!-----------------------------------------------------------------------
! Purpose:
!   MG cloud scheme
!   Provides the CAM interface to the prognostic cloud water and ice parameterization
! Author:
!   SHI Xiangjun 2008/11/08 add to gamil
!-----------------------------------------------------------------------

module MG

    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physconst,     only: gravit, latvap, latice
    use chemistry,     only: chem_is
    use phys_buffer,   only: pbuf_add, pbuf_times
    use constituents,  only: pcnst, cnst_add, nonadvec

    ! debug modules
    use pmgrid,        only: masterproc

    implicit none

    private

    save

    public stratiform_register
    public stratiform_init_cnst
    public stratiform_implements_cnst
    public stratiform_init
    public stratiform_tend

    ! TAG-1
    ! integrate the choosing MG scheme into namelist
#include<RK_or_MG.h> ! sxj-2008-11-13

    integer, parameter :: ncnstmax=4                      ! number of constituents
    character(8), parameter :: cnst_names(ncnstmax) = &   ! constituent names
        (/'CLDLIQ','CLDICE','NUMLIQ','NUMICE'/)
    integer ncnst                                         ! number of constituents (can vary)

    logical use_shfrc     ! local copy of flag from convect_shallow_use_shfrc

    integer ixcldliq      ! cloud liquid amount index
    integer ixcldice      ! cloud ice amount index
    integer ixnumliq      ! cloud liquid number index
    integer ixnumice      ! cloud ice water index
    integer qcwat_idx     ! qcwat index in physics buffer
    integer lcwat_idx     ! lcwat index in physics buffer
    integer tcwat_idx     ! tcwat index in physics buffer
    integer cld_idx       ! cld index in physics buffer
    integer concld_idx    ! concld index in physics buffer

    integer ndx_totevapr  ! pbuf index for nevapr+evapsnow
    integer ndx_totprecp  ! pbuf index for prain+prodsnow

contains

    ! register the constituents (cloud liquid and cloud ice) and the fields in the physics buffer

    subroutine  stratiform_register

        use constituents, only: cnst_add
        use physconst,    only: mwdry, cpair
        use phys_buffer,  only: pbuf_add
        use tracers,      only: ixcldw ! added by SHI Xiangjun

        integer idx

        ! TAG-1
        ! put the choosing of stratiform cloud scheme in namelist
        if ( RK_or_MG .eq. 'MG' ) then
            ncnst = 4
        else
            write(6,*) "RK_or_MG=",RK_or_MG,"--MG.F90"
            write(6,*) "RK_or_MG shoud be MG when call this subroutine"
            call endrun
        end if

        ! register cloud water and determine index
        call cnst_add(cnst_names(1), nonadvec, mwdry, cpair, 0._r8, ixcldliq, &
            longname='Grid box averaged cloud liquid amount')
        ! DONG Li: The index of cloud water is duplicate
        ixcldw = ixcldliq
        call cnst_add(cnst_names(2), nonadvec, mwdry, cpair, 0._r8, ixcldice, &
            longname='Grid box averaged cloud ice amount')
        call cnst_add(cnst_names(3), nonadvec, mwdry, cpair, 0._r8, ixnumliq, &
            longname='Grid box averaged cloud liquid number')
        call cnst_add(cnst_names(4), nonadvec, mwdry, cpair, 0._r8, ixnumice, &
            longname='Grid box averaged cloud ice number')
        ! request physics buffer space for fields that persist across timesteps.
        call pbuf_add('QCWAT',  'global', 1, pver, pbuf_times, qcwat_idx)
        call pbuf_add('LCWAT',  'global', 1, pver, pbuf_times, lcwat_idx)
        call pbuf_add('TCWAT',  'global', 1, pver, pbuf_times, tcwat_idx)
        call pbuf_add('CLD',    'global', 1, pver, pbuf_times, cld_idx)
        call pbuf_add('CONCLD', 'global', 1, pver, pbuf_times, concld_idx)
        call pbuf_add('QINI',   'physpkg', 1, pver, pcnst, idx)
        call pbuf_add('TINI',   'physpkg', 1, pver, 1, idx)
        call pbuf_add('QME',    'physpkg', 1, pver, 1, idx)
        call pbuf_add('PRAIN',  'physpkg', 1, pver, 1, idx)
        call pbuf_add('NEVAPR', 'physpkg', 1, pver, 1, idx)
        if (chem_is('trop_mozart')) then
            call pbuf_add('TOTEVAPR', 'physpkg', 1, pver, 1, ndx_totevapr)
            call pbuf_add('TOTPRECP', 'physpkg', 1, pver, 1, ndx_totprecp)
        end if
        call pbuf_add('REI',    'physpkg', 1, pver, 1, idx)
        call pbuf_add('REL',    'physpkg', 1, pver, 1, idx)
        call pbuf_add('REL_FN', 'physpkg', 1, pver, 1, idx)  ! REL at fixed number for indirect rad forcing

    end subroutine stratiform_register

    ! return true if specified constituent is implemented by this package

    function stratiform_implements_cnst(name)
        character(*), intent(in) :: name
        logical stratiform_implements_cnst
        integer m

        stratiform_implements_cnst = .false.

        do m = 1, ncnst
            if (name == cnst_names(m)) then
                stratiform_implements_cnst = .true.
                return
            end if
        end do

    end function stratiform_implements_cnst

    ! initialize the cloud water mixing ratios (liquid and ice), if they are
    ! not read from the initial file

    subroutine stratiform_init_cnst(name, q)
        character(*), intent(in)  :: name
        real(r8), intent(out) :: q(:,:,:) ! mass mixing ratio

        if ( name == 'CLDLIQ' ) then
            q = 0.0_r8
        else if ( name == 'CLDICE' ) then
            q = 0.0_r8
        else if ( name == 'NUMLIQ' ) then
            q = 0.0_r8
        else if ( name == 'NUMICE' ) then
            q = 0.0_r8
        end if

        return
    end subroutine stratiform_init_cnst

    ! initialize the cloud water parameterization

    subroutine stratiform_init
        use cldwat,        only: inimc
        use cldwat2m,      only: inimmc
        use constituents,  only: cnst_get_ind, cnst_name, cnst_longname
        use history,       only: addfld, add_default, phys_decomp
        use tracers,       only: sflxnam
        use physconst,     only: tmelt, rh2o, rhodair
        use dycore,        only: dycore_is

        integer :: m, mm
        ! TAG-2
        !character(len=16) :: deep_scheme
        ! initialization routine for prognostic cloud water
        !deep_scheme='ZM'  ! **

        if (RK_or_MG=='MG') then
            call inimmc
        else
            ! TAG-3
            write(6,*) "RK_or_MG=",RK_or_MG,"--MG.F90"
            write(6,*) "RK_or_MG shoud be MG when call this subroutine"
            call endrun
        end if

        ! find out whether shfrc from convect_shallow will be used in cldfrc
        ! shfrc!  cloud fraction from shallow convection scheme  ! gamil does not have this part so --sxj--
        ! TAG-2
        use_shfrc = .false.
        ! register history variables
        ! ---sxj--some variables addfld in history.F90-bldfld
        !do m = 1,ncnst
        !call cnst_get_ind(cnst_names(m), mm)
        !call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
        !call addfld (sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
        !call add_default (cnst_name(mm), 1, ' ')
        !call add_default (sflxnam(mm),   1, ' ')
        !end do
        !call addfld ('FWAUT   ','fraction',pver, 'A','Relative importance of liquid autoconversion' ,phys_decomp)
        !call addfld ('FSAUT   ','fraction',pver, 'A','Relative importance of ice autoconversion'    ,phys_decomp)
        !call addfld ('FRACW   ','fraction',pver, 'A','Relative  importance of rain accreting liquid',phys_decomp)
        !call addfld ('FSACW   ','fraction',pver, 'A','Relative  importance of snow accreting liquid',phys_decomp)
        !call addfld ('FSACI   ','fraction',pver, 'A','Relative  importance of snow accreting ice'   ,phys_decomp)
        !call addfld ('CME     ','kg/kg/s ',pver, 'A','Rate of cond-evap within the cloud'           ,phys_decomp)
        ! sxj add aer_act
        call addfld ('aer_act ','cm-3     ',pver, 'A','aerosol activate number                  '    ,phys_decomp)
        call add_default ('aer_act ', 1, ' ')
        call addfld ('CMEICE  ','kg/kg/s ',pver, 'A','Rate of cond-evap of ice within the cloud'    ,phys_decomp)
        call addfld ('CMELIQ  ','kg/kg/s ',pver, 'A','Rate of cond-evap of liq within the cloud'    ,phys_decomp)
        call addfld ('ICE2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of ice to precip'          ,phys_decomp)
        call addfld ('LIQ2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of liq to precip'          ,phys_decomp)
        call addfld ('ZMDLF   ','kg/kg/s ',pver, 'A','Detrained liquid water from ZM convection'    ,phys_decomp)
        call addfld ('PRODPREC','kg/kg/s ',pver, 'A','Rate of conversion of condensate to precip'   ,phys_decomp)
        call addfld ('EVAPPREC','kg/kg/s ',pver, 'A','Rate of evaporation of falling precip'        ,phys_decomp)
        call addfld ('EVAPSNOW','kg/kg/s ',pver, 'A','Rate of evaporation of falling snow'          ,phys_decomp)
        call addfld ('HPROGCLD','W/kg'    ,pver, 'A','Heating from prognostic clouds'               ,phys_decomp)
        call addfld ('HCME    ','W/kg'    ,pver, 'A','Heating from cond-evap within the cloud'      ,phys_decomp)
        call addfld ('HEVAP   ','W/kg'    ,pver, 'A','Heating from evaporation of falling precip'   ,phys_decomp)
        call addfld ('HFREEZ  ','W/kg'    ,pver, 'A','Heating rate due to freezing of precip'       ,phys_decomp)
        call addfld ('HMELT   ','W/kg'    ,pver, 'A','Heating from snow melt'                       ,phys_decomp)
        call addfld ('HREPART ','W/kg'    ,pver, 'A','Heating from cloud ice/liquid repartitioning' ,phys_decomp)
        call addfld ('REPARTICE','kg/kg/s',pver, 'A','Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
        call addfld ('REPARTLIQ','kg/kg/s',pver, 'A','Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)
        !call addfld ('FICE    ','fraction',pver, 'A','Fractional ice content within cloud'          ,phys_decomp)
        !call addfld ('ICWMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud water mixing ratio'       ,phys_decomp)
        !call addfld ('ICIMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud ice mixing ratio'         ,phys_decomp)
        call addfld ('PCSNOW  ','m/s'     ,1   , 'A','Snow fall from prognostic clouds'             ,phys_decomp)
        call addfld ('DQSED   ','kg/kg/s' ,pver, 'A','Water vapor tendency from cloud sedimentation',phys_decomp)
        call addfld ('DLSED   ','kg/kg/s' ,pver, 'A','Cloud liquid tendency from sedimentation'     ,phys_decomp)
        call addfld ('DISED   ','kg/kg/s' ,pver, 'A','Cloud ice tendency from sedimentation'        ,phys_decomp)
        call addfld ('HSED    ','W/kg'    ,pver, 'A','Heating from cloud sediment evaporation'      ,phys_decomp)
        call addfld ('SNOWSED ','m/s'     ,1   , 'A','Snow from cloud ice sedimentation'            ,phys_decomp)
        call addfld ('RAINSED ','m/s'     ,1   , 'A','Rain from cloud liquid sedimentation'         ,phys_decomp)
        call addfld ('PRECSED ','m/s'     ,1   , 'A','Precipitation from cloud sedimentation'       ,phys_decomp)
        !call add_default ('FICE    ', 1, ' ')
        !call addfld ('CNVCLD  ','fraction',1,    'A','Vertically integrated convective cloud amount',phys_decomp)
        !call addfld ('CLDST   ','fraction',pver, 'A','Stratus cloud fraction',phys_decomp)
        !call addfld ('CONCLD  ','fraction',pver, 'A','Convective cloud cover',phys_decomp)
        !call add_default ('CONCLD', 1, ' ')
        ! history variables for new microphysics
        call addfld ('MPDT   ','K/s     ',pver, 'A','T tendency - Morrison microphysics',phys_decomp)
        call addfld ('MPDQ   ','kg/kg/s     ',pver, 'A','Q tendency - Morrison microphysics',phys_decomp)
        call addfld ('ICWNC   ','m-3     ',pver, 'A','Prognostic in-cloud water number conc'       ,phys_decomp)
        call addfld ('ICINC   ','m-3     ',pver, 'A','Prognostic in-cloud ice number conc'         ,phys_decomp)
        call addfld ('EFFLIQ  ','Micron  ',pver, 'A','Prognostic droplet effective radius'       ,phys_decomp)
        call addfld ('EFFLIQ_IND','Micron  ',pver, 'A','Prognostic droplet effective radius (indirect effect)'       ,phys_decomp)
        call addfld ('EFFICE  ','Micron  ',pver, 'A','Prognostic ice effeictive radius'         ,phys_decomp)
        call addfld ('WSUB  ','m/s     ',pver, 'A','Diagnostic sub-grid vertical velocity'         ,phys_decomp)
        call addfld ('CDNUMC  ','#/m2    ',1, 'A','Vertically-integrated droplet concentration'         ,phys_decomp)
        call add_default('ICWNC   ', 1, ' ')  !!sxj
        call add_default('ICINC   ', 1, ' ')  !!sxj
        call add_default('EFFLIQ  ', 1, ' ')  !!sxj
        call add_default('EFFICE  ', 1, ' ')  !!sxj
        call add_default('CDNUMC  ', 1, ' ')  !!sxj
        call add_default('WSUB    ', 1, ' ')  !!sxj
        call addfld ('CCN1    ','#/cm3   ',pver, 'A','CCN concentration at S=0.02%',phys_decomp)
        call addfld ('CCN2    ','#/cm3   ',pver, 'A','CCN concentration at S=0.05%',phys_decomp)
        call addfld ('CCN3    ','#/cm3   ',pver, 'A','CCN concentration at S=0.1%',phys_decomp)
        call add_default('CCN3    ', 1, ' ')  !!sxj
        call addfld ('CCN4    ','#/cm3   ',pver, 'A','CCN concentration at S=0.2%',phys_decomp)
        call addfld ('CCN5    ','#/cm3   ',pver, 'A','CCN concentration at S=0.5%',phys_decomp)
        call addfld ('CCN6    ','#/cm3   ',pver, 'A','CCN concentration at S=1.0%',phys_decomp)
        ! diagnostic precip
        call addfld ('QRAIN   ','kg/kg   ',pver, 'A','Diagnostic grid-mean rain mixing ratio'         ,phys_decomp)
        call addfld ('QSNOW   ','kg/kg   ',pver, 'A','Diagnostic grid-mean snow mixing ratio'         ,phys_decomp)
        call addfld ('NRAIN   ','m-3     ',pver, 'A','Diagnostic grid-mean rain number conc'         ,phys_decomp)
        call addfld ('NSNOW   ','m-3     ',pver, 'A','Diagnostic grid-mean snow number conc'         ,phys_decomp)
        ! averaging for cloud particle number and size
        call addfld ('AWNC   ','m-3     ',pver, 'A','Average cloud water number conc'         ,phys_decomp)
        call addfld ('AWNI   ','m-3     ',pver, 'A','Average cloud ice number conc'         ,phys_decomp)
        call addfld ('AREL  ','Micron  ',pver, 'A','Average droplet effective radius'       ,phys_decomp)
        call addfld ('AREI  ','Micron  ',pver, 'A','Average ice effective radius'       ,phys_decomp)
        call add_default('AWNC   ', 1, ' ')  !!sxj
        call add_default('AWNI   ', 1, ' ')  !!sxj
        call add_default('AREL   ', 1, ' ')  !!sxj
        call add_default('AREI   ', 1, ' ')  !!sxj
        !sxj
        call addfld ('aer_data  ','  kg    ',pver, 'A','read in aerosol data for test '       ,phys_decomp)
        call add_default('aer_data  ', 1, ' ')  !!sxj
        ! frequency arrays for above
        call addfld ('FREQL  ','fraction  ',pver, 'A','Fractional occurance of liquid'       ,phys_decomp)
        call addfld ('FREQI  ','fraction  ',pver, 'A','Fractional occurance of ice'       ,phys_decomp)
        ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
        call addfld ('AQRAIN   ','kg/kg   ',pver, 'A','Average rain mixing ratio'         ,phys_decomp)
        call addfld ('AQSNOW   ','kg/kg   ',pver, 'A','Average snow mixing ratio'         ,phys_decomp)
        call addfld ('ANRAIN   ','m-3     ',pver, 'A','Average rain number conc'         ,phys_decomp)
        call addfld ('ANSNOW   ','m-3     ',pver, 'A','Average snow number conc'         ,phys_decomp)
        call addfld ('ADRAIN   ','Micron  ',pver, 'A','Average rain effective Diameter'         ,phys_decomp)
        call addfld ('ADSNOW   ','Micron  ',pver, 'A','Average snow effective Diameter'         ,phys_decomp)
        call addfld ('FREQR  ','fraction  ',pver, 'A','Fractional occurance of rain'       ,phys_decomp)
        call addfld ('FREQS  ','fraction  ',pver, 'A','Fractional occurance of snow'       ,phys_decomp)
        ! Average cloud top particle size and number (liq, ice) and frequency
        call addfld ('ACTREL  ','Micron  ',1, 'A','Average Cloud Top droplet effective radius'       ,phys_decomp)
        call addfld ('ACTREI  ','Micron  ',1, 'A','Average Cloud Top ice effective radius'       ,phys_decomp)
        call addfld ('ACTNL  ','1/m3  ',1, 'A','Average Cloud Top droplet number'       ,phys_decomp)
        call addfld ('ACTNI  ','1/m3  ',1, 'A','Average Cloud Top ice number'       ,phys_decomp)
        call addfld ('FCTL  ','fraction  ',1, 'A','Fractional occurance of cloud top liquid'       ,phys_decomp)
        call addfld ('FCTI  ','fraction  ',1, 'A','Fractional occurance of cloud top ice'       ,phys_decomp)
        call add_default('ACTREL  ', 1, ' ')  !!sxj
        call add_default('ACTREI  ', 1, ' ')  !!sxj
        call add_default('ACTNL   ', 1, ' ')  !!sxj
        call add_default('ACTNI   ', 1, ' ')  !!sxj

        if (masterproc) then
            write(6, "('Notice: MG::stratiform_init: finished')")
        end if

        return
    end subroutine stratiform_init

    ! Interface to sedimentation, detrain, cloud fraction and  microphysics subroutines
    ! DONG Li: sst is added by SHI Xiangjun, is it duplicate with ts?
    subroutine stratiform_tend(state, ptend_all, dtime, icefrac, landfrac, ocnfrac, &
                               snowh, dlf, cmfmc, cmfmc2, ts, sst, zdu, &
                               prec_str, snow_str, prec_sed, snow_sed, prec_pcw, snow_pcw, &
                               qtend, ttend, ltend, fice, pbuf)
        use shr_kind_mod,     only: r8 => shr_kind_r8
        use ppgrid
        use physics_types,    only: physics_state, physics_ptend, physics_tend
        use physics_types,    only: physics_ptend_init, physics_update, physics_tend_init
        use physics_types,    only: physics_ptend_sum,  physics_state_copy
        use history,          only: outfld
        use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
        use constituents,     only: cnst_get_ind
        use cldwat,           only: pcond, cldwat_fice
        use cldwat2m,         only: mmicro_pcond
        use aerosol_index,    only: naer_all,naer
        use aerosol_mass_interface, only: collect_aer_masses, aer_mass
        use time_manager,     only: is_first_step
        use cloud_fraction_MG,   only: cldfrc_MG

        implicit none

        real(r8), parameter :: pnot = 1.0e5_r8          ! reference pressure

        type(physics_state), intent(in)  :: state       ! state variables
        type(physics_ptend), intent(out) :: ptend_all   ! package tendencies
        type(pbuf_fld), intent(inout) :: pbuf(pbuf_size_max)
        real(r8), intent(in) :: dtime                   ! timestep
        real(r8), intent(in) :: icefrac (pcols)         ! sea ice fraction (fraction)                   ! DONG Li: These variables are from "comsrf",
        real(r8), intent(in) :: landfrac(pcols)         ! land fraction (fraction)                      !          can we use "comsrf" directly?
        real(r8), intent(in) :: ocnfrac (pcols)         ! ocean fraction (fraction)                     !
        real(r8), intent(in) :: snowh(pcols)            ! Snow depth over land, water equivalent (m)    !
        real(r8), intent(in) :: dlf(pcols,pver)         ! detrained water from ZM
        real(r8), intent(in) :: cmfmc(pcols,pverp)      ! convective mass flux--m sub c
        real(r8), intent(in) :: cmfmc2(pcols,pver)      ! shallow convective mass flux--m sub c
        real(r8), intent(in) :: ts(pcols)               ! surface temperature
        real(r8), intent(in) :: sst(pcols)

        real(r8), intent(in) :: qtend(pcols,pver)       !
        real(r8), intent(in) :: ttend(pcols,pver)       !  made as input arguments by SHI Xiangjun 2009/03/19
        real(r8), intent(in) :: ltend(pcols,pver)       !

        real(r8), intent(in) :: zdu(pcols,pver)         ! detrainment rate from deep convection
        real(r8), intent(out) :: fice(pcols,pver)       ! fractional ice content within cloud, added by SHI Xiangjun
        real(r8), intent(out) :: prec_str(pcols)        ! [Total] sfc flux of precip from stratiform (m/s)
        real(r8), intent(out) :: snow_str(pcols)        ! [Total] sfc flux of snow from stratiform   (m/s)
        real(r8), intent(out) :: prec_sed(pcols)        ! surface flux of total cloud water from sedimentation
        real(r8), intent(out) :: snow_sed(pcols)        ! surface flux of cloud ice from sedimentation
        real(r8), intent(out) :: prec_pcw(pcols)        ! sfc flux of precip from microphysics(m/s)
        real(r8), intent(out) :: snow_pcw(pcols)        ! sfc flux of snow from microphysics (m/s)

        type(physics_state) state1                      ! local copy of the state variable
        type(physics_tend ) tend                        ! physics tendencies (empty, needed for physics_update call)
        type(physics_ptend) ptend_loc                   ! package tendencies
        integer i, k, m
        integer lchnk                  ! chunk identifier
        integer ncol                   ! number of atmospheric columns

        ! physics buffer fields
        integer itim, ifld
        real(r8), pointer :: cld(:,:)    ! cloud fraction
        real(r8), pointer :: concld(:,:) ! convective cloud fraction
        real(r8), pointer :: qme(:,:)
        real(r8), pointer :: prain(:,:)
        real(r8), pointer :: nevapr(:,:)
        real(r8), pointer :: rel(:,:)     ! liquid effective drop radius (microns)
        real(r8), pointer :: rei(:,:)     ! ice effective drop size (microns)
        real(r8) shfrc(pcols,pver)                ! cloud fraction from shallow convection scheme
        real(r8), pointer :: rel_fn(:,:)  ! ice effective drop size at fixed number (indirect effect) (microns)
        real(r8), pointer :: kkvh(:,:)    ! heat flux for cldwat  !
        !! local variables for stratiform_sediment
        real(r8) :: rain(pcols)                       ! surface flux of cloud liquid
        real(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
        real(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)
        !! local variables for cldfrc
        real(r8)  cldst(pcols,pver)     ! cloud fraction
        real(r8)  clc(pcols)            ! column convective cloud amount
        real(r8) rhdfda(pcols,pver)     ! d_RH/d_cloud_fraction    ====wlin
        real(r8) rhu00(pcols,pver)      ! RH limit, U00             ====wlin
        real(r8) relhum(pcols,pver)       ! RH, output to determine drh/da
        real(r8) rhu002(pcols,pver)       ! same as rhu00 but for perturbed rh
        real(r8) cld2(pcols,pver)         ! same as cld but for perturbed rh
        real(r8) concld2(pcols,pver)      ! same as concld but for perturbed rh
        real(r8) cldst2(pcols,pver)       ! same as cldst but for perturbed rh
        real(r8) relhum2(pcols,pver)      ! RH after  perturbation
        real(r8) :: pmid(pcols,pver)      ! midpoint pressures
        real(r8) :: t(pcols,pver)         ! temperature
        real(r8) :: q(pcols,pver)         ! specific humidity
        real(r8) :: omga(pcols,pver)      ! vertical pressure velocity
        real(r8) :: phis(pcols)           ! surface geopotential
        real(r8) :: pdel(pcols,pver)      ! pressure depth of layer
        real(r8) :: ps(pcols)             ! surface pressure
        !! local variables for microphysics
        real(r8) :: rdtime                        ! 1./dtime
        ! real(r8) :: qtend (pcols,pver)            ! moisture tendencies
        !real(r8) :: ttend (pcols,pver)            ! temperature tendencies
        !real(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
        real(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
        real(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
        real(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
        real(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
        real(r8) :: cmeheat (pcols,pver)          ! heating rate due to phase change of precip
        real(r8) :: prodsnow(pcols,pver)          ! local production of snow
        real(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
        !real(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
        real(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
        real(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
        real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
        real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
        real(r8) fwaut(pcols,pver)
        real(r8) fsaut(pcols,pver)
        real(r8) fracw(pcols,pver)
        real(r8) fsacw(pcols,pver)
        real(r8) fsaci(pcols,pver)
        real(r8) cmeice(pcols,pver)   ! Rate of cond-evap of ice within the cloud
        real(r8) cmeliq(pcols,pver)   ! Rate of cond-evap of liq within the cloud
        real(r8) ice2pr(pcols,pver)   ! rate of conversion of ice to precip
        real(r8) liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
        real(r8) liq2snow(pcols,pver) ! rate of conversion of liquid to snow
        real(r8) temp(pcols)
        real(r8) res(pcols,pver)
        !! variables for morrison microphysics
        real(r8) :: dum1,dum2
        real(r8) :: qc(pcols,pver)
        real(r8) :: qi(pcols,pver)
        real(r8) :: nc(pcols,pver)
        real(r8) :: ni(pcols,pver)
        real(r8) :: icinc(pcols,pver)             ! in cloud ice number conc
        real(r8) :: cdnumc(pcols)                 ! vertically-integrated droplet concentration
        real(r8) :: icwnc(pcols,pver)             ! in cloud water number conc
        real(r8) :: effliq(pcols,pver)            ! in cloud liq eff rad
        real(r8) :: effice(pcols,pver)            ! in cloud ice eff rad
        real(r8) :: effliq_fn(pcols,pver)         ! in cloud liq eff rad at fixed number concentration
        real(r8) :: wsub(pcols,pver)              ! sub-grid vertical velocity (m/s)
        !! output from mmicro_pcond
        real(r8) :: tlat(pcols,pver)
        real(r8) :: qvlat(pcols,pver)
        real(r8) :: qcten(pcols,pver)
        real(r8) :: qiten(pcols,pver)
        real(r8) :: ncten(pcols,pver)
        real(r8) :: niten(pcols,pver)
        real(r8) :: effc(pcols,pver)
        real(r8) :: effc_fn(pcols,pver)     ! liquid effective radius at fixed number (for indirect calc)
        real(r8) :: effi(pcols,pver)
        real(r8) :: prect(pcols)
        real(r8) :: preci(pcols)
        !! averaging arrays for effective radius and number....
        real(r8) :: efiout(pcols,pver)
        real(r8) :: efcout(pcols,pver)
        real(r8) :: ncout(pcols,pver)
        real(r8) :: niout(pcols,pver)
        real(r8) :: freqi(pcols,pver)
        real(r8) :: freql(pcols,pver)
        !! average cloud top radius & number
        real(r8) :: ctrel(pcols)
        real(r8) :: ctrei(pcols)
        real(r8) :: ctnl(pcols)
        real(r8) :: ctni(pcols)
        real(r8) :: fcti(pcols)
        real(r8) :: fctl(pcols)
        real(r8) mass_to_mmr(pcols,pver) ! conversion of layer mass to mass mixing ratio
        real(r8) aer_mmr(pcols,pver,naer_all) ! aerosol mass mixing ratio
        !! WATER TRACERS
        integer itrip            ! counter of water tracer triplets
        integer ixwtvap, ixwtliq, ixwtice   ! constituent indicies for vap, liq, ice of triplet
        real(r8) dum(pcols,pver)        ! dummy argument
        !real(r8) :: wtprec_sed(pcols,pwspc) ! surface flux of total cloud water from sedimentation
        !real(r8) :: wtsnow_sed(pcols,pwspc) ! surface flux of cloud ice from sedimentation
        !real(r8) :: wtrain_sed(pcols,pwspc) ! surface flux of cloud liquid

        if (RK_or_MG.ne.'MG') then
            write(6,*) "RK_or_MG=",RK_or_MG,"--MG.F90_stratiform_tend"
            write(6,*) "RK_or_MG shoud be MG when call this subroutine"
            call endrun
        end if
        !
        lchnk = state%lchnk
        ncol  = state%ncol
        call physics_state_copy(state,state1)   ! copy state to local state1.
        call physics_ptend_init(ptend_loc)  ! initialize local ptend type
        call physics_ptend_init(ptend_all)  ! initialize output ptend type
        call physics_tend_init(tend)        ! tend here is just a null place holder
        !! Associate pointers with physics buffer fields
        itim = pbuf_old_tim_idx()
        if (itim>1) then
            write(*,*) "itim always be 1 now gamil",itim
            call endrun
        endif
        !ifld = pbuf_get_fld_idx('QCWAT')
        !qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
        !ifld = pbuf_get_fld_idx('TCWAT')
        !tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
        !ifld = pbuf_get_fld_idx('LCWAT')
        !lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
        ifld = pbuf_get_fld_idx('CLD')
        cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
        ifld = pbuf_get_fld_idx('CONCLD')
        concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
        ifld = pbuf_get_fld_idx('QME')
        qme  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('PRAIN')
        prain  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('NEVAPR')
        nevapr  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('REL')
        rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('REI')
        rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('REL_FN')
        rel_fn  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
        ifld = pbuf_get_fld_idx('KVH')
        kkvh => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
        !! if first timestep, initialize heatflux....in pbuf at all time levels.
        if (is_first_step()) then
            kkvh(:,:)= 0.0_r8
        end if
        !
        ! ++detrain ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Put the detraining cloud water from convection into the cloud and environment.
        ptend_loc%name  = 'pcwdetrain'
        ptend_loc%lq(ixcldliq) = .TRUE.
        ptend_loc%lq(ixcldice) = .TRUE.
        ptend_loc%lq(ixnumliq) = .TRUE.
        ptend_loc%lq(ixnumice) = .TRUE.
        ptend_loc%ls           = .TRUE.
        ! Put all of the detraining cloud water from convection into the large scale cloud.
        ! put detraining cloud water into liq and ice based on temperature partition
        do k = 1,pver
            do i = 1,state1%ncol
                if (state1%t(i,k) > 268.15_r8) then
                    dum1=0.0_r8
                else if (state1%t(i,k) < 238.15_r8) then
                    dum1=1.0_r8
                    if (state1%t(i,k) < 50._r8) then
                        write(*,*) "state1%t(i,k) < 150._r8"
                        !                call endrun
                    endif
                else
                    dum1=(268.15_r8-state1%t(i,k))/30._r8
                end if
                ptend_loc%q(i,k,ixcldliq) = dlf(i,k)*(1._r8-dum1)
                ptend_loc%q(i,k,ixcldice) = dlf(i,k)*dum1
                ! calculate detrainment of Nc  ! 'dif'------detrain water ---sxj---
                ! assume detrained cloud water has mean volume radius of 8 micron
                dum2 = dlf(i,k)*(1._r8-dum1)
                ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8)*0.125  !20101114ljli
                dum2 = dlf(i,k)*dum1
                if (state1%t(i,k) < 233.15_r8) then !!! maybe wrong   233.15_r8 should accord with the T 238.15 above ---sxj---
                    ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)
                endif
                ! account for latent heat release during freezing
                ptend_loc%s(i,k) = dlf(i,k)*dum1*latice
            end do
        end do
        call outfld('ZMDLF' ,dlf  , pcols,state1%lchnk)
        ! add tendency from this process to tend from other processes here
        call physics_ptend_sum(ptend_loc, ptend_all, state)
        call physics_update(state1, tend, ptend_loc, dtime)
        call physics_ptend_init(ptend_loc)
        !
        !++++ cldfrc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! cloud fraction after transport and convection,
        ! derive the relationship between rh and cld from
        ! the employed cloud scheme
        pmid(:ncol,:pver) = state1%pmid(:ncol,:pver)
        t(:ncol,:pver) = state1%t(:ncol,:pver)
        q(:ncol,:pver) = state1%q(:ncol,:pver,1)
        omga(:ncol,:pver) = state1%omega(:ncol,:pver)
        pdel(:ncol,:pver) =  state1%pdel(:ncol,:pver)
        ps(:ncol) = state1%pint(:ncol,pverp)
        phis(:ncol) = state1%phis(:ncol)
        call t_startf("cldfrc_MG")
        call cldfrc_MG(state1, cmfmc, cmfmc2, zdu,&
                       landfrac, ocnfrac,         &
                       ts, sst, ps, snowh,        &
                       cld, clc, concld, cldst,   &
                       relhum, rhu00, 0)
        ! in order to calculate rhdfda so call cldfrc_MG again with index=1
        call cldfrc_MG(state1, cmfmc, cmfmc2, zdu,&
                       landfrac, ocnfrac,         &
                       ts, sst, ps, snowh,        &
                       cld2, clc, concld2, cldst2,&
                       relhum2, rhu002, 1)
        call t_stopf("cldfrc_MG")
        ! cldfrc does not define layer cloud for model layer at k=1
        ! so set rhu00(k=1)=2.0 to not calculate cme for this layer
        rhu00(:ncol,1) = 2.0_r8
        ! Add following to estimate rhdfda
        do k=1,pver
            do i=1,ncol
                if (relhum(i,k) < rhu00(i,k) ) then
                    rhdfda(i,k)=0.0_r8
                else if (relhum(i,k) >= 1.0_r8 ) then
                    rhdfda(i,k)=0.0_r8
                else
                    if((cld2(i,k) - cld(i,k) ) < 1.e-4_r8 ) then
                        rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8   !instead of 0.0
                    else
                        rhdfda(i,k)=0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
                    endif
                endif
            enddo
        enddo

        !+ mp +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! cloud water and ice parameterizations
        call t_startf('stratiform_microphys')
        ! Initialize chunk id and size
        lchnk = state1%lchnk
        ncol  = state1%ncol
        rdtime = 1._r8/dtime
        ! Define fractional amount of cloud condensate in ice phase
        call cldwat_fice(ncol, state1%t, fice, fsnow)
        !tmax_fice  = tmelt - 10._r8       ! max temperature for cloud ice formation
        !tmin_fice  = tmax_fice - 30._r8   ! min temperature for cloud ice formation
        !tmax_fsnow = tmelt                ! max temperature for transition to convective snow
        !tmin_fsnow = tmelt-5._r8          ! min temperature for transition to convective snow
        ptend_loc%name         = 'cldwat_MG'
        ptend_loc%ls           = .true.
        ptend_loc%lq(1)        = .true.
        ptend_loc%lq(ixcldice) = .true.
        ptend_loc%lq(ixcldliq) = .true.
        ptend_loc%lq(ixnumliq) = .true.
        ptend_loc%lq(ixnumice) = .true.

        qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
        qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
        nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
        ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)

        !write(*,*) ixcldliq,ixcldice,ixnumliq,ixnumice
        !stop 22

        !qtend(:ncol,:pver) = (state1%q(:ncol,:pver,1) - qcwat(:ncol,:pver))*rdtime
        !ttend(:ncol,:pver) = (state1%t(:ncol,:pver)   - tcwat(:ncol,:pver))*rdtime
        !ltend(:ncol,:pver) = (state1%q(:ncol,:pver,ixcldliq)+state1%q(:ncol,:pver,ixcldice) - &
        !      lcwat(:ncol,:pver))*rdtime


        !load aer_mass in aerosol_radiation_interface.(CAM3.5.8)
        ! gamil aer_mass in aerosol_mass_interface
        call collect_aer_masses(state1)                                     !!!sxj--2009-0309
        !   aer_mass(pcols, pver, naer_all, begchunk:endchunk)   read from MATCH

        call outfld('aer_data',aer_mass(1:pcols,1:pver,1,lchnk), pcols,lchnk)    !!sxj
        !      call outfld('AREL',efcout,    pcols,lchnk)
        mass_to_mmr(:ncol,:)=gravit*state1%rpdel(:ncol,:)
        !do m=1,naer_all
        do m=1,naer_all
            aer_mmr(:ncol,:,m)=aer_mass(:ncol,:,m,lchnk)*mass_to_mmr(:ncol,:)
        enddo
        !   aer_mmr kg/kg

        !write(6,*)  aer_mmr(ncol,5,5)
        !write(*,*) "MG.F90--line631"
        !call endrun

        call t_startf('mmicro_pcond')

        call mmicro_pcond (lchnk,   ncol, dtime,&
            state1%t  , ttend      , state1%q(1,1,1), qtend   , ltend,&
                                !qc,qi,nc,ni, state1%pmid , state1%pdel , cld, aer_mmr, rhdfda, rhu00, fice, &
                                !qc,qi,nc,ni, state1%pmid , state1%pdel , cld,           rhdfda, rhu00, fice, & !!!sxj--20081113
            qc,qi,nc,ni, state1%pmid , state1%pdel , cld, aer_mmr, rhdfda, rhu00, fice, &  !!sxj--20090311 take aerosol into account
            tlat,qvlat,qcten,qiten,ncten,niten,effc,effc_fn,effi, &
            prect,preci,kkvh,wsub,nevapr,evapsnow,prain,prodsnow,qme,landfrac,snowh)

        call t_stopf('mmicro_pcond')

        !-----------
        if (chem_is('trop_mozart')) then
            pbuf(ndx_totevapr)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = nevapr(:ncol,:pver) + evapsnow(:ncol,:pver)
            pbuf(ndx_totprecp)%fld_ptr(1,1:ncol,1:pver,lchnk,1) = prain(:ncol,:pver) + prodsnow(:ncol,:pver)
        endif

        do i = 1,ncol
            do k = 1,pver
                ptend_loc%s(i,k)          =tlat(i,k)
                ptend_loc%q(i,k,1)        =qvlat(i,k)
                ptend_loc%q(i,k,ixcldliq) =qcten(i,k)
                ptend_loc%q(i,k,ixcldice) =qiten(i,k)
                ptend_loc%q(i,k,ixnumliq) =ncten(i,k)
                ptend_loc%q(i,k,ixnumice) =niten(i,k)
                !write(*,'(6f13.5)') tlat(i,k),qvlat(i,k),qcten(i,k),qiten(i,k),ncten(i,k),niten(i,k)
            end do
        end do

        prec_pcw(:ncol) = prect(:ncol)
        snow_pcw(:ncol) = preci(:ncol)
        prec_sed(:ncol) =0._r8
        snow_sed(:ncol) =0._r8
        !prec_str(:ncol) = prec_pcw(:ncol)+prec_sed(:ncol)-rliq(:ncol)  !!!!  sxj ????????--- *******%%%%%%%%%
        prec_str(:ncol) = prec_pcw(:ncol)+prec_sed(:ncol)
        snow_str(:ncol) = snow_pcw(:ncol)+snow_sed(:ncol)


        call physics_ptend_sum(ptend_loc,ptend_all, state)     !!!
        !write(*,*) "MG.F90 ptend_sum before****",lchnk
        !call endrun



        ! Set the name of the final package tendencies. Note that this
        ! is a special case in physics_update, so a change here must be
        ! matched there.
        ptend_all%name = 'stratiform'
        ! used below
        call physics_update (state1, tend, ptend_loc, dtime)



        call physics_ptend_init(ptend_loc)

        !      accumulate prec and snow
        ! Save off q and t after cloud water
        do k=1,pver
            !        qcwat(:ncol,k) = state1%q(:ncol,k,1)
            !        tcwat(:ncol,k) = state1%t(:ncol,k)
            !        lcwat(:ncol,k) = state1%q(:ncol,k,ixcldliq) + state1%q(:ncol,k,ixcldice)
            rel(:ncol,k) = effc(:ncol,k)
            rel_fn(:ncol,k) = effc_fn(:ncol,k)
            rei(:ncol,k) = effi(:ncol,k)
        end do

        ! Compute in cloud ice and liquid mixing ratios (output only)
        do k=1,pver
            do i = 1,ncol
                icimr(i,k) = (state1%q(i,k,ixcldice) + dtime*ptend_loc%q(i,k,ixcldice)) /&
                    max(0.01_r8,cld(i,k))
                icwmr(i,k) = (state1%q(i,k,ixcldliq) + dtime*ptend_loc%q(i,k,ixcldliq)) /&
                    max(0.01_r8,cld(i,k))
                icinc(i,k) = (state1%q(i,k,ixnumice) + dtime*ptend_loc%q(i,k,ixnumice)) /&
                    max(0.01_r8,cld(i,k))*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
                icwnc(i,k) = (state1%q(i,k,ixnumliq) + dtime*ptend_loc%q(i,k,ixnumliq)) /&
                    max(0.01_r8,cld(i,k))*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
                effliq(i,k) = effc(i,k)
                effliq_fn(i,k) = effc_fn(i,k)
                effice(i,k) = effi(i,k)
            end do
        end do

        !--------------------- OUTPUT FIELDS FOR HISTORY ---------------------

        ! Column droplet concentration
        do i = 1,ncol
            cdnumc(i)=0._r8
            do k=1,pver
                cdnumc(i)=cdnumc(i)+state1%q(i,k,ixnumliq)*state1%pdel(i,k)/9.816_r8
            end do
        end do


        ! Averaging for new output fields
        efcout(:,:)=0._r8
        efiout(:,:)=0._r8
        ncout(:,:)=0._r8
        niout(:,:)=0._r8
        freql(:,:)=0._r8
        freqi(:,:)=0._r8
        do i = 1,ncol
            do k=1,pver
                !if (cld(i,k).gt.0.01.and.icwmr(i,k).gt.1.e-4_r8) then  !sxj ---?????
                if (cld(i,k).gt.0.01.and.icwmr(i,k).gt.1.e-6_r8) then  !sxj ---?????
                    efcout(i,k)=effc(i,k)
                    ncout(i,k)=icwnc(i,k)
                    freql(i,k)=1._r8
                endif
                if (cld(i,k).gt.0.01.and.icimr(i,k).gt.1.e-6_r8) then
                    efiout(i,k)=effi(i,k)
                    niout(i,k)=icinc(i,k)
                    freqi(i,k)=1._r8
                endif
            end do
        end do

        !Add averaged out fields.
        call outfld('AREL',efcout,    pcols,lchnk)
        call outfld('AREI',efiout,    pcols,lchnk)
        call outfld('AWNC',ncout,    pcols,lchnk)
        call outfld('AWNI',niout,    pcols,lchnk)
        call outfld('FREQL',freql,    pcols,lchnk)
        call outfld('FREQI',freqi,    pcols,lchnk)

        !Cloud top effective radius and number....
        fcti(:)=0._r8
        fctl(:)=0._r8
        ctrel(:)=0._r8
        ctrei(:)=0._r8
        ctnl(:)=0._r8
        ctni(:)=0._r8
        do i = 1,ncol
            do k=1,pver
                if (cld(i,k).gt.0.01.and.icwmr(i,k).gt.1.e-7_r8) then
                    ctrel(i)=effc(i,k)
                    ctnl(i)=icwnc(i,k)
                    fctl(i)=1._r8
                    exit
                endif

                if (cld(i,k).gt.0.01.and.icimr(i,k).gt.1.e-7_r8) then
                    ctrei(i)=effi(i,k)
                    ctni(i)=icinc(i,k)
                    fcti(i)=1._r8
                    exit
                endif

            enddo
        enddo

        call outfld('ACTREL',ctrel,    pcols,lchnk)
        call outfld('ACTREI',ctrei,    pcols,lchnk)
        call outfld('ACTNL',ctnl,    pcols,lchnk)
        call outfld('ACTNI',ctni,    pcols,lchnk)
        call outfld('FCTL',fctl,    pcols,lchnk)
        call outfld('FCTI',fcti,    pcols,lchnk)

        !write(*,*) "MG.F90--789"
        !call endrun

        ! microphysics variables to output fields
        call outfld('MPDT',tlat,    pcols,lchnk)
        call outfld('MPDQ',qvlat,    pcols,lchnk)
        call outfld('ICINC',icinc,    pcols,lchnk)
        call outfld('ICWNC',icwnc,    pcols,lchnk)
        call outfld('EFFLIQ',effliq,    pcols,lchnk)
        call outfld('EFFLIQ_IND',effliq_fn,    pcols,lchnk)
        call outfld('EFFICE',effice,    pcols,lchnk)
        call outfld('WSUB',wsub,    pcols,lchnk)
        call outfld('CDNUMC',cdnumc,    pcols,lchnk)

        call outfld('ICIMR',icimr,    pcols,lchnk)
        call outfld('ICWMR',icwmr,    pcols,lchnk)
        call outfld('CME'  ,qme     , pcols,lchnk)
        call outfld('PRODPREC',prain, pcols,lchnk)
        call outfld('EVAPPREC',nevapr, pcols,lchnk)
        !   call outfld('EVAPSNOW',evapsnow, pcols,lchnk)

        call t_stopf('stratiform_microphys')

        call t_startf("cldfrc_MG")
        call cldfrc_MG(state1, cmfmc, cmfmc2, zdu, &
                       landfrac, ocnfrac,         &
                       ts, sst, ps, snowh,        &
                       cld, clc, concld, cldst,   &
                       relhum, rhu00, 0)
        call t_stopf("cldfrc_MG")

        call outfld('CONCLD  ', concld, pcols, lchnk)
        call outfld('CLDST   ', cldst,  pcols, lchnk)
        call outfld('CNVCLD  ', clc,    pcols, lchnk)

    endsubroutine stratiform_tend

end module MG
