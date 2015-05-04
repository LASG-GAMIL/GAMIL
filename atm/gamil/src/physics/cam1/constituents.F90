#include <misc.h>
#include <params.h>

module constituents

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Contains data and functions for manipulating advected and non-advected constituents
    !
    ! Public functions/subroutines:
    !   initindx, add_cnst
    ! 
    ! Author: B.A. Boville
    ! 
    !-----------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst,    only: r_universal

    implicit none

    private

    save

    public cnst_add      ! add a constituent to the list of advected (or nonadvected) constituents
    public cnst_get_ind  ! get the index of a constituent
    public cnst_chk_dim  ! check that number of constituents added equals dimensions (pcnst, pnats)

    integer, parameter, public :: pcnst  = PCNST       ! number of advected constituents (including water vapor)
    integer, parameter, public :: pnats  = PNATS       ! number of non-advected constituents
    integer, parameter, public :: ppcnst = pcnst+pnats ! total number of constituents

    integer, parameter, public :: advected = 0         ! type value for constituents which are advected
    integer, parameter, public :: nonadvec = 1         ! type value for constituents which are not advected

    character(8), public :: cnst_name(ppcnst)          ! constituent names (including any non-advected)
    character(128), public :: cnst_longname(ppcnst)    ! long name of constituents

    ! Constants for each tracer
    real(r8), public :: cnst_cp  (ppcnst)             ! specific heat at constant pressure (J/kg/K)
    real(r8), public :: cnst_cv  (ppcnst)             ! specific heat at constant volume (J/kg/K)
    real(r8), public :: cnst_mw  (ppcnst)             ! molecular weight (kg/kmole)
    real(r8), public :: cnst_rgas(ppcnst)             ! gas constant ()
    real(r8), public :: qmin     (ppcnst)             ! minimum permitted constituent concentration (kg/kg)
    real(r8), public :: qmincg   (ppcnst)             ! for backward compatibility only

    ! Volume mixing ratios for trace gases at the Earth's surface
    real(r8), public :: co2vmr                        ! co2   volume mixing ratio 
    real(r8), public :: n2ovmr                        ! n2o   volume mixing ratio 
    real(r8), public :: ch4vmr                        ! ch4   volume mixing ratio 
    real(r8), public :: f11vmr                        ! cfc11 volume mixing ratio 
    real(r8), public :: f12vmr                        ! cfc12 volume mixing ratio 

    integer, private :: padv = 0                      ! index pointer to last advected tracer
    integer, private :: pnad = pcnst                  ! index pointer to last non-advected tracer

CONTAINS

    subroutine cnst_add(name, type, mwc, cpc, qminc, ind, longname)
        use pmgrid, only: iam
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: Register a constituent for treament by physical parameterizations and 
        !          transport (if type=advected)
        ! 
        ! Method: 
        ! <Describe the algorithm(s) used in the routine.> 
        ! <Also include any applicable external references.> 
        ! 
        ! Author:  B.A. Boville
        ! 
        !-----------------------------Arguments---------------------------------
        !
        character(*), intent(in) :: name
        character(*), intent(in), optional :: longname

        integer,  intent(in) :: type   ! flag indicating advected or nonadvected
        real(r8), intent(in) :: mwc    ! constituent molecular weight (kg/kmol)
        real(r8), intent(in) :: cpc    ! constituent specific heat at constant pressure (J/kg/K)
        real(r8), intent(in) :: qminc  ! minimum value of mass mixing ratio (kg/kg)
                                       ! normally 0., except water 1.E-12, for radiation.

        integer, intent(out) :: ind    ! global constituent index (in q array)

        if (type == advected) then
            ! set tracer index and check validity, advected tracer
            padv = padv+1
            ind  = padv
            if (padv > pcnst) then
                write(6, "('Error: constitutes::cnst_add: advected tracer index greater than pcnst ', i)") pcnst
                call endrun
            end if
            if (iam .eq. 0) then
                write(6, "('Notice: constitutes::cnst_add: add advected tracer ', a, ' at index ', i)") name, ind
            end if
        else if (type == nonadvec) then
            ! set tracer index and check validity, non-advected tracer
            pnad = pnad+1
            ind  = pnad
            if (pnad > ppcnst) then
                write(6, "('Error: constitutes::cnst_add: non-advected tracer index greater than pcnst+pnats ', i)") ppcnst
                call endrun
            end if
            if (iam .eq. 0) then
                write(6, "('Notice: constitutes::cnst_add: add non-advected tracer ', a, ' at index ', i)") name, ind
            end if
        else
            ! unrecognized type value
            write(6, "('Error: constitutes::cnst_add: input type flag ""', a, '"" invalid')") type
            call endrun
        end if

        ! set tracer name and constants
        cnst_name(ind) = name
        if (present(longname)) then
            cnst_longname(ind) = longname
        else
            cnst_longname(ind) = name
        end if

        cnst_cp(ind) = cpc
        cnst_mw(ind) = mwc
        qmin   (ind) = qminc
        qmincg (ind) = qminc
        if (ind == 1) qmincg = 0.  ! This crap is replicate what was there before ****

        cnst_rgas(ind) = r_universal*mwc
        cnst_cv  (ind) = cpc-cnst_rgas(ind)

        return
    end subroutine cnst_add

    subroutine cnst_get_ind (name, ind)
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: Get the index of a constituent 
        ! 
        ! Method: 
        ! <Describe the algorithm(s) used in the routine.> 
        ! <Also include any applicable external references.> 
        ! 
        ! Author:  B.A. Boville
        ! 
        !-----------------------------Arguments---------------------------------
        !
        character(len=*), intent(in) :: name ! constituent name

        integer, intent(out)   :: ind    ! global constituent index (in q array)

        !---------------------------Local workspace-----------------------------
        integer :: m                                   ! tracer index

        !-----------------------------------------------------------------------

        ! Find tracer name in list
        do m = 1, ppcnst
            if (name == cnst_name(m)) then
                ind  = m
                return
            end if
        end do

        ! Unrecognized name
        write(6,*) 'CNST_GET_IND, name:', name,  ' not found in list:', cnst_name(:)
        call endrun

    end subroutine cnst_get_ind

    !==============================================================================
    subroutine cnst_chk_dim
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: Check that the number of registered constituents of each type is the
        !          same as the dimension
        ! 
        ! Method: 
        ! <Describe the algorithm(s) used in the routine.> 
        ! <Also include any applicable external references.> 
        ! 
        ! Author:  B.A. Boville
        ! 
        !-----------------------------------------------------------------------
        !
        if (padv /= pcnst) then
            write(6,*)'CNST_CHK_DIM: number of advected tracer ',padv, ' not equal to pcnst = ',pcnst
            call endrun ()
        endif
        if (pnad /= ppcnst) then
            write(6,*)'CNST_CHK_DIM: number of non-advected tracers ',pnad, ' not equal to pcnst+pnats = ', &
                ppcnst
            call endrun ()
        endif

    end subroutine cnst_chk_dim
end module constituents
