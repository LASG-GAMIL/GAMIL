#include <misc.h>
#include <params.h>

subroutine camice(srf_state,srfflx)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! CAM sea ice surface fluxes.
!
! Method: 
! 
! Author:
! 
!-----------------------------------------------------------------------
!
! $Id: camice.F90,v 1.1.4.3 2002/06/15 13:50:09 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use pspect
  use ice_data
  use comsrf, only: surface_state,srfflx_parm,icefrac,snowhice,sicthk, &
	tsice,asdirice,asdifice,aldirice,aldifice
  use phys_grid, only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use time_manager, only: get_nstep, get_step_size, get_curr_calday
  use ice_dh, only:prognostic_icesnow
  implicit none
!
! Input/Output arguments
!
   type(surface_state), intent(inout), dimension(begchunk:endchunk) :: srf_state
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx

#include <comctl.h>

!---------------------------Local variables-----------------------------
  integer :: nstep          ! current timestep number
  integer :: dtime          ! timestep size [seconds]
  real(r8) rtime            ! calendar day for next timestep
  real(r8) lats(pcols)            ! 
  real(r8) lons(pcols)            ! 
  real(r8) cdaynext         ! calendar day for next timestep
  real(r8) cosznext(pcols)  ! cosine solar zenith angle next timestep
  integer ncol              ! number of columns in chunk
  integer c             ! chunk index
  integer idum1,idum2,idum3,idum4,i ! temporary variables
  real(r8) snowfall(pcols,begchunk:endchunk)  ! total snowfall rate
!-----------------------------------------------------------------------
!
! Calendar day for next time step
!
  call t_startf ('camice_st')
  nstep = get_nstep()
  dtime = get_step_size()
  rtime=dtime
  cdaynext = get_curr_calday(offset=dtime)
!
! set up snowfall here so it doesn't have to be private in the omp call
!
  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
     do i = 1,ncol
	if (prognostic_icesnow) then
           snowfall(i,c)=srf_state(c)%precsc(i)+srf_state(c)%precsl(i)
        else
           snowfall(i,c)=0.	
        end if
     end do
  end do
  call t_stopf ('camice_st')

!$OMP PARALLEL DO PRIVATE (C, NCOL, LATS, LONS, COSZNEXT,I)

  do c=begchunk,endchunk
     ncol = get_ncols_p(c)

! Sea ice surface fluxes and temperatures

     call seaice (c, ncol, rtime, icefrac(1,c), tsice(1,c), &
                  sicthk(1,c), snowhice(1,c), srf_state(c)%ubot, &
                     srf_state(c)%vbot, srf_state(c)%tbot, &
                  srf_state(c)%qbot, srf_state(c)%thbot, srf_state(c)%zbot, &
                     srf_state(c)%pbot ,srf_state(c)%flwds, &
                  srf_state(c)%sols, srf_state(c)%soll, srf_state(c)%solsd, &
                     srf_state(c)%solld, asdirice(1,c), &
                  aldirice(1,c), asdifice(1,c), aldifice(1,c), &
    	             snowfall(1,c), srf_state(c)%tssub, &
                  srfflx(c)%cflx, srfflx(c)%wsx, srfflx(c)%wsy, &
                     srfflx(c)%ts, srfflx(c)%shf, &
         	  srfflx(c)%lhf, srfflx(c)%lwup, srfflx(c)%tref)
		 
!
! Albedos for next time step 
!
! Note the total albedo here that is returned to the atmosphere 
! model is based on a weighted sum of the albedo over ice and ocean
! using fractional areas from this time step. The absorbed shortwave over
! sea ice in the next step uses ice albedos that are saved at there
! present value but with a NEW fractional area that is input prior to 
! the next time through the sea ice model.  Hence
! there is a time step mismatch in the absorbed solar over sea ice. 
! CCSM would not allow such a thing, but here we are specifying sst, 
! over the ocean fraction anyway so it doesn't really matter. 

     call get_rlat_all_p(c, ncol, lats)
     call get_rlon_all_p(c, ncol, lons)
     call zenith (cdaynext, lats, lons, cosznext, ncol)
     call albice(c,ncol, &
                 srf_state(c)%tbot,snowhice(1,c),cosznext, &
                 srfflx(c)%asdir, srfflx(c)%aldir, &
                 srfflx(c)%asdif, srfflx(c)%aldif)
!
! save off ice albedos for sea ice routine per email Bitz.
! I should note that I made one change to the "physics" from John's
! original fracice implementation. John had the absorbed solar by the
! sea ice equal to the gridcell average.  This is pretty far off when
! the sea ice fraction is small. I realize that it is standard practise
! in many models, but it doesn't have to be.  Therefore I have compute a
! special srfrad over ice and I send the ice albedos to the restart
! file.
!
     do i = 1,ncol
        asdirice(i,c)=srfflx(c)%asdir(i)
        aldirice(i,c)=srfflx(c)%aldir(i)
        asdifice(i,c)=srfflx(c)%asdif(i)
        aldifice(i,c)=srfflx(c)%aldif(i)
     end do

  end do

  return
end subroutine camice
