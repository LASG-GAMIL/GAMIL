#include <misc.h>
#include <params.h>

module cloudsimulator
!-----------------------------------------------------------------------
!Purpose: CAM interface to
!         Name:		ISCCP Simulator ICARUS/SCOPS
!         What:		Simulate ISCCP cloud products from GCM inputs
!         Version:	3.4
!         Authors:	Steve Klein (sak@gfdl.noaa.gov)
!         		Mark Webb (Mark.Webb@MetOffice.com)
!
!Author:  W. Lin,       Original version
!         J. Rosinski,  Modifications, March 2003
!         B. Eaton,     Update to simulator version 3.4, March 2004
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, pver
   use history,      only: addfld, add_default, phys_decomp, outfld, fillvalue
   use icarus_scops, only: ntau, npres

   implicit none
   private
   save

!Public functions/subroutines

   public :: &
      cloudsimulator_init,    &
      cloudsimulator_run

   logical, public :: doisccp = .false.    ! whether to do ISCCP calcs and I/O     

CONTAINS

subroutine cloudsimulator_init

   use icarus_scops, only: icarus_scops_init

   !------------------------------------------------------------------------------

   call addfld ('FISCCP1 ','mixed   ',npres*ntau,'A', &
                'grid box fraction covered by each ISCCP D level cloud types',phys_decomp)
   call addfld ('TCLDAREA','fraction',1,'A','Total cloud area',phys_decomp)
   call addfld ('MEANPTOP','mb      ',1,'A','Mean cloud top pressure',phys_decomp)
   call addfld ('MEANTAU ','unitless',1,'A','Mean optical thickness',phys_decomp)
   call addfld ('MEANTTOP','K       ',1,'A','Mean cloud top temperature',phys_decomp)
   call addfld ('CLOUDY  ','fraction',1,'A','Binary flag for TCLDAREA > cldmin',phys_decomp)

   call add_default ('FISCCP1 ', 1, ' ')
   call add_default ('TCLDAREA', 1, ' ')
   call add_default ('MEANPTOP', 1, ' ')
   call add_default ('MEANTAU ', 1, ' ')
   call add_default ('MEANTTOP', 1, ' ')
   call add_default ('CLOUDY  ', 1, ' ')

   call icarus_scops_init

end subroutine cloudsimulator_init

subroutine cloudsimulator_run(state, ts, concld, cld, cliqwp, &
                              cicewp, rel, rei, emis, coszrs  ) 

   use physics_types,   only: physics_state
   use icarus_scops,    only: isccp_cloud_types

   type(physics_state), intent(in) :: state

   real(r8), intent(in) :: ts(pcols)           ! skin temperature
   real(r8), intent(in) :: concld(pcols,pver)  ! convective cloud cover
   real(r8), intent(in) :: cld(pcols,pver)     ! cloud cover
   real(r8), intent(in) :: cliqwp(pcols,pver)  ! in-cloud liquid water path
   real(r8), intent(in) :: cicewp(pcols,pver)  ! in-cloud ice water path
   real(r8), intent(in) :: rel(pcols,pver)     ! Liquid cloud particle effective radius
   real(r8), intent(in) :: rei(pcols,pver)     ! Ice effective drop size (microns)
   real(r8), intent(in) :: emis(pcols,pver)    ! Cloud longwave emissivity
   real(r8), intent(in) :: coszrs(pcols)       ! cosine solar zenith angle (to tell if day or night)

! Local variables:

   integer, parameter :: debug = 0     ! set to non-zero value to print out inputs
                                       ! with step debug
   integer, parameter :: debugcol = 0  ! set to non-zero value to print out column
                                       ! decomposition with step debugcol

   integer, parameter :: top_height = 1 ! 1 = adjust top height using both a computed
                                        ! infrared brightness temperature and the visible
					! optical depth to adjust cloud top pressure. Note
					! that this calculation is most appropriate to compare
					! to ISCCP data during sunlit hours.
   integer, parameter :: overlap = 3    ! 3=max/rand

   integer,  parameter :: nsubcol = 50      ! # of columns each grid box is subdivided into
   real(r8), parameter :: emsfc_lw = 0.99   ! longwave emissivity of surface at 10.5 microns

! constants for optical depth calculation from radcswmx.F90
   real(r8), parameter :: abarl = 2.817e-02    ! A coefficient for extinction optical depth
   real(r8), parameter :: bbarl = 1.305        ! b coefficient for extinction optical depth
   real(r8), parameter :: abari = 3.448e-03    ! A coefficient for extinction optical depth
   real(r8), parameter :: bbari = 2.431        ! b coefficient for extinction optical depth
   real(r8), parameter :: cldmin = 1.0e-80_r8  ! note: cldmin much less than cldmin from cldnrh
   real(r8), parameter :: cldeps = 0.0_r8      ! 

   integer :: lchnk                         ! chunk identifier
   integer :: ncol                          ! number of atmospheric columns

   integer  :: sunlit(pcols)                ! 1 for day points, 0 for night time
   integer  :: seed1(pcols)                  ! seed values for random number generator
   integer  :: seed2(pcols)                  ! seed values for random number generator
   integer  :: seed3(pcols)                  ! seed values for random number generator
   integer  :: seed4(pcols)                  ! seed values for random number generator

   real(r8) :: tau(pcols,pver)              ! optical depth

   real(r8) :: fq_isccp(pcols,ntau,npres)   ! the fraction of the model grid box covered by
                                            ! each of the 49 ISCCP D level cloud types
   real(r8) :: totalcldarea(pcols)          ! the fraction of model grid box columns
                                            ! with cloud somewhere in them.  This should
					    ! equal the sum over all entries of fq_isccp
   real(r8) :: meanptop(pcols)              ! mean cloud top pressure (mb) - linear averaging
                                            ! in cloud top pressure.
   real(r8) :: meanttop(pcols)              ! mean cloud top temp (k) - linear averaging

   real(r8) :: meantaucld(pcols)            ! mean optical thickness 
                                            ! linear averaging in albedo performed.

   real(r8) :: boxtau(pcols,nsubcol)        ! optical thickness in each column
      
   real(r8) :: boxptop(pcols,nsubcol)       ! cloud top pressure (mb) in each column

   real(r8) :: fq_isccp_s1(pcols,ntau*npres)! accumulated fq_isccp

   integer :: i, k, it, ip, itaunpres       ! indices

   real(r8) :: cloudy(pcols)   ! cloudy flag, which may be used to derived mean values under cloudy
        		       ! conditions. The cloud flag itself is freuency of cloudy condition
	 		       ! when average over an accumulation period.
!------------------------------------------------------------------------------------------------

   call t_startf ('cloudsimulator_run')

   lchnk = state%lchnk
   ncol  = state%ncol

   ! compute the optical depth (code from radcswmx)
   do k=1,pver
      do i=1,ncol
         if(cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
            tau(i,k) = (abarl + bbarl/rel(i,k)) * cliqwp(i,k) + &
                       (abari + bbari/rei(i,k)) * cicewp(i,k)
         else
            tau(i,k) = 0.0
         endif
      end do
   end do

   sunlit = 0
   do i=1,ncol
      if (coszrs(i) > 0.) sunlit(i) = 1
      seed1(i) = (state%pmid(i,pver) - int(state%pmid(i,pver)))      * 1000000000
      seed2(i) = (state%pmid(i,pver-1) - int(state%pmid(i,pver-1)))  * 1000000000
      seed3(i) = (state%pmid(i,pver-2) - int(state%pmid(i,pver-2)))  * 1000000000
      seed4(i) = (state%pmid(i,pver-3) - int(state%pmid(i,pver-3)))  * 1000000000
   end do

   call ISCCP_CLOUD_TYPES( &
      debug,               &
      debugcol,            &
      ncol,                &
      sunlit(:ncol),       &
      pver,                &
      nsubcol,             &
      seed1(:ncol),        &
      seed2(:ncol),        &
      seed3(:ncol),        &
      seed4(:ncol),        &
      state%pmid(:ncol,:), &
      state%pint(:ncol,:), &
      state%q(:ncol,:,1),  &
      cld(:ncol,:),        &
      concld(:ncol,:),     &
      tau(:ncol,:),        &
      tau(:ncol,:),        &
      top_height,          &
      overlap,             &
      ts(:ncol),           &
      emsfc_lw,            &
      state%t(:ncol,:),    &
      emis(:ncol,:),       &
      emis(:ncol,:),       &
      fq_isccp(:ncol,:,:), &
      totalcldarea(:ncol), &
      meanptop(:ncol),     &
      meanttop(:ncol),     &
      meantaucld(:ncol),   &
      boxtau(:ncol,:),     &
      boxptop(:ncol,:)     )
            
   do i=1,ncol
      if (coszrs(i) > 0.) then
         if (totalcldarea(i) >= cldmin) then
            cloudy(i) = 1.0

            ! save standard ISCCP type of 7x7 clouds
            do ip=1,npres
               do it=1,ntau
                  itaunpres = (ip-1)*ntau+it
                  fq_isccp_s1(i,itaunpres) = fq_isccp(i,it,ip)
               end do
            end do
            
         else                             !cloud free in the (daytime) grid box
            fq_isccp_s1(i,:) = 0.
            totalcldarea(i)  = 0.
            meanptop(i)      = fillvalue
            meantaucld(i)    = fillvalue
            meanttop(i)      = fillvalue
            cloudy(i)        = 0.
         endif
      else                                ! nighttime
         fq_isccp_s1(i,:) = fillvalue
         totalcldarea(i)  = fillvalue
         meanptop(i)      = fillvalue
         meantaucld(i)    = fillvalue
         meanttop(i)      = fillvalue
         cloudy(i)        = fillvalue
      end if
   end do
!
! dont need to call outfld if all points are nighttime
!
   if (any(coszrs(:ncol) > 0.)) then
      call outfld('FISCCP1 ',fq_isccp_s1, pcols,lchnk)
      call outfld('TCLDAREA',totalcldarea,pcols,lchnk)
      call outfld('MEANPTOP',meanptop    ,pcols,lchnk)
      call outfld('MEANTAU ',meantaucld  ,pcols,lchnk)
      call outfld('MEANTTOP',meanttop    ,pcols,lchnk)
      call outfld('CLOUDY  ',cloudy      ,pcols,lchnk)
   end if

   call t_stopf ('cloudsimulator_run')

end subroutine cloudsimulator_run

!#######################################################################

end module cloudsimulator
