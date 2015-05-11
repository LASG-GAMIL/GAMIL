
#include <misc.h> !
#include <params.h> !

module aerosol_mass_interface
!----------------------------------------------------------------------- 
! 
! Purpose: collect aerosol mass 
!  sxj  2009-03-09
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pcols, pver, pverp, begchunk, endchunk
use prescribed_aerosols, only: aerosol_mass_get, scenario_carbon_scale
use aerosol_index, only: idxCARBONfirst, numCARBON, idxDUSTfirst, &
   numDUST, idxSUL, idxSSLTA, idxSSLTC, idxVOLC, idxOCPHO, idxOCPHI, idxBCPHO, &
   idxBCPHI, idxBG, idxVOLC, idxSSLTfirst, numSSLT, &
   naer_all, num_aer_groups, group, &
   idxCARg, idxSULg, idxSSLTg, idxDUSTg, idxVOLCg, idxBGg
!use aerosol_intr, only: set_aerosol_from_prognostics
  USE pmgrid,      ONLY: masterproc,iam     !! for debug-test 
#ifdef SPMD               !  for debug sxj--
   use mpishorthand     
#endif                     
implicit none
                                                                               
private
save
! portion of each species group to use in computation
! of aerosol forcing in driving the climate
  real(r8) :: sulscl  = 1._r8
  real(r8) :: carscl  = 1._r8
  real(r8) :: ssltscl = 1._r8
  real(r8) :: dustscl = 1._r8
  real(r8) :: volcscl = 1._r8


  real(r8), public, allocatable :: aer_mass(:, :, :, :) 


! Procedure to get masses of aerosols in daylit 
! columns for shortwave radiative forcing computation
  public collect_aer_masses

  public aerosol_mass_init



!===============================================================================
CONTAINS
!===============================================================================

subroutine aerosol_mass_init()

  allocate (aer_mass(pcols, pver, naer_all, begchunk:endchunk) )

end subroutine aerosol_mass_init



subroutine collect_aer_masses(state)

   use physics_types,   only: physics_state

   type(physics_state), intent(in) :: state

   integer lchnk !  chunk index
   real(r8) :: int_scales(naer_all)    ! scaling factors for aerosols

   lchnk = state%lchnk

   call get_int_scales(int_scales)

   !write(*,*) "aerosol_mass_interface.F90 line66"
   !call endrun

   call aerosol_mass_get(lchnk, state%pint, aer_mass(:,:,:,lchnk), int_scales)

!#ifdef SPMD               !  for debug sxj--
!   CALL mpibarrier (mpicom)  
!   write(*,*) "aerosol_mass_interface.F90_line73"                     
!   write(*,*) iam,lchnk,"aer_mass(:,24,1,lchnk)"
!   write(*,*) aer_mass(:,24,1,lchnk)
  ! call endrun
!#endif

! overwrite with prognostics aerosols
   !call set_aerosol_from_prognostics (state, aer_mass(:,:,:,lchnk), int_scales)

end subroutine collect_aer_masses


                                                                                
subroutine get_int_scales(scales)
  real(r8), intent(out)::scales(naer_all)  ! scale each aerosol by this amount
  integer i                                  ! index through species
                                                                                
  scales(idxBG) = 1._r8
  scales(idxSUL) = sulscl
  do i = idxSSLTfirst, idxSSLTfirst+numSSLT-1
    scales(i) = ssltscl
  enddo

  do i = idxCARBONfirst, idxCARBONfirst+numCARBON-1
    scales(i) = carscl
  enddo
                                                                                
  do i = idxDUSTfirst, idxDUSTfirst+numDUST-1
    scales(i) = dustscl
  enddo

  scales(idxVOLC) = volcscl
                                                                                
  return
end subroutine get_int_scales


end module aerosol_mass_interface
