module aerosol_index
!----------------------------------------------------------------------- 
! 
! Purposes: 
!   Store a uniform set of indices for aerosols across
!     radiation, aerosol_radiation_interface, aerosol_optics,
!     prescribed_aerosols, prognostic_aerosols
!-----------------------------------------------------------------------

   implicit none

! naer_all is total number of species
! naer is number of species in climatology
! naer_all = naer + 1 (background "species") + 1 (volcanic)
   integer, public, parameter :: naer_all = 13
   integer, public, parameter :: naer = 11

   real, public, parameter :: wgt_sscm = 6.0 / 7.0 ! Fraction of total seasalt mass in coarse mode

! indices to aerosol array (species portion)
   integer, public, parameter :: &
      idxSUL   =  1, &
      idxSSLTA =  2, & ! accumulation mode
      idxSSLTC =  3, & ! coarse mode
      idxOCPHO =  8, &
      idxBCPHO =  9, &
      idxOCPHI =  10, &
      idxBCPHI = 11, &
      idxBG    = 12, &
      idxVOLC  = 13

! indices to sections of array that represent 
! groups of aerosols
   integer, public, parameter :: &
      idxSSLTfirst    = 2, numSSLT  = 2, &
      idxDUSTfirst    = 4, &
      numDUST         = 4, &
      idxCARBONfirst = 8, &
      numCARBON      = 4

! names of aerosols are they are represented in
! the climatology file and for the purposes
! of outputting masses.  Appended '_V' indicates field has been vertically summed.
   character(len=8), public, parameter :: aerosol_name(naer_all) =  &
     (/"MSUL_V  "&
      ,"MSSLTA_V"&
      ,"MSSLTC_V"& 
      ,"MDUST1_V"&
      ,"MDUST2_V"&
      ,"MDUST3_V"&
      ,"MDUST4_V"&
      ,"MOCPHO_V"&
      ,"MBCPHO_V"&
      ,"MOCPHI_V"&
      ,"MBCPHI_V"&
      ,"Backgrnd"&
      ,"MVOLC_V "/)

! number of different "groups" of aerosols
   integer, parameter, public :: num_aer_groups=6

! which group does each bin belong to?
   integer, dimension(naer_all), parameter, public ::  &
      group =(/1,2,2,3,3,3,3,4,4,4,4,5,6/)

! name of each group
   character(len=10), dimension(num_aer_groups), parameter :: &
      aerosol_names = (/'sul  ','sslt ','dust ','car  ','bg   ','volc '/)
                                                                                
! name of each group
   integer, public, parameter :: &
       idxSULg  = 1, &
       idxSSLTg = 2, &
       idxDUSTg = 3, &
       idxCARg  = 4, &
       idxBGg   = 5, &
       idxVOLCg = 6
      
! which band is the "visible" band for output to diagnostics
   integer, parameter, public :: idxVIS = 8    ! index to visible band

  private
  save


end module aerosol_index
