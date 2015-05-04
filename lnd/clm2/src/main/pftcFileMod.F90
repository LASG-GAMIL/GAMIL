#include <misc.h>
#include <preproc.h>

module pftcFileMod

!=======================================================================
CONTAINS
!=======================================================================

  subroutine pftconrd

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read and initialize vegetation (PFT) constants 
! 
! Method: 
! 
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id: pftcFileMod.F90,v 1.3.6.6 2002/06/15 13:50:35 erik Exp $
!-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar       !parameters
    use clm_varctl       !run control variables
    use pft_varcon       !vegetation type constants
    use fileutils, only : opnfil, getfil, relavu, getavu
    use spmdMod          !spmd variables and routines
    implicit none

! ------------------------ local variables ------------------------
    integer :: i,n              !loop indices
    character(len=256) :: locfn !local file name
    integer :: ier              !error code
! -----------------------------------------------------------------

! Set specific vegetation type values

    ncorn  = 15
    nwheat = 16
    ntree = 8    !value for last type of tree
    noveg = 0    !value for non-vegetated

! Assign unit number to file. Get local file. Open file and read PFT's.
! Close and release file.

    if (masterproc) then
       write (6,*) 'Attempting to read PFT physiological data .....'
       n = getavu()
       call getfil (fpftcon, locfn, 0)
       call opnfil (locfn, n, 'f')
       do i = 1, numpft
          read (n,*)  pftname(i),              &
                      z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
                      vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
                      rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
                      taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
                      roota_par(i), rootb_par(i) 
       end do
       call relavu (n)
    endif

! Define PFT zero to be bare ground

    pftname(noveg) = 'not_vegetated' 
    z0mr(noveg) = 0.
    displar(noveg) = 0.
    dleaf(noveg) = 0.
    c3psn(noveg) = 1.
    vcmx25(noveg) = 0.
    mp(noveg) = 9.
    qe25(noveg) = 0.
    rhol(noveg,1) = 0.
    rhol(noveg,2) = 0.
    rhos(noveg,1) = 0.
    rhos(noveg,2) = 0.
    taul(noveg,1) = 0.
    taul(noveg,2) = 0.
    taus(noveg,1) = 0.
    taus(noveg,2) = 0.
    xl(noveg) = 0.
    roota_par(noveg) = 0.
    rootb_par(noveg) = 0.

#if ( defined SPMD )
    call mpi_bcast (z0mr, size(z0mr), mpir8, 0, mpicom, ier)
    call mpi_bcast (displar, size(displar), mpir8, 0, mpicom, ier)
    call mpi_bcast (dleaf, size(dleaf), mpir8, 0, mpicom, ier)
    call mpi_bcast (c3psn, size(c3psn), mpir8, 0, mpicom, ier)
    call mpi_bcast (vcmx25, size(vcmx25), mpir8, 0, mpicom, ier)
    call mpi_bcast (mp, size(mp), mpir8, 0, mpicom, ier)
    call mpi_bcast (qe25, size(qe25), mpir8, 0, mpicom, ier)
    call mpi_bcast (rhol, size(rhol), mpir8, 0, mpicom, ier)
    call mpi_bcast (rhos, size(rhos), mpir8, 0, mpicom, ier)
    call mpi_bcast (taul, size(taul), mpir8, 0, mpicom, ier)
    call mpi_bcast (taus, size(taus), mpir8, 0, mpicom, ier)
    call mpi_bcast (xl, size(xl), mpir8, 0, mpicom, ier)
    call mpi_bcast (roota_par, size(roota_par), mpir8, 0, mpicom, ier)
    call mpi_bcast (rootb_par, size(rootb_par), mpir8, 0, mpicom, ier)
#endif

    if (masterproc) then
       write (6,*) 'Successfully read PFT physiological data'
       write (6,*)
    endif

    return 
  end subroutine pftconrd

end module pftcFileMod
