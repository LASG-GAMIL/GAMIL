#include <misc.h>
#include <params.h>

module buffer

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !   Definition and initialization of time-cycling physics arrays
    !
    ! Author: 
    ! 
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use tracers,      only: pcnst, pnats
    use ppgrid,       only: pcols, pver, begchunk, endchunk
    use prognostics,  only: ptimelevels
    use infnan

    real(r8), allocatable :: qrs(:,:,:)      ! shortwave radiative heating rate 
    real(r8), allocatable :: qrl(:,:,:)      ! longwave  radiative heating rate 

    real(r8), allocatable :: qcwat(:,:,:,:)
    real(r8), allocatable :: tcwat(:,:,:,:)
    real(r8), allocatable :: cld(:,:,:,:)
    real(r8), allocatable :: lcwat(:,:,:,:)

    real(r8), allocatable :: pblht(:,:)      ! PBL height
    real(r8), allocatable :: tpert(:,:)      ! temperature perturbation (PBL)
    real(r8), allocatable :: qpert(:,:,:)    ! moisture/constituent perturb.(PBL)

contains

    subroutine initialize_buffer
        !
        ! Allocate memory
        !
        allocate (qcwat(pcols,pver,begchunk:endchunk,ptimelevels))
        allocate (tcwat(pcols,pver,begchunk:endchunk,ptimelevels))
        allocate (cld  (pcols,pver,begchunk:endchunk,ptimelevels))
        allocate (lcwat(pcols,pver,begchunk:endchunk,ptimelevels))

        allocate (qrs  (pcols,pver,begchunk:endchunk))     
        allocate (qrl  (pcols,pver,begchunk:endchunk))     
        allocate (pblht(pcols,begchunk:endchunk))       
        allocate (tpert(pcols,begchunk:endchunk))       
        allocate (qpert(pcols,pcnst+pnats,begchunk:endchunk)) 
        !
        ! Initialize to NaN or Inf
        !
        qcwat (:,:,:,:) = inf
        tcwat (:,:,:,:) = inf
        cld   (:,:,:,:) = inf
        lcwat (:,:,:,:) = inf

        qrs   (:,:,:) = inf
        qrl   (:,:,:) = inf
        pblht (:,:) = inf
        tpert (:,:) = inf
        qpert (:,:,:) = inf

        return
    end subroutine initialize_buffer

end module buffer
