module comfm1

    !! (wh 2003.07.09)
    !! (wh 2003.10.23)  change the arrays into allocatable ones
    !! (wh 2003.12.01)  ghi added

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid, only: plond, beglatexdyn, endlatexdyn, plev, plevp
    use infnan

    implicit none

    integer  :: itime

    real(r8), allocatable :: su (:,:,:)
    real(r8), allocatable :: sv (:,:,:)
    real(r8), allocatable :: st (:,:,:)
    real(r8), allocatable :: sq (:,:,:)

    real(r8), allocatable :: u  (:,:,:)
    real(r8), allocatable :: v  (:,:,:)
    real(r8), allocatable :: t  (:,:,:)
    real(r8), allocatable :: q  (:,:,:)
    real(r8), allocatable :: ws (:,:,:)
    real(r8), allocatable :: wpa(:,:,:)
    real(r8), allocatable :: ghi(:,:,:)
    real(r8), allocatable :: pes(:,:)       ! ps - pmtop (for fm2003)
    real(r8), allocatable :: ghs(:,:)       

contains

    subroutine initialize_comfm1
        !
        ! Purpose:  Allocate and initialize the comfm1 arrays.
        !
        allocate (su  (plond, beglatexdyn:endlatexdyn, plev))
        allocate (sv  (plond, beglatexdyn:endlatexdyn, plev))
        allocate (st  (plond, beglatexdyn:endlatexdyn, plev))
        allocate (sq  (plond, beglatexdyn:endlatexdyn, plev))

        allocate (u   (plond, beglatexdyn:endlatexdyn, plev))
        allocate (v   (plond, beglatexdyn:endlatexdyn, plev))
        allocate (t   (plond, beglatexdyn:endlatexdyn, plev))
        allocate (q   (plond, beglatexdyn:endlatexdyn, plev))
        allocate (ws  (plond, beglatexdyn:endlatexdyn, plevp))
        allocate (wpa (plond, beglatexdyn:endlatexdyn, plev))
        allocate (ghi (plond, beglatexdyn:endlatexdyn, plevp))
        allocate (pes (plond, beglatexdyn:endlatexdyn))
        allocate (ghs (plond, beglatexdyn:endlatexdyn))

        su (:,:,:) = 0.0
        sv (:,:,:) = 0.0
        st (:,:,:) = 0.0
        sq (:,:,:) = 0.0

        u  (:,:,:) = inf
        v  (:,:,:) = inf 
        t  (:,:,:) = inf 
        q  (:,:,:) = inf
        ws (:,:,:) = inf
        wpa(:,:,:) = inf
        ghi(:,:,:) = inf
        pes(:,:)   = inf
        ghs(:,:)   = inf

        return
    end subroutine initialize_comfm1

end module comfm1
