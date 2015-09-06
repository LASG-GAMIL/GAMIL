module comfm1

    !! (wh 2003.07.09)
    !! (wh 2003.10.23)  change the arrays into allocatable ones
    !! (wh 2003.12.01)  ghi added

    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid, only: beglatexdyn, endlatexdyn, plev, plevp
    use mpi_gamil, only: register_comm_array, beglonex, endlonex
    use infnan

    implicit none

    integer  :: itime

    real(r8), allocatable :: su (:,:,:)
    real(r8), allocatable :: sv (:,:,:)
    real(r8), allocatable :: st (:,:,:)
    real(r8), allocatable :: sq (:,:,:)
    real(r8), allocatable :: sut (:,:,:)
    real(r8), allocatable :: svt (:,:,:)
    real(r8), allocatable :: stt (:,:,:)

    real(r8), allocatable :: u  (:,:,:)
    real(r8), allocatable :: v  (:,:,:)
    real(r8), allocatable :: t  (:,:,:)
    real(r8), allocatable :: q  (:,:,:)
    real(r8), allocatable :: ws (:,:,:)
    real(r8), allocatable :: wpa(:,:,:)
    real(r8), allocatable :: ghi(:,:,:)
    real(r8), allocatable :: pes(:,:)       ! ps - pmtop (for fm2003)
    real(r8), allocatable :: ghs(:,:)    
    real(r8), allocatable :: u0  (:,:,:)
    real(r8), allocatable :: v0  (:,:,:)
    real(r8), allocatable :: ws0  (:,:,:)
    real(r8), allocatable :: qt  (:,:,:)
    real(r8), allocatable :: dp  (:,:)
    real(r8), allocatable :: pp  (:,:)
    real(r8), allocatable :: fac  (:,:,:)
    real(r8), allocatable :: fbc  (:,:,:)
    real(r8), allocatable :: uuk  (:,:)
    real(r8), allocatable :: hhk  (:,:)
    real(r8), allocatable :: dus  (:,:)
    real(r8), allocatable :: dps2  (:,:)

    real(r8), allocatable :: ply  (:,:,:)
    real(r8), allocatable :: tb (:,:,:)
    real(r8), allocatable :: qk (:,:,:)
    real(r8), allocatable :: cb (:,:,:)
    real(r8), allocatable :: dcb (:,:,:)
    real(r8), allocatable :: cb0 (:,:,:)
    real(r8), allocatable :: p (:,:)
    real(r8), allocatable :: c0 (:,:)
    integer, allocatable :: nigw(:)
    integer, allocatable :: nigw_2D(:,:)
    
    real(r8), allocatable :: uu (:,:,:)
    real(r8), allocatable :: vv (:,:,:)
    real(r8), allocatable :: tt (:,:,:)
    real(r8), allocatable :: dps0 (:,:)
    real(r8), allocatable :: dps1 (:,:)
    real(r8), allocatable :: hps (:,:)
    real(r8), allocatable :: hh (:,:,:)
    real(r8), allocatable :: ttz (:,:,:)
    real(r8), allocatable :: uz (:,:,:)
    real(r8), allocatable :: vz (:,:,:)
    real(r8), allocatable :: ttv (:,:,:)
    real(r8), allocatable :: dps (:,:)


    real(r8), allocatable :: du (:,:,:)
    real(r8), allocatable :: dv (:,:,:)
    real(r8), allocatable :: dtt (:,:,:)
    real(r8), allocatable :: du0 (:,:,:)
    real(r8), allocatable :: dv0 (:,:,:)
    real(r8), allocatable :: dtt0 (:,:,:)
    real(r8), allocatable :: du1 (:,:,:)
    real(r8), allocatable :: dv1 (:,:,:)
    real(r8), allocatable :: dtt1 (:,:,:)

    real(r8), allocatable :: uk (:,:,:)
    real(r8), allocatable :: vk (:,:,:)
    real(r8), allocatable :: ttk (:,:,:)
    real(r8), allocatable :: psk (:,:)

    real(r8), allocatable :: dq (:,:,:)
    real(r8), allocatable :: uq (:,:,:)
    real(r8), allocatable :: vq (:,:,:)
    real(r8), allocatable :: wq (:,:,:)
    real(r8), allocatable :: pq (:,:)

    real(r8), allocatable :: hu (:,:)
    real(r8), allocatable :: hv (:,:)
    real(r8), allocatable :: cu (:,:)
    real(r8), allocatable :: cv (:,:)
    real(r8), allocatable :: h (:,:)

    real(r8), allocatable :: tuu (:,:,:)
    real(r8), allocatable :: tvv (:,:,:)
    real(r8), allocatable :: A (:,:,:)
    real(r8), allocatable :: QH (:,:,:)
    real(r8), allocatable :: QHSTAR (:,:,:)
    real(r8), allocatable :: USTAR (:,:,:)
    real(r8), allocatable :: VSTAR (:,:,:)
    real(r8), allocatable :: BETA (:,:,:)
    real(r8), allocatable :: FX (:,:,:)
    real(r8), allocatable :: FY (:,:,:)

    real(r8), allocatable :: QC (:,:,:)
    real(r8), allocatable :: QR (:,:,:)
    real(r8), allocatable :: QI (:,:,:)
    real(r8), allocatable :: QS (:,:,:)
    real(r8), allocatable :: QG (:,:,:)
    real(r8), allocatable :: D (:,:,:)
    real(r8), allocatable :: DT (:,:,:)
    real(r8), allocatable :: DS (:,:,:)
    real(r8), allocatable :: DA (:,:,:)
    real(r8), allocatable :: DB (:,:,:)
    real(r8), allocatable :: VR (:,:,:)
    real(r8), allocatable :: HQK (:,:,:)
    real(r8), allocatable :: TK (:,:,:)
    real(r8), allocatable :: HVK (:,:,:)
    real(r8), allocatable :: HUK (:,:,:)
    real(r8), allocatable :: ROT (:,:,:)
    real(r8), allocatable :: RLNT (:,:,:)
    real(r8), allocatable :: RDLN (:,:,:)
    real(r8), allocatable :: RDLT (:,:,:)
    real(r8), allocatable :: TW (:,:,:)
    
    real(r8), allocatable :: psb (:,:)
    real(r8), allocatable :: tsb (:,:)


contains

    subroutine initialize_comfm1
        !
        ! Purpose:  Allocate and initialize the comfm1 arrays.
        !
        allocate (su  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (sv  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (st  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (sq  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (sut  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (svt  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (stt  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))

        allocate (u   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (v   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (t   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (q   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (ws  (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (wpa (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (ghi (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (pes (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (ghs (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (u0   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (v0   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (ws0   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qt   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dp   (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (p   (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (fac   (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (fbc   (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (uuk   (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (hhk   (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (dus   (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (dps2   (beglonex:endlonex, beglatexdyn:endlatexdyn))

        allocate (ply   (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (tb   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qk   (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (cb  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dcb  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (cb0  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (c0  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (nigw  (beglatexdyn:endlatexdyn))
        allocate (nigw_2D  (beglonex:endlonex, beglatexdyn:endlatexdyn))

        allocate (uu  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (vv  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (tt  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))

        
        allocate (dps0  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (dps1  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (hps  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (hh  (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp))
        allocate (ttz  (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp+1))
        allocate (uz  (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp+1))
        allocate (vz  (beglonex:endlonex, beglatexdyn:endlatexdyn, plevp+1))
        allocate (ttv  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dps  (beglonex:endlonex, beglatexdyn:endlatexdyn))


        allocate (du  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dv  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dtt  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (du0  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dv0  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dtt0  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (du1  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dv1  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (dtt1  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))

        allocate (uk  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (vk  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (ttk  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (psk  (beglonex:endlonex, beglatexdyn:endlatexdyn))

        allocate (dq  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (uq  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (vq  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (wq  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (pq  (beglonex:endlonex, beglatexdyn:endlatexdyn))

        allocate (hu  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (hv  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (cu  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (cv  (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (h  (beglonex:endlonex, beglatexdyn:endlatexdyn))

        allocate (tuu  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (tvv  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (A  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (QH  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (QHSTAR  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (USTAR  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (VSTAR  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (BETA  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (FX  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (FY  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))


        allocate (qc  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qr  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qi  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qs  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (qg  (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (D (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (DT (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (DS (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (DA (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (DB (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (VR (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (HQK (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (TK (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (HVK (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (HUK (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (ROT (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (RLNT (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (RDLN (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (RDLT (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))
        allocate (TW (beglonex:endlonex, beglatexdyn:endlatexdyn, plev))

        allocate (psb (beglonex:endlonex, beglatexdyn:endlatexdyn))
        allocate (tsb (beglonex:endlonex, beglatexdyn:endlatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,su(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,sv(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,st(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,sq(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,sut(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,svt(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,stt(:,beglatexdyn,1))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,u(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,v(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,t(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,q(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,ws(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,wpa(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,ghi(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,pes(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,ghs(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,u0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,v0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,ws0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qt(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dp(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,p(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,fac(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,fbc(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,uuk(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,hhk(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dus(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dps2(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,ply(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,tb(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qk(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,cb(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dcb(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,cb0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,c0(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,uu(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,vv(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,tt(:,beglatexdyn,1))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dps0(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dps1(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,hps(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp,1,1,hh(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp+1,1,1,ttz(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp+1,1,1,uz(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plevp+1,1,1,vz(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,ttv(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,dps(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,du(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dv(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dtt(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,du0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dv0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dtt0(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,du1(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dv1(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dtt1(:,beglatexdyn,1))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,uk(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,vk(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,ttk(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,psk(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,dq(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,uq(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,vq(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,wq(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,pq(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,hu(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,hv(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,cu(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,cv(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,h(:,beglatexdyn))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,tuu(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,tvv(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,A(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,QH(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,QHSTAR(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,USTAR(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,VSTAR(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,BETA(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,FX(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,FY(:,beglatexdyn,1))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qc(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qr(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qi(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qs(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,qg(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,D(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,DT(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,DS(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,DA(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,DB(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,VR(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,HQK(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,TK(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,HVK(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,HUK(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,ROT(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,RLNT(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,RDLN(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,RDLT(:,beglatexdyn,1))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,plev,1,1,TW(:,beglatexdyn,1))

        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,tsb(:,beglatexdyn))
        call register_comm_array(beglonex,endlonex,beglatexdyn,endlatexdyn,1,1,1,1,psb(:,beglatexdyn))

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

        psb (:,:) = inf
        tsb (:,:) = inf

        return
    end subroutine initialize_comfm1

end module comfm1
