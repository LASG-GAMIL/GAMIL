!***************************************************************
!  This is a source file of GAMIL, which registers all variables 
!  into C-Coupler library. This file was initially finished by
!  Dr. Li Liu. If you have any problem, please contact Dr. Li 
!  Liu via liuli-cess@tsinghua.edu.cn
!***************************************************************


module diag_all_variables

    use c_coupler_interface_mod
    use comfm1
    use pmgrid, only: beglatexdyn, endlatexdyn, plev, plevp
    use mpi_gamil, only: register_comm_array, beglonex, endlonex, ibeg0, iend0, jbeg0, jend0,ibeg1,iend1


   interface diag_one_dyn_variable ; module procedure &
        diag_one_dyn_variable_2D, &
        diag_one_dyn_variable_3D
   end interface



contains

    subroutine diag_one_dyn_variable_2D(var, hint)
       implicit none
       character(len=*), intent(in) :: hint      
       real                         :: var(beglonex:endlonex, beglatexdyn:endlatexdyn)
       real                         :: tmp_array(beglonex:endlonex, beglatexdyn:endlatexdyn)
       integer                      :: i,j

        tmp_array=0.0
        do j=jbeg0,jend0
        do i=ibeg1,iend1
           tmp_array(i,j)=var(i,j)
        enddo
        enddo

        call c_coupler_check_sum_for_external_data(tmp_array,(endlatexdyn-beglatexdyn+1)*(endlonex-beglonex+1), hint)
        

    end subroutine diag_one_dyn_variable_2D


    subroutine diag_one_dyn_variable_3D(var, hint)
       implicit none
       character(len=*), intent(in) :: hint    
       real                         :: var(beglonex:endlonex,beglatexdyn:endlatexdyn,plev)
       real                         :: tmp_array(beglonex:endlonex,beglatexdyn:endlatexdyn,plev)
       integer                      :: i,j,k

        tmp_array=0.0
        do k=1,plev
        do j=jbeg0,jend0
        do i=ibeg1,iend1
           tmp_array(i,j,k)=var(i,j,k)
        enddo
        enddo
        enddo

        call c_coupler_check_sum_for_external_data(tmp_array,(endlatexdyn-beglatexdyn+1)*(endlonex-beglonex+1)*plev, hint)
        

    end subroutine diag_one_dyn_variable_3D



    subroutine diag_dyn_variables
       implicit none
       
       call diag_one_dyn_variable(du, "okokok du")
       call diag_one_dyn_variable(dv, "okokok dv")
       call diag_one_dyn_variable(dtt, "okokok dtt")
       call diag_one_dyn_variable(du0, "okokok du0")
       call diag_one_dyn_variable(dv0, "okokok dv0")
       call diag_one_dyn_variable(dtt0, "okokok dtt0")
       call diag_one_dyn_variable(dps0, "okokok dps0")
       call diag_one_dyn_variable(du1, "okokok du1")
       call diag_one_dyn_variable(dv1, "okokok dv1")
       call diag_one_dyn_variable(dtt1, "okokok dtt1")
       call diag_one_dyn_variable(dps1, "okokok dps1")
       call diag_one_dyn_variable(uu, "okokok uu")
       call diag_one_dyn_variable(vv, "okokok vv")
       call diag_one_dyn_variable(tt, "okokok tt")
       call diag_one_dyn_variable(p, "okokok p")
       call diag_one_dyn_variable(ply, "okokok ply")
       call diag_one_dyn_variable(dps, "okokok dps")
       call diag_one_dyn_variable(uk, "okokok uk")
       call diag_one_dyn_variable(vk, "okokok vk")
       call diag_one_dyn_variable(ttk, "okokok ttk")
       call diag_one_dyn_variable(psk, "okokok psk")
       call diag_one_dyn_variable(tb, "okokok tb")
       call diag_one_dyn_variable(cb, "okokok cb")
       call diag_one_dyn_variable(dcb, "okokok dcb")
       call diag_one_dyn_variable(hps, "okokok hps")
       call diag_one_dyn_variable(c0, "okokok c0")
       call diag_one_dyn_variable(cb0, "okokok cb0")

       call diag_one_dyn_variable(su, "okokok su")
       call diag_one_dyn_variable(sv, "okokok sv")
       call diag_one_dyn_variable(st, "okokok st")
       call diag_one_dyn_variable(sq, "okokok sq")
       call diag_one_dyn_variable(u, "okokok u")
       call diag_one_dyn_variable(v, "okokok v")
       call diag_one_dyn_variable(t, "okokok t")
       call diag_one_dyn_variable(q, "okokok q")
       call diag_one_dyn_variable(ws, "okokok ws")
       call diag_one_dyn_variable(wpa, "okokok wpa")
       call diag_one_dyn_variable(ghi, "okokok ghi")
       call diag_one_dyn_variable(pes, "okokok pes")
       call diag_one_dyn_variable(ghs, "okokok ghs")

    end subroutine diag_dyn_variables



    subroutine diag_phys_variables(phys_state, phys_state0, phys_tend)
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_types,  only: physics_state, physics_tend
    use phys_grid,      only: get_ncols_p
    implicit none
    real, allocatable :: chunk_array_2D(:,:)        
    real, allocatable :: chunk_array_3D(:,:,:) 
    integer               :: lchnk, ncols, i, k
    type(physics_state),intent(in)   :: phys_state(begchunk:endchunk)
    type(physics_state),intent(in)   :: phys_state0(begchunk:endchunk)
    type(physics_tend ),intent(in)   :: phys_tend(begchunk:endchunk)

    allocate(chunk_array_2D(pcols,begchunk:endchunk))
    allocate(chunk_array_3D(pcols,pver,begchunk:endchunk))

    chunk_array_2D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do i = 1, ncols
          chunk_array_2D(i,lchnk)     =  phys_state(lchnk)%ps(i)
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_2D,pcols*(endchunk-begchunk+1), "okokok phys_state ps")

    chunk_array_2D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do i = 1, ncols
          chunk_array_2D(i,lchnk)     =  phys_state(lchnk)%phis(i)
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_2D,pcols*(endchunk-begchunk+1), "okokok phys_state phis")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%t(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state t")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%t1(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state t1")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%t2(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state t2")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%u(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state u")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%v(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state v")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%s(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state s")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%omega(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state omega")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%pmid(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state pmid")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%pdel(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state pdel")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%rpdel(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state rpdel")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%lnpmid(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state lnpmid")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%exner(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state exner")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%zm(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state zm")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%pint(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state pint")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%lnpint(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state lnpint")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%zi(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state zi")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%q(i,k,1)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state q")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%q1(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state q1")


    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state(lchnk)%q2(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state q2")

    chunk_array_2D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do i = 1, ncols
          chunk_array_2D(i,lchnk)     =  phys_state0(lchnk)%ps(i)
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_2D,pcols*(endchunk-begchunk+1), "okokok phys_state0 ps")

    chunk_array_2D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do i = 1, ncols
          chunk_array_2D(i,lchnk)     =  phys_state0(lchnk)%phis(i)
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_2D,pcols*(endchunk-begchunk+1), "okokok phys_state0 phis")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%t(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 t")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%t1(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 t1")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%t2(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 t2")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%u(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 u")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%v(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 v")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%s(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 s")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%omega(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 omega")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%pmid(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 pmid")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%pdel(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 pdel")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%rpdel(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 rpdel")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%lnpmid(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 lnpmid")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%exner(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 exner")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%zm(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 zm")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%pint(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 pint")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%lnpint(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 lnpint")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%zi(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 zi")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%q(i,k,1)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 q")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%q1(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 q1")


    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_state0(lchnk)%q2(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_state0 q2")


    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_tend(lchnk)%dtdt(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_tend dtdt")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_tend(lchnk)%dudt(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_tend dudt")

    chunk_array_3D = 0.0
    do lchnk = begchunk, endchunk
       ncols = get_ncols_p(lchnk)
       do k = 1, pver
       do i = 1, ncols
          chunk_array_3D(i,k,lchnk)     =  phys_tend(lchnk)%dvdt(i,k)
       end do
       end do
    end do
    call c_coupler_check_sum_for_external_data(chunk_array_3D,pcols*(endchunk-begchunk+1)*pver, "okokok phys_tend dvdt")

    deallocate(chunk_array_2D)
    deallocate(chunk_array_3D)

    end subroutine diag_phys_variables

end module diag_all_variables

