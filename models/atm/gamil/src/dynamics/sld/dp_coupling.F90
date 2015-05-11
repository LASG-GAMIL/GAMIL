!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,        only: pcols, pver
   use rgrid,         only: nlon
   use pmgrid
   use phys_grid
   use physics_types, only: physics_state, physics_tend
   use constituents,  only: pcnst, pnats

   implicit none

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(ps, t3, u3, v3, q3, omga, phis, phys_state)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
!------------------------------------------------------------------------------
    use physconst,     only: cappa
!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (plond, beglat:endlat)            ! surface pressure
    real(r8), intent(in) :: t3  (plond, plev, beglatex:beglatex+numlats-1)  ! temperature
    real(r8), intent(in) :: u3  (plond, plev, beglatex:beglatex+numlats-1)  ! u-wind component
    real(r8), intent(in) :: v3  (plond, plev, beglatex:beglatex+numlats-1)  ! v-wind component
    real(r8), intent(in) :: q3  (plond, plev, pcnst+pnats, beglatex:beglatex+numlats-1) ! constituents
    real(r8), intent(in) :: omga(plond, plev, beglat:endlat)      ! vertical velocity
    real(r8), intent(in) :: phis(plond, beglat:endlat)            ! Surface geopotential

    type(physics_state), intent(out), dimension(begchunk:endchunk) :: phys_state
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,k,j,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Copy data from dynamics data structure to physics data structure
!-----------------------------------------------------------------------
    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, ncol, lons)
          call get_lat_all_p(lchnk, ncol, lats)
          phys_state(lchnk)%ncol  = ncol
          phys_state(lchnk)%lchnk = lchnk

          do i=1,ncol
             phys_state(lchnk)%ps   (i)     = ps  (lons(i),lats(i))
             phys_state(lchnk)%phis (i)     = phis(lons(i),lats(i))
          end do

          do k=1,plev
             do i=1,ncol
                phys_state(lchnk)%t    (i,k)   = t3  (lons(i),k,lats(i))
                phys_state(lchnk)%u    (i,k)   = u3  (lons(i),k,lats(i))
                phys_state(lchnk)%v    (i,k)   = v3  (lons(i),k,lats(i))
                phys_state(lchnk)%omega(i,k)   = omga(lons(i),k,lats(i))
             end do
          end do

          do m=1,pcnst+pnats
             do k=1,plev
                do i=1,ncol
                   phys_state(lchnk)%q    (i,k,m) = q3  (lons(i),k,m,lats(i))
                end do
             end do
          end do

       end do

    else

       tsize = 4 + (pcnst+pnats)

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

       do j=beglat,endlat

          call block_to_chunk_send_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

             bbuffer(bpter(i,0))   = ps  (i,j)
             bbuffer(bpter(i,0)+1) = phis(i,j)

             do k=1,plev

                bbuffer(bpter(i,k))   = t3  (i,k,j)
                bbuffer(bpter(i,k)+1) = u3  (i,k,j)
                bbuffer(bpter(i,k)+2) = v3  (i,k,j)
                bbuffer(bpter(i,k)+3) = omga(i,k,j)

                do m=1,pcnst+pnats
                   bbuffer(bpter(i,k)+3+m) = q3  (i,k,m,j)
                end do

             end do

          end do

       end do

       call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)
          phys_state(lchnk)%ncol  = ncol
          phys_state(lchnk)%lchnk = lchnk

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncol

             phys_state(lchnk)%ps   (i)     = cbuffer(cpter(i,0))
             phys_state(lchnk)%phis (i)     = cbuffer(cpter(i,0)+1)

             do k=1,plev

                phys_state(lchnk)%t    (i,k)   = cbuffer(cpter(i,k))
                phys_state(lchnk)%u    (i,k)   = cbuffer(cpter(i,k)+1)
                phys_state(lchnk)%v    (i,k)   = cbuffer(cpter(i,k)+2)
                phys_state(lchnk)%omega (i,k)   = cbuffer(cpter(i,k)+3)

                do m=1,pcnst+pnats
                   phys_state(lchnk)%q (i,k,m) = cbuffer(cpter(i,k)+3+m)
                end do

             end do

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)

    endif

!-----------------------------------------------------------------------
! Fill auxilliary arrays in physics data structure
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)

! pressure arrays
       call plevs0(ncol, pcols, pver, &
                   phys_state(lchnk)%ps,   phys_state(lchnk)%pint,    &
                   phys_state(lchnk)%pmid, phys_state(lchnk)%pdel)

! log(pressure) arrays and Exner function
       do k=1,pver+1
          do i=1,ncol
             phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
          end do
       end do
       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

    end do

    return
  end subroutine d_p_coupling

!===============================================================================
  subroutine p_d_coupling(phys_state, phys_tend, t2, fu, fv, flx_net, q3)
!------------------------------------------------------------------------------
! Coupler for converting physics output variables into dynamics input variables
!------------------------------Arguments--------------------------------
    type(physics_state),intent(in), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend), intent(in), dimension(begchunk:endchunk) :: phys_tend

    real(r8), intent(out) :: t2(plond, plev, beglat:endlat)        ! temp tendency
    real(r8), intent(out) :: fu(plond, plev, beglat:endlat)        ! u wind tendency
    real(r8), intent(out) :: fv(plond, plev, beglat:endlat)        ! v wind tendency
    real(r8), intent(out) :: flx_net(plond,beglat:endlat)          ! net flux
    real(r8), intent(out) :: q3(plond, plev, pcnst+pnats, beglatex:beglatex+numlats-1) ! constituents
!
!---------------------------Local workspace-----------------------------
    real(r8), allocatable, dimension(:) :: &
       bbuffer, cbuffer              ! transpose buffers

    integer :: i,j,k,m,lchnk         ! indices
    integer :: ncol                  ! number of columns in current chunk
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(plon,0:plev)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
!-----------------------------------------------------------------------

    if (local_dp_map) then

!$OMP PARALLEL DO PRIVATE (LCHNK, NCOL, I, K, M, LONS, LATS)

      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         call get_lon_all_p(lchnk, ncol, lons)
         call get_lat_all_p(lchnk, ncol, lats)

         do k=1,plev
            do i=1,ncol
               t2(lons(i),k,lats(i))   = phys_tend(lchnk)%dTdt (i,k)
               fu(lons(i),k,lats(i))   = phys_tend(lchnk)%dudt (i,k)
               fv(lons(i),k,lats(i))   = phys_tend(lchnk)%dvdt (i,k)
            end do
         end do

         do i=1,ncol
            flx_net(lons(i),lats(i))   = phys_tend(lchnk)%flx_net(i)
         end do

         do m=1,pcnst+pnats
            do k=1,plev
               do i=1,ncol
                  q3(lons(i),k,m,lats(i)) = phys_state(lchnk)%q(i,k,m)
               end do
            end do
         end do
      end do

    else

       tsize = 3 + (pcnst+pnats)

       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

       do lchnk = begchunk,endchunk
          ncol = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,plev+1,tsize,cpter)

          do i=1,ncol

             cbuffer(cpter(i,0)) = phys_tend(lchnk)%flx_net(i)

             do k=1,plev

                cbuffer(cpter(i,k))   = phys_tend(lchnk)%dTdt (i,k)
                cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt (i,k)
                cbuffer(cpter(i,k)+2) = phys_tend(lchnk)%dvdt (i,k)

                do m=1,pcnst+pnats
                   cbuffer(cpter(i,k)+2+m) = phys_state(lchnk)%q(i,k,m)
                end do

             end do

          end do

       end do

       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)

       do j=beglat,endlat

          call chunk_to_block_recv_pters(j,plon,plev+1,tsize,bpter)

          do i=1,nlon(j)

             flx_net(i,j) = bbuffer(bpter(i,0))

             do k=1,plev

                t2(i,k,j) = bbuffer(bpter(i,k))
                fu(i,k,j) = bbuffer(bpter(i,k)+1)
                fv(i,k,j) = bbuffer(bpter(i,k)+2)

                do m=1,pcnst+pnats
                   q3(i,k,m,j) = bbuffer(bpter(i,k)+2+m)
                end do

             end do

          end do

       end do

       deallocate(bbuffer)
       deallocate(cbuffer)

    endif

    return
  end subroutine p_d_coupling
end module dp_coupling
