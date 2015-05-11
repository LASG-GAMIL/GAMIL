#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  dynpkg --- Driver for the NASA finite-volume dynamical core
!
! !INTERFACE:
subroutine dynpkg( im,      jm,      km,      u,       v,      &
                   ps,      delp,    pe,      pk,      pkz,    &
                   peln,    ptop,    jfirst,  jlast,           &
                   kfirst,  klast,   klastp,  phis,            &
                   omga,    cp,      nq,      q3,      ndt,    &
                   om,      rair,    ae,      ns0,     nc,     &
                   pt,      convt,   iord,    jord,   &
                   kord,    te0,     uxy,     vxy,     psxy,   &
                   delpxy,  pexy,    pkxy,    pkzxy,   pelnxy, &
                   phisxy,  omgaxy,  q3xy,    ptxy,      &
                   ifirstxy,ilastxy, jfirstxy,jlastxy, dcaf,   &
                   rayf,    ideal )

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use dynamics_vars, only : coslon, sinlon, cosl5, sinl5, acosp,  &
                       acap, rcap, sine, cosp, sinp, cose,   &
                       ns, icd, jcd, ng_c, ng_d, ng_s, q3t, yzt, xyt
   use pmgrid, only: twod_decomp, myid_z, npr_z, npr_y, iam
   use physconst, only: cappa, gravit

#if defined( SPMD )
   use spmd_dyn, only : comm_z, inter_ijk, inter_ikjp, inter_q3
   use mod_comm, only : mp_send3d, mp_recv3d
   use redistributemodule, only : redistributetype, redistributestart,  &
                  redistributefinish
#endif

   implicit none

! !INPUT PARAMETERS:
   integer, intent(in):: im        ! dimension in east-west
   integer, intent(in):: jm        ! dimension in North-South
   integer, intent(in):: km        ! number of Lagrangian layers
   integer, intent(in):: jfirst    ! starting latitude index for MPI
   integer, intent(in):: jlast     ! ending latitude index for MPI
   integer, intent(in):: kfirst    ! starting vertical index for MPI
   integer, intent(in):: klast     ! ending vertical index for MPI
   integer, intent(in):: klastp    ! klast, except km+1 when klast=km
   integer, intent(in):: nq        ! total # of tracers to be advected
   integer, intent(in):: nc        ! declared dimension of q3
   integer, intent(in):: ndt       ! the large time step in seconds
                                   ! Also the mapping time step in this setup
   integer, intent(in):: iord      ! parameter controlling monotonicity in E-W
                                   ! recommendation: iord=4
   integer, intent(in):: jord      ! parameter controlling monotonicity in N-S 
                                   ! recommendation: jord=4
   integer, intent(in):: kord      ! parameter controlling monotonicity in mapping
   integer, intent(in):: ifirstxy, ilastxy, jfirstxy, jlastxy  ! xy decomposition
!                                    index ranges
                                   ! recommendation: kord=4
   real(r8), intent(in):: phis(im,jfirst:jlast)   ! surface geopotential (grav*zs)
   real(r8), intent(in):: om       ! angular velocity of earth's rotation  
   real(r8), intent(in):: cp       ! heat capacity of air at constant pressure
   real(r8), intent(in):: ae       ! radius of the earth (m)
   real(r8), intent(in):: rair     ! Gas constant of the air
   real(r8), intent(in):: ptop     ! Pressure at model top (interface pres)
   logical,  intent(in):: convt    ! Flag to control pt output
                                   ! If true: output pt, the virtual temperature
                                   ! If false: pt is updated
   logical,  intent(in):: dcaf     ! Dry convection flag (use only in ideal_phys case)
   logical,  intent(in):: rayf     ! Rayleigh friction flag (off by default)
   logical,  intent(in):: ideal    ! True for ideal physics

   real(r8), intent(in):: phisxy(ifirstxy:ilastxy,jfirstxy:jlastxy)  ! xy-decomposed version of phis

! !INPUT/OUTPUT PARAMETERS:
   integer, intent(inout):: ns0    ! total number of splits for Lagrangian
                                   ! dynamics; a suitable value will be automatically
                                   ! determined if ns0 = 0 (i.e., if you don't know what value
                                   ! to be used try ns0=0
   real(r8), intent(inout):: pt(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) 
                           ! scaled (virtual) potential temperature
   real(r8), intent(inout):: ps(im,jfirst:jlast)                ! Surface pressure (pa) 
   real(r8), intent(inout):: delp(im,jfirst:jlast,kfirst:klast) ! Pressure thickness
   real(r8), intent(inout):: u(im,jfirst-ng_d:jlast+ng_s,kfirst:klast)    ! u wind velocities, staggered grid
   real(r8), intent(inout):: v(im,jfirst-ng_s:jlast+ng_d,kfirst:klast)    ! v wind velocities, staggered grid
   real(r8), intent(inout):: q3(im,jfirst-ng_d:jlast+ng_d,kfirst:klast,nc) ! Tracers

!--------------------------------------------------------------------------------------
! The following three arrays must be pre-computed as input to benergy(). They are NOT
! needed if consv=.F.; updated on output (to be used by physdrv)
! Please refer to routine pkez on the algorithm for computing pkz
! from pe and pk
!--------------------------------------------------------------------------------------

   real(r8), intent(inout):: pe(im,kfirst:klast+1,jfirst:jlast)    ! Pres at layer edges 
   real(r8), intent(inout):: pk(im,jfirst:jlast,kfirst:klast+1)    ! pe**cappa
   real(r8), intent(inout):: pkz(im,jfirst:jlast,kfirst:klast)     ! finite-volume mean of pk

!--------------------------------------------------------------------------------------
! !OUTPUT (input values are not used):
!--------------------------------------------------------------------------------------
   real(r8), intent(out):: te0                   ! Total energy before dynamics
   real(r8), intent(out):: peln(im,kfirst:klast+1,jfirst:jlast) ! log pressure (pe) at layer edges
   real(r8), intent(out):: omga(im,kfirst:klast,jfirst:jlast)   ! vertical pressure velocity (pa/sec)

! The following variables are xy-decomposed instanciations:
   real(r8), intent(out):: uxy(ifirstxy:ilastxy,jfirstxy:jlastxy+1,km)  ! ghosted N1
   real(r8), intent(out):: vxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
   real(r8), intent(out):: psxy(ifirstxy:ilastxy,jfirstxy:jlastxy)
   real(r8), intent(out):: delpxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
   real(r8), intent(out):: pexy(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy)
   real(r8), intent(out):: pkxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1)
   real(r8), intent(out):: pkzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
   real(r8), intent(out):: pelnxy(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy)
   real(r8), intent(out):: omgaxy(ifirstxy:ilastxy,km,jfirstxy:jlastxy)
   real(r8), intent(out):: q3xy(ifirstxy:ilastxy,jfirstxy:jlastxy,km,nc)
   real(r8), intent(out):: ptxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)

! !DESCRIPTION:
!
! Developer: Shian-Jiann Lin, NASA/GSFC; email: lin@dao.gsfc.nasa.gov
!
! Top view of D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!               u(i,j+1)
!                 |
!      v(i,j)---delp(i,j)---v(i+1,j)
!                 |
!               u(i,j)
!
! External routine required: the user needs to supply a subroutine to set up
!                            "Eulerian vertical coordinate" for remapping purpose.
!                             Currently this routine is named as set_eta()
!                             In principle any terrian following vertical
!                             coordinate can be used. The input to fvcore
!                             need not be on the same vertical coordinate
!                             as the output.
!                             If SPMD is defined the Pilgrim communication
!                             library developed by Will Sawyer will be needed.
!
! Remarks: values at poles for both u and v need not be defined; but values for
!          all other scalars needed to be defined at both poles (as polar cap mean
!          quantities). Tracer advection is done "off-line" using the
!          large time step. Consistency is maintained by using the time accumulated
!          Courant numbers and horizontal mass fluxes for the FFSL algorithm.
!          The input "pt" can be either dry potential temperature
!          defined as T/pkz (adiabatic case) or virtual potential temperature
!          defined as T*/pkz (full phys case). IF convt is true, pt is not updated.
!          Instead, virtual temperature is ouput.
!          ipt is updated if convt is false.
!          The user may set the value of nx to optimize the SMP performance
!          The optimal valuse of nx depends on the total number of available
!          shared memory CPUs per node (NS). Assuming the maximm MPI 
!          decomposition is used in the y-direction, set nx=1 if the
!          NS <=4; nx=4 if NS=16.
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Initial SMP version delivered to Will Sawyer
!   WS  99.10.03:  1D MPI completed and tested; 
!   WS  99.10.11:  Additional documentation
!   WS  99.10.28:  benergy and te_map added; arrays pruned
!   SJL 00.01.01:  SMP and MPI enhancements; documentation
!   WS  00.07.13:  Changed PILGRIM API
!   WS  00.08.28:  SPMD instead of MPI_ON
!   AAM 00.08.10:  Add kfirst:klast
!   WS  00.12.19:  phis now distr., LLNL2DModule initialized here
!   WS  01.02.02:  bug fix: parsplit only called for FIRST time
!   WS  01.04.09:  Added initialization of ghost regions
!   WS  01.06.10:  Removed if(first) section; use module
!   AAM 01.06.27:  Extract te_map call into separate routine
!   AAM 01.07.13:  Get rid of dynpkg2; recombine te_map;
!                  perform forward transposes for 2D decomposition
!   WS  01.12.10:  Ghosted PT (changes benergy, cd_core, te_map, hswf)
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local variables

   integer i, j, k, iq          ! Loop indicies
   real(r8) umax                ! Maximum winds, m/s
   parameter (umax = 300.0)

   logical consv                ! Flag to force conservation of toal energy
   parameter (consv = .true.)

   integer    nx          ! # of split pieces in x-direction; for performance, the
   parameter (nx = 4)     ! user may set nx=1 if there is NO shared memory multitasking
   integer ipe, it
   integer nsplit, n, n2
   integer incount, outcount

! Geometric arrays

! Move the following 3D arrays to an initialization routine?
   real(r8), allocatable :: worka(:,:,:),dp0(:,:,:),cx(:,:,:),cy(:,:,:)
   real(r8), allocatable :: mfx(:,:,:), mfy(:,:,:)
   real(r8), allocatable :: delpf(:,:,:), uc(:,:,:), vc(:,:,:)
   real(r8), allocatable :: dwz(:,:,:), pkc(:,:,:), wz(:,:,:)
   real(r8), allocatable :: dpt(:,:,:)
   real(r8), allocatable :: pkcc(:,:,:), wzc(:,:,:)
! The following variables are work arrays for xy=>yz transpose
   real(r8), allocatable :: pkkp(:,:,:), wzkp(:,:,:), pekp(:,:,:)
! The following variables are xy instanciations
   real(r8), allocatable :: mfxxy(:,:,:), dp0xy(:,:,:), wzxy(:,:,:)
! ps3 and psxy3 are dummy 3d variants of ps and psxy, resp.
   real(r8), allocatable :: ps3(:,:,:), psxy3(:,:,:)


   double precision zamda, zam5
   logical first, fill

   integer imh
   real(r8) dt
   real(r8) bdt

   data first /.true./
   data fill  /.true./              ! perform a simple filling algorithm
                                        ! in case negatives were found
! Allocate temporary work arrays
! Change later to use pointers for SMP performance???
! (prime candidates: uc, vc, delpf)

      allocate( worka(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   dp0(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   mfx(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   mfy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate(    cx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    cy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate( delpf(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    uc(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    vc(im,jfirst-2:   jlast+2,   kfirst:klast) )
      allocate(   dpt(im,jfirst-1:   jlast+1,   kfirst:klast) )
      allocate(   dwz(im,jfirst-1:    jlast,    kfirst:klast+1) )
      allocate(   pkc(im,jfirst-1:   jlast+1,   kfirst:klast+1) ) 
      allocate(    wz(im,jfirst-1:   jlast+1,   kfirst:klast+1) )
      allocate(  pkcc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(   wzc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(pkkp(im,jfirst:jlast,kfirst:klastp))
      allocate(wzkp(im,jfirst:jlast,kfirst:klastp))
      allocate(pekp(im,kfirst:klastp,jfirst:jlast))
      allocate(wzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1))
      allocate(mfxxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(dp0xy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(ps3(im,jfirst:jlast,kfirst:klast))
      allocate(psxy3(ifirstxy:ilastxy,jfirstxy:jlastxy,km))

! First touch pkc and wz??? (bufferpack is multitask in vertical but geop
! computations are parallel in j-loop)

   if ( km > 1 ) then         ! not shallow water equations

      if( consv ) then
! Compute globally integrated Total Energy (te0)

         call t_startf('benergy')
         call benergy(im,     jm,     km,     u,      v,         &
                      pt,     ng_d,   ng_s,   delp,   pe,        &
                      pk,     pkz,    phis,   cp,     te0,       &
                      mfx,    dp0,    jfirst, jlast,  kfirst,    &
                      klast,  klastp)
         call t_stopf('benergy') 

      endif

   endif

! Determine splitting
      bdt = ndt

! Second level splitting

      n2 = max ( 1, ns/4 )
      nsplit = (ns+n2-1) / n2

      dt = bdt / float(nsplit*n2)

  do 2000 n=1, n2

   if( nq > 0 ) then

!$omp parallel do private(i, j, k)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=1,im
! Save initial delp field before the small-time-step
! Initialize the CFL number accumulators: (cx, cy)
! Initialize total mass fluxes: (mfx, mfy)
               dp0(i,j,k) = delp(i,j,k)
                cx(i,j,k) = 0.
                cy(i,j,k) = 0.
               mfx(i,j,k) = 0.
               mfy(i,j,k) = 0.
            enddo
         enddo
      enddo

   endif

   do it=1, nsplit

      if(it == nsplit .and. n == n2) then
         ipe = 1                     ! end of fvcore; output pe for te_map
      elseif(it == 1 .and. n == 1) then
         ipe = -1                    ! start of cd_core
      else
         ipe = 0
      endif

! Call the Lagrangian dynamical core using small tme step

      call t_startf('cd_core')

      call cd_core(im,     jm,    km,     nq,     nx,              &
                   jfirst, jlast, kfirst, klast,                   &
                   klastp,  u,    v,      pt,     delp,            &
                   pe,      pk,   dt,     ptop,   umax,            &
                   ae,      rcap, cp,     cappa,  icd,             &
                   jcd,     iord, jord,   ng_c,   ng_d,            &
                   ng_s,    ipe,  om,     phis,                    &
                   cx,      cy,   mfx,    mfy,    delpf,           &
                   uc,      vc,   pkz,    dpt,    worka,           &
                   dwz, pkc, wz, phisxy, ptxy, pkxy,               &
                   pexy, pkcc, wzc, wzxy, delpxy, pkkp, wzkp,      &
                   pekp, ifirstxy, ilastxy, jfirstxy, jlastxy)

      call t_stopf('cd_core')

   enddo

   if( nq .ne. 0 ) then

! Perform large-tme-step scalar transport using the accumulated CFL and
! mass fluxes
 
      call t_startf('trac2d')
      call trac2d( dp0,    q3,     nc,     nq,     cx,             &
                   cy,     mfx,    mfy,    iord,   jord,           &
                   ng_d,   fill,   im,     jm,     km,             &
                   jfirst, jlast,  kfirst, klast,  pkz,            &
                   worka  )
      call t_stopf('trac2d')

   endif

2000  continue


#if defined (SPMD)
   if (twod_decomp .eq. 1) then
!
! Transpose ps, u, v, and q3 from yz to xy decomposition
!
! Note: pt, pe and pk will have already been transposed through
! call to geopk in cd_core. geopk does not actually require
! secondary xy decomposition; direct 16-byte technique works just
! as well, perhaps better. However, transpose method is used on last
! call to avoid having to compute these three transposes now.
!
      call t_startf('transpose_fwd')

! Embed ps in 3D array, per requirement of Pilgrim

!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               mfx(i,j,k) = ps(i,j)
            enddo
         enddo
      enddo
      call redistributestart (inter_ijk, .true., mfx)
!
! TEMPORARY!
!
!$omp parallel do private(i,j,k,iq)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               yzt(i,j,k) = u(i,j,k)
            enddo
         enddo
      enddo
      call redistributefinish(inter_ijk, .true., mfxxy)

!$omp parallel do private(i,j)
      do j = jfirstxy,jlastxy
         do i = ifirstxy,ilastxy
            psxy(i,j) = mfxxy(i,j,1)
         enddo
      enddo

      call redistributestart (inter_ijk, .true., yzt) ! send U
!
! TEMPORARY!
!
!$omp parallel do private(i,j,k,iq)
      do iq = 1,nc
         do k = kfirst,klast
            do j = jfirst,jlast
               do i = 1,im
                  q3t(i,j,k,iq) = q3(i,j,k,iq)
               enddo
            enddo
         enddo
      enddo
      call redistributefinish(inter_ijk, .true., xyt) ! recv UXY

      call redistributestart (inter_q3, .true., q3t)
!
! TEMPORARY!
!
!$omp parallel do private(i,j,k)
      do k = 1,km
         do j = jfirstxy,jlastxy
            do i = ifirstxy,ilastxy
               uxy(i,j,k) = xyt(i,j,k)
            enddo
         enddo
      enddo
!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               yzt(i,j,k) = v(i,j,k)
            enddo
         enddo
      enddo
      call redistributefinish(inter_q3, .true., q3xy)

      call redistributestart (inter_ijk, .true., yzt)  ! send V
      call redistributefinish(inter_ijk, .true., vxy)  ! recv VXY

      call t_stopf('transpose_fwd')
    endif
#endif

    if ( km > 1 ) then           ! not shallow water equations

! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.

      call t_startf('te_map')
      if (twod_decomp .eq. 1) then
! 
! te_map requires uxy, vxy, psxy, pexy, pkxy, phisxy, q3xy, and ptxy
!
         call te_map(consv,  convt,  psxy,  omgaxy, pexy,            &
                     delpxy, pkzxy,  pkxy,  ndt,    im,              &
                     jm,     km,     nx,    jfirstxy, jlastxy,       &
                     0,      0,      1,     0,        0,             &
                     ifirstxy, ilastxy,              &
                     nq,     uxy,    vxy,   ptxy,   q3xy,            &
                     phisxy, cp,     cappa, kord,   pelnxy,          &
                     te0,    mfxxy,  dp0xy, nc )
!
! te_map computes uxy, vxy, psxy, delpxy, pexy, pkxy, pkzxy,
! pelnxy, omgaxy, q3xy and ptxy.
!
      else
         call te_map(consv,  convt,  ps,    omga,   pe,              &
                     delp,   pkz,    pk,    ndt,    im,              &
                     jm,     km,     nx,    jfirst, jlast,           &
                     ng_d,   ng_d,   ng_s,  ng_s,   ng_d,            &
                     1,      im,                      &
                     nq,     u,      v,     pt,     q3,              &
                     phis,   cp,     cappa, kord,   peln,            &
                     te0,    mfx,    dp0,   nc )
      endif
      call t_stopf('te_map')

    endif

#if defined( SPMD )
    if (twod_decomp .eq. 1) then

       call t_startf('transpose_bck1')

       if ( .not. convt ) then
!
! Transpose delpxy to delp for simplified physics (for full_phys,
! delp is recomputed after physics advance)
!
          call redistributestart (inter_ijk, .false., delpxy)
          call redistributefinish(inter_ijk, .false., delp)
       endif

!
! Transpose pexy into pekp, then embed in pe and perform boundary update
! (pexy is not needed for physics update)
!

       call redistributestart (inter_ikjp, .false., pexy)
       call redistributefinish(inter_ikjp, .false., pekp)

!$omp parallel do private(i,j,k)
       do j = jfirst,jlast
          do k = kfirst,klastp
             do i = 1,im
                pe(i,k,j) = pekp(i,k,j)
             enddo
          enddo
       enddo

       if (npr_z > 1) then
! Send one level to the PE below
          call mp_send3d( iam-npr_y, iam+npr_y, im, km+1, jm,            &
                          1, im, kfirst, klast+1, jfirst, jlast,         &
                          1, im, kfirst, kfirst, jfirst, jlast, pe )
          call mp_recv3d( iam+npr_y, im, km+1, jm,                       &
                          1, im, kfirst, klast+1, jfirst, jlast,         &
                          1, im, klast+1, klast+1, jfirst, jlast, pe )
       endif
!
! Transpose psxy into ps, using 3D temporary arrays
! (psxy is not needed for physics update)
!
       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                psxy3(i,j,k) = psxy(i,j)
             enddo
          enddo
       enddo

       call redistributestart (inter_ijk, .false., psxy3)
       call redistributefinish(inter_ijk, .false., ps3)

       do j=jfirst,jlast
          do i=1,im
             ps(i,j) = ps3(i,j,kfirst)
          enddo
       enddo

       call t_stopf('transpose_bck1')

    endif
#endif

    deallocate( mfy )
    deallocate( mfx )
    deallocate(  cy )
    deallocate(  cx )
    deallocate( dp0 )
    deallocate( delpf )
    deallocate( uc    )
    deallocate( vc    )
    deallocate( dpt   )
    deallocate( pkc   )
    deallocate( dwz   )
    deallocate(  wz   )
    deallocate( worka )
    deallocate( pkcc )
    deallocate( wzc )
    deallocate( pkkp )
    deallocate( wzkp )
    deallocate( pekp )
    deallocate( wzxy )
    deallocate( mfxxy )
    deallocate( dp0xy )
    deallocate( ps3 )
    deallocate( psxy3 )

!----------------------------------------------------------
! Idealized physics: do Held-Suarez-Williamson-Lin forcing.
! Since actual variable names depend on whether we are using
! 2D decomposition, branching is required.
!----------------------------------------------------------
    if (ideal) then
       call t_startf('ideal_phys')
       if (twod_decomp .eq. 1) then
!--------------------------------------------------------------------------
! For 2D decomposition, hswf requires u3sxy, v3sxy, ptxy, pexy and 
! pkzxy, and computes u3sxy, v3sxy and ptxy.
!--------------------------------------------------------------------------
          call hswf( im,   jm,   km,   jfirstxy, jlastxy,         &
                     ifirstxy,         ilastxy,                   &
                     uxy,  vxy,  ptxy, 0, 0, 1, 0, 0,             &
                     pexy,     pkzxy,           &
                     ndt,  cappa,      gravit,   rair,    dcaf,   &
                    .true.,      rayf, sinp,     cosp,    sine,   &
                     cose, coslon,     sinlon )
       else
          call hswf( im,   jm,   km,   jfirst,   jlast,           &
                     1,    im,   u,    v,        pt,              &
                     ng_d, ng_d, ng_s, ng_s,     ng_d,            &
                     pe,   pkz,                                   &
                     ndt,  cappa,      gravit,   rair,    dcaf,   &
                    .true.,      rayf, sinp,     cosp,    sine,   &
                     cose, coslon,     sinlon )
       endif
       call t_stopf('ideal_phys')
    endif

!EOC
end subroutine dynpkg
!-----------------------------------------------------------------------
