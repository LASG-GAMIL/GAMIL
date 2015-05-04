#include <misc.h>
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: cd_core --- Dynamical core for both C- and D-grid Lagrangian
!                       dynamics
!
! !INTERFACE:
 subroutine cd_core(im,   jm,   km,    nq,     nx,                 &
                    jfirst, jlast , kfirst, klast, klastp, u,      &
                    v,   pt,   delp,   pe,     pk,                 &
                    dt,   ptopin, umax,   ae,     rcap,            &
                    cp,   akap, iord_c, jord_c, iord_d, jord_d,    &
                    ng_c, ng_d, ng_s, ipe,    om,     hs,          &
                    cx3  ,  cy3, mfx, mfy,                         &
                    delpf, uc, vc, ptc, dpt, ptk,                  &
                    wz3, pkc, wz,  hsxy, ptxy, pkxy,               &
                    pexy, pkcc, wzc, wzxy, delpxy, pkkp, wzkp,     &
                    pekp, ifirstxy, ilastxy, jfirstxy, jlastxy )

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use sw_core, only : c_sw, d_sw
   use pft_module, only : fftfax, pft2d, pft_cf
   use dynamics_vars, only : sinp, cosp, cose, acosp,              &
                       sinlon, coslon, cosl5, sinl5,               &
                       ak, bk, ptop, ns, yzt
   use pmgrid, only: twod_decomp, myid_y, myid_z, npr_y, npr_z, iam

#if defined( SPMD )
   use redistributemodule, only: redistributestart, redistributefinish
   use spmd_dyn, only: comm_y, comm_z, inter_ijk, inter_ijkp
   use mod_comm, only: mp_send4d_ns, mp_recv4d_ns,                 &
                       mp_send2_ns, mp_recv2_ns,                   &
                       mp_send3d_2, mp_recv3d_2,                   &
                       mp_send3d, mp_recv3d

#define CPP_PRT_PREFIX  if(iam==0)
#else
#define CPP_PRT_PREFIX
#endif

   implicit none

! !INPUT PARAMETERS:

! Input paraterers:
  integer im, jm, km
  integer nx                 ! # of split pieces in longitude direction
  integer jfirst
  integer jlast
  integer kfirst
  integer klast
  integer klastp             ! klast, except km+1 when klast=km
  integer ifirstxy,ilastxy   ! xy-decomp. longitude ranges
  integer jfirstxy,jlastxy   ! xy-decomp. latitude ranges
  integer ipe                ! ipe=1:  end of cd_core()
                             ! ipe=-1: start of cd_core()
                             ! ipe=0 :
  integer nq                 ! # of tracers to be advected by trac2d
  integer iord_c, jord_c     ! scheme order on C grid in X and Y dir.
  integer iord_d, jord_d     ! scheme order on D grid in X and Y dir.
  integer ng_c               ! ghost latitudes on C grid
  integer ng_d               ! ghost lats on D (Max NS dependencies, ng_d >= ng_c)
  integer ng_s               ! max(ng_c+1,ng_d) significant if ng_c = ng_d

  real(r8) ae                ! Radius of the Earth (m)
  real(r8) om                ! rotation rate
  real(r8) ptopin
  real(r8) umax
  real(r8) dt            !small time step in seconds
  real(r8) rcap
  real(r8) cp
  real(r8) akap

! Input time independent arrays:
  real(r8) hs(im,jfirst:jlast)         !surface geopotential
  real(r8) hsxy(ifirstxy:ilastxy,jfirstxy:jlastxy) !surface geopotential XY-decomp.

! !INPUT/OUTPUT PARAMETERS:
! Prognostic variables:
  real(r8)   u(im,jfirst-ng_d:jlast+ng_s,kfirst:klast) ! u-Wind (m/s)
  real(r8)   v(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) ! v-Wind (m/s)
  real(r8) delp(im,jfirst:jlast,kfirst:klast)          ! Delta pressure (pascal)
  real(r8)   pt(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)! Scaled-Potential temperature

! Input/output: accumulated winds & mass fluxes on c-grid for large-
!               time-step transport
  real(r8) cx3(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)! Accumulated Courant number in X
  real(r8) cy3(im,jfirst:jlast+1,kfirst:klast)        ! Accumulated Courant number in Y
  real(r8) mfx(im,jfirst:jlast,kfirst:klast)          ! Mass flux in X  (unghosted)
  real(r8) mfy(im,jfirst:jlast+1,kfirst:klast)        ! Mass flux in Y

! !OUTPUT PARAMETERS:
  real(r8) pe(im,kfirst:klast+1,jfirst:jlast)         ! Edge pressure (pascal)
  real(r8) pk(im,jfirst:jlast,kfirst:klast+1)         ! Pressure to the kappa
  real(r8) ptxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) ! Potential temperature XY decomp
  real(r8) pkxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) ! P-to-the-kappa XY decomp
  real(r8) pexy(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy) ! Edge pressure XY decomp

! Input work arrays:
  real(r8) delpf(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)  ! filtered delp
  real(r8)   uc(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)   ! u-Winds on C-grid
  real(r8)   vc(im,jfirst-2:   jlast+2,   kfirst:klast)   ! v-Winds on C-grid

  real(r8) ptc(im,jfirst:jlast,kfirst:klast)
  real(r8) ptk(im,jfirst:jlast,kfirst:klast)
  real(r8) dpt(im,jfirst-1:jlast+1,kfirst:klast)
  real(r8) wz3(im,jfirst-1:jlast  ,kfirst:klast+1)
  real(r8) pkc(im,jfirst-1:jlast+1,kfirst:klast+1) 
  real(r8)  wz(im,jfirst-1:jlast+1,kfirst:klast+1)
  real(r8) pkcc(im,jfirst:jlast,kfirst:klast+1) 
  real(r8) wzc(im,jfirst:jlast,kfirst:klast+1) 
  real(r8) wzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1)
  real(r8) delpxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
  real(r8) pkkp(im,jfirst:jlast,kfirst:klastp)
  real(r8) wzkp(im,jfirst:jlast,kfirst:klastp)
  real(r8) pekp(im,kfirst:klastp,jfirst:jlast)

#if defined ( HIGH_P )
  real(r8) wzz(im,jfirst:jlast+1,kfirst:klast+1)
#endif

! ! !DESCRIPTION:
!    Perform a dynamical update for one small time step; the small
!    time step is limitted by the fastest wave within the Lagrangian control-
!    volume 
!
! !REVISION HISTORY:
!     SJL  99.01.01:   Original SMP version
!     WS   99.04.13:   Added jfirst:jlast concept
!     SJL  99.07.15:   Merged c_core and d_core to this routine
!     WS   99.09.07:   Restructuring, cleaning, documentation
!     WS   99.10.18:   Walkthrough corrections; frozen for 1.0.7
!     WS   99.11.23:   Pruning of some 2-D arrays
!     SJL  99.12.23:   More comments; general optimization; reduction
!                      of redundant computation & communication
!     WS   00.05.14:   Modified ghost indices per Kevin's definition
!     WS   00.07.13:   Changed PILGRIM API
!     WS   00.08.28:   Cosmetic changes: removed old loop limit comments
!     AAM  00.08.30:   Introduced kfirst,klast
!     WS   00.12.01:   Replaced MPI_ON with SPMD; hs now distributed
!     WS   01.04.11:   PILGRIM optimizations for begin/endtransfer
!     WS   01.05.08:   Optimizations in the call of c_sw and d_sw
!     AAM  01.06.27:   Reinstituted 2D decomposition for use in ccm
!     WS   01.12.10:   Ghosted PT, code now uses mod_comm primitives
!     WS   01.12.31:   Removed vorticity damping, ghosted U,V,PT
!     WS   02.01.15:   Completed transition to mod_comm
!     WS   02.07.04:   Fixed 2D decomposition bug dest/src for mp_send3d
!     WS   02.09.04:   Integrated fvgcm-1_3_71 zero diff. changes by Lin
!
!EOP
!---------------------------------------------------------------------
!BOC
! Local 2D arrays:
      real(r8)   wk(im,jfirst:  jlast+2)
      real(r8)  wk1(im,jfirst-1:jlast+1)
      real(r8)  wk2(im,jfirst-ng_d:jlast+ng_d)
      real(r8)  wk3(im,jfirst-1:jlast+1)

      real(r8) p1d(im)

! Local scalars
      real(r8) ptmp, pint, press
      real(r8)  rat, ycrit
      real(r8)  dt0, dt5
      real(r8)  dl, dp

      integer i, j, k
      integer ks
      integer js1g1, js2g0, js2g1, js2gc
      integer jn2g0, jn1g1, jn1gc
      integer iord , jord
      integer ktot, ktotp

      real(r8)  pi
      real(r8)  zt_c   
      real(r8)  zt_d  
      real(r8)  tau, fac, pk4
      real(r8)  tiny
      parameter (tiny = 1.e-10)

! Declare permanent local arrays
      integer  ifax(13)                      !ECMWF fft
      real(r8), allocatable, save :: trigs(:)
      real(r8), allocatable, save :: fc(:), f0(:)
      real(r8), allocatable, save :: dc(:,:), de(:,:), sc(:), se(:)
      real(r8), allocatable, save :: cdx(:,:), cdy(:,:)
      real(r8), allocatable, save :: dtdx(:), dtdxe(:), txe5(:), dtxe5(:)
      real(r8), allocatable, save :: dyce(:),   dx(:) ,  rdx(:),    cy(:)
      real(r8), allocatable, save :: dtdx2(:), dtdx4(:),  dxdt(:), dxe(:)
      real(r8), allocatable, save :: cye(:),    dycp(:),  rdxe(:)

      real(r8) rdy, dtdy, dydt, dtdy5, tdy5
      data  dt0 / 0./

      save ifax
      save dtdy, dydt, dtdy5, tdy5, rdy
      save dl, dp
      save zt_c, zt_d

#if defined( SPMD )
      integer dest, src
#endif

!******************************************************************
!******************************************************************
!
! IMPORTANT CODE OPTIONS - SEE BELOW
!
!******************************************************************
!******************************************************************

! Option for which version of geopk to use with yz decomposition.
! Version geopk16 (geopk16byte =true) uses 16-byte reals to preserve accuracy
!   through round-off (order of summation varies with z decomposition).
!   This version avoids transposes and instead involves semi-global
!   communication in Z.
! If geopk16byte=false, variables are transposed to/from xy decomposition
!   for use in geopk.
! On last small timestep (ipe=1) for D-grid, the version of geopk that uses transposes
!   is called regardless, as some transposed quantities are required for
!   the te_map phase.

      logical geopk16byte
!     data geopk16byte / .true. /
      data geopk16byte / .false. /
      logical geopkc16, geopkd16

      geopkc16 = .false.
      geopkd16 = .false.
      if (geopk16byte) then
        geopkc16 = .true.
        if (ipe /= 1) geopkd16 = .true.
      endif

!******************************************************************

      ktot = klast - kfirst + 1
      ktotp = ktot + 1

#if defined( SPMD )
      call t_startf('send_uv')
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,       &
                         kfirst, klast, ng_s, ng_d, u )
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,       &
                         kfirst, klast, ng_d, ng_s, v )
      call t_stopf('send_uv')
#endif

! Set general loop limits
! jfirst >= 1; jlast <= jm
      js1g1  = max(1,jfirst-1)
      js2g0  = max(2,jfirst)
      js2g1  = max(2,jfirst-1)
      jn2g0  = min(jm-1,jlast)
      jn1g1  = min(jm,jlast+1)

! Construct C-grid dependent loop limits
      js2gc  = max(2,jfirst-ng_c)     ! NG latitudes on S (starting at 2)
!
! WS: 00.04.13 : An ugly hack to handle JORD>1 JCD=1
!
      if ( ng_c == 1 .AND. ng_d > 1 ) THEN
        js2gc  = max(2,jfirst-2) 
      endif
      jn1gc  = min(jm,jlast+ng_c)     ! ng_c latitudes on N (ending at jm)

      if( abs(dt0-dt) > 0.1 ) then
        if ( .not. allocated( dtdx ) ) then
          allocate(dtdx(jm),dtdx2(jm), dtdx4(jm), dtdxe(jm), dxdt(jm), &
                   dxe(jm),  cye(jm),dycp(jm),rdxe(jm),                &
                   txe5(jm), dtxe5(jm),dyce(jm),                       &
                   dx(jm),rdx(jm),cy(jm) )

          allocate( trigs(3*im/2+1) )
          allocate( sc(js2g0:jn2g0),    se(js2g0:jn1g1)    )
          allocate( dc(im,js2g0:jn2g0), de(im,js2g0:jn1g1) )

          call fftfax(im, ifax, trigs)

! Determine ycrit such that effective DX >= DY
          pi  = 4.d0 * datan(1.d0)
          rat = float(im)/float(2*(jm-1))
          ycrit = acos( min(0.81, rat) ) * (180./pi)

          call pft_cf(im, jm, js2g0, jn2g0, jn1g1, sc, se, dc, de,  &
                      cosp, cose, ycrit)

          allocate( cdx(js2g0:jn1g1,kfirst:klast) )
          allocate( cdy(js2g0:jn1g1,kfirst:klast) )

          allocate( f0(jfirst-ng_s:jlast+ng_d) ) ! 000304 bug fix: ng_s not ng_d
          allocate( fc(js2gc:jn1gc) )

          do j=max(1,jfirst-ng_s),min(jm,jlast+ng_d)        ! 000304 bug fix
             f0(j) = (om+om)*sinp(j)
          enddo

! Compute coriolis parameter at cell corners.
          do j=js2gc, jn1gc                    ! Not the issue with ng_c = ng_d 
             fc(j) = 0.5*(f0(j) + f0(j-1))
          enddo

        endif

        dt0 = dt
        dt5 = 0.5*dt

        pi  = 4.d0 * datan(1.d0)
        dl  = (pi+pi)/im
        dp  = pi/(jm-1)

        rdy   = 1./(ae*dp)
        dtdy  = dt *rdy
        dtdy5 = dt5*rdy
        dydt  = (ae*dp) / dt
        tdy5  = 0.5/dtdy

        do j=2,jm-1
          dx(j)    = dl*ae*cosp(j)
          rdx(j)   = 1. / dx(j)
          dtdx(j)  = dt / dx(j)
          dxdt(j)  = dx(j) / dt
          dtdx2(j) = 0.5*dtdx(j)
          dtdx4(j) = 0.5*dtdx2(j)
          dycp(j)  = ae*dp/cosp(j)
          cy(j)    =  rdy * acosp(j)
        enddo

        do j=2,jm
          dxe(j)   = ae*dl*cose(j)
          rdxe(j)  = 1. / dxe(j)
          dtdxe(j) = dt / dxe(j)
          dtxe5(j) = 0.5*dtdxe(j)
          txe5(j)  = 0.5/dtdxe(j)
           cye(j)  =  1. / (ae*cose(j)*dp)
          dyce(j)  = ae*dp/cose(j)
        enddo

! C-grid
        zt_c = abs(umax*dt5) / (dl*ae)
        CPP_PRT_PREFIX write(6,*) 'c-core: ', (180./pi)*acos(zt_c)
! D-grid
        zt_d = abs(umax*dt) / (dl*ae)
        CPP_PRT_PREFIX write(6,*) 'd-coret: ', (180./pi)*acos(zt_d)

        if ( ptopin /= ptop) then
             write(6,*) 'PTOP as input to cd_core != ptop from dynamics_vars'
             stop
        endif

!-----------------------------------------
! Divergence damping coeff. dx(2)*dy/(TAU)
!-----------------------------------------
          do k=kfirst,klast
             press = 0.5 * ( ak(k)+ak(k+1) + (bk(k)+bk(k+1))*1.E5 )
             tau = 8. * (1.+ tanh(1.0*log(ptop/press)) )
             tau = max(1., tau) / (128.*abs(dt))
            do j=js2g0,jn1g1
               fac = tau * ae / cose(j)
               cdx(j,k) = fac*dp
               cdy(j,k) = fac*dl
            enddo
          enddo
      endif

      if ( ipe == -1 .or. ns == 1 ) then          ! starting cd_core
!$omp parallel do private(i, j, k, wk, wk2)
         do k=kfirst,klast
            do j=jfirst,jlast
               do i=1,im
                  delpf(i,j,k) = delp(i,j,k)
               enddo
            enddo
            call pft2d( delpf(1,js2g0,k), sc(js2g0), dc(1,js2g0),  &
                        im, jn2g0-js2g0+1, ifax, trigs, wk, wk2 )
         enddo
      endif

#if defined( SPMD )
      call t_startf('recv_uv')
      call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,       &
                         kfirst, klast, ng_s, ng_d, u )
      call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,       &
                         kfirst, klast, ng_d, ng_s, v )
      call t_stopf('recv_uv')

      call t_startf('ghost_pt_delpf')
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,        &
                         kfirst, klast, ng_d, ng_d, pt )
      if ( ipe == -1 .or. ns == 1 ) then          ! starting cd_core
         call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,       &
                            kfirst, klast, ng_d, ng_d, delpf )
      endif                         ! end if ipe = -1 check
      call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,          &
                         kfirst, klast, ng_d, ng_d, pt )
      if ( ipe == -1 .or. ns == 1 ) then          ! starting cd_core
         call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,       &
                            kfirst, klast, ng_d, ng_d, delpf )
      endif                          ! end if ipe = -1 check
      call t_stopf('ghost_pt_delpf')
#endif

  call t_startf('c_core')

!$omp parallel do private(i, j, k, iord, jord)    

  do  k=kfirst,klast     ! This is the main parallel loop.

     if ( k <= km/8 ) then
         iord = 1
         jord = 1
     else
         iord = iord_c
         jord = jord_c
     endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the C-grid
!-----------------------------------------------------------------

        call c_sw( u(1,jfirst-ng_d,k),   v(1,jfirst-ng_s,k),         &
                   pt(1,jfirst-ng_d,k),  delp(1,jfirst,k),           &
                   uc(1,jfirst-ng_d,k),  vc(1,jfirst-2,k),           &
                   ptc(1,jfirst,k),      delpf(1,jfirst-ng_d,k),     &
                   ptk(1,jfirst,k),                                  &
                   cosp,   acosp,   cose,   coslon,   sinlon,        &
                   dxdt,   dxe,     dtdx2,  dtdx4,    dtxe5,   rdxe, &
                   dycp,   dydt,    dtdy5,  cye,      fc,            &
                   ifax,   trigs,   dc(1,js2g0),      sc,            &
                   zt_c,   tiny,    rcap,   im,       jm,            &     
                   jfirst, jlast,   ng_c,   ng_d,     ng_s,          &
                   js2g0,  jn2g0,   js2gc,  jn1gc,                   &
                   iord,   jord,    cosl5,  sinl5 )
  enddo

  call t_stopf('c_core')


! MPI note: uc, vc, ptk, and ptc computed within the above k-look from jfirst to jlast
! Needed by D-core: uc(jfirst-ng_d:jlast+ng_d), vc(jfirst:jlast+1)


  call t_startf('c_geop')
#if defined( SPMD )

      if (twod_decomp == 1) then

       if (geopkc16) then

!
! Stay in yz space and use semi-global z communications and 16-byte reals
!

        call geopk16(ptop, pe, ptk, pkcc, wzc, hs, ptc, 0, im, jm, km,    &
                     jfirst, jlast, 1, im, cp, akap, kfirst, klast)

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
        do k = kfirst, klast+1
          do j = jfirst, jlast
            do i = 1, im
              pkc(i,j,k) = pkcc(i,j,k)
              wz(i,j,k) = wzc(i,j,k)
            enddo
          enddo
        enddo

       else

!
! Transpose to xy decomposition
!

        call redistributestart(inter_ijk, .true., ptk)
        call redistributefinish(inter_ijk, .true., delpxy)
        call redistributestart(inter_ijk, .true., ptc)
        call redistributefinish(inter_ijk, .true.,  ptxy)

        call geopk(ptop, pexy, delpxy, pkxy, wzxy, hsxy, ptxy, 0, im, jm, km,  &
                   jfirstxy, jlastxy, ifirstxy, ilastxy,                  &
                   cp, akap, nx, 0)

!
! Transpose back to yz decomposition.
! delpxy, ptxy and pexy are not output quantities on this call.
! pkkp and wzkp are holding arrays, whose specific z-dimensions
!    are required by Pilgrim.
!

        call redistributestart(inter_ijkp, .false., pkxy)
        call redistributefinish(inter_ijkp, .false., pkkp)
        call redistributestart(inter_ijkp, .false., wzxy)
        call redistributefinish(inter_ijkp, .false., wzkp)

!$omp parallel do private(i, j, k)
        do k = kfirst, klastp
           do j = jfirst, jlast
              do i = 1, im
                 pkc(i,j,k) = pkkp(i,j,k)
                 wz(i,j,k) = wzkp(i,j,k)
              enddo
           enddo
        enddo

        if (npr_z > 1) then
!
! Fill in klast+1
!
         call mp_send3d_2( iam-npr_y, iam+npr_y, im, jm, km+1,             &
                           1, im, jfirst-1, jlast+1, kfirst, klast+1,      &
                           1, im, jfirst, jlast, kfirst, kfirst, pkc, wz )
         call mp_recv3d_2( iam+npr_y, im, jm, km+1,                        &
                           1, im, jfirst-1, jlast+1, kfirst, klast+1,      &
                           1, im, jfirst, jlast, klast+1, klast+1, pkc, wz)
        endif       ! npr_z > 1

       endif       ! geopkc16

      else
#endif

        if (geopkc16) then

!
! Use 16-byte reals (for compatibility with 2-D decomposition)
!

          call geopk16(ptop, pe, ptk, pkcc, wzc, hs, ptc, 0, im, jm, km,         &
                       jfirst, jlast, 1, im, cp, akap, 1, km)

        else

!
! Use 8-byte reals (standard)
!

          call geopk(ptop, pe, ptk, pkcc, wzc, hs, ptc, 0, im, jm, km,           &
                     jfirst, jlast, 1, im,                                    &
                     cp, akap, nx, 0)

        endif

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
       do k = kfirst, klast+1
          do j = jfirst, jlast
            do i = 1, im
               pkc(i,j,k) = pkcc(i,j,k)
               wz(i,j,k) = wzc(i,j,k)
            enddo
          enddo
       enddo

#if defined( SPMD )
      endif

#endif
   call t_stopf('c_geop')


! Upon exit from geopk, the quantities pe, pkc and wz will have been
! updated at klast+1


#if defined( SPMD )
      call t_startf('send_pkc_wz')
!
! pkc & wz need to be ghosted only at jfirst-1
!
      dest = iam+1
      src  = iam-1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      if ( mod(iam,npr_y) == 0 ) src = -1
      call mp_send3d_2( dest, src, im, jm, km+1,                   &
                        1, im, jfirst-1, jlast+1, kfirst, klast+1,    &
                        1, im, jlast, jlast, kfirst, klast+1, pkc, wz)
      call t_stopf('send_pkc_wz')
#endif


#if defined( HIGH_P )
      call t_startf('highp2')
      call highp2(pkc,  wz,  wz3,    wzz,   dpt,          &
                  im,   jm,   km,  jfirst, jlast,         &
                  kfirst, klast, klastp, nx )
      call t_stopf('highp2')

     call t_startf('c_u_pgrad')
!$omp parallel do private(i, j, k, p1d, wk, wk2)
     do k=kfirst,klast
         do j=js2g0,jn2g0
            do i=1,im
               p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
            enddo

            uc(1,j,k) = uc(1,j,k) + dtdx2(j) / (p1d(1)+p1d(im)) *         &
                (dpt(im,j,k)-dpt(1,j,k)-wz3(1,j,k)+wz3(1,j,k+1))
            do i=2,im
               uc(i,j,k) = uc(i,j,k) + dtdx2(j) / (p1d(i)+p1d(i-1)) *     &
                 (dpt(i-1,j,k)-dpt(i,j,k)-wz3(i,j,k)+wz3(i,j,k+1))
            enddo
         enddo          ! j-loop

         call pft2d(uc(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,        &
                    jn2g0-js2g0+1, ifax, trigs, wk, wk2 )
     enddo
     call t_stopf('c_u_pgrad')
#else

   call t_startf('c_u_loop')
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, p1d, wk, wk2)

   do k=kfirst,klast
      do j=js2g0,jn2g0
         do i=1,im
            p1d(i) = pkc(i,j,k+1) - pkc(i,j,k)
         enddo

         uc(1,j,k) = uc(1,j,k) + dtdx2(j) * (                        &
                (wz(im,j,k+1)-wz(1,j,k))*(pkc(1,j,k+1)-pkc(im,j,k))  &
              + (wz(im,j,k)-wz(1,j,k+1))*(pkc(im,j,k+1)-pkc(1,j,k))) &
                                 / (p1d(1)+p1d(im))
         do i=2,im
            uc(i,j,k) = uc(i,j,k) + dtdx2(j) * (                     &
              (wz(i-1,j,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i-1,j,k))  &
            + (wz(i-1,j,k)-wz(i,j,k+1))*(pkc(i-1,j,k+1)-pkc(i,j,k))) &
                                 / (p1d(i)+p1d(i-1))
         enddo

      enddo
      call pft2d(uc(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,     &
                 jn2g0-js2g0+1, ifax, trigs, wk, wk2 )
   enddo 
   call t_stopf('c_u_loop')
#endif

#if defined( SPMD )
      call t_startf('recv_pkc_wz')
      call mp_recv3d_2( src, im, jm, km+1,                          &
                        1, im, jfirst-1, jlast+1, kfirst, klast+1,    &
                        1, im, jfirst-1, jfirst-1, kfirst, klast+1, pkc, wz)
      call t_stopf('recv_pkc_wz')

      call t_startf('send_uc')
      call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,                &
                         kfirst, klast, ng_d, ng_d, uc )
      call t_stopf('send_uc')
#endif

  call t_startf('c_v_pgrad')
!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, wk, wk1 )

! pkc and wz need only to be ghosted jfirst-1

  do k=kfirst,klast
     do j=js1g1,jlast
        do i=1,im
           wk1(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
        enddo
     enddo

     do j=js2g0,jlast
        do i=1,im
           vc(i,j,k) = vc(i,j,k) + dtdy5/(wk1(i,j)+wk1(i,j-1)) *  &
#if defined ( HIGH_P )
           ( dpt(i,j-1,k)-dpt(i,j,k)-wzz(i,j,k)+wzz(i,j,k+1) )
#else
           ( (wz(i,j-1,k+1)-wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j-1,k))  &
          +  (wz(i,j-1,k)-wz(i,j,k+1))*(pkc(i,j-1,k+1)-pkc(i,j,k)) )
#endif
       enddo
     enddo

     call pft2d(vc(1,js2g0,k), se(js2g0), de(1,js2g0), im,    &
                jlast-js2g0+1, ifax, trigs, wk, wk1 )
  enddo

  call t_stopf('c_v_pgrad')

#if defined( SPMD )
      call t_startf('recv_uc')
      call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,               &
                         kfirst, klast, ng_d, ng_d, uc )
      call t_stopf('recv_uc')

      call t_startf('ghost_vc')
! vc only needs to be ghosted at jlast+1
      dest = iam-1
      src  = iam+1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( dest, src, im, jm, km,                      &
                      1, im, jfirst-2, jlast+2, kfirst, klast,       &
                      1, im, jfirst, jfirst, kfirst, klast, vc )            
      call mp_recv3d( src, im, jm, km,                             &
                      1, im, jfirst-2, jlast+2, kfirst, klast,       &
                      1, im, jlast+1, jlast+1, kfirst, klast, vc )
      call t_stopf('ghost_vc')
#endif

      call t_startf('d_core')

!$omp parallel do private(i, j, k, iord, jord)    
      do k=kfirst,klast

        if( k <= km/8 ) then
          if( k == 1 ) then
            iord = 1
            jord = 1
          else
            iord = min(2, iord_d)
            jord = min(2, jord_d)
          endif
        else
          iord = iord_d
          jord = jord_d
        endif

!-----------------------------------------------------------------
! Call the vertical independent part of the dynamics on the D-grid
!-----------------------------------------------------------------

     call d_sw( u(1,jfirst-ng_d,k),      v(1,jfirst-ng_s,k),     &
                uc(1,jfirst-ng_d,k),    vc(1,jfirst-2,k),        &
                pt(1,jfirst-ng_d,k),   delp(1,jfirst,k),         &
                delpf(1,jfirst-ng_d,k), cx3(1,jfirst-ng_d,k),    &
                cy3(1,jfirst,k),        mfx(1,jfirst,k),         &
                mfy(1,jfirst,k), cdx(js2g0,k),  cdy(js2g0,k),    &
                dtdx,   dtdxe,  dtxe5,  txe5,  dyce,  rdx,  cy,  &
                dx,  f0(jfirst-ng_d), js2g0,  jn1g1, im,  jm,    &
                jfirst, jlast,  ng_d,  ng_s,   nq,    iord,      &
                jord,   zt_d,   rcap,  tiny,   dtdy,             &
                dtdy5,  tdy5,   rdy,    cosp,  acosp, cose,      &
                coslon, sinlon, cosl5, sinl5 )

  enddo

  call t_stopf('d_core')

  call t_startf('d_geop')

#if defined( SPMD )
      if (twod_decomp == 1) then

       if (geopkd16) then
!

!
! Stay in yz space and use semi-global z communications and 16-byte reals

         call geopk16(ptop, pe, delp, pkcc, wzc, hs, pt, ng_d, im, jm, km,      &
                      jfirst, jlast, 1, im, cp, akap, kfirst, klast)

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
        do k = kfirst, klast+1
          do j = jfirst, jlast
            do i = 1, im
              pkc(i,j,k) = pkcc(i,j,k)
              wz(i,j,k) = wzc(i,j,k)
            enddo
          enddo
        enddo

       else
!

!
! Transpose to xy decomposition

        call redistributestart(inter_ijk, .true., delp)
        call redistributefinish(inter_ijk, .true., delpxy)
!
! Temporary solution to redistribute the unghosted version of PT
! Art: can we phase out the redistribute of PT, by using on GEOPK16?
!
!$omp parallel do private(i,j,k)
        do k=kfirst,klast
          do j=jfirst,jlast
            do i=1,im
              yzt(i,j,k) = pt(i,j,k)
            enddo
          enddo
        enddo
        call redistributestart(inter_ijk, .true., yzt)
        call redistributefinish(inter_ijk, .true., ptxy)

        call geopk(ptop, pexy, delpxy, pkxy, wzxy, hsxy, ptxy, 0, im, jm, km,    &
                   jfirstxy, jlastxy, ifirstxy, ilastxy,                    &
                   cp, akap, nx, ipe)

!
! Transpose back to yz decomposition
!

        call redistributestart(inter_ijkp, .false., pkxy)
        call redistributefinish(inter_ijkp, .false., pkkp)
        call redistributestart(inter_ijkp, .false., wzxy)
        call redistributefinish(inter_ijkp, .false., wzkp)

!$omp parallel do private(i, j, k)
        do k = kfirst, klastp
          do j = jfirst, jlast
            do i = 1, im
              pkc(i,j,k) = pkkp(i,j,k)
              wz(i,j,k) = wzkp(i,j,k)
            enddo
          enddo
        enddo

!
! pexy is an output quantity from geopk if ipe=1
! Save pexy and pkxy to avoid main transpose; pe not needed.
!

         if (npr_z > 1) then
!
! Fill in klast+1
!
           call mp_send3d_2( iam-npr_y, iam+npr_y, im, jm, km+1,             &
                             1, im, jfirst-1, jlast+1, kfirst, klast+1,      &
                             1, im, jfirst, jlast, kfirst, kfirst, pkc, wz )
           call mp_recv3d_2( iam+npr_y, im, jm, km+1,                        &
                             1, im, jfirst-1, jlast+1, kfirst, klast+1,      &
                             1, im, jfirst, jlast, klast+1, klast+1, pkc, wz)
         endif       ! npr_z > 1

       endif       ! geopkd16

      else
#endif

        if (geopkd16) then

!
! Use 16-byte reals (for compatibility with 2-D decomposition)
!

          call geopk16(ptop, pe, delp, pkcc, wzc, hs, pt, ng_d, im, jm, km,          &
                       jfirst, jlast, 1, im, cp, akap, 1, km)

        else

!
! Use 8-byte reals (standard)
!

          call geopk(ptop, pe, delp, pkcc, wzc, hs, pt, ng_d, im, jm, km,            &
                     jfirst, jlast, 1, im,                                     &
                     cp, akap, nx, ipe)

        endif

!
! Geopk does not need j ghost zones of pkc and wz
!

!$omp parallel do private(i, j, k)
        do k = kfirst, klast+1
          do j = jfirst, jlast
            do i = 1, im
              pkc(i,j,k) = pkcc(i,j,k)
              wz(i,j,k) = wzc(i,j,k)
            enddo
          enddo
        enddo

#if defined( SPMD )
      endif
#endif
!
! Upon exit from geopk, the quantities pe, pkc and wz will have been
!      updated at klast+1

      call t_stopf('d_geop')


#if defined( SPMD )
! Exchange boundary regions on north and south for pkc and wz
      call t_startf('send_pkc_wz')
      call mp_send2_ns( im, jm, km+1, jfirst, jlast,           &
                        kfirst, klast+1, 1, pkc, wz)
      call t_stopf('send_pkc_wz')
#endif

   if ( ipe /= 1 ) then          !  not the last call
!
! Perform some work while sending data on the way
!

!$omp parallel do private(i, j, k, wk, wk2)

    do k=kfirst,klast
         do j=jfirst,jlast
            do i=1,im
               delpf(i,j,k) = delp(i,j,k)
             enddo
          enddo
          call pft2d( delpf(1,js2g0,k), sc(js2g0), dc(1,js2g0),  &
                      im, jn2g0-js2g0+1, ifax, trigs, wk, wk2 )
    enddo

   else
! Last call
!$omp parallel do private(i, j, k)
    do k=kfirst,klast+1
        do j=jfirst,jlast
           do i=1,im
             pk(i,j,k) = pkc(i,j,k)
           enddo
         enddo
    enddo
   endif

#if defined( SPMD )
      call t_startf('recv_pkc_wz')
      call mp_recv2_ns( im, jm, km+1, jfirst, jlast,          &
                        kfirst, klast+1, 1, pkc, wz)
      call t_stopf('recv_pkc_wz')
      if ( ipe /= 1 ) then          !  not the last call
         call t_startf('send_delpf')
         call mp_send4d_ns( im, jm, km, 1, jfirst, jlast,         &
                            kfirst, klast, ng_d, ng_d, delpf )
         call t_stopf('send_delpf')
      endif
#endif


#if defined ( HIGH_P )

      call t_startf('highp')

      call highp2(pkc,  wz,  wz3,    wzz,   dpt,          &
                  im,   jm,   km,  jfirst, jlast,         &
                  kfirst, klast, klastp, nx )

      call t_stopf('highp')

#else

!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k)

   do k=kfirst,klast
      do j=js1g1,jn1g1                  ! dpt needed NS
         do i=1,im                       ! wz, pkc ghosted NS
            dpt(i,j,k)=(wz(i,j,k+1)+wz(i,j,k))*(pkc(i,j,k+1)-pkc(i,j,k))
         enddo
      enddo
   enddo

#endif

!  GHOSTING:   wz (input) NS ; pkc (input) NS

      call t_startf('d-pgrad')

!$omp parallel do private(i, j, k, wk3, wk1)

   do 4500 k=kfirst,klast+1

#if ( !defined HIGH_P )
      if (k == 1) then
        do j=js2g0,jlast
          do i=1,im
            wz3(i,j,1) = 0.
            wz(i,j,1) = 0.
          enddo
        enddo
        pk4 = 4.*ptop**akap
        do j=js2g0,jn1g1
          do i=1,im
            pkc(i,j,1) = pk4
          enddo
        enddo
        goto 4500
      endif
#endif

      do j=js2g1,jn2g0                             ! wk3 needed S
#if defined ( HIGH_P )
          do i=1,im
             wk3(i,j) = wz3(i,j,k) 
          enddo
#else
            wk3(1,j) = (wz(1,j,k)+wz(im,j,k)) *       &
                     (pkc(1,j,k)-pkc(im,j,k))
          do i=2,im
            wk3(i,j) = (wz(i,j,k)+wz(i-1,j,k)) *      & 
                       (pkc(i,j,k)-pkc(i-1,j,k))
          enddo
#endif
      enddo

      do j=js2g1,jn2g0                               
         do i=1,im-1
            wk1(i,j) = wk3(i,j) + wk3(i+1,j)        
         enddo
            wk1(im,j) = wk3(im,j) + wk3(1,j)      ! wk3 ghosted S
      enddo

      if ( jfirst == 1 ) then
          do i=1,im
            wk1(i, 1) = 0.
          enddo
      endif

      if ( jlast == jm ) then
          do i=1,im
            wk1(i,jm) = 0.
          enddo
      endif

      do j=js2g0,jlast                          ! wk1 ghosted S
         do i=1,im
            wz3(i,j,k) = wk1(i,j) + wk1(i,j-1)
         enddo
      enddo

! N-S walls

      do j=js2g0,jn1g1                     ! wk1 needed N

#if defined ( HIGH_P )
          do i=1,im                          ! wz, pkc ghosted NS
             wk1(i,j) = wzz(i,j,k)
          enddo
#else
          do i=1,im                         ! wz, pkc ghosted NS
             wk1(i,j) = (wz(i,j,k)+wz(i,j-1,k))*(pkc(i,j,k)-pkc(i,j-1,k))
          enddo
#endif
      enddo

      do j=js2g0,jn1g1                    ! wk3 needed N
            wk3(1,j) = wk1(1,j) + wk1(im,j)      ! wk1 ghosted N
          do i=2,im
            wk3(i,j) = wk1(i,j) + wk1(i-1,j)   ! wk1 ghosted N
          enddo
      enddo

      do j=js2g0,jn2g0
         do i=1,im
            wz(i,j,k) = wk3(i,j) + wk3(i,j+1)  ! wk3 ghosted N
         enddo
      enddo

      do j=js1g1,jn1g1
         wk1(1,j) = pkc(1,j,k) + pkc(im,j,k)
         do i=2,im
            wk1(i,j) = pkc(i,j,k) + pkc(i-1,j,k)
         enddo
      enddo
 
      do j=js2g0,jn1g1
         do i=1,im
            pkc(i,j,k) = wk1(i,j) + wk1(i,j-1)
         enddo
      enddo
4500  continue

!     GHOSTING:   dpt (loop 4000) NS ; pkc (loop 4500) N
!
! Beware k+1 references directly below (AAM)
!
!$omp parallel do private(i, j, k, wk, wk1, wk2, wk3)

   do 6000 k=kfirst,klast

      do j=js1g1,jn1g1
         wk1(1,j) = dpt(1,j,k) + dpt(im,j,k)
         do i=2,im
            wk1(i,j) = dpt(i,j,k) + dpt(i-1,j,k)
         enddo
      enddo
 
      do j=js2g0,jn1g1
         do i=1,im
            wk2(i,j) = wk1(i,j) + wk1(i,j-1)
             wk(i,j) = pkc(i,j,k+1) - pkc(i,j,k)
         enddo
      enddo

      do j=js2g0,jlast
         do i=1,im-1
            wk3(i,j) = uc(i,j,k)  +  dtdxe(j)/(wk(i,j) + wk(i+1,j))      &
                       * (wk2(i,j)-wk2(i+1,j)+wz3(i,j,k+1)-wz3(i,j,k))
         enddo
          wk3(im,j) = uc(im,j,k)  +  dtdxe(j)/(wk(im,j) + wk(1,j))       &
                      * (wk2(im,j)-wk2(1,j)+wz3(im,j,k+1)-wz3(im,j,k))
      enddo

      do j=js2g0,jn2g0                  ! Assumes wk2 ghosted on N
         do i=1,im
            wk1(i,j) = vc(i,j,k)  +  dtdy/(wk(i,j)+wk(i,j+1)) *          &
                       (wk2(i,j)-wk2(i,j+1)+wz(i,j,k+1)-wz(i,j,k))
         enddo
      enddo

#if ( !defined ALT_PFT )
      call pft2d( wk3(1,js2g0), se(js2g0), de(1,js2g0), im,           &
                  jlast-js2g0+1,  ifax, trigs, wk, wk2 )
      call pft2d( wk1(1,js2g0), sc(js2g0), dc(1,js2g0), im,           &
                  jn2g0-js2g0+1,  ifax, trigs, wk, wk2 )
#endif

      do j=js2g0,jn2g0
         do i=1,im
            v(i,j,k) = v(i,j,k) + wk1(i,j)
            u(i,j,k) = u(i,j,k) + wk3(i,j)
         enddo
      enddo
 
      if ( jlast == jm ) then
         do i=1,im
            u(i,jlast,k) = u(i,jlast,k) + wk3(i,jlast)
         enddo
      endif

#if defined ( ALT_PFT )
      call pft2d( u(1,js2g0,k), se(js2g0), de(1,js2g0), im,    &
                  jlast-js2g0+1,  ifax, trigs, wk, wk2 )
      call pft2d( v(1,js2g0,k), sc(js2g0), dc(1,js2g0), im,    &
                  jn2g0-js2g0+1,  ifax, trigs, wk, wk2 )
#endif

6000  continue
      call t_stopf('d-pgrad')

#if defined( SPMD )
      if ( ipe /= 1 ) then
         call t_startf('recv_delpf')
         call mp_recv4d_ns( im, jm, km, 1, jfirst, jlast,       &
                            kfirst, klast, ng_d, ng_d, delpf )
         call t_stopf('recv_delpf')
      endif
#endif

      return
!EOC
 end
!-----------------------------------------------------------------------
