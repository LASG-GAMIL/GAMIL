! $Id: ice_tstm.F,v 1.9.2.3 2002/06/20 01:49:08 eaton Exp $
!=======================================================================
!---! Energy-conserving sea ice model
!---! Routines to solve heat-equation using linear method
!---!
!---! author C. M. Bitz
!---!
!---! See Bitz, C.M., and W.H. Lipscomb, 1999: 
!---! An energy-conserving thermodynamic model of sea ice,
!---! J. Geophys. Res., 104, 15,669-15,677. 
!---!     
!---! See Bitz, C.M., M.M. Holland, A.J. Weaver, M. Eby:
!---! Simulating the ice-thickness distribution in a climate model
!---! accepted to J. Geophys. Res.
!---! (prepublications copy available upon request)
!---!
!---! The author grants permission to the public to copy and use this
!---! software without charge, provided that this Notice and any statement
!---! of authorship are reproduced on all copies and any publications that
!---! result from the use of this software must (1) refer to the publications 
!---! listed above and (2) acknowledge the origin and author of the model.
!---! This software is without warranty, expressed or implied, and the
!---! author assumes no liability or responsibility for its use. 
!=======================================================================

      module ice_tstm

      use ice_constants

      implicit none

!=======================================================================

      contains

!=======================================================================

      subroutine tstm(   dtau, tmz1d, sal1d,   Tf1
     $                ,  area,     hi,    hs,  dswr
     $                , dswrv, dswrn,   flwd, dflwup
     $                ,  dflh, dfsh ,  asnow
     $                ,  tbot,    Ib,     F, condb
     $                ,    ts,    ti,  flwup,   flh
     $                ,   fsh,linpts)

!---!-------------------------------------------------------------------
!---! This routine calculates the evolution of the ice interior and 
!---! surface temperature from the heat equation and surface energy 
!---! balance
!---!
!---! The albedo is fixed for this calculation
!---!
!---! Solves the heat equation which is non-linear for saline ice
!---! by linearizing and then iterating a set of equations
!---! Scheme is backwards Euler giving a tridiagonal
!---! set of equations implicit in temperature
!---! (Tried crank-nicholson but it behaves poorly when
!---! thin ice has a weird initial temperature profile. 
!---!
!---!    The solution to the heat equation ignores the insulating effect 
!---!    of snow if it less than hsmin (1 cm) thick, but I do not like to 
!---!    "kill" it when it is that thick because sometimes the snowfall 
!---!    rate is really small...
!---! Must have a sufficient amount of snow to solve heat equation in snow
!---! hsmin is the minimum depth of snow in order to solve for ti(0)
!---! if snow thickness < hsmin then do not change ti(0)
!---!
!---! The number of equations that must be solved by the tridiagonal solver
!---! depends on whether the surface is melting and whether there is snow. 
!---! Four cases are possible:
!---! 1 = freezing w/ snow, 2 = freezing w/ no snow,
!---! 3 = melting w/ snow, and 4 = melting w/ no snow
!---!
!---! NOTE:
!---!
!---! If you are using an sst/ice dataset that was interpolated from a 
!---! different resolution the temperature profiles may not be in equilibrium. 
!---! The model will try to adjust the profiles to allow for convergence.
!---! The user should see the following message on the first time step of
!---! a run as the profile is adjusted:
!---!
!---! WARNING: ice_tstm ::profile reset xx points at chunck nn see ice_tstm.F 
!---! for more information'
!---!
!---! Occasionally the ice dataset may be so far out of equilibrium that the
!---! model blows up on the first timestep.  In this circumstance you might
!---! try setting the namelist variable 'reset_csim_iceprops = .T'. and 
!---! restarting the run.  This will force all ice profiles to an initial 
!---! an initial state of freezing and allow the model to spin up to
!---! equilibrium. If this warning is generated after the first few time
!---! steps it may indicate a problem with the ice model.
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &   dtau    ! timestep                                     (s)          
     &,  tmz1d(ni) ! melting temp of each layer                   (C)
     &,  sal1d(ni+1) ! salinity of each layer                       (ppt)
     &,  area    ! area of the ice/snow
     &,  hi,hs    ! ice and snow thickness                       (m)
     &,  asnow      ! 1 - patchy snow frac
     &,  dswr    ! above srfc net dnwd shortwave, positive down (W/m**2)
     &,  dswrv   ! dswr in vis (wvlngth < 700nm)                (W/m**2)
     &,  dswrn   ! dswr in nir (wvlngth > 700nm)                (W/m**2)
     &,  flwd     ! dnwd longwave flux (always positive)         (W/m**2)
     &,  dfsh    ! derriv wrt ts of dnwd sensible flux          (W/m**2)
     &,  dflh    ! derriv wrt ts of dnwd latent flux            (W/m**2)
     &,  dflwup   ! derriv wrt ts of upwd longwave flux          (W/m**2)
     &,  Tf1      ! freezing temp of water below ice             (C)

      real (kind=dbl_kind), intent(out) ::
     &   tbot    ! bottom ice temp                              (C)
     &,  condb   ! conductive flux at bottom srfc                (W/m**2)
     &,  F       ! net flux at top srfc including conductive flux (W/m**2)
     &,  Ib      ! solar passing through the bottom ice srfc    (W/m**2)

      real (kind=dbl_kind), intent(inout) ::
     &   ts       ! srfc temp of snow or ice  (C)
     &,  ti(0:ni) ! snow(0) and ice(1:ni) interior temp (C)
     &,  fsh     ! dnwd sensible flux (W/m**2)
     &,  flh     ! dnwd latent flux (always negative) (W/m**2)
     &,  flwup    ! upwd longwave flux (always negative) (W/m**2)

      integer (kind=int_kind), intent(inout) ::
     &  linpts     ! counter for number of ice points reset to linear profile

      ! local variables
      real (kind=dbl_kind) ::
     $   ti_old(0:ni)  ! initial temperature in the ice and snow
     &,  ki(0:ni+1)    ! layer conductivity divided by layer thickness
     &,  zeta(0:ni)    ! the terms in heat equation that are independent of ti
     &,  eta(0:ni)  ! time step / ice layer thickness / fresh ice heat capacity
     &,  tinext(-1:ni) ! incremented ts and ti 
     &,  dti(-1:ni)    ! incremental changes to ts and ti
         ! a,b,c are vectors that describe the diagonal and off-diagonal 
         ! elements of the matrix [A], such that [A] ti = r
     &,  a(-1:ni),b(-1:ni),c(-1:ni),d(-1:ni),r(-1:ni)
     &,  cpiz(ni)     ! saline ice heat capacity 
     &,  Iabs(ni+1) ! solar absorbed in each layer

      real (kind=dbl_kind) ::
     &   ts_old        ! initial temperature of the ice/snow srfc
     &,  Fo      ! net flux at top srfc excluding conductive flux in ice/snow
     &,  dFo_dt  ! derivative of Fo wrt temperature
     &,  condt   ! conductive flux at top srfc
     &,  melts   ! the srfc melting temp (either TSMELT or TIMELT)
     &,  alph, bet   ! parameters for maintaining 2nd order accurate diff at boundar
     &,  z       ! vertical cordinate
     &,  tinterf
         ! a,b,c are vectors that describe the diagonal and off-diagonal 
         ! elements of the matrix [A], such that [A] ti = r
     &,  dhi            ! ice layer thickness
     &,  dt_dhi,dt_hs   ! time step / ice or snow thickness
     &,  errit         ! the absolute value of the maximum dti
     &,  absorb     ! sum of Iabs
     &,  Io      ! solar penetrating top srfc, positive down    (W/m**2)
     &,  Iovis    ! solar penetrating top in vis (wvlngth < 700nm)
     &,  Ionir    ! solar penetrating top in nir (wvlngth > 700nm)

      real (kind=dbl_kind), parameter ::
     &   hsmin = 0.01_dbl_kind  ! minimum allowable snow thickness for heat eq.
     &,  T_errmax = 5.0e-4_dbl_kind     ! error tolerance for temp diff (C) 

      integer (kind=int_kind) ::
     &   N         ! number of equations solved by tridiag solver
     &,  layer      ! counter for ice layers
     &,  ipc        ! counter for iterations of dti

      logical (kind=log_kind) ::
     &   convrg       ! flag that is true if temperature converges
     &,  complt       ! flag that is true if a layer melts internally

      logical (kind=log_kind) :: verbos
!      verbos = .true.
      verbos = .false.

      !-----------------------------------------------------------------
      ! allocate local arrays
      !-----------------------------------------------------------------
!      allocate(ti_old(0:ni))
!      allocate(tinext(-1:ni))
!      allocate(dti(-1:ni))
!      allocate(a(-1:ni))
!      allocate(b(-1:ni))
!      allocate(c(-1:ni))
!      allocate(d(-1:1))
!      allocate(r(-1:ni))
!      allocate(ki(0:ni+1))
!      allocate(zeta(0:ni))
!      allocate(eta(0:ni))
!      allocate(cpiz (ni))
!      allocate(Iabs(ni+1))

      !-----------------------------------------------------------------
      ! setup helpful parameters
      !-----------------------------------------------------------------
      dhi    = hi/real(ni)
      dt_dhi = dtau/dhi
      dt_hs = c0
      if (hs.gt.hsmin) dt_hs=dtau/hs

      ts_old = ts
      do layer = 0,ni
         ti_old(layer) = ti(layer)
      enddo
      tbot = Tf1

      call conductiv(sal1d,hsmin,ki,ti,tbot,hs,dhi)

      ! the solar radiation absorbed internally
      z = c0

! changing io from 0.3 to the two band version here
! allows a lot more of the light that penetrates
!  the upper surface to pass through to the ocean
      ! transmit 0.7 of visible and none of the near ir
      Iovis = dswrv * 0.7_dbl_kind * asnow
      Ionir = dswrn * 0.0_dbl_kind * asnow
      Io = Iovis + Ionir
!      print*,Io,Iovis,Ionir

      Iabs(1) = Io
      do layer = 1,ni
        z = z+dhi
        if (z.lt.p1) then ! no absorbtion in the first 10 cm
          Iabs(layer+1) = Io
        else
          Iabs(layer+1) = Iovis*exp(-kappav*(z-p1))
     $                  + Ionir*exp(-kappan*(z-p1))
        endif
        Iabs(layer) = (Iabs(layer)-Iabs(layer+1))
      enddo

      ! all fluxes are positive down
      ! even flwup, which has an explicit UP in its variable name,
      ! is ALWAYS NEGATIVE
      ! flhn (latent heat flux) takes heat away from ice/snow
      ! from sublimation when it is negative 
      dFo_dt = dflh + dfsh + dflwup
      Fo     = dswr - Io + flwd + flh + fsh + flwup

      if (hs.gt.hsmin) then
         alph  = c2*(c2*hs + dhi)/(hs+dhi)
         bet   = -c2*hs*hs/(c2*hs+dhi)/(hs+dhi)
         melts = Tsmelt
      else
         alph  =  c3
         bet   =  -p33
         melts = Timelt
      endif

      ts_old = ts
      ts     = min(ts,melts-p01)   ! absolutely necessary

      !-----------------------------------------------------------------
      ! beginning of iterative proceedure
      !-----------------------------------------------------------------

 500  continue
      ipc = 1
 1000 continue

      !-----------------------------------------------------------------
      ! setup terms that depend on the iterated temperature
      !-----------------------------------------------------------------
      do layer = 1,ni
        cpiz(layer) = rcpi+rLfidepressT*sal1d(layer)/Ti_old(layer)/
     $       min(Ti(layer),tmz1d(layer))
      enddo

      eta(0) = c0
      if (hs.gt.hsmin) eta(0) = dtau/hs
      do layer = 1,ni
         eta(layer) = dt_dhi/cpiz(layer)
      enddo
      
      zeta(0) = rcps*Ti_old(0)
      do layer = 1,ni
        zeta(layer) = Ti_old(layer)+eta(layer)*Iabs(layer)
      enddo

      if (ts.lt.melts-puny) then
         ! solve heat equation for ice/snow using flux BC

         if (hs.gt.hsmin) then
      !-----------------------------------------------------------------
      !     case = 1      ! case of freezing with snow layer
      !-----------------------------------------------------------------

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,1)

            a(0) = -eta(0)*ki(0)*(alph+bet)
            c(0) = eta(0)*(bet*ki(0)-ki(1))
            b(0) = rcps+eta(0)*(alph*ki(0)+ki(1))
            r(0) = zeta(0)

            a(-1) = c0
            c(-1) = ki(0)*alph
            d(-1) = ki(0)*bet
            b(-1) = dFo_dt-c(-1)-d(-1)
            r(-1) = -Fo+dFo_dt*ts_old

            ! row operation to get rid of d(-1)
            b(-1) = c(0)*b(-1)-d(-1)*a(0)
            c(-1) = c(0)*c(-1)-d(-1)*b(0)
            r(-1) = c(0)*r(-1)-d(-1)*r(0)

            N = ni+2
            call tridag(a(-1),b(-1),c(-1),r(-1),tinext(-1),N)
            dti(-1) = tinext(-1)-ts
            ts = tinext(-1)
            do layer = 0,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo

            if (verbos)
     $           write(*,2000) 'I. iterations',real(ipc),ts
     $           ,(ti(layer),layer = 0,4),(dti(layer),layer = -1,4)

         else
      !-----------------------------------------------------------------
      !     case = 2      ! case of freezing with no snow layer
      !-----------------------------------------------------------------
           
            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,2)

            a(1) = -eta(1)*ki(1)*(alph+bet)
            c(1) = -eta(1)*(ki(2)-bet*ki(1))
            b(1) = c1+eta(1)*(ki(2)+alph*ki(1))
            r(1) = zeta(1)

            a(0) = c0
            c(0) = ki(1)*alph
            d(0) = ki(1)*bet
            b(0) = dFo_dt-c(0)-d(0)
            r(0) = -Fo+dFo_dt*ts_old

            ! row operation to get rid of d(0)
            b(0) = c(1)*b(0)-d(0)*a(1)
            c(0) = c(1)*c(0)-d(0)*b(1)
            r(0) = c(1)*r(0)-d(0)*r(1)

            N = ni+1
            call tridag(a(0),b(0),c(0),r(0),tinext(0),N)

            dti(-1) = c0
            dti(0)  = tinext(0)-ts
            ts = tinext(0)
            do layer = 1,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo
            ti(0) = min(ts,melts)

            if (verbos)
     $           write(*,2000) 'II. iterations', real(ipc),ts,
     $           (ti(layer),layer = 0,4),(dti(layer),layer = 0,4)

         endif
         ts = min(ts,melts)
      else
         ts = melts
         if (hs.gt.hsmin) then
      !-----------------------------------------------------------------
      !     case = 3      ! case of melting with snow layer
      !-----------------------------------------------------------------
            
            a(0) = -eta(0)*ki(0)*(alph+bet)
            c(0) = eta(0)*(bet*ki(0)-ki(1))
            b(0) = rcps+eta(0)*(alph*ki(0)+ki(1))
            r(0) = zeta(0)-a(0)*ts
            a(0) = c0

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,1)

            N = ni+1
            call tridag(a(0),b(0),c(0),r(0),tinext(0),N)

            dti(-1) = c0
            do layer = 0,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo

            if (verbos)
     $           write(*,2000)'III. T =  ', real(ipc)
     $           ,(ti(layer),layer = 0,4),(dti(layer),layer = 1,4)

         else           
      !-----------------------------------------------------------------
      !     case = 4      ! case of melting with no snow layer
      !-----------------------------------------------------------------
            
            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,2)

            a(1) = -eta(1)*ki(1)*(alph+bet)
            c(1) = -eta(1)*(ki(2)-bet*ki(1))
            b(1) = c1+eta(1)*(ki(2)+alph*ki(1))
            r(1) = zeta(1)-a(1)*Ts
            a(1) = c0
           
            N = ni
            call tridag(a(1),b(1),c(1),r(1),tinext(1),N)

            dti(-1) = c0
            dti(0) = c0
            do layer = 1,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo
            ti(0) = ts

            if (verbos)
     $           write(*,2000)'IV. T =  ', real(ipc)
     $           ,(ti(layer),layer = 0,4),(dti(layer),layer = 1,4)
         endif
      endif

      !-----------------------------------------------------------------
      ! end iterative proceedure, see if need to reiterate
      !-----------------------------------------------------------------
      errit = c0
      do layer = -1,ni
         errit = max(errit,abs(dti(layer)))
      enddo
      if ((errit.lt.T_errmax).or.(ipc.gt.20)) go to 6000 ! done iterating

      ipc = ipc+1
      go to 1000    ! to beginning of iterative routine

      !-----------------------------------------------------------------
      ! continue from here when done iterating
      ! begin with error checking
      !-----------------------------------------------------------------
 6000 continue

      complt = .false.
      do layer = 1,ni
         if (ti(layer).gt.tmz1d(layer)) complt = .true.
      enddo

      convrg = .true.
      if (errit.gt.T_errmax) then
        if (errit.gt.0.005_dbl_kind) then
          write(6,*) 'WARNING NO CONVERGENCE in icemodel',errit
          convrg = .false.
        else
          write(6,*) 'WARNING POOR CONVERGENCE in icemodel',errit
        endif
      endif

      if (complt .or. .not.convrg) then
        if (complt .and. convrg) then
          write(6,*) 'WARNING CONVERGES to ti > melting in icemodel'
        endif
        if (complt .and. .not.convrg) 
     $       write(6,*) '& final iteration has ti > melting'
!        write(6,*) '***************************************************'
!        write(6,*) 'HINT:If the ice model blows up due to this problem'
!        write(6,*) 'you might try setting the namelist variable'
!        write(6,*) 'reset_csim_iceprops = .T.'
!        write(6,*) 'and restarting the run'
!        write(6,*) '***************************************************'
        if (verbos) then
           write(6,*) 'diagnostics useful for debugging'
           write(6,*) area,hi,hs,dswr,flwd
           write(6,*) 'tstmnew computed T(z):'
           write(6,*) ts,(ti(layer),layer = 0,ni)
           write(6,*) 'tstmnew started with T(z):'
           write(6,*) ts_old,(ti_old(layer),layer = 0,ni),tbot
        end if        
        ts = ts_old
        tinterf  =  (hi*ksno*Ts+hs*kice*tbot)/(hi*ksno+hs*kice)
        Ti(0) = p5*(Ts+tinterf)
        do layer = 1,ni
          Ti(layer) = tinterf+(real(layer)-p5)*(tbot-tinterf)/real(ni)
        enddo
!        write(6,*) 'setting the temperature profile to be linear T(z):'
!        write(6,*) ts,(ti(layer),layer = 0,ni)
        linpts=linpts+1
      elseif (ti(0).gt.(Tsmelt+puny)) then
        if (hs.gt.hsmin) then
          write(6,*) 'WARNING snow temperature > melting in ice model'
          write(6,*) '  likely no problem if area is puny'
          write(6,*) '  area=',area
          write(6,*) 'setting snow temperature equal to melting'
        endif
        ti(0) = Tsmelt
      endif

      !-----------------------------------------------------------------
      ! finish up by updating the sfc fluxes, etc.
      !-----------------------------------------------------------------

      if (hs.gt.hsmin) then
         condt = ki(0)*(alph*(ti(0)-ts)+bet*(ti(1)-ts))
      else
         condt = ki(1)*(alph*(ti(1)-ts)+bet*(ti(2)-ts))
      endif
      condb = ki(ni+1)*( c3*(ti(ni)-tbot) - (ti(ni-1)-tbot)/c3  ) 

      fsh  = fsh+dfsh*(Ts-Ts_old)
      flh  = flh+dflh*(Ts-Ts_old)
      flwup = flwup+dflwup*(Ts-Ts_old)
      Fo   = Fo+dFo_dt*(Ts-Ts_old)
      F    = Fo + condt

      ! in the rare event that F<0 set it equal to 0 and adjust sensible heat 
      if (F .lt. 0.) then 
        Fo   = -condt
        F    = c0
        fsh = Fo-dswr+Io-flwd-flh-flwup
      endif

      absorb = c0
      do layer = 1,ni
         absorb = absorb+Iabs(layer)
      enddo
      Ib = Io-absorb  ! solar passing through the bottom of the ice

!      print*, '(tstmnew)', F, condt,condb,(-condt+condb)*dtau

      !-----------------------------------------------------------------
      ! deallocate arrays
      !-----------------------------------------------------------------
!      deallocate(ti_old)
!      deallocate(tinext)
!      deallocate(dti)
!      deallocate(a)
!      deallocate(b)
!      deallocate(c)
!      deallocate(d)
!      deallocate(r)
!      deallocate(ki)
!      deallocate(zeta)
!      deallocate(eta)
!      deallocate(cpiz)
!      deallocate(Iabs)

 2000 format(A15, 50(f8.3))

      end subroutine tstm

!=======================================================================

      subroutine getabc(a,b,c,r,ti,tbot,zeta,k,eta,lfirst)

!---!-------------------------------------------------------------------
!---!     compute elements of tridiagonal matrix
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: 
     $        lfirst          ! start with this layer 

      real (kind=dbl_kind), intent(in) ::
     &        ti   (0:ni)   ! temperature of ice-snow layers
     $       ,tbot            ! temperature of ice bottom srfc
     $       ,zeta (0:ni)   ! 
     $       ,k    (0:ni+1) ! ice-snow conductivitiy
     $       ,eta  (0:ni)   !

      real (kind=dbl_kind), intent(out) ::
     &        a    (-1:ni)  ! sub-diagonal elements
     $       ,b    (-1:ni)  ! diagonal elements
     $       ,c    (-1:ni)  ! super-diagonal elements
     $       ,r    (-1:ni)  ! constants (indep. of ti)

      ! parameters for maintaining 2nd order accurate diff at bottom boundary
      real (kind=dbl_kind), parameter ::
     &   alph = c3
     &,  bet  = -p33

      ! local variable
      integer (kind=int_kind) :: layer

      ! if there is snow lfirst = 1 otherwise it is 2
      do layer = lfirst,(ni-1)
        a(layer) = -eta(layer)*k(layer)
        c(layer) = -eta(layer)*k(layer+1)
        b(layer) = c1-c(layer)-a(layer)
        r(layer) = zeta(layer)
      enddo
      a(ni) = -eta(ni)*(k(ni)-bet*k(ni+1))
      c(ni) = c0
      b(ni) = c1+eta(ni)*(k(ni)+alph*k(ni+1))
      r(ni) = zeta(ni)+eta(ni)*(alph+bet)*k(ni+1)*tbot

      end subroutine getabc

!=======================================================================

      subroutine tridag(a,b,c,r,u,nrows)

      integer (kind=int_kind), intent(in) :: nrows  ! number of rows
      integer, parameter :: nimax = ni+2
      real (kind=dbl_kind), intent(in) ::
     &        a      (nrows)  ! sub-diagonal elements
     $       ,b      (nrows)  ! diagonal elements
     $       ,c      (nrows)  ! super-diagonal elements
     $       ,r      (nrows)  ! constants (indep. of ti)

      real (kind=dbl_kind), intent(out) ::
     &        u      (nrows)  ! solution

      integer (kind=int_kind) :: layer
      real (kind=dbl_kind) ::
     &        bet             ! work variable

      real (kind=dbl_kind), dimension (nimax) ::
     $ gam   ! work array
      
!      allocate(gam(nrows))

      bet = b(1)
      u(1) = r(1)/bet
      do layer = 2,nrows
        gam(layer) = c(layer-1)/bet
        bet = b(layer)-a(layer)*gam(layer)
        u(layer) = (r(layer)-a(layer)*u(layer-1))/bet
      enddo
      do layer = nrows-1,1,-1
        u(layer) = u(layer)-gam(layer+1)*u(layer+1)
      enddo

!      deallocate(gam)

      end subroutine tridag

!=======================================================================

      subroutine conductiv(sal1d,hsmin,ki,ti,tbot,hs,dhi)

      real (kind=dbl_kind), intent(in) ::
     &   hsmin         ! minimum allowable snow thickness for heat eq.
     &,  sal1d(ni+1)  ! salinity of each layer
     &,  hs,ti(0:ni),tbot,dhi
      real (kind=dbl_kind), intent(out) ::
     &   ki(0:ni+1)

      integer (kind=int_kind) :: layer
      real (kind=dbl_kind), parameter ::
     &   tmax = -p1
     &,  beta  = 0.1172_dbl_kind   ! param for conductivity        (W/m)
     &,  kimin = 0.1000_dbl_kind   ! minimum conductivity in ice   (W/m/deg)

      do layer = 2,ni
         ki(layer) = kice+beta*(sal1d(layer)+sal1d(layer+1))/c2
     $        /min(tmax,(ti(layer-1)+ti(layer))/c2)
         ki(layer) = max(ki(layer),kimin)
         ki(layer) = ki(layer)/dhi
      enddo
      ki(ni+1) = kice+beta*sal1d(ni+1)/tbot
      ki(ni+1) = max(ki(ni+1),kimin)
      ki(ni+1) = ki(ni+1)/dhi
      ki(1) = kice+beta*sal1d(1)/min(tmax,ti(1))
      ki(1) = max(ki(1),kimin)
      ki(1) = ki(1)/dhi
!      ki(1) = kice/dhi
      if (hs.gt.hsmin) then
         ki(0) = ksno/hs
         ki(1) = c2*ki(1)*ki(0)/(ki(1)+ki(0))
      else
         ki(0) = c0
      endif
              
      end subroutine conductiv

!=======================================================================

      end module ice_tstm

!=======================================================================
