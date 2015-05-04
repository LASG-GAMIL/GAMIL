c=======================================================================
!---! Energy-conserving sea ice model
!---! Routines to grow/melt ice and adjust temperature profile
!---!
!---! author C. M. Bitz
!---!
!---! See Bitz, C.M., and W.H. Lipscomb, 1999: 
!---! An energy-conserving thermodynamic model of sea ice,
!---! J. Geophys. Res., 104, 15,669-15,677. 
!---!     
!---! The author grants permission to the public to copy and use this
!---! software without charge, provided that this Notice and any statement
!---! of authorship are reproduced on all copies and any publications that
!---! result from the use of this software must (1) refer to the publications 
!---! listed above and (2) acknowledge the origin and author of the model.
!---! This software is without warranty, expressed or implied, and the
!---! author assumes no liability or responsibility for its use. 
c=======================================================================

      module ice_dh

      use ice_kinds_mod
      use ice_constants
!      use ice_itd

      implicit none

      integer (kind=int_kind) :: errflag

! ADD FLAG FOR ATMOSPHERE MODEL VERSION OF THIS ROUTINE
! Skips the last half of this routine
! Maybe someday this model will be used with a mixed layer ocean
! and then this could be set to false
      logical (kind=log_kind), parameter :: fixice = .true.
!
! Add flag if prognostic snow on ice, and flag to reset csim iceprops.
!
      logical (kind=log_kind) :: prognostic_icesnow,reset_csim_iceprops

c=======================================================================

      contains

c=======================================================================

      logical function icemodel_is( name )
!
! Input arguments
!
      character(len=*) :: name
      
      if (name == 'csim' .or. name == 'CSIM' ) then
         icemodel_is = .true.
      else
         icemodel_is = .false.
      end if
      
      return

      end function icemodel_is

c=======================================================================

      subroutine dh(  dtsub,    sal1d,    tiz
     $              ,  tbot,       hi,       hs,   fbot
     $              ,  fnet,    condb,      flh
     $              ,  dhib,     dhit,      dhs,   subi
     $              ,  subs,     dhif,     dhsf,     qi
     $              ,  focn, i,j )

!---!-------------------------------------------------------------------
!---! Computes the thickness changes at the top and bottom
!---! and adjusts layer energy of melt
!---! does not allow h<0  
!---! Focn= actual flux of heat from the ocean layer under sea ice 
!---! (equal to fbot unless all the ice melts away)
!---! compensates for rare case of melting entire slab through 
!---!-------------------------------------------------------------------

!      use ice_state
!      use ice_diagnostics

      real (kind=dbl_kind), intent(in) :: 
     &   dtsub                ! timestep
     &,  sal1d   (plevmx+1)     ! ice salinity                           (ppt)
     &,  tiz   (0:plevmx)       ! snow/ice internal temp                   (C)
     &,  Tbot                 ! ice bottom in                            (C)
     &,  hi                   ! initial ice thickness                    (m)
     &,  hs                   ! initial snow thickness                   (m)
     &,  fbot                 ! flx from ocean, potent.             (W/m**2)
     &,  fnet                 ! net flx at top srf incl. cond. flx  (W/m**2)
     &,  condb                ! cond. flx at bot.                   (W/m**2)
     &,  flh                  ! latent heat flx                     (w/m**2)

      integer (kind=int_kind), intent(in) :: 
     &   i,j ! grid location for debugging

      ! thickness changes from grow/melt (default) or sublimate/flooding      
      real (kind=dbl_kind), intent(out) :: 
     &   dhib                 ! ice bot, dhib<0 if melt                  (m)
     &,  dhit                 ! ice top, dhit<=0                         (m)
     &,  dhs                  ! snow top, dhit<=0                        (m)
     &,  subi                 ! ice top, subi<0 if sublimating           (m)
     &,  subs                 ! snow, subs<0 if sublimating              (m)
     &,  dhif                 ! ice top from flooding, dhif>0            (m)
     &,  dhsf                 ! snow from flooding, dhsf<0               (m)

     &,  qi(plevmx)             ! energy of melt of ice per unit vol. (J/m**3)
     &,  focn                 ! actual flx of heat used from ocn    (w/m**2)

      real (kind=dbl_kind) :: 
     &   delti(plevmx)          ! evolving ice layer thickness 
     &,  delts                ! evolving snow thickness
     &,  sumr                 ! dummy for adjusting delti
     &,  rmvi                 ! another dummy for adjusting delti
     &,  dtop                 ! dhit+subi for adjust
      
     &,  qigrow               ! energy of melt of ice that grows 
     &,  ebot                 ! heat available to grow/melt at bottom
     &,  etop                 ! heat available to melt at top

     &,  qs                   ! energy of melt of snow per unit vol. (J/m**3)
     &,  enet                 ! sum of energy of melt of ice and snow (J/m**2)
     &,  evnet                ! sum of energy of melt and vapor. of ice and snow (J/m**2)

     &,  qiflood              ! energy of melt of flooded ice      (W/m**2)

      real (kind=dbl_kind) ::     hi_tw          ! ice thickness   (m)

      integer (kind=int_kind) :: layer

      logical (kind=log_kind) :: verbos
c      verbos = .true.
      verbos = .false.

      errflag = 0

      dhib  = c0
      dhit  = c0
      dhs   = c0
      dhif  = c0
      dhsf  = c0
      subi  = c0
      subs  = c0
      qiflood = c0

      ! thickness of snow and each ice layer 
      delti(1) = hi/ni
      do layer = 2,ni
        delti(layer) = delti(1)
      enddo
      delts = hs

      ! energy of melt per unit vol for snow and ice for each layer and sum
      qs = -rLfs
      enet = c0
      do layer = 1,ni
        qi(layer) = energ(tiz(layer),sal1d(layer))
        enet = enet + qi(layer)
      enddo
      enet = enet*delti(1) + qs*hs

      !-----------------------------------------------------------------
      ! sublimate/condense 
      !-----------------------------------------------------------------
      etop  = -flh*dtsub           !  etop>0
      evnet = -rLvs*hs-rLvi*hi+enet+etop
      if ( evnet .ge. 0. ) then    !  should never happen
        subi = -hi
        subs = -hs
        focn = condb + evnet/dtsub   
        write(6,*)  flh,dtsub,-rLvs*hs,-rLvi*hi,enet,etop
        write(6,*) 'sublimate away all sea ice'
        write(6,*) 'something is probably seriously wrong'
!        call print_state('ice state at dh stop',i,j)
        stop
      endif
!      call t_startf ('srfsub')
      call srfsub(    qi,   qs, delti, delts,
     $              subi, subs,  etop,  enet )
!      call t_stopf ('srfsub')

      ! adjust the layer thickness to reflect subl/cond changes
      delts = hs + subs
      rmvi = subi
      do layer = 1,ni
        sumr = max( -delti(layer), rmvi )
        rmvi = rmvi-sumr
        delti(layer) = delti(layer) + sumr
      enddo

      !-----------------------------------------------------------------
      ! melt at top srf melt      
      ! may melt when Tsf < melting, but alway a neglible amount
      !-----------------------------------------------------------------
      if ( fnet .gt. 0. ) then
        etop = fnet * dtsub
        enet = enet + etop
        if ( enet .ge. 0. ) then          !    remotely possible
          dhit = -( hi + subi )
          dhs  = -( hs + subs )
          focn = condb + enet/dtsub
          return
        endif
!        call t_startf ('srfmelt')
        call srfmelt(   qi,   qs, delti, delts,
     $                dhit,  dhs,  etop  )
!        call t_stopf ('srfmelt')
        ! adjust the layer thickness to reflect melt changes
        delts = delts + dhs
        rmvi = dhit
        do layer = 1,ni
          sumr = max( -delti(layer), rmvi )
          rmvi = rmvi - sumr
          delti(layer) = delti(layer) + sumr
        enddo
      endif

      if (fixice) return

      dtop = dhit + subi

      !-----------------------------------------------------------------
      ! melt/grow at bot srf
      !-----------------------------------------------------------------
      focn = fbot
      ebot = dtsub * ( condb - fbot )
      if (ebot .le. 0. ) then            !    grow at bottom
        qigrow = energ( tbot, salnew ) 
        dhib = ebot/qigrow
      else                               !    melt at bottom
        qigrow = c0                      !    on purpose
        if ( (enet+ebot) .ge. 0. ) then  !    remotely possible
          dhib = -( hi + dtop )
          dhs  = -( hs + subs )
          focn = condb + enet / dtsub
          return
        endif
!        call t_startf ('botmelt')
        call botmelt(   qi,   qs, delti, delts, 
     $                dhib,  dhs,  ebot  )
!        call t_stopf ('botmelt')
      endif

      !-----------------------------------------------------------------
      ! stop if error occurred in srfsub, srfmelt or botmelt
      !-----------------------------------------------------------------
      if (errflag.ne.0) then
!         call print_state('state at dh error   ',i,j)
         stop
      endif

      !-----------------------------------------------------------------
      ! check to see if there is any ice left after top/bottom 
      ! melt/growth
      !-----------------------------------------------------------------
      hi_tw = hi+dhib+dtop

      if ( hi_tw .le. 0. ) then
         dhib = -(hi+dtop)
         ! convert any snow that might be left into sea ice
         dhsf = min(-(hs+dhs),c0)
         dhif = -rhos / rhoi * dhsf
         do layer = 1,ni
           qi(layer) = qs * rhoi / rhos
         enddo
      else

      !-----------------------------------------------------------------
      ! flooding
      !-----------------------------------------------------------------
!        call t_startf ('freeboard')
        call freeboard(hs,hi,dhs,qs,dhsf,dhif,qiflood)
!        call t_stopf ('freeboard')
        if (verbos) write(6,*) 'fld',dhib,dtop,dhif,dhsf

      !-----------------------------------------------------------------
      ! adjust layers
      !-----------------------------------------------------------------
!        call t_startf ('adjust')
        call adjust(hi,dhib,dtop,dhif,dhsf,qiflood,qigrow,qi) 
!        call t_stopf ('adjust')

      endif

      end subroutine dh

c=======================================================================

      subroutine freeboard(hs,hi,dhs,qs,dhsf,dhif,qiflood)

!---!-------------------------------------------------------------------
!---! freeboard adjustment due to flooding ... snow-ice formation
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &    hi  ! initial ice thickness                    (m)
     &,   hs  ! initial snow thickness                   (m)
     &,   dhs ! snow top, dhit<=0                        (m)
     &,   qs  ! energy of melt of snow per unit vol. (J/m**3)

      real (kind=dbl_kind), intent(out) ::
     &    dhif ! ice top from flooding, dhif>0            (m)
     &,   dhsf ! snow from flooding, dhsf<0               (m)
     &,   qiflood ! energy of melt of flooded ice    (W/m**2)

      real (kind=dbl_kind) :: zintfc  ! height of snow/ice interf. wrt ocn (m)

      zintfc  =  hi - (rhos*hs + rhoi*hi)/rhow
      if (( zintfc .lt. 0. ).and.(hs+dhs.gt.0.)) then
        dhsf    =  rhoi / rhos * zintfc
        dhsf    =  max( dhsf, -(hs+dhs) )
        dhif    =  -rhos / rhoi * dhsf
        qiflood =  qs * rhoi / rhos
      endif

      end subroutine freeboard

c=======================================================================

      subroutine srfsub(  qi,   qs, delti, delts,
     $                   subi, subs,  etop,  enet )

!---!-------------------------------------------------------------------
!---! compute the sea ice and snow thickness changes from  
!---! sublimation/condensation
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &   qi (1:plevmx), qs      ! energy of melt of ice/snow per vol  (J/m**3)
     &,  delti(plevmx), delts   ! thickness of ice/snow layer              (m)

      real (kind=dbl_kind), intent(out) :: 
     &   subi, subs  ! subl/cond. amount for ice/snow           (m)

      real (kind=dbl_kind), intent(inout) :: 
     &   etop                 ! energy avail to sub/cond ice/snow   (J/m**2)
     &,  enet                 ! energy needed to melt all ice/snow  (J/m**2)

      real (kind=dbl_kind) :: 
     &   ue                   ! energy of melt + vapor for ice/snow (J/m**3)
     &,  subil                ! subl from this layer

      integer (kind=int_kind) :: layer

      subs = c0
      subi = c0
      if ( delts .gt. 0. ) then
        ! convert etop into snow thickness
        ! positive if cond/negative if subl
        ue = qs - rLvs               
        subs = etop/ue     
        if ( (delts + subs) .ge. 0. ) then
          ! either all condensation becomes snow
          ! or sublimate some snow but no ice
          etop = c0
          enet = enet + subs * qs
          return
        else
          ! sublimate all of the snow and some ice too
          subs = -delts
          etop = etop + delts * ue
          enet = enet + subs  * qs
        endif
      endif
      do layer = 1,ni
        ! convert etop into ice thickness
        ! positive if cond/negative if subl
        ue = qi(layer) - rLvi
        subil = etop/ue
        if ( (delti(layer)+subil) .ge. 0. ) then
          ! either all condensation becomes ice 
          ! or sublimate some ice from this layer
          subi = subi + subil
          etop = c0
          enet = enet + subil * qi(layer)
          return
        else
          ! sublimate all of the layer 
          subi = subi - delti(layer)  
          etop = etop + delti(layer) * ue
          enet = enet - delti(layer) * qi(layer)
        endif
      enddo
      
      if (iam.eq.0) then
      write(6,*)  'ERROR in srfsub',etop
      end if
      errflag = 1

      end subroutine srfsub
      
c=======================================================================

      subroutine srfmelt(   qi,   qs, delti, delts,
     $                    dhit,  dhs,  etop  ) 

!---!-------------------------------------------------------------------
!---! melt ice/snow from the top srf
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &   qi (1:plevmx), qs      ! energy of melt of ice/snow per vol (J/m**3)
     &,  delti(plevmx), delts   ! thickness of ice/snow layer             (m)

      real (kind=dbl_kind), intent(out) :: dhit, dhs ! ice/snow thickness change               (m)

      real (kind=dbl_kind), intent(inout) :: etop    ! energy avail to melt ice and snow  (J/m**2)

      real (kind=dbl_kind) :: dhitl  ! melt from this layer          (m)
      integer (kind=int_kind) :: layer

      dhit = c0
      dhs  = c0

      if ( delts .gt. 0. ) then
        ! convert etop into snow thickness
        dhs  =  etop/qs
        if ( (delts + dhs) .ge. 0. ) then
          ! melt only some of the snow
          etop = c0
          return
        else
          ! melt all of the snow and some ice too
          dhs  = -delts
          etop = etop+qs*delts
        endif
      endif

      do layer = 1,ni
        ! convert etop into ice thickness
        dhitl = etop/qi(layer)
        if ( (delti(layer)+dhitl) .ge. 0. ) then
          ! melt some ice from this layer
          dhit = dhit + dhitl
          etop = c0
          return
        else
          ! melt all of the ice in this layer
          dhit = dhit - delti(layer)   
          etop = etop + delti(layer) * qi(layer)
        endif
      enddo
      
      write(6,*)  'ERROR in srfmelt',etop
      errflag = 1

      end subroutine srfmelt
      
c=======================================================================

      subroutine botmelt(   qi,   qs, delti, delts,
     $                    dhib,  dhs,  ebot  )

!---!-------------------------------------------------------------------
!---! melt from bottom
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &   qi (1:plevmx), qs      ! energy of melt of ice/snow per vol (J/m**3)
     &,  delti(plevmx), delts   ! thickness of ice/snow layer             (m)

      real (kind=dbl_kind), intent(out) :: dhib  ! ice thickness change (m)

      real (kind=dbl_kind), intent(inout) :: 
     &   ebot                 ! energy avail to melt ice and snow  (J/m**2)
     &,  dhs                  ! snow thickness change                   (m)

      real (kind=dbl_kind) :: 
     &   dhibl                ! melt from this ice layer                (m)
     &,  dhsl                 ! melt from this snow layer               (m)
      integer (kind=int_kind) :: layer

      dhib = c0

      do layer = ni,1,-1
         dhibl = ebot/qi(layer)
         if ( (delti(layer)+dhibl) .ge. 0. ) then
            dhib = dhib + dhibl
            ebot = c0
            return
         else
            dhib = dhib - delti(layer)   
            ebot = ebot + delti(layer) * qi(layer)
         endif
      enddo

      ! finally melt snow if necessary
      dhsl = ebot/qs
      if ( (delts + dhsl) .ge. 0. ) then
         dhs = dhs + dhsl
         ebot = c0
         return
      endif

      write(6,*)  'ERROR in botmelt',ebot
      errflag = 1

      end subroutine botmelt

c=======================================================================

      subroutine adjust(hi0,dhib,dhit,dhif,dhsf,qiflood,qigrow,qi_tw)

!---!-------------------------------------------------------------------
!---! Adjusts temperature profile to account for changing
!---! the layer spacing due to growth/melt (incl. subl/cond, flooding)
!---! At start the energy of melting was computed
!---! after updating tiz from the heat equation
!---! hi is the thickness prior to changes from dhib and dhit
!---! hi_tw is the thickness after making these changes
!---! dhib<0 if there is melt at the bottom
!---! dhit<0 if there is melt at the top
!---! generally _tw is a suffix to label the adjusted variables
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: 
     &   hi0                  ! initial ice thickness                    (m)
     &,  dhib                 ! ice bot, dhib<0 if melt                  (m)
     &,  dhit                 ! ice top, dhit<=0                         (m)
     &,  dhif                 ! ice top from flooding, dhif>0            (m)
     &,  dhsf                 ! snow from flooding, dhsf<0               (m)

     &,  qiflood              ! qi for flooded ice                  (J/m**3)
     &,  qigrow               ! qi for ice growing on bot           (J/m**3)

      real (kind=dbl_kind), intent(inout) :: 
     &   qi_tw(plevmx)          ! energy of melt of ice per unit vol. (J/m**3)


      ! the following pairs are before & after adjusting
      real (kind=dbl_kind) ::     hi, hi_tw          ! ice thickness   (m)
     &,      z(plevmx+2)        ! vertical layer position        (m)
     &,      z_tw(plevmx+1)     ! vertical layer position        (m)
     &,      D, D_tw          ! layer thickness                (m)
      integer (kind=int_kind) ::  k, k_tw          ! index of layer

      real (kind=dbl_kind) :: qi(plevmx)  ! qi of initial layers  (J/m**3)
     &,  fract(plevmx,plevmx)     ! fract of layer k that makes up layer k_tw
      integer (kind=int_kind) :: layer             

      ! first check to see if there is any ice left after top/bottom 
      ! melt/growth
      hi = hi0
      hi_tw = hi0+dhib+dhit

      ! there is ice left so allow snow-ice conversion from flooding
      hi_tw = hi_tw+dhif

      if ( (abs(dhib) .gt. puny) .or. (abs(dhit) .gt. puny)
     $     .or. (dhif.gt.0)) then
         do layer = 1,ni
            qi(layer) = qi_tw(layer)
         enddo
      
         ! layer thickness      
         D = hi/ni 
         D_tw = hi_tw/ni 

         ! z is positive down and zero is relative to the top 
         ! of the ice from the old time step
         z(1) = -dhit          ! necessary
         z_tw(1) = -dhit-dhif
         do layer = 2,ni
           z(layer) = D*(layer-1)
           z_tw(layer) = z_tw(1)+D_tw*(layer-1) 
         enddo
         z(ni+1) = hi + min(dhib,c0)
         z_tw(ni+1) = z_tw(1)+hi_tw

         do k_tw = 1,ni
           qi_tw(k_tw) = c0
           do k = 1,ni
             fract(k,k_tw) = ( min( z_tw(k_tw+1), z(k+1)) -
     $                         max( z_tw(k_tw)  , z(k)  )   )
             if (fract(k,k_tw).gt.0.)
     $            qi_tw(k_tw) = qi_tw(k_tw)+fract(k,k_tw)*qi(k)
           enddo
         enddo

         if (dhif.gt.D_tw) then
           do k = 1,ni
             qi_tw(k) = qi_tw(k) + qiflood*max(c0,
     $                   min(D_tw,dhif-D_tw*(k-1)))
           enddo
         else
           qi_tw(1) = qi_tw(1) + dhif*qiflood
         endif

         if (dhib.gt.D_tw) then
           z(ni+2) = hi + dhib
           do k = 1,ni
             qi_tw(k) = qi_tw(k) +  qigrow*max(c0,
     $                   (min(z(ni+2),Z_tw(k+1))-max(hi,Z_tw(k))))
           enddo
         else
           qi_tw(ni) = qi_tw(ni) + max(dhib,c0)*qigrow
         endif

         do k_tw = 1,ni
           qi_tw(k_tw) = qi_tw(k_tw)/D_tw
         enddo

      endif

      end subroutine adjust

c=======================================================================

      real function energ(Tmp ,sal)

!---!-------------------------------------------------------------------
!---! compute the energy of melting per unit volume  (J/m**3)
!---! relative to melting (negative quantity)
!---!-------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) ::
     &    Tmp          ! midpt temperature of ice layer       (C)
     &,   sal          ! midpt salinity of ice layer        (ppt)  

      energ = -rLfi - rcpi*(-depressT*sal-Tmp)-rLfidepressT*sal/Tmp

      end function energ

c=======================================================================

      end module ice_dh

c=======================================================================


