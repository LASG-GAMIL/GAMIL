module tp_core
!BOP
!
! !MODULE: tp_core --- Utilities for the transport core
!
!
! !PUBLIC MEMBER FUNCTIONS:
      public tp2c, tp2d, xtp, fxppm, xmist, steepx, lmppm
      public huynh, ytp, ymist, fyppm, tpcc, ycc
!
! !DESCRIPTION:
!
!      This module provides 
!
!      \begin{tabular}{|l|l|} \hline \hline
!       tp2c  &   \\ \hline
!       tp2d  &   \\ \hline 
!       xtp  &   \\ \hline 
!       fxppm  &   \\ \hline 
!       xmist  &   \\ \hline 
!       steepx  &   \\ \hline 
!       lmppm  &   \\ \hline 
!       huynh  &   \\ \hline 
!       ytp  &   \\ \hline 
!       ymist  &   \\ \hline 
!       fyppm  &   \\ \hline 
!       tpcc  &   \\ \hline 
!       ycc  &   \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   01.01.15   Lin        Routines coalesced into this module
!   01.03.26   Sawyer     Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tp2c --- Perform transport on a C grid
!
! !INTERFACE: 
 subroutine tp2c(dh, va, h, crx, cry, im, jm,                            &
                 iord, jord, ng, fx, fy, ffsl,                           &
                 rcap, acosp, xfx, yfx, cosp, id, jfirst, jlast)
!-----------------------------------------------------------------------

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
   integer im, jm                  ! Dimensions
   integer jfirst, jlast           ! Latitude strip
   integer iord, jord              ! Interpolation order in x,y
   integer ng                      ! Max. NS dependencies
   integer id                      ! density (0)  (mfx = C)
   real (r8) rcap                  ! Ask S.-J. (polar constant?)
   real (r8) acosp(jm)             ! Ask S.-J. (difference to cosp??)
   logical ffsl(jm)                ! Use flux-form semi-Lagrangian trans.?
                                      ! (N*NG S*NG)
   real (r8) cosp(jm)                   ! Critical angle
   real (r8) va(im,jfirst:jlast)        ! Courant  (unghosted)
   real (r8) h(im,jfirst-ng:jlast+ng)   ! Pressure ( N*NG S*NG )
   real (r8) crx(im,jfirst-ng:jlast+ng) ! Ask S.-J. ( N*NG S*NG )
   real (r8) cry(im,jfirst:jlast+1)     ! Ask S.-J. ( N like FY )
   real (r8) xfx(im,jfirst:jlast)       ! Ask S.-J. ( unghosted like FX )
   real (r8) yfx(im,jfirst:jlast+1)     ! Ask S.-J. ( N like FY )

! !OUTPUT PARAMETERS:
   real (r8) dh(im,jfirst:jlast)        ! Ask S.-J. ( unghosted )
   real (r8) fx(im,jfirst:jlast)        ! Flux in x ( unghosted )
   real (r8) fy(im,jfirst:jlast+1)      ! Flux in y ( N, see tp2c )

! !DESCRIPTION:
!     Perform transport on a C grid.   The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!     Ask S.-J. how exactly this differs from TP2C.
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   integer i, j, js2g0, jn2g0
   real (r8)  sum1

   js2g0  = max(2,jfirst)          !  No ghosting
   jn2g0  = min(jm-1,jlast)        !  No ghosting

   call tp2d(va, h, crx, cry, im, jm, iord, jord, ng,fx, fy, ffsl,    &
             xfx, yfx, cosp, id, jfirst, jlast)

   do j=js2g0,jn2g0
      do i=1,im-1
         dh(i,j) = fx(i,j) - fx(i+1,j) + (fy(i,j)-fy(i,j+1))*acosp(j)
      enddo
   enddo

   do j=js2g0,jn2g0
      dh(im,j) = fx(im,j) - fx(1,j) + (fy(im,j)-fy(im,j+1))*acosp(j)
   enddo

! Poles
   if ( jfirst ==  1 ) then
!       sum1 = - SUM( fy(1:im, 2) ) * rcap
        sum1 = 0.
        do i=1,im
          sum1 = sum1 + fy(i,2)
        enddo
          sum1 = -sum1*rcap
        do i=1,im
          dh(i, 1) = sum1
        enddo
   endif
   
   if ( jlast == jm ) then
!       sum1 = SUM( fy(1:im,jm) ) * rcap
        sum1 = 0.
        do i=1,im
          sum1 = sum1 + fy(i,jm)
        enddo
          sum1 = sum1*rcap
        do i=1,im
          dh(i,jm) = sum1
        enddo
   endif
   return
!EOC
 end subroutine tp2c
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tp2d --- Perform transport on a D grid
!
! !INTERFACE: 
 subroutine tp2d(va, q, crx, cry, im, jm, iord, jord, ng, fx, fy,        &
                 ffsl, xfx, yfx, cosp, id, jfirst, jlast)
!-----------------------------------------------------------------------
! !USES:

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
   integer im, jm                    ! Dimensions
   integer jfirst, jlast             ! Latitude strip
   integer iord, jord                ! Interpolation order in x,y
   integer ng                        ! Max. NS dependencies
   integer id                        ! density (0)  (mfx = C)
                                     ! mixing ratio (1) (mfx = mass flux)
   logical ffsl(jm)                  ! Use flux-form semi-Lagrangian trans.?
                                     ! ghosted N*ng S*ng
   real (r8) cosp(jm)                     ! Critical angle
   real (r8) va(im,jfirst:jlast)          ! Courant  (unghosted)
   real (r8) q(im,jfirst-ng:jlast+ng)     ! transported scalar ( N*NG S*NG )
   real (r8) crx(im,jfirst-ng:jlast+ng)   ! Ask S.-J. ( N*NG S*NG )
   real (r8) cry(im,jfirst:jlast+1)       ! Ask S.-J. ( N like FY )
   real (r8) xfx(im,jfirst:jlast)         ! Ask S.-J. ( unghosted like FX )
   real (r8) yfx(im,jfirst:jlast+1)       ! Ask S.-J. ( N like FY )

! !OUTPUT PARAMETERS:
   real (r8) fx(im,jfirst:jlast)          ! Flux in x ( unghosted )
   real (r8) fy(im,jfirst:jlast+1)        ! Flux in y ( N, see tp2c )

! !DESCRIPTION:
!     Perform transport on a D grid.   The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!
!
! !REVISION HISTORY:
!   WS  99.04.13:  Added jfirst:jlast concept
!       99.04.21:  Removed j1 and j2 (j1=2, j2=jm-1 consistently)
!       99.04.27:  Removed dc, wk2 as arguments (local to YTP)
!       99.04.27:  Removed adx as arguments (local here)
!   SJL 99.07.26:  ffsl flag added
!   WS  99.09.07:  Restructuring, cleaning, documentation
!   WS  99.10.22:  NG now argument; arrays pruned
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
   integer i, j, iad, jp, js2g0, js2gng, jn2g0, jn2gng
   real (r8) adx(im,jfirst-ng:jlast+ng)
   real (r8) wk1(im)
   real (r8)   dm(-im/3:im+im/3)
   real (r8) qtmp(-im/3:im+im/3)
   real (r8)   al(-im/3:im+im/3)
   real (r8)   ar(-im/3:im+im/3)
   real (r8)   a6(-im/3:im+im/3)

! Number of ghost latitudes
    js2g0  = max(2,jfirst)          !  No ghosting
    js2gng = max(2,jfirst-ng)       !  Number needed on S
    jn2g0  = min(jm-1,jlast)        !  No ghosting
    jn2gng = min(jm-1,jlast+ng)     !  Number needed on N
    iad = 1

!!!      do j=2,jm-1
    do j=js2gng,jn2gng               !  adx needed on N*ng S*ng

       call xtp(im,  ffsl(j), wk1, q(1,j),                &
                crx(1,j), iad, crx(1,j), cosp(j), 0,      &
                dm, qtmp, al, ar, a6)

       do i=1,im-1
          adx(i,j) = q(i,j) + 0.5 *                       &
                     (wk1(i)-wk1(i+1) + q(i,j)*(crx(i+1,j)-crx(i,j)))
       enddo
          adx(im,j) = q(im,j) + 0.5 *                     &
                      (wk1(im)-wk1(1) + q(im,j)*(crx(1,j)-crx(im,j)))
    enddo

! WS 99.09.07 : Split up north and south pole

     if ( jfirst-ng <= 1 ) then
        do i=1,im 
          adx(i, 1) = q(i,1)
        enddo
     endif 
     if ( jlast+ng >= jm ) then
        do i=1,im 
          adx(i,jm) = q(i,jm)
        enddo
     endif

     call ytp(im,jm,fy, adx,cry,yfx,ng,jord,0,jfirst,jlast)

!!!   do j=2,jm-1
      do j=js2g0,jn2g0
        do i=1,im
           jp = j-va(i,j)
           wk1(i) = q(i,j) +0.5*va(i,j)*(q(i,jp)-q(i,jp+1))
        enddo

        call xtp(im,  ffsl(j), fx(1,j), wk1,                  &
                 crx(1,j), iord, xfx(1,j), cosp(j), id,       &
                 dm, qtmp, al, ar, a6)
      enddo
    return
!EOC
 end subroutine tp2d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xtp
!
! !INTERFACE: 
 subroutine xtp(im, ffsl,  fx,  q,  c,  iord,  mfx,            &
                cosa, id, dm, qtmp, al, ar, a6)
!-----------------------------------------------------------------------

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none
 
! !INPUT PARAMETERS:
   integer id               ! ID = 0: density (mfx = C)
                            ! ID = 1: mixing ratio (mfx is mass flux)

   integer im               ! Total longitudes
   real (r8) c(im)          ! Courant numbers
   real (r8) q(im)
   real (r8) mfx(im)
   logical ffsl
   integer iord
   real (r8) cosa

! !INPUT/OUTPUT PARAMETERS:
   real (r8) qtmp(-im/3:im+im/3)   ! Input work arrays:
   real (r8)   dm(-im/3:im+im/3)
   real (r8)   al(-im/3:im+im/3)
   real (r8)   ar(-im/3:im+im/3)
   real (r8)   a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
   real (r8) fx(im)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! Local:
   real (r8)       cos_upw               !critical cosine for upwind
   real (r8)       cos_van               !critical cosine for van Leer
   real (r8)       cos_ppm               !critical cosine for ppm

   parameter (cos_upw = 0.05)       !roughly at 87 deg.
   parameter (cos_van = 0.25)       !roughly at 75 deg.
   parameter (cos_ppm = 0.25)

   integer i, imp
   real (r8) qmax, qmin
   real (r8) rut, tmp
   integer iu, itmp, ist
   integer isave(im)
   integer iuw, iue

   imp = im + 1

   do i=1,im
      qtmp(i) = q(i)
   enddo

   if( ffsl ) then
! Flux-Form Semi-Lagrangian transport

! Figure out ghost zone for the western edge:
      iuw =  -c(1)
      iuw = min(0, iuw)
 
      do i=iuw, 0
         qtmp(i) = q(im+i)
      enddo 

! Figure out ghost zone for the eastern edge:
      iue = im - c(im)
      iue = max(imp, iue)

      do i=imp, iue
         qtmp(i) = q(i-im)
      enddo

      if( iord == 1 .or. cosa < cos_upw) then
      do i=1,im
        iu = c(i)
      if(c(i) .le. 0.) then
        itmp = i - iu
        isave(i) = itmp - 1
      else
        itmp = i - iu - 1
        isave(i) = itmp + 1
      endif
        fx(i) = (c(i)-iu) * qtmp(itmp)
      enddo
      else

      do i=1,im
! 2nd order slope
         tmp = 0.25*(qtmp(i+1) - qtmp(i-1))
         qmax = max(qtmp(i-1), qtmp(i), qtmp(i+1)) - qtmp(i)
         qmin = qtmp(i) - min(qtmp(i-1), qtmp(i), qtmp(i+1))
         dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
      enddo

 
      do i=iuw, 0
         dm(i) = dm(im+i)
      enddo 

      do i=imp, iue
         dm(i) = dm(i-im)
      enddo

      if(iord .ge. 3 .and. cosa .gt. cos_ppm) then
         call fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6,         &
                    iuw, iue, ffsl, isave)
      else
      do i=1,im
            iu  = c(i)
            rut = c(i) - iu
         if(c(i) .le. 0.) then
            itmp = i - iu
            isave(i) = itmp - 1
            fx(i) = rut*(qtmp(itmp)-dm(itmp)*(1.+rut))
         else
            itmp = i - iu - 1
            isave(i) = itmp + 1
            fx(i) = rut*(qtmp(itmp)+dm(itmp)*(1.-rut))
         endif
      enddo
      endif

      endif

      do i=1,im
      if(c(i) .ge. 1.) then
        do ist = isave(i),i-1
           fx(i) = fx(i) + qtmp(ist)
        enddo
      elseif(c(i) .le. -1.) then
        do ist = i,isave(i)
           fx(i) = fx(i) - qtmp(ist)
        enddo
      endif
      enddo

      if(id .ne. 0) then
         do i=1,im
            fx(i) =  fx(i)*mfx(i)
         enddo
      endif

   else
! Regular PPM (Eulerian without FFSL extension)

      qtmp(imp) = q(1)
      qtmp(  0) = q(im)

      if(iord == 1 .or. cosa < cos_upw) then
         do i=1,im
            iu = float(i) - c(i)
            fx(i) = mfx(i)*qtmp(iu)
         enddo
      else

         qtmp(-1)    = q(im-1)
         qtmp(imp+1) = q(2)

         if(iord > 0 .or. cosa < cos_van) then
            call xmist(im, qtmp, dm, 2)
         else
            call xmist(im, qtmp, dm, iord)
         endif

         dm(0) = dm(im)

         if( abs(iord).eq.2 .or. cosa .lt. cos_van ) then
            do i=1,im
               iu = float(i) - c(i)
               fx(i) =  mfx(i)*(qtmp(iu)+dm(iu)*(sign(1.,c(i))-c(i)))

!              if(c(i) .le. 0.) then
!                 fx(i) = qtmp(i) - dm(i)*(1.+c(i))
!              else
!                 fx(i) = qtmp(i-1) + dm(i-1)*(1.-c(i))
!              endif
!                 fx(i) = fx(i)*mfx(i)

            enddo
         else
            call fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6,       &
                       iuw, iue, ffsl, isave)
         endif
      endif

   endif
   return
!EOC
 end subroutine xtp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: xmist
!
! !INTERFACE: 
 subroutine xmist(im,  q,  dm,  id)
!-----------------------------------------------------------------------

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im                   ! Total number of longitudes
 integer id                   ! ID = 0: density (mfx = C)
                              ! ID = 1: mixing ratio (mfx is mass flux)
 real(r8)  q(-im/3:im+im/3)   ! Input latitude 

! !OUTPUT PARAMETERS:
 real(r8) dm(-im/3:im+im/3)   ! 

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC

 real(r8) r24
 parameter( r24 = 1./24.)

 integer i
 real(r8) qmin, qmax

    if(id .le. 2) then
       do i=1,im
          dm(i) = r24*(8.*(q(i+1) - q(i-1)) + q(i-2) - q(i+2))
       enddo
    else
       do i=1,im
          dm(i) = 0.25*(q(i+1) - q(i-1))
       enddo
    endif

    if( id < 0 ) return

! Apply monotonicity constraint (Lin et al. 1994, MWR)
      do i=1,im
         qmax = max( q(i-1), q(i), q(i+1) ) - q(i)
         qmin = q(i) - min( q(i-1), q(i), q(i+1) )
         dm(i) = sign( min(abs(dm(i)), qmax, qmin), dm(i) )
      enddo
  return
!EOC
 end subroutine xmist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fxppm
!
! !INTERFACE: 
 subroutine fxppm(im, c, mfx,  p, dm, fx, iord, al, ar, a6,        &
                  iuw, iue, ffsl, isave)
!-----------------------------------------------------------------------
!
! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im, iord
 real (r8)  c(im)
 real (r8) p(-im/3:im+im/3)
 real (r8) dm(-im/3:im+im/3)
 real (r8) mfx(im)
 integer iuw, iue
 logical ffsl
 integer isave(im)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) al(-im/3:im+im/3)
 real (r8) ar(-im/3:im+im/3)
 real (r8) a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
 real (r8) fx(im)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r3, r23
 parameter ( r3 = 1./3., r23 = 2./3. )

 integer i, lmt
 integer iu, itmp
 real (r8) ru
 logical steep

  if( iord == 6 ) then
      steep = .true.
  else
      steep = .false.
  endif

  do i=1,im
     al(i) = 0.5*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
  enddo

  if( steep ) call steepx( im, p, al(1), dm )

     do i=1,im-1
        ar(i) = al(i+1)
     enddo
        ar(im) = al(1)

  if(iord == 7) then
     call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
  else
     if(iord .eq. 3 .or. iord .eq. 5) then
         do i=1,im
            a6(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
         enddo
     endif
     lmt = iord - 3
     call lmppm( dm(1), a6(1), ar(1), al(1), p(1), im, lmt )
  endif

  if( ffsl ) then

      do i=iuw, 0
         al(i) = al(im+i)
         ar(i) = ar(im+i)
         a6(i) = a6(im+i)
      enddo

      do i=im+1, iue
         al(i) = al(i-im)
         ar(i) = ar(i-im)
         a6(i) = a6(i-im)
      enddo

      do i=1,im
            iu = c(i)
            ru = c(i) - iu
         if(c(i) .gt. 0.) then
            itmp = i - iu - 1
            isave(i) = itmp + 1
            fx(i) = ru*(ar(itmp)+0.5*ru*(al(itmp)-ar(itmp) +     &
                        a6(itmp)*(1.-r23*ru)) )
         else
            itmp = i - iu
            isave(i) = itmp - 1
            fx(i) = ru*(al(itmp)-0.5*ru*(ar(itmp)-al(itmp) +     &
                        a6(itmp)*(1.+r23*ru)) )
         endif
      enddo

  else
         al(0) = al(im)
         ar(0) = ar(im)
         a6(0) = a6(im)
      do i=1,im
         if(c(i) .gt. 0.) then
            fx(i) = ar(i-1) + 0.5*c(i)*(al(i-1) - ar(i-1) +   &
                    a6(i-1)*(1.-r23*c(i)) )
      else
            fx(i) = al(i) - 0.5*c(i)*(ar(i) - al(i) +         &
                    a6(i)*(1.+r23*c(i)))
      endif
            fx(i) = mfx(i) * fx(i)
      enddo
  endif
  return
!EOC
 end subroutine fxppm
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: steepx
!
! !INTERFACE: 
 subroutine  steepx(im, p, al, dm)
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im
 real (r8)  p(-im/3:im+im/3)
 real (r8) dm(-im/3:im+im/3)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) al(im)

! !DESCRIPTION:
!   
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i
 real (r8) r3
 parameter ( r3 = 1./3. )

 real (r8) dh(0:im)
 real (r8) d2(0:im+1)
 real (r8) eta(0:im)
 real (r8) xxx, bbb, ccc

   do i=0,im
      dh(i) = p(i+1) - p(i)
   enddo

! Needs dh(0:im)
   do i=1,im
      d2(i) = dh(i) - dh(i-1)
   enddo
   d2(0) = d2(im)
   d2(im+1) = d2(1)

! needs p(-1:im+2), d2(0:im+1)
   do i=1,im
      if( d2(i+1)*d2(i-1).lt.0. .and. p(i+1).ne.p(i-1) ) then
          xxx    = 1. - 0.5 * ( p(i+2) - p(i-2) ) / ( p(i+1) - p(i-1) )
          eta(i) = max(0., min(xxx, 0.5) )
      else
          eta(i) = 0.
      endif
    enddo

    eta(0) = eta(im)

! needs eta(0:im), dh(0:im-1), dm(0:im)
   do i=1,im
      bbb = ( 2.*eta(i  ) - eta(i-1) ) * dm(i-1) 
      ccc = ( 2.*eta(i-1) - eta(i  ) ) * dm(i  ) 
      al(i) = al(i) + 0.5*( eta(i-1) - eta(i)) * dh(i-1) + (bbb - ccc) * r3
   enddo
   return
!EOC
 end subroutine steepx
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: lmppm
!
! !INTERFACE: 
 subroutine lmppm(dm, a6, ar, al, p, im, lmt)
!-----------------------------------------------------------------------

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im   ! Total longitudes
 integer lmt  ! LMT = 0: full monotonicity
              ! LMT = 1: Improved and simplified full monotonic constraint
              ! LMT = 2: positive-definite constraint
              ! LMT = 3: Quasi-monotone constraint
 real(r8) p(im)
 real(r8) dm(im)

! !OUTPUT PARAMETERS:
 real(r8) a6(im)
 real(r8) ar(im)
 real(r8) al(im)

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 real (r8) r12
 parameter ( r12 = 1./12. )

 real (r8) da1, da2, fmin, a6da
 real (r8) dr, dl

 integer i

! LMT = 0: full monotonicity
! LMT = 1: Improved and simplified full monotonic constraint
! LMT = 2: positive-definite constraint
! LMT = 3: Quasi-monotone constraint

  if( lmt == 0 ) then

! Full constraint
  do i=1,im
     if(dm(i) .eq. 0.) then
         ar(i) = p(i)
         al(i) = p(i)
         a6(i) = 0.
     else
         da1  = ar(i) - al(i)
         da2  = da1**2
         a6da = a6(i)*da1
         if(a6da .lt. -da2) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
         elseif(a6da .gt. da2) then
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
         endif
     endif
  enddo

  elseif( lmt == 1 ) then

! Improved (Lin 2001?) full constraint
      do i=1,im
           da1 = dm(i) + dm(i)
            dl = sign(min(abs(da1),abs(al(i)-p(i))), da1)
            dr = sign(min(abs(da1),abs(ar(i)-p(i))), da1)
         ar(i) = p(i) + dr
         al(i) = p(i) - dl
         a6(i) = 3.*(dl-dr)
      enddo

  elseif( lmt == 2 ) then
! Positive definite constraint
      do 250 i=1,im
      if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 250
      fmin = p(i) + 0.25*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
      if(fmin.ge.0.) go to 250
      if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
            ar(i) = p(i)
            al(i) = p(i)
            a6(i) = 0.
      elseif(ar(i) .gt. al(i)) then
            a6(i) = 3.*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
      else
            a6(i) = 3.*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
      endif
250   continue

  elseif(lmt .eq. 3) then
! Quasi-monotone constraint
      do i=1,im
         da1 = 4.*dm(i)
          dl = sign(min(abs(da1),abs(al(i)-p(i))), da1)
          dr = sign(min(abs(da1),abs(ar(i)-p(i))), da1)
         ar(i) = p(i) + dr
         al(i) = p(i) - dl
         a6(i) = 3.*(dl-dr)
      enddo
  endif
  return
!EOC
 end subroutine lmppm
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: huynh --- Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !INTERFACE: 
 subroutine huynh(im, ar, al, p, d2, d1)
!-----------------------------------------------------------------------

! !USES:

 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im
 real(r8)  p(im)

! !OUTPUT PARAMETERS:
 real(r8) ar(im)
 real(r8) al(im)
 real(r8) d2(im)
 real(r8) d1(im)

! !DESCRIPTION:
!
!   Enforce Huynh's 2nd constraint in 1D periodic domain
!
! !REVISION HISTORY:
!   99.01.01 Lin      Creation
!   01.03.27 Sawyer   Additional ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer  i
 real(r8) pmp
 real(r8) lac
 real(r8) pmin
 real(r8) pmax

! Compute d1 and d2
      d1(1) = p(1) - p(im)
      do i=2,im
         d1(i) = p(i) - p(i-1)
      enddo

      do i=1,im-1
         d2(i) = d1(i+1) - d1(i)
      enddo
      d2(im) = d1(1) - d1(im)

! Constraint for AR
!            i = 1
         pmp   = p(1) + 2.0 * d1(1)
         lac   = p(1) + 0.5 * (d1(1)+d2(im)) + d2(im) 
         pmin  = min(p(1), pmp, lac)
         pmax  = max(p(1), pmp, lac)
         ar(1) = min(pmax, max(ar(1), pmin))

      do i=2, im
         pmp   = p(i) + 2.0*d1(i)
         lac   = p(i) + 0.5*(d1(i)+d2(i-1)) + d2(i-1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         ar(i) = min(pmax, max(ar(i), pmin))
      enddo

! Constraint for AL
      do i=1, im-1
         pmp   = p(i) - 2.0*d1(i+1)
         lac   = p(i) + 0.5*(d2(i+1)-d1(i+1)) + d2(i+1)
         pmin  = min(p(i), pmp, lac)
         pmax  = max(p(i), pmp, lac)
         al(i) = min(pmax, max(al(i), pmin))
      enddo

! i=im
         i = im
         pmp    = p(im) - 2.0*d1(1)
         lac    = p(im) + 0.5*(d2(1)-d1(1)) + d2(1)
         pmin   = min(p(im), pmp, lac)
         pmax   = max(p(im), pmp, lac)
         al(im) = min(pmax, max(al(im), pmin))

! compute A6 (d2)
      do i=1, im
         d2(i) = 3.*(p(i)+p(i)  - (al(i)+ar(i)))
      enddo
    return
!EOC
 end subroutine huynh
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ytp
!
! !INTERFACE: 
 subroutine ytp(im, jm, fy, q, c, yfx, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  Max. NS dependencies
 integer jord                        !  order of subgrid dist
 integer iv                          !  Scalar=0, Vector=1
 real (r8) q(im,jfirst-ng:jlast+ng)       !  advected scalar N*jord S*jord
 real (r8) c(im,jfirst:jlast+1)           !  Courant   N (like FY)
 real (r8) yfx(im,jfirst:jlast+1)         !  Backgrond mass flux

! !OUTPUT PARAMETERS:
 real (r8) fy(im,jfirst:jlast+1)          !  Flux      N (see tp2c)

! !DESCRIPTION:
!     This routine calculates the flux FX.  The method chosen
!     depends on the order of the calculation JORD (currently
!     1, 2 or 3).  
!
! !CALLED FROM:
!     cd_core
!     tp2d
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.13:  Added jfirst:jlast concept
!  WS  99.04.21:  Removed j1 and j2 (j1=2, j2=jm-1 consistently)
!                 removed a6,ar,al from argument list
!  WS  99.04.27:  DM made local to this routine
!  WS  99.09.09:  Documentation; indentation; cleaning
!  WS  99.10.22:  Added NG as argument; pruned arrays
!  SJL 99.12.24:  Revised documentation; optimized for better cache usage
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
 integer i, j, jt
 integer js2g0, jn1g1

! work arrays (should pass in eventually for performance enhancement):
 real (r8) dm(im,jfirst-ng:jlast+ng)

!     real (r8) ar(im,jfirst-1:jlast+1)  ! AR needs to be ghosted on NS
!     real (r8) al(im,jfirst-1:jlast+2)  ! AL needs to be ghosted on N2S
!     real (r8) a6(im,jfirst-1:jlast+1)  ! A6 needs to be ghosted on NS

   js2g0  = max(2,jfirst)       ! No ghosting
   jn1g1  = min(jm,jlast+1)     ! Ghost N*1
     
   if(jord == 1) then
!!!        do j=2,jm
        do j=js2g0,jn1g1
          do i=1,im
            jt = float(j) - c(i,j)
            fy(i,j) = q(i,jt)
          enddo
        enddo
   else

!
! YMIST requires q on NS;  Only call to YMIST here
!
        call ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)

        if( abs(jord) .ge. 3 ) then
 
          call fyppm(c,q,dm,fy,im,jm,ng,jord,iv,jfirst,jlast)

        else
!
! JORD can either have the value 2 or -2 at this point
!
!!!          do j=2,jm
          do j=js2g0,jn1g1
            do i=1,im
              jt = float(j) - c(i,j)
              fy(i,j) = q(i,jt) + (sign(1.,c(i,j))-c(i,j))*dm(i,jt)
            enddo
          enddo
        endif
   endif

!!!      do j=2,jm
      do j=js2g0,jn1g1
        do i=1,im
          fy(i,j) = fy(i,j)*yfx(i,j)
        enddo
      enddo
    return
!EOC
 end subroutine ytp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ymist
!
! !INTERFACE: 
 subroutine ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  NS dependencies
 integer jord                        !  order of subgrid distribution
 integer iv                          !  Scalar (==0) Vector (==1)
 real (r8) q(im,jfirst-ng:jlast+ng)  !  transported scalar  N*ng S*ng

! !OUTPUT PARAMETERS:
 real (r8) dm(im,jfirst-ng:jlast+ng)      !  Slope only N*(ng-1) S*(ng-1) used

! !DESCRIPTION:
!     Calculate the slope of the pressure.  The number of ghost
!     latitudes (NG) depends on what method (JORD) will be used
!     subsequentally.    NG is equal to MIN(ABS(JORD),3).
!
! !CALLED FROM:
!     ytp
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.13:  Added jfirst:jlast concept
!  WS  99.09.09:  Documentation; indentation; cleaning
!  SJL 00.01.06:  Documentation
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local variables

 integer i, j, jm1, im2, js2gng1, jn2gng1
 real (r8) qmax, qmin, tmp

    js2gng1 = max(2,   jfirst-ng+1)     !  Number needed on S
    jn2gng1 = min(jm-1,jlast+ng-1)      !  Number needed on N

    jm1 = jm - 1
    im2 = im / 2

!!!      do j=2,jm1
      do j=js2gng1,jn2gng1
        do i=1,im
           dm(i,j) = 0.25*(q(i,j+1) - q(i,j-1))
        enddo
      enddo

   if( iv == 0 ) then

        if ( jfirst-ng <= 1 ) then
! S pole
          do i=1,im2
            tmp = 0.25*(q(i,2)-q(i+im2,2))
            qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) =  - dm(i-im2, 1)
          enddo
        endif

        if ( jlast+ng >= jm ) then
! N pole
          do i=1,im2
            tmp = 0.25*(q(i+im2,jm1)-q(i,jm1))
            qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) =  - dm(i-im2,jm)
          enddo
        endif

   else

        if ( jfirst-ng <= 1 ) then
! South
          do i=1,im2
            tmp  = 0.25*(q(i,2)+q(i+im2,2))
            qmax = max(q(i,2),q(i,1), -q(i+im2,2)) - q(i,1)
            qmin = q(i,1) - min(q(i,2),q(i,1),-q(i+im2,2))
            dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i, 1) = dm(i-im2, 1)
          enddo
        endif

        if ( jlast+ng >= jm ) then
! North
          do i=1,im2
            tmp  = -0.25*(q(i+im2,jm1)+q(i,jm1))
            qmax = max(-q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
            qmin = q(i,jm) - min(-q(i+im2,jm1),q(i,jm), q(i,jm1))
            dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
          enddo

          do i=im2+1,im
            dm(i,jm) = dm(i-im2,jm)
          enddo
        endif

   endif

   if( jord > 0 ) then
!
! Applies monotonic slope constraint (off if jord less than zero)
!
!!!        do j=2,jm1
        do j=js2gng1,jn2gng1
          do i=1,im
            qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
            qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
            dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
          enddo
        enddo
   endif
    return
!EOC
 end subroutine ymist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: fyppm
!
! !INTERFACE: 
 subroutine fyppm(c,  q,  dm, flux, im, jm, ng, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer ng                          !  Max. NS dependencies
 integer jord                        !  Approximation order
 integer iv                          !  Scalar=0, Vector=1
 real (r8)  q(im,jfirst-ng:jlast+ng) !  mean value needed only N*2 S*2
 real (r8) dm(im,jfirst-ng:jlast+ng) !  Slope     needed only N*2 S*2
 real (r8)  c(im,jfirst:jlast+1)     !  Courant   N (like FLUX)

! !INPUT/OUTPUT PARAMETERS:
 real (r8) ar(im,jfirst-1:jlast+1)   ! AR needs to be ghosted on NS
 real (r8) al(im,jfirst-1:jlast+2)   ! AL needs to be ghosted on N2S
 real (r8) a6(im,jfirst-1:jlast+1)   ! A6 needs to be ghosted on NS

! !OUTPUT PARAMETERS:
 real (r8) flux(im,jfirst:jlast+1)   !  Flux      N (see tp2c)

! !DESCRIPTION:
!
!   NG is passed from YTP for convenience -- it is actually 1 more in NS
!   than the actual number of latitudes needed here.  But in the shared-memory 
!   case it becomes 0, which is much cleaner.
!
! !CALLED FROM:
!      ytp
!
! !REVISION HISTORY:
!
!  SJL 99.04.13:  Delivery
!  WS  99.04.19:  Added jfirst:jlast concept; FYPPM only called from YTP
!  WS  99.04.21:  Removed j1, j2  (j1=2, j2=jm-1 consistently)
!                 removed a6,ar,al from argument list
!  WS  99.09.09:  Documentation; indentation; cleaning
!  WS  99.10.22:  Added ng as argument; Pruned arrays
!  WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC
 real (r8)   r3, r23
 parameter ( r3 = 1./3., r23 = 2./3. )
 integer i, j, imh, jm1, lmt
 integer js1g1, js2g0, js2g1, jn1g2, jn1g1, jn2g1
!     logical steep

!     if(jord .eq. 6) then
!        steep = .true.
!     else
!        steep = .false.
!     endif

      imh = im / 2
      jm1 = jm - 1

      js1g1  = max(1,jfirst-1)         ! Ghost S*1
      js2g0  = max(2,jfirst)           ! No ghosting
      js2g1  = max(2,jfirst-1)         ! Ghost S*1
      jn1g1  = min(jm,jlast+1)         ! Ghost N*1
      jn1g2  = min(jm,jlast+2)         ! Ghost N*2
      jn2g1  = min(jm-1,jlast+1)       ! Ghost N*1

!!!      do j=2,jm
      do j=js2g1,jn1g2                 ! AL needed N2S
        do i=1,im                      ! P, dm ghosted N2S2 (at least)
          al(i,j) = 0.5*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
      enddo

! Yeh's steepening procedure; to be implemented
!     if(steep) call steepy(im,   jm,   jfirst,   jlast,       &
!                           ng,    q,       al,   dm )

!!!      do j=1,jm1
      do j=js1g1,jn2g1                 ! AR needed NS
        do i=1,im
          ar(i,j) = al(i,j+1)          ! AL ghosted N2S
        enddo
      enddo

! WS 990726 :  Added condition to decide if poles are on this processor

! Poles:

   if( iv == 0 ) then

        if ( jfirst == 1 ) then
          do i=1,imh
            al(i,    1) = al(i+imh,2)
            al(i+imh,1) = al(i,    2)
          enddo
        endif

        if ( jlast == jm ) then
          do i=1,imh
            ar(i,    jm) = ar(i+imh,jm1)
            ar(i+imh,jm) = ar(i,    jm1)
          enddo
        endif

   else

        if ( jfirst == 1 ) then
          do i=1,imh
            al(i,    1) = -al(i+imh,2)
            al(i+imh,1) = -al(i,    2)
          enddo
        endif

        if ( jlast == jm ) then
          do i=1,imh
            ar(i,    jm) = -ar(i+imh,jm1)
            ar(i+imh,jm) = -ar(i,    jm1)
          enddo
        endif

   endif

   if( jord == 3 .or. jord == 5 ) then
!!!      do j=1,jm
      do j=js1g1,jn1g1               ! A6 needed NS
        do i=1,im
          a6(i,j) = 3.*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
        enddo
      enddo
   endif

      lmt = jord - 3

!!!        do j=1,jm
!       do j=js1g1,jn1g1             !  A6, AR, AL needed NS
!         call lmppm(dm(1,j),a6(1,j),ar(1,j),al(1,j),q(1,j),im,lmt)
!       enddo

        call lmppm(dm(1,js1g1), a6(1,js1g1), ar(1,js1g1),               &
                   al(1,js1g1),  q(1,js1g1), im*(jn1g1-js1g1+1), lmt)

!!!      do j=2,jm
      do j=js2g0,jn1g1                 ! flux needed N
        do i=1,im
          if(c(i,j).gt.0.) then
            flux(i,j) = ar(i,j-1) + 0.5*c(i,j)*(al(i,j-1) - ar(i,j-1) +  &
                        a6(i,j-1)*(1.-r23*c(i,j)) )
          else
            flux(i,j) = al(i,j) - 0.5*c(i,j)*(ar(i,j) - al(i,j) +        &
                        a6(i,j)*(1.+r23*c(i,j)))
          endif
        enddo
      enddo
    return
!EOC
 end subroutine fyppm 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: tpcc
!
! !INTERFACE: 
 subroutine tpcc(va,   ymass,  q,   crx,  cry,  im,   jm,  ng_c, ng_d,   &
                 iord, jord,   fx,  fy,   ffsl, cose, jfirst, jlast,     &
                 dm,   qtmp,   al,  ar,   a6 )       
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
  integer im, jm                    ! Dimensions
  integer ng_c                      ! 
  integer ng_d                      ! 
  integer jfirst, jlast             ! Latitude strip
  integer iord, jord                ! Interpolation order in x,y
  logical ffsl(jm)                  ! Flux-form semi-Lagrangian transport?
  real (r8) cose(jm)                ! Critical cosine  (replicated)
  real (r8) va(im,jfirst:jlast)     ! Courant (unghosted like FX)
  real (r8) q(im,jfirst-ng_d:jlast+ng_d) !
  real (r8) crx(im,jfirst-ng_c:jlast+ng_c)
  real (r8) cry(im,jfirst:jlast)    ! Courant # (ghosted like FY)
  real (r8) ymass(im,jfirst:jlast)  ! Background y-mass-flux (ghosted like FY)

! Input 1D work arrays:
  real (r8)   dm(-im/3:im+im/3)
  real (r8) qtmp(-im/3:im+im/3)
  real (r8)   al(-im/3:im+im/3)
  real (r8)   ar(-im/3:im+im/3)
  real (r8)   a6(-im/3:im+im/3)

! !OUTPUT PARAMETERS:
  real (r8) fx(im,jfirst:jlast)     ! Flux in x (unghosted)
  real (r8) fy(im,jfirst:jlast)     ! Flux in y (unghosted since iv==0)

! !DESCRIPTION:
!     In this routine the number 
!     of north ghosted latitude min(abs(jord),2), and south ghosted
!     latitudes is XXXX
!
! !CALLED FROM:
!     cd_core
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Delivery
!   WS  99.04.13:  Added jfirst:jlast concept
!   WS  99.05.10:  Replaced JNP with JM, JMR with JM-1, IMR with IM
!   WS  99.05.10:  Removed fvcore.h and JNP, IMH, IML definitions
!   WS  99.10.20:  Pruned arrays
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

  real (r8) adx(im,jfirst-1:jlast+2)
  integer north, south
  integer i, j, jp, im2, js2g0, js2gs, jn2g0, jn1g0, jn1gn
  real (r8) wk1(im)
  real (r8) fx1(im)

    im2 = im/2
    north = min(2,abs(jord))         ! north == 1 or 2
    south = north-1                  ! south == 0 or 1
    js2g0 = max(2,jfirst)
    js2gs = max(2,jfirst-south)
    jn2g0 = min(jm-1,jlast)
    jn1gn = min(jm,jlast+north)
    jn1g0 = min(jm,jlast)

! This loop must be ghosted N*NG, S*NG

!!!      do j=2,jm
      do j=js2gs,jn1gn
        call xtp( im, ffsl(j), wk1, q(1,j),           &
                  crx(1,j), 1, crx(1,j), cose(j), 0,  &
                  dm, qtmp, al, ar, a6 )

        do i=1,im-1
          adx(i,j) = q(i,j) + 0.5 *                   &
                     (wk1(i)-wk1(i+1) + q(i,j)*(crx(i+1,j)-crx(i,j)))
        enddo

        adx(im,j) = q(im,j) + 0.5 *                   &
                    (wk1(im)-wk1(1) + q(im,j)*(crx(1,j)-crx(im,j)))
      enddo

      call ycc(im, jm, fy, adx, cry, ymass, jord, 0,jfirst,jlast)

! For Scalar only!!!
      if ( jfirst-ng_d <= 1 ) then
        do i=1,im2
          q(i,1) = q(i+im2,  2)
        enddo
        do i=im2+1,im
           q(i,1) = q(i-im2,  2)
        enddo
      endif

      if ( jlast == jm ) then
        do i=1,im2
          fx1(i) = q(i+im2,jm)
        enddo
        do i=im2+1,im
           fx1(i) = q(i-im2,jm)
        enddo

        do i=1,im
          if(va(i,jm) .gt. 0.) then
            adx(i,jm) = q(i,jm) + 0.5*va(i,jm)*(q(i,jm-1)-q(i,jm))
          else
            adx(i,jm) = q(i,jm) + 0.5*va(i,jm)*(q(i,jm)-fx1(i))
          endif
        enddo
      endif

!!!      do j=2,jm-1
      do j=js2g0,jn2g0
        do i=1,im
          jp = j-va(i,j)
! jp = j     if va < 0
! jp = j -1  if va < 0
! q needed max(1, jfirst-1)
          adx(i,j) = q(i,j) + 0.5*va(i,j)*(q(i,jp)-q(i,jp+1))
        enddo
      enddo

!!!      do j=2,jm
      do j=js2g0,jn1g0
        call xtp( im, ffsl(j), fx(1,j), adx(1,j),                 &
                  crx(1,j), iord, crx(1,j), cose(j), 0,           &
                  dm, qtmp, al, ar, a6 )
      enddo
    return
!EOC
 end subroutine tpcc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ycc
!
! !INTERFACE: 
 subroutine ycc(im, jm, fy, q, vc, ymass, jord, iv, jfirst, jlast)
!-----------------------------------------------------------------------

! !USES:
 use shr_kind_mod, only : r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer im, jm                      !  Dimensions
 integer jfirst, jlast               !  Latitude strip
 integer jord                        !  Approximation order
 integer iv                          !  Scalar=0, Vector=1
 real (r8) q(im,jfirst-1-iv:jlast+2)      !  Field (N*2 S*(iv+1))
 real (r8) vc(im,jfirst-iv:jlast)         !  Courant  (like FY)
 real (r8) ymass(im,jfirst-iv:jlast)      !  background mass flux

! !OUTPUT PARAMETERS:
 real (r8) fy(im,jfirst-iv:jlast)         !  Flux (S if iv=1)

! !DESCRIPTION:
!     Will Sawyer's note: In this routine the number 
!     of ghosted latitudes NG is min(abs(jord),2).  The scalar/vector
!     flag determines whether the flux FY needs to be ghosted on the
!     south.  If called from CD\_CORE (iv==1) then it does, if called
!     from TPCC (iv==0) it does not.  
!
! !CALLED FROM:
!     cd_core
!     tpcc
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.04.19:  Added jfirst:jlast concept
!   WS  99.04.27:  DC removed as argument (local to this routine); DC on N
!   WS  99.05.10:  Replaced JNP with JM, JMR with JM-1, IMR with IM
!   WS  99.05.10:  Removed fvcore.h
!   WS  99.07.27:  Built in tests for SP or NP
!   WS  99.09.09:  Documentation; indentation; cleaning; pole treatment
!   WS  99.09.14:  Loop limits
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!
!EOP
!---------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
  real (r8) dc(im,jfirst-iv:jlast+1)
  real (r8) qmax, qmin
  integer i, j, jt, im2, js2giv, js3giv, jn2g1, jn2g0

   im2 = im/2

   js2giv = max(2,jfirst-iv)
   js3giv = max(3,jfirst-iv)
   jn2g1  = min(jm-1,jlast+1)
   jn2g0  = min(jm-1,jlast)
      
   if(jord == 1) then
!!!        do j=2,jm-1
        do j=js2giv,jn2g0                      ! FY needed on S*iv
          do i=1,im
! jt=j if vc > 0; jt=j+1 if vc <=0
            jt = float(j+1)  - vc(i,j)         ! VC ghosted like fy
            fy(i,j) = q(i,jt)*ymass(i,j)       ! ymass ghosted like fy
          enddo                                ! q ghosted N*1, S*iv
        enddo

   else

!!!        do j=3,jm-1
        do j=js3giv,jn2g1                      ! dc needed N*1, S*iv
          do i=1,im
            dc(i,j) = 0.25*(q(i,j+1)-q(i,j-1)) ! q ghosted N*2, S*(iv+1)
          enddo
        enddo

        if(iv.eq.0) then
! Scalar.

! WS 99.07.27 : Split loops in SP and NP regions, added SP/NP tests

          if ( jfirst-iv <= 2 ) then
            do i=1,im2
              dc(i, 2) = 0.25 * ( q(i,3) - q(i+im2,2) )
            enddo

            do i=im2+1,im
              dc(i, 2) = 0.25 * ( q(i,3) - q(i-im2,2) )
            enddo
          endif

          if ( jlast == jm ) then
            do i=1,im2
              dc(i,jm) = 0.25 * ( q(i+im2,jm) - q(i,jm-1) )
            enddo

            do i=im2+1,im
              dc(i,jm) = 0.25 * ( q(i-im2,jm) - q(i,jm-1) )
            enddo
          endif

        else
! Vector winds

! WS 99.07.27 : Split loops in SP and NP regions, added SP/NP tests

          if ( jfirst-iv <= 2 ) then
            do i=1,im2
              dc(i, 2) =  0.25 * ( q(i,3) + q(i+im2,2) )
            enddo

            do i=im2+1,im
              dc(i, 2) =  0.25 * ( q(i,3) + q(i-im2,2) )
            enddo
          endif

          if ( jlast == jm ) then
            do i=1,im2
              dc(i,jm) = -0.25 * ( q(i,jm-1) + q(i+im2,jm) )
            enddo

            do i=im2+1,im
              dc(i,jm) = -0.25 * ( q(i,jm-1) + q(i-im2,jm) )
            enddo
          endif

        endif

        if( jord > 0 ) then
! Monotonic constraint
!!!          do j=3,jm-1
          do j=js3giv,jn2g1            ! DC needed N*1, S*iv
            do i=1,im                  ! P ghosted N*2, S*(iv+1)
              qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
              qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
              dc(i,j) = sign(min(abs(dc(i,j)),qmin,qmax),dc(i,j))
            enddo
          enddo
!
! WS 99.08.03 : Following loop split into SP and NP part
!
          if ( jfirst-iv <= 2 ) then
            do i=1,im
              dc(i, 2) = 0.
            enddo
          endif
          if ( jlast == jm ) then
            do i=1,im
              dc(i,jm) = 0.
            enddo
          endif
        endif

!!!        do j=2,jm-1
       do j=js2giv,jn2g0                   ! fy needed S*iv
         do i=1,im                       
           jt = float(j+1)  - vc(i,j)      ! vc, ymass ghosted like fy
           fy(i,j) = (q(i,jt)+(sign(1.,vc(i,j))-vc(i,j))*dc(i,jt))*ymass(i,j)
         enddo
       enddo
    endif
    return
!EOC
 end subroutine ycc
!-----------------------------------------------------------------------

end module tp_core
