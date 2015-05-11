module mapz_module

 public map1_ppm, mapn_ppm

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_ppm --- Piecewise parabolic mapping, variant 1
!
! !INTERFACE:
  subroutine map1_ppm( km,   pe1,   q1,  kn,   pe2,   q2,                &
                       ng_s, ng_n, itot, i1, i2,                         &
                       j, jfirst, jlast, iv, kord)

! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real(r8), intent(in) ::  pe1(itot,km+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) ::  pe2(itot,kn+1)  ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(r8), intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km) ! Field input

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    incorporated latest FVGCM version
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real(r8)       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)
      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, k, l, ll, k0
      real(r8)    pl, pr, qsum, delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

      return
!EOC
 end subroutine map1_ppm
!----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  mapn_ppm --- Piecewise parabolic mapping, multiple tracers
!
! !INTERFACE:
 subroutine mapn_ppm(km,   pe1,   q1, nq,                  &
                     kn,   pe2,   q2, ng_s, ng_n,          &
                     itot, i1, i2, j,                      &
                     jfirst, jlast, iv, kord)

! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      integer, intent(in) :: itot              ! Total latitudes
      integer, intent(in) :: iv                ! Mode: 0 ==  constituents  1 == ???
      integer, intent(in) :: kord              ! Method order
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: jfirst            ! Starting latitude
      integer, intent(in) :: jlast             ! Finishing latitude
      integer, intent(in) :: ng_s              ! Ghosted latitudes south
      integer, intent(in) :: ng_n              ! Ghosted latitudes north
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension
      integer, intent(in) :: nq                ! Number of tracers

      real(r8), intent(in) :: pe1(itot,km+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real(r8), intent(in) :: pe2(itot,kn+1)   ! pressure at layer edges 
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate
      real(r8), intent(in) ::  q1(itot,jfirst-ng_s:jlast+ng_n,km,nq) ! Field input
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout)::  q2(itot,jfirst-ng_s:jlast+ng_n,kn,nq) ! Field output

! !DESCRIPTION:
!
!     Perform piecewise parabolic method on a given latitude    
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY: 
!    02.04.04   Sawyer    incorporated latest FVGCM version, ProTeX
!    02.06.20   Sawyer    made Q2 inout since the args for Q1/Q2 same
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real(r8)       r3, r23
      parameter (r3 = 1./3., r23 = 2./3.)

      real(r8)   dp1(i1:i2,km)
      real(r8)  q4(4,i1:i2,km)

      integer i, k, l, ll, k0, iq
      real(r8)    pl, pr, qsum, delp, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

    do 2000 iq=1,nq
      do k=1,km
         do i=i1,i2
            q4(1,i,k) = q1(i,j,k,iq)
         enddo
      enddo

! Compute vertical subgrid distribution
      call ppm2m( q4, dp1, km, i1, i2, iv, kord )

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k,iq) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                          *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+     &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*             &
                     (r3*(1.+pl*(1.+pl))))
              do ll=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,ll+1) ) then
! Whole layer..
                 qsum = qsum + dp1(i,ll)*q4(1,i,ll)
                 else
                 delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp / dp1(i,ll)
                 qsum = qsum + delp*(q4(2,i,ll)+0.5*esl*               &
                       (q4(3,i,ll)-q4(2,i,ll)+q4(4,i,ll)*(1.-r23*esl)))
                 k0 = ll
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,j,k,iq) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue
2000  continue

      return
!EOC
 end subroutine mapn_ppm
!----------------------------------------------------------------------- 


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppm2m --- Piecewise parabolic method for fields
!
! !INTERFACE:
 subroutine ppm2m(a4, delp, km, i1, i2, iv, kord)

! !USES:
 use shr_kind_mod, only: r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
 integer, intent(in):: iv      ! iv =-1: winds
                               ! iv = 0: positive definite scalars
                               ! iv = 1: others
 integer, intent(in):: i1      ! Starting longitude
 integer, intent(in):: i2      ! Finishing longitude
 integer, intent(in):: km      ! vertical dimension
 integer, intent(in):: kord    ! Order (or more accurately method no.):
                               ! 
 real (r8), intent(in):: delp(i1:i2,km)     ! layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real (r8), intent(inout):: a4(4,i1:i2,km)  ! Interpolated values

! !DESCRIPTION:
!
!   Perform the piecewise parabolic method 
! 
! !REVISION HISTORY: 
!   ??.??.??    Lin        Creation
!   02.04.04    Sawyer     Newest release from FVGCM
!   02.04.23    Sawyer     Incorporated minor algorithmic change to 
!                          maintain CAM zero diffs (see comments inline)
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! local arrays:
      real(r8)  dc(i1:i2,km)
      real(r8)  h2(i1:i2,km)
      real(r8) delq(i1:i2,km)
      real(r8) df2(i1:i2,km)
      real(r8) d4(i1:i2,km)

! local scalars:
      real(r8) fac
      real(r8) a1, a2, c1, c2, c3, d1, d2
      real(r8) qmax, qmin, cmax, cmin
      real(r8) qm, dq, tmp

      integer i, k, km1, lmt
      real(r8) qmp, pmp
      real(r8) lac
      integer it

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
            c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
            c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k)
            qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))
             dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
            df2(i,k) = tmp
         enddo
      enddo

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do k=3,km1
      do i=i1,i2
      c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
      a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
      a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
      a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                                delp(i,k-1)*a1*dc(i,k  ) )
      enddo
      enddo

      call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
      dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
      c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
      c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1**2)
      a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
      a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
      dc(i,1) =  a4(1,i,1) - a4(2,i,1)
! No over- and undershoot condition
      cmax = max(a4(1,i,1), a4(1,i,2))
      cmin = min(a4(1,i,1), a4(1,i,2))
      a4(2,i,2) = max(cmin,a4(2,i,2))
      a4(2,i,2) = min(cmax,a4(2,i,2))
      enddo

      if( iv == 0 ) then
         do i=i1,i2
!
! WS: 02.04.23  Algorithmic difference with FVGCM.  FVGCM does this:
!
!!!            a4(2,i,1) = a4(1,i,1)
!!!            a4(3,i,1) = a4(1,i,1)
!
!     CAM does this:
!
            a4(2,i,1) = max(0.,a4(2,i,1))
            a4(2,i,2) = max(0.,a4(2,i,2))
         enddo
      elseif ( iv == -1 ) then
! Winds:
        if( km > 32 ) then
          do i=i1,i2
! More dampping: top layer as the sponge
             a4(2,i,1) = a4(1,i,1)
             a4(3,i,1) = a4(1,i,1)
          enddo
        else
          do i=i1,i2
             if( a4(1,i,1)*a4(2,i,1) <=  0. ) then
                 a4(2,i,1) = 0.
             else
                 a4(2,i,1) = sign(min(abs(a4(1,i,1)),    &
                                      abs(a4(2,i,1))),   &
                                          a4(1,i,1)  )
            endif
          enddo
        endif
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
         dc(i,km) = a4(3,i,km) -  a4(1,i,km)
! No over- and under-shoot condition
         cmax = max(a4(1,i,km), a4(1,i,km1))
         cmin = min(a4(1,i,km), a4(1,i,km1))
         a4(2,i,km) = max(cmin,a4(2,i,km))
         a4(2,i,km) = min(cmax,a4(2,i,km))
      enddo

! Enforce constraint at the surface

      if ( iv == 0 ) then
! Positive definite scalars:
           do i=i1,i2
              a4(3,i,km) = max(0., a4(3,i,km))
           enddo
      elseif ( iv == -1 ) then
! Winds:
           do i=i1,i2
              if( a4(1,i,km)*a4(3,i,km) <=  0. ) then
                  a4(3,i,km) = 0.
              else
                  a4(3,i,km) = sign( min(abs(a4(1,i,km)),   &
                                         abs(a4(3,i,km))),  &
                                             a4(1,i,km)  )
              endif
           enddo
      endif

      do k=1,km1
         do i=i1,i2
            a4(3,i,k) = a4(2,i,k+1)
         enddo
      enddo
 
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
 
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!****6***0*********0*********0*********0*********0*********0**********72
! Huynh's 2nd constraint
!****6***0*********0*********0*********0*********0*********0**********72
      do k=2, km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2
!           h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))
!    &               / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )
!    &               * delp(i,k)**2
! Method#3
            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      if( kord == 7 ) then
         fac = 1.5           ! original quasi-monotone
      else
         fac = 0.125         ! full monotone
      endif

      do k=3, km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(3,i,k) = min(max(a4(3,i,k), qmin), qmax)
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         qmin  = min(a4(1,i,k), qmp, lac)
         qmax  = max(a4(1,i,k), qmp, lac)
         a4(2,i,k) = min(max(a4(2,i,k), qmin), qmax)
! Recompute A6
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to prevent negatives when kord=7
         if (iv == 0 .and. kord == 7) then
             call kmppm(dc(i1,k), a4(1,i1,k), it, 2)
         endif
      enddo

      else
 
         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

      do k=3, km-2
      if( kord /= 4) then
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
      endif
         call kmppm(dc(i1,k), a4(1,i1,k), it, lmt)
      enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call kmppm(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      return
!EOC
 end subroutine ppm2m
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  kmppm --- Perform piecewise parabolic method in vertical
!
! !INTERFACE:
 subroutine kmppm(dm, a4, itot, lmt)

! !USES:
 use shr_kind_mod, only: r8 => shr_kind_r8
 implicit none

! !INPUT PARAMETERS:
      real(r8), intent(in)::     dm(*)     ! ??????
      integer, intent(in) ::     itot      ! Total Longitudes
      integer, intent(in) ::     lmt       ! 0: Standard PPM constraint
                                           ! 1: Improved full monotonicity constraint (Lin)
                                           ! 2: Positive definite constraint
                                           ! 3: do nothing (return immediately)

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: a4(4,*)   ! ???????
                                           ! AA <-- a4(1,i)
                                           ! AL <-- a4(2,i)
                                           ! AR <-- a4(3,i)
                                           ! A6 <-- a4(4,i)

! !DESCRIPTION:
!
!    Writes a standard set of data to the history buffer. 
!
! !REVISION HISTORY: 
!    00.04.24   Lin       Last modification
!    01.03.26   Sawyer    Added ProTeX documentation
!    02.04.04   Sawyer    Incorporated newest FVGCM version
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      real(r8)       r12
      parameter (r12 = 1./12.)

      real(r8) qmp
      integer i
      real(r8) da1, da2, a6da
      real(r8) fmin

! Developer: S.-J. Lin, NASA-GSFC
! Last modified: Apr 24, 2000

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

      return
!EOC
 end subroutine kmppm
!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  steepz --- Calculate attributes for PPM
!
! !INTERFACE:
 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)

! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: km                   ! Total levels
      integer, intent(in) :: i1                   ! Starting longitude
      integer, intent(in) :: i2                   ! Finishing longitude
      real(r8), intent(in) ::  dp(i1:i2,km)       ! grid size
      real(r8), intent(in) ::  dq(i1:i2,km)       ! backward diff of q
      real(r8), intent(in) ::  d4(i1:i2,km)       ! backward sum:  dp(k)+ dp(k-1) 
      real(r8), intent(in) :: df2(i1:i2,km)       ! first guess mismatch
      real(r8), intent(in) ::  dm(i1:i2,km)       ! monotonic mismatch

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) ::  a4(4,i1:i2,km)  ! first guess/steepened

!
! !DESCRIPTION:
!   This is complicated stuff related to the Piecewise Parabolic Method
!   and I need to read the Collela/Woodward paper before documenting
!   thoroughly.
!
! !REVISION HISTORY: 
!   ??.??.??    Lin?       Creation
!   01.03.26    Sawyer     Added ProTeX documentation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, k
      real(r8) alfa(i1:i2,km)
      real(r8)    f(i1:i2,km)
      real(r8)  rat(i1:i2,km)
      real(r8)  dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) = (rat(i,k+1) - rat(i,k))                             &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k))) 
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

      return
!EOC
 end subroutine steepz
!----------------------------------------------------------------------- 

end module mapz_module
