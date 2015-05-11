!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  dry_adj --- Dry atmosphere adjustment (???)
!
! !INTERFACE:
subroutine dry_adj (itot,     km,     rdt,     pc,              &
                    pl1,      fu,     fv,      pt,              &
                    u,        v,      dp )
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none

! !INPUT PARAMETERS:
   integer itot               ! Horizontal (longitude) points
   integer km                 ! Vertical levels
   real(r8) pc                ! Critical pressure
   real(r8) rdt               ! Delta time?
   real(r8) pl1(km)           ! Edge pressure
   real(r8) dp(itot,km)       ! Change in pressure

! !INPUT/OUTPUT PARAMETERS:
   real(r8) u(itot,km)        ! Winds in X
   real(r8) v(itot,km)        ! Winds in Y
   real(r8) pt(itot,km)       ! Potential temperatur
   real(r8) fu(itot,km)       ! ???
   real(r8) fv(itot,km)       ! ???

! !DESCRIPTION:
!   Calculate the dry atmosphere mean adjustment (???) for a given
!   latitude.
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.03.26   Sawyer  Added ProTeX documentation
!   01.06.13   Mirin   2-D decomposition
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
   integer i, k
   integer npt
   real(r8) tiny
   parameter (tiny = 0.01)
   real(r8) ut(itot,km)
   real(r8) vt(itot,km)
   integer klow
!-----------------------------------------------------------------------

   klow = km
   do k=1, km
      if( pc < pl1(k) ) then
          klow = k
          go to 123
      endif
   enddo

123   continue

   do k=1,klow
      do i=1,itot
         ut(i,k) = u(i,k)
         vt(i,k) = v(i,k)
      enddo
   enddo

! Top down

   npt = 0
   do k=1, klow-1
      do i=1,itot
         if((pt(i,k)+tiny) .lt. pt(i,k+1)) then
             pt(i,k) = (pt(i,k)*dp(i,k)+pt(i,k+1)*dp(i,k+1))       &
                     / (dp(i,k) + dp(i,k+1))
             pt(i,k+1) = pt(i,k)
             u(i,k) = (u(i,k)*dp(i,k)+u(i,k+1)*dp(i,k+1))          &
                    / (dp(i,k) + dp(i,k+1))
             u(i,k+1) = u(i,k)
             v(i,k) = (v(i,k)*dp(i,k)+v(i,k+1)*dp(i,k+1))          &
                    / (dp(i,k) + dp(i,k+1))
             v(i,k+1) = v(i,k)
             npt = npt + 1
         endif
      enddo
   enddo

! From Bottom up

    if(npt .ne. 0) then
       do k=klow, 2
          do i=1,itot
             if(pt(i,k) .gt. (pt(i,k-1)+tiny)) then
                pt(i,k) = (pt(i,k)*dp(i,k)+pt(i,k-1)*dp(i,k-1))     &
                        / (dp(i,k) + dp(i,k-1))
                pt(i,k-1) = pt(i,k)
                u(i,k) = (u(i,k)*dp(i,k)+u(i,k-1)*dp(i,k-1))        &
                        / (dp(i,k) + dp(i,k-1))
                u(i,k-1) = u(i,k)
                v(i,k) = (v(i,k)*dp(i,k)+v(i,k-1)*dp(i,k-1))        &
                        / (dp(i,k) + dp(i,k-1))
                v(i,k-1) = v(i,k)
                npt = npt + 1
             endif
          enddo
       enddo

       do k=1,klow
          do i=1,itot
             fu(i,k) = fu(i,k) + (u(i,k) - ut(i,k)) * rdt
             fv(i,k) = fv(i,k) + (v(i,k) - vt(i,k)) * rdt
          enddo
       enddo
    endif
!EOC
!-----------------------------------------------------------------------
end subroutine dry_adj
