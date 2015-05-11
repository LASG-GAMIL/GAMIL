#include <misc.h>
#include <params.h>

subroutine tfilt_massfix (ztodt   ,lat     ,u3m1    ,v3m1    ,t3m1    , &
                          q3      ,q3m1    ,ps      ,cwava   ,alpha   , &
                          etamid  ,qfcst   ,div     ,phis    ,omga    , &
                          dpsl    ,dpsm    ,nlon    ,t3      ,beta )
!-----------------------------------------------------------------------
!
! Purpose:
! Atmosphere and constituent mass fixer
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use commap
  use history, only: outfld
  use constituents, only: pcnst, pnats, qmin
  use tracers, only: tottnam

  implicit none

#include <comctl.h>
#include <comlun.h>
#include <comqfl.h>
#include <comtfc.h>
!
!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                          ! timestep
  integer , intent(in)   :: lat                            ! latitude index
  real(r8), intent(in)   :: u3m1    (plond,plev)             ! u
  real(r8), intent(in)   :: v3m1    (plond,plev)             ! v
  real(r8), intent(inout):: t3m1    (plond,plev)             ! T
  real(r8), intent(in)   :: q3      (plond,plev,pcnst+pnats) ! q + consts (time level n-1)
  real(r8), intent(inout):: q3m1    (plond,plev,pcnst+pnats) ! q + consts (time level n  )
  real(r8), intent(inout):: ps    (plond)                  ! Ps
  real(r8), intent(in)   :: cwava                          ! weight for global integrals
  real(r8), intent(in)   :: alpha (pcnst)                  ! slt fixer coefficient
  real(r8), intent(in)   :: etamid(plev)                   ! vertical coords at midpoints 
  real(r8), intent(in)   :: qfcst (plond,plev,pcnst)       ! slt moisture forecast
  real(r8), intent(in)   :: div   (plond,plev)             ! divergence
  real(r8), intent(in)   :: phis  (plond)                  ! Geopotential field
  real(r8), intent(out)  :: omga  (plond,plev)             ! vertical motion
  real(r8), intent(in)   :: dpsl  (plond)                  ! long comp of grad ln(ps)
  real(r8), intent(in)   :: dpsm  (plond)                  ! lat comp  of grad ln(ps)
  integer , intent(in)   :: nlon                           ! number of longitudes
  real(r8), intent(in)   :: t3(plond,plev) ! temperature
  real(r8), intent(in)   :: beta                           ! energy fixer coefficient
!
!---------------------------Local workspace-----------------------------
!
  integer ifcnt                     ! Counter
  real(r8) qfcst1(plond,plev,pcnst) ! slt moisture forecast temporary
  real(r8) rpmid (plond,plev)       ! 1./pmid
  real(r8) pdel  (plond,plev)       ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) pint  (plond,plevp)      ! pressure at model interfaces (n  )
  real(r8) pmid  (plond,plev)       ! pressure at model levels (time n)
! real(r8) utend (plond,plev)       ! du/dt
! real(r8) vtend (plond,plev)       ! dv/dt
! real(r8) ttend (plond,plev)       ! dT/dt
! real(r8) qtend (plond,plev,pcnst) ! dq/dt
! real(r8) pstend(plond)            ! d(ps)/dt
  real(r8) psl   (plond)            ! sea level pressure
  real(r8) corm                     ! fixer limit
  real(r8) wm                       ! accumulator 
  real(r8) absf                     ! absolute value of fixer
  real(r8) worst                    ! largest fixer contribution at each model level
  logical lfixlim                   ! flag to turn on fixer limiter

  real(r8) ta    (plond,plev,pcnst) ! total advection of constituents
  real(r8) dqfx3 (plond,plev,pcnst) ! q tendency due to mass adjustment
  real(r8) coslat                   ! cosine(latitude)
  real(r8) rcoslat(plond)           ! 1./cosine(latitude)
! real(r8) engt                     ! Thermal   energy integral
! real(r8) engk                     ! Kinetic   energy integral
! real(r8) engp                     ! Potential energy integral

  integer i, k, m                   ! indices
!
!-----------------------------------------------------------------------
  qfcst1(1:nlon,:,:) = qfcst(i1:nlon+i1-1,:,:)
!
  coslat  = cos(clat(lat))
  do i=1,nlon
     rcoslat(i) = 1./coslat
  enddo
  lfixlim = .true.
  corm    = 0.1
!
! Fix temperature to counter-act energy imbalance created by dynamics only
!
  do k=1,plev
    do i=1,nlon
      t3m1(i,k) = t3m1(i,k) + beta*abs(t3m1(i,k) - t3(i,k))
    end do
  end do
!
! Set average dry mass to specified constant preserving horizontal
! gradients of ln(ps). Proportionality factor was calculated in STEPON
! for nstep=0 or SCAN2 otherwise from integrals calculated in INIDAT
! and SCAN2 respectively.
! Set p*.
!
  do i=1,nlon
     ps(i) = ps(i)*fixmas
  end do
!
! Set current time pressure arrays for model levels etc.
!
  call plevs0(nlon    ,plond   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
  rpmid(:nlon,:plev) = 1./pmid(:nlon,:plev)
!
! Compute q tendency due to mass adjustment
! If LFIXLIM = .T., then:
! Check to see if fixer is exceeding a desired fractional limit of the
! constituent mixing ratio ("corm").  If so, then limit the fixer to
! that specified limit.
!
  do m=1,pcnst
     do k=1,plev
        do i=1,nlon
           dqfx3(i,k,m) = alpha(m)*etamid(k)*abs(qfcst1(i,k,m) - q3(i,k,m))
        end do

        if (lfixlim) then
           ifcnt = 0
           worst = 0.
           wm    = 0.
           do i = 1,nlon
              absf = abs(dqfx3(i,k,m))
              if (absf.gt.corm) then
                 ifcnt = ifcnt + 1
                 worst = max(absf,worst)
                 wm = wm + absf
                 dqfx3(i,k,m) = sign(corm,dqfx3(i,k,m))
              endif
           end do
           if (ifcnt.gt.0) then
              wm = wm/float(ifcnt)
              write (6,1000) m,corm,ifcnt,k,lat,wm,worst
           endif
        endif

        do i=1,nlon
           dqfx3(i,k,m) = qfcst1(i,k,m)*dqfx3(i,k,m)/ztodt
#ifdef HADVTEST
           q3m1(i,k,m) = qfcst1(i,k,m)
#else
           q3m1(i,k,m) = qfcst1(i,k,m) + ztodt*dqfx3(i,k,m)
#endif
           ta(i,k,m) = (q3m1(i,k,m) - q3(i,k,m))/ztodt
        end do
     end do
  end do

!
! Check for and correct invalid constituents
!
  call qneg3('TFILT_MASSFIX',lat   ,nlon    ,plond   ,plev    , &
             pcnst+pnats,qmin ,q3m1(1,1,1))
!
! Send slt tendencies to the history tape
!
  do m=1,pcnst
     call outfld(tottnam(m),ta(1,1,m)   ,plond   ,lat     )
  end do
!
! Calculate vertical motion field
!
  call omcalc(rcoslat ,div     ,u3m1    ,v3m1    ,dpsl    , &
              dpsm    ,pmid    ,pdel    ,rpmid   ,pint(1,plevp), &
              omga    ,nlon    )

  call plevs0(nlon    ,plond   ,plev    ,ps      ,pint    ,pmid    ,pdel)
!
! Compute time tendencies:comment out since currently not on h-t
!
!      do k=1,plev
!        do i=1,nlon
!          ttend(i,k) = (t3m1(i,k)-tm2(i,k))/ztodt
!          utend(i,k) = (u3m1(i,k)-um2(i,k))/ztodt
!          vtend(i,k) = (v3m1(i,k)-vm2(i,k))/ztodt
!        end do  
!      end do  
!      do m=1,pcnst
!        do k=1,plev
!          do i=1,nlon
!            qtend(i,k,m) = (q3m1(i,k,m) - qm2(i,k,m))/ztodt
!          end do
!        end do
!      end do
!      do i=1,nlon
!        pstend(i) = (ps(i) - psm2(i))/ztodt
!      end do

!
! do m=1,pcnst
!    call outfld(tendnam(m),qtend(1,1,m),plond   ,lat     )
! end do
! call outfld('UTEND   ',utend   ,plond   ,lat     )
! call outfld('VTEND   ',vtend   ,plond   ,lat     )
! call outfld('TTEND   ',ttend   ,plond   ,lat     )
! call outfld('LPSTEN  ',pstend  ,plond   ,lat     )

  call plevs0(nlon    ,plond   ,plev    ,ps      ,pint    ,pmid    ,pdel)

  return
1000 format(' TIMEFILTER: WARNING: fixer for tracer ',i3,' exceeded ', &
            f8.5,' for ',i5,' points at k,lat = ',2i4,' Avg/Worst = ',1p2e10.2)
end subroutine tfilt_massfix
