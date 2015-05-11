#include <misc.h>
#include <params.h>

! (wanhui 2003.04.03)
! -------------------
	subroutine vpar (  pmtop, p0, sig, sigl, dsig )

        implicit none

#include <PARADYN>

        real*8  p0          ! model bottom pressure,     intent( in)
        real*8  pmtop       ! model top pressure,        intent(out)
        real*8  sig (nz)    ! sigma at model interfaces, intent(out)
        real*8  sigl(nl)    ! sigma at model levels,     intent(out)
	real*8  dsig(nl)    ! vertical stepsize,         intent(out)

	real*8  ps0         ! referrence surface pressre of cam2
        real*8  hyai(nz)    ! parameter 'a' in cam2 pressure formula
	real*8  hybi(nz)    ! parameter 'b' in cam2 pressure formula
	integer k           ! layer index

        data ps0  /1.0d3/   ! unit: hPa

        data hyai /0.00219406700000001d0, 0.00489520900000001d0,
     _             0.009882418d0,         0.01805201d0,
     _             0.02983724d0,          0.0446233400000002d0,
     _             0.0616058700000002d0,  0.0785124300000004d0,
     _             0.0773127100000002d0,  0.0759013100000003d0,
     _             0.0742408600000002d0,  0.0722874400000002d0,
     _             0.0699893299999998d0,  0.06728574d0,
     _             0.06410509d0,          0.0603632200000002d0,
     _             0.0559611100000001d0,  0.0507822500000001d0,
     _             0.0446896000000001d0,  0.0375219099999999d0,
     _             0.0290894900000001d0,  0.02084739d0,
     _             0.01334443d0,          0.00708499000000001d0,
     _             0.00252136d0,          0.0d0,
     _             0.0d0/

        data hybi /0.0d0,                 0.0d0,
     _             0.0d0,                 0.0d0,
     _             0.0d0,                 0.0d0,
     _             0.0d0,                 0.0d0,
     _             0.01505309d0,          0.03276228d0,
     _             0.05359622d0,          0.0781062700000006d0,
     _             0.1069411d0,           0.140863700000001d0,
     _             0.180772d0,            0.227722d0,
     _             0.282956200000001d0,   0.347936400000002d0,
     _             0.4243822d0,           0.514316800000003d0,
     _             0.620120200000002d0,   0.723535500000004d0,
     _             0.817676800000001d0,   0.896215300000001d0,
     _             0.953476100000003d0,   0.9851122d0,
     _             1.0d0/

!--------------------------------------------------------------

	pmtop = hyai(1)*ps0

        sig(1)= 0.0d0
        do k=1,nl
           sig (k+1)= (hyai(k+1)*ps0+hybi(k+1)*p0-pmtop)/(p0-pmtop)
	   sigl(k  )= 0.5d0*( sig(k)+sig(k+1) )
           dsig(k  )= sig(k+1)-sig(k)
        enddo
!	write(*,*) 'pmtop:',pmtop
!	write(*,'(1x,i3,3f20.16)')  (k,sig(k),sigl(k),dsig(k),k=1,nl)
!	write(*,'(1x,i3,f20.16)')   nz,sig(k)

	return
        end

