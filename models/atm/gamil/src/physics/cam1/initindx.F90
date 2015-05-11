
subroutine initindx
!-----------------------------------------------------------------------
!
! Purpose: Registers constituents and determines index values.
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, ppcnst, cnst_add, advected, nonadvec, cnst_chk_dim, cnst_name
  use chemistry,    only: chem_register_cnst
  use physconst,    only: mwdry, cpair, mwh2o, cph2o
  use tracers
  use phys_buffer,  only: pbuf_init  !sxj-2008-11-10
  use MG,           only: stratiform_register  !sxj-2008-11-10
  use vertical_diffusion, only: vd_register  !sxj
  USE pmgrid,  ONLY: masterproc,iam !--sxj--debug-test

  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
#include <RK_or_MG.h> !sxj
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index

  logical, parameter :: cldw_adv=.false.  ! true => cloud water is treated as advected tracer

  character*3 trnum   ! Advected species number (Character)

!
!-----------------------------------------------------------------------
! Register constituents and tracers and set starting indexes
! Note that public starting indexes should not be required, these should be
! private to the package (module) implementing the tracer and may be obtained
! by calling cnst_get_ind
!-----------------------------------------------------------------------
!
  ! Initialize physics buffer
  call pbuf_init() !sxj-2008-11-10

! Register water vapor
  call cnst_add('Q', advected, mwh2o, cph2o, 1.0E-12_r8, ixmoist, &
                longname='Specific humidity')
!
! Register advected tracers and determine starting index
  do m = 1, nusr_adv
     write(unit=trnum,fmt='(i3)') m+100
     call cnst_add('ADV'//trnum(2:3), advected, mwdry, cpair, 0._r8, mm, &
                  longname='Advected tracer no. '//trnum(2:3))
     if (m == 1) ixuadv = mm
  end do
!
! Register non-advected tracers and determine starting index (shouldn't need
  do m = 1, nusr_nad
     write(unit=trnum,fmt='(i3)') m+100
     call cnst_add('NAD'//trnum(2:3), nonadvec, mwdry, cpair, 0._r8, mm, &
                  longname='Non-advected tracer no. '//trnum(2:3))
     if (m == 1) ixunad = mm
  end do
!
!
! Register chemical constituents
  if (trace_gas) then
     call chem_register_cnst
  endif


  if (RK_or_MG=='MG') then
! Register MG stratiform constituents ! sxj-2008-11-10
! cloud water ice  /two-moment
     call stratiform_register
	 call vd_register
	 if (masterproc) write(6,*) "stratiform_register***sxj********"
  elseif (RK_or_MG=='RK') then
! Register cloud water and determine index (either advected or non-adv).
     if (cldw_adv) then  !(cldw_adv false)  pcnst pnats
         call cnst_add('CWAT', advected, mwdry, cpair, 0._r8, ixcldw, &
                    longname='Total Grid box averaged Condensate Amount (liquid + ice)')
     else
         call cnst_add('CWAT', nonadvec, mwdry, cpair, 0._r8, ixcldw, &
                    longname='Total Grid box averaged Condensate Amount (liquid + ice)')
     endif
! add 3 nonadvec constituent   make pnad(last non-advect tracer number equal pcnst<1>+pnats<4>)
     call cnst_add('TEST1', nonadvec, mwdry, cpair, 0._r8, mm)
	 call cnst_add('TEST2', nonadvec, mwdry, cpair, 0._r8, mm)
	 call cnst_add('TEST3', nonadvec, mwdry, cpair, 0._r8, mm)
  endif
!
! Register advected test tracers and determine starting index
  if (trace_test1 .or. trace_test2 .or. trace_test3) &
       call cnst_add('TEST1', advected, mwdry, cpair, 0._r8, ixtrct)
  if (trace_test2 .or. trace_test3) &
       call cnst_add('TEST2', advected, mwdry, cpair, 0._r8, mm)
  if (trace_test3) &
       call cnst_add('TEST3', advected, mwdry, cpair, 0._r8, mm)
!
! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim
!
! Set default names for non-water advected and non-advected tracers
! Set names of advected and non-advected tracer diagnostics
!
  do m=1,ppcnst
     dcconnam(m) = 'DC'//cnst_name(m)
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do
  do m=1,pcnst
     hadvnam(m)  = 'HA'//cnst_name(m)
     vadvnam(m)  = 'VA'//cnst_name(m)
     fixcnam(m)  = 'DF'//cnst_name(m)
     tendnam(m)  = 'TE'//cnst_name(m)
     tottnam(m)  = 'TA'//cnst_name(m)
  end do

  return
end subroutine initindx
