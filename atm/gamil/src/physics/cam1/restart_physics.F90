#include <misc.h>
#include <params.h>

!!  (wanhui 2003.07.10)
!!------------------------

module restart_physics

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid
    use phys_grid,    only: read_chunk_from_field, write_field_from_chunk, get_ncols_p
    use pmgrid,       only: masterproc
    use prognostics,  only: ptimelevels, n3, n3m2
    use buffer
    use radae,        only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer
    use comsrf
    use ioFileMod
    use phys_buffer
#ifdef COUP_CSM
    use ccsm_msg,     only: initialize_ccsm_msg, write_restart_ccsm, read_restart_ccsm
#endif

    implicit none

    private
    !
    ! Public interfaces
    !
    public write_restart_physics
    public read_restart_physics
    public get_abs_restart_filepath

    !
    ! Private data
    !
    character(len=256) :: pname  ! Full abs-ems restart filepath
    !
    ! Filename specifier for restart abs-ems file
    ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
    !
    character(len=256) :: rafilename_spec = '%c.cam2.ra.%y-%m-%d-%s'   ! abs-ems restart


CONTAINS

    subroutine write_restart_physics (nrg, nrg2)
        use filenames, only: mss_irt, mss_wpass, get_archivedir, interpret_filename_spec
        ! for nlend and aeres
#include <comctl.h>
        !
        ! Input arguments
        !
        integer :: nrg
        integer :: nrg2
        integer :: kvh_idx
        !
        ! Local workspace
        !
        real(r8) tmpfield(pcols,begchunk:endchunk)
        real(r8) tmpfield3d(pcols,plevmx,begchunk:endchunk)
        real(r8) tmpfield3d2(pcols,pverp,begchunk:endchunk)
        integer i                 ! loop index
        integer n3tmp             ! timestep index
        character(len=256) fname  ! abs-ems restart filename
        integer ioerr             ! I/O status
        integer  :: ncol          ! number of vertical columns
        !
        ! Buffer module variables
        !
        call write_field_from_chunk(nrg,1,1,1,pblht)
        call write_field_from_chunk(nrg,1,1,1,tpert)
        call write_field_from_chunk(nrg,1,pver,1,qrs)
        call write_field_from_chunk(nrg,1,pver,1,qrl)
        call write_field_from_chunk(nrg,1,pcnst+pnats,1,qpert)
        !
        ! cld, qcwat, and tcwat are physics things, but have dynamics time levels
        !
        call write_field_from_chunk(nrg,1,pver,1,cld(1,1,begchunk,n3  ))
        call write_field_from_chunk(nrg,1,pver,1,cld(1,1,begchunk,n3m2))

        call write_field_from_chunk(nrg,1,pver,1,qcwat(1,1,begchunk,n3  ))
        call write_field_from_chunk(nrg,1,pver,1,qcwat(1,1,begchunk,n3m2))

        call write_field_from_chunk(nrg,1,pver,1,tcwat(1,1,begchunk,n3  ))
        call write_field_from_chunk(nrg,1,pver,1,tcwat(1,1,begchunk,n3m2))

        call write_field_from_chunk(nrg,1,pver,1,lcwat(1,1,begchunk,n3  ))
        call write_field_from_chunk(nrg,1,pver,1,lcwat(1,1,begchunk,n3m2))
        !!
        !
        ! Comsrf module variables
        !
#if (! defined COUP_CSM)
        call write_field_from_chunk(nrg,1,1,1,fsnt)
#endif
        call write_field_from_chunk(nrg,1,1,1,fsns)
#if (! defined COUP_CSM)
        call write_field_from_chunk(nrg,1,1,1,flnt)
        call write_field_from_chunk(nrg,1,1,1,flns)
#endif
        do i = begchunk, endchunk
            tmpfield(:,i) = srfflx_state2d(i)%asdir(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i = begchunk, endchunk
            tmpfield(:,i) = srfflx_state2d(i)%asdif(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i = begchunk,endchunk
            tmpfield(:,i) = srfflx_state2d(i)%aldir(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i = begchunk, endchunk
            tmpfield(:,i) = srfflx_state2d(i)%aldif(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

#if (! defined COUP_CSM)
        call write_field_from_chunk(nrg,1,1,1,asdirice)
        call write_field_from_chunk(nrg,1,1,1,asdifice)
        call write_field_from_chunk(nrg,1,1,1,aldirice)
        call write_field_from_chunk(nrg,1,1,1,aldifice)
        call write_field_from_chunk(nrg,1,1,1,tsice)
#endif

        do i=begchunk,endchunk
            tmpfield(:,i) = srfflx_state2d(i)%lwup(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        call write_field_from_chunk(nrg,1,1,1,landfrac)
        call write_field_from_chunk(nrg,1,1,1,landm)
        call write_field_from_chunk(nrg,1,1,1,sgh)
        do i=begchunk,endchunk
            tmpfield(:,i) = srfflx_state2d(i)%ts(:)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            tmpfield3d(:,:,i) = surface_state2d(i)%tssub(:,:)
        end do
        call write_field_from_chunk(nrg,1,plevmx,1,tmpfield3d)
        call write_field_from_chunk(nrg,1,1,1,sicthk)
        call write_field_from_chunk(nrg,1,1,1,snowhland)
#if (! defined COUP_CSM)
        call write_field_from_chunk(nrg,1,1,1,snowhice)
#else
        snowhice = 0.
#endif
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%flwds(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%sols(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%soll(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%solsd(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%solld(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)
        call write_field_from_chunk(nrg,1,1,1,trefmxav)
        call write_field_from_chunk(nrg,1,1,1,trefmnav)
        call write_field_from_chunk(nrg,1,1,1,icefrac)
        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%zbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%ubot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%vbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%thbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%qbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%pbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

        do i=begchunk,endchunk
            ncol = get_ncols_p(i)
            tmpfield(:ncol,i) = surface_state2d(i)%tbot(:ncol)
        end do
        call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%ts(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%asdir(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%aldir(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%asdif(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%aldif(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%wsx(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%wsy(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%lhf(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%shf(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%lwup(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
        tmpfield(:ncol,i) = srfflx_parm2d(i)%cflx(:,1)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield(:ncol,i) = srfflx_parm2d(i)%tref(:ncol)
      end do
      call write_field_from_chunk(nrg,1,1,1,tmpfield)

      kvh_idx = pbuf_get_fld_idx('KVH')
      do i=begchunk,endchunk
       ncol = get_ncols_p(i)
       tmpfield3d2(:ncol,:pverp,i) = pbuf(kvh_idx)%fld_ptr(1,1:ncol,1:pverp,i,1)
      end do
      call write_field_from_chunk(nrg,1,pverp,1,tmpfield3d2)
#ifdef COUP_CSM
        call write_restart_ccsm ()
#endif
        !
        !-----------------------------------------------------------------------
        ! Write the abs/ems restart dataset if necessary
        !-----------------------------------------------------------------------
        !
        if (aeres) then
            if (masterproc) then
                fname = interpret_filename_spec( rafilename_spec )
                pname = trim(get_archivedir('rest'))//fname
                call opnfil(fname, nrg2, 'u')
                write(nrg,iostat=ioerr) pname
                if (ioerr /= 0 ) then
                    write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
                    call endrun
                end if
            endif

            call write_field_from_chunk(nrg2, 1, pverp*pverp,1, abstot_3d(1,1,1,begchunk))
            call write_field_from_chunk(nrg2, 1, pver*4,     1, absnxt_3d(1,1,1,begchunk))
            call write_field_from_chunk(nrg2, 1, pverp,      1, emstot_3d(1,1,begchunk))

            if (masterproc) then
                close(nrg2)
                call putfil (fname, pname, mss_wpass, mss_irt, (.not. nlend) )
            end if
        end if

        return
    end subroutine write_restart_physics

    !#######################################################################

    subroutine read_restart_physics (nrg, nrg2, aeres )
        !
        ! Arguments
        !
        integer, intent(in) :: nrg
        integer, intent(in) :: nrg2
        integer :: kvh_idx

        logical, intent(in) :: aeres
        !
        ! Local workspace
        !
        real(r8) tmpfield(pcols,begchunk:endchunk)
        real(r8) tmpfield3d(pcols,plevmx,begchunk:endchunk)
        real(r8) tmpfield3d2(pcols,pverp,begchunk:endchunk)
        integer i                 ! loop index
        integer n3tmp             ! timestep index
        character*80  locfn       ! Local filename
        integer ioerr             ! I/O status
        !
        ! Buffer module variables
        !
        call initialize_buffer ()

        call read_chunk_from_field(nrg,1,1,1,pblht)
        call read_chunk_from_field(nrg,1,1,1,tpert)
        call read_chunk_from_field(nrg,1,pver,1,qrs)
        call read_chunk_from_field(nrg,1,pver,1,qrl)
        call read_chunk_from_field(nrg,1,pcnst+pnats,1,qpert)
        !
        ! cld, qcwat, and tcwat are physics things, but have dynamics time levels
        !
        !!
        call read_chunk_from_field(nrg,1,pver,1,cld(1,1,begchunk,n3  ))
        call read_chunk_from_field(nrg,1,pver,1,cld(1,1,begchunk,n3m2))

        call read_chunk_from_field(nrg,1,pver,1,qcwat(1,1,begchunk,n3  ))
        call read_chunk_from_field(nrg,1,pver,1,qcwat(1,1,begchunk,n3m2))

        call read_chunk_from_field(nrg,1,pver,1,tcwat(1,1,begchunk,n3  ))
        call read_chunk_from_field(nrg,1,pver,1,tcwat(1,1,begchunk,n3m2))

        call read_chunk_from_field(nrg,1,pver,1,lcwat(1,1,begchunk,n3  ))
        call read_chunk_from_field(nrg,1,pver,1,lcwat(1,1,begchunk,n3m2))
        !!
        !
        ! Comsrf module variables
        !
        call initialize_comsrf
#if (! defined COUP_CSM)
        call read_chunk_from_field(nrg,1,1,1,fsnt)
#endif
        call read_chunk_from_field(nrg,1,1,1,fsns)
#if (! defined COUP_CSM)
        call read_chunk_from_field(nrg,1,1,1,flnt)
        call read_chunk_from_field(nrg,1,1,1,flns)
#endif
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%asdir(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%asdif(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%aldir(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%aldif(:) = tmpfield(:,i)
        end do

#if (! defined COUP_CSM)
        call read_chunk_from_field(nrg,1,1,1,asdirice)
        call read_chunk_from_field(nrg,1,1,1,asdifice)
        call read_chunk_from_field(nrg,1,1,1,aldirice)
        call read_chunk_from_field(nrg,1,1,1,aldifice)
        call read_chunk_from_field(nrg,1,1,1,tsice)
#endif

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%lwup(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,landfrac)
        call read_chunk_from_field(nrg,1,1,1,landm)
        call read_chunk_from_field(nrg,1,1,1,sgh)
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            srfflx_state2d(i)%ts(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,plevmx,1,tmpfield3d)
        do i=begchunk,endchunk
            surface_state2d(i)%tssub(:,:) = tmpfield3d(:,:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,sicthk)
#if (! defined COUP_CSM)
        call read_chunk_from_field(nrg,1,1,1,snowhland)
#endif
        call read_chunk_from_field(nrg,1,1,1,snowhice)
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%flwds(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%sols(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%soll(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%solsd(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%solld(:) = tmpfield(:,i)
        end do
        call read_chunk_from_field(nrg,1,1,1,trefmxav)
        call read_chunk_from_field(nrg,1,1,1,trefmnav)
        call read_chunk_from_field(nrg,1,1,1,icefrac)

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%zbot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%ubot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%vbot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%thbot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%qbot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%pbot(:) = tmpfield(:,i)
        end do

        call read_chunk_from_field(nrg,1,1,1,tmpfield)
        do i=begchunk,endchunk
            surface_state2d(i)%tbot(:) = tmpfield(:,i)
        end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%ts(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%asdir(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%aldir(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%asdif(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%aldif(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%wsx(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%wsy(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%lhf(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%shf(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%lwup(:) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
       srfflx_parm2d(i)%cflx(:,1) = tmpfield(:,i)
      end do

      call read_chunk_from_field(nrg,1,1,1,tmpfield)
      do i=begchunk,endchunk
	srfflx_parm2d(i)%tref(:) = tmpfield(:,i)
      end do
      call read_chunk_from_field(nrg,1,pverp,1,fld_kvh)

#ifdef COUP_CSM
        call initialize_ccsm_msg ()
        call read_restart_ccsm ()
#endif
        !
        !-----------------------------------------------------------------------
        ! Read the abs/ems restart dataset if necessary
        !-----------------------------------------------------------------------
        !
        call initialize_radbuffer ()
        if (aeres) then
            if (masterproc) then
                read(nrg,iostat=ioerr) pname
                if (ioerr /= 0 ) then
                    write (6,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
                    call endrun
                end if
                call getfil (pname, locfn)
                call opnfil (locfn, nrg2, 'u')
            endif

            call read_chunk_from_field(nrg2, 1, pverp*pverp,1,abstot_3d(1,1,1,begchunk))
            call read_chunk_from_field(nrg2, 1, pver*4,     1,absnxt_3d(1,1,1,begchunk))
            call read_chunk_from_field(nrg2, 1, pverp,      1,emstot_3d(1,1,begchunk))

            if (masterproc) close(nrg2)
        end if

        return
    end subroutine read_restart_physics

    character(len=256) function get_abs_restart_filepath ( )
        !
        ! Return the full filepath to the abs-ems restart file
        !
        get_abs_restart_filepath = pname
    end function get_abs_restart_filepath

end module restart_physics
