MODULE mo_convect_tables

  !-----------------------------------------------------------------------
  ! *mo_convect_tables* - tables for convective adjustment code
  !
  !           d.salmond  cray (uk)   12/8/91
  !
  !
  ! table lookups replaced
  !
  !           A. Rhodin mpi 12/98
  !
  ! When replacing the table lookups the following code has been used:
  !
  !   tlucua :  c2es*EXP(MERGE(c3les,c3ies,lo)*(tt-tmelt) &
  !           /      (tt-MERGE(c4les,c4ies,lo)))
  !
  !   tlucub :  MERGE(c5alvcp,c5alscp,lo) / (tt-MERGE(c4les,c4ies,lo))**2
  !
  !   tlucuc :  MERGE(alvdcp, alsdcp, lo)
  !
  !   tlucuaw:  c2es*EXP(c3les*(tt-tmelt)*(1./(tt-c4les)))
  !
  !   with:     lo = tt > tmelt
  !
  ! compile with option -DNOLOOKUP in order to replace lookup tables
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  IMPLICIT NONE

  SAVE

  !----------------
  ! Public entities
  !----------------
  PRIVATE
  !-----------------------------
  ! lookup tables -- obsolescent
  !-----------------------------
  PUBLIC :: jptlucu1          ! lookup table lower bound
  PUBLIC :: jptlucu2          ! lookup table upper bound
  PUBLIC :: tlucua            ! table -- e_s*Rd/Rv
  PUBLIC :: tlucub            ! table -- for derivative calculation: d es/ d t
  PUBLIC :: tlucuc            ! table -- l/cp
  PUBLIC :: tlucuaw           ! table
  PUBLIC :: set_lookup_tables ! initialization routine
  PUBLIC :: lookupoverflow    ! lookup table overflow flag
  PUBLIC :: lookuperror       ! error handling routine 
  !-----------------
  ! Module variables
  !-----------------

  INTEGER, PARAMETER :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER, PARAMETER :: jptlucu2 = 370000  ! lookup table upper bound

  LOGICAL :: lookupoverflow

  REAL(r8) :: tlucua(jptlucu1:jptlucu2)        ! table - e_s*Rd/Rv
  REAL(r8) :: tlucub(jptlucu1:jptlucu2)        ! table - for derivative calculation
  REAL(r8) :: tlucuc(jptlucu1:jptlucu2)        ! table - l/cp
  REAL(r8) :: tlucuaw(jptlucu1:jptlucu2)       ! table

CONTAINS

  SUBROUTINE set_lookup_tables
    !----------------------------
    ! -- Initialise lookup tables
    !      called from 'setphys'
    !----------------------------
    USE mo_constants, ONLY: alsdcp, alvdcp, c2es, c3ies, c3les, c4ies, c4les, &
                            c5alscp, c5alvcp, tmelt, rd, rv

    REAL(r8), PARAMETER :: clavm1 = -6096.9385
    REAL(r8), PARAMETER :: clavm2 =    21.2409642
    REAL(r8), PARAMETER :: clavm3 =    -2.711193
    REAL(r8), PARAMETER :: clavm4 =     1.673952
    REAL(r8), PARAMETER :: clavm5 =     2.433502

    REAL(r8), PARAMETER :: ciavm1 = -6024.5282
    REAL(r8), PARAMETER :: ciavm2 =    29.32707
    REAL(r8), PARAMETER :: ciavm3 =     1.0613868
    REAL(r8), PARAMETER :: ciavm4 =    -1.3198825
    REAL(r8), PARAMETER :: ciavm5 =    -0.49382577

    REAL(r8)    :: tt
    INTEGER :: it  
    LOGICAL :: lo
    LOGICAL :: lamip2
    
    lamip2 = .true.
    tt = jptlucu1 * 0.001
    DO it = jptlucu1, jptlucu2
      lo = tt > tmelt
      
      IF (lamip2) THEN
        tlucua(it) = EXP((MERGE(clavm1, ciavm1, lo)/tt           &
                        + MERGE(clavm2, ciavm2, lo)              &
                        + MERGE(clavm3, ciavm3, lo)*tt/100.      &
                        + MERGE(clavm4, ciavm4, lo)*tt*tt/1.e5   &
                        + MERGE(clavm5, ciavm5, lo)*LOG(tt)))*rd/rv
      ELSE
        tlucua(it) = c2es*EXP(MERGE(c3les,  c3ies,  lo)*(tt-tmelt) &
                   *  (1./(tt-MERGE(c4les,  c4ies,  lo))))
      END IF
      tlucub(it)  =          MERGE(c5alvcp,c5alscp,lo) &
                  *  (1./(tt-MERGE(c4les,  c4ies,  lo)))**2
      tlucuc(it)  =          MERGE(alvdcp, alsdcp, lo)

      tlucuaw(it) = c2es*EXP(c3les*(tt-tmelt)*(1./(tt-c4les)))

      tt = tt + 0.001
    END DO

  END SUBROUTINE set_lookup_tables

  SUBROUTINE lookuperror (name)

!   USE mo_exception,  ONLY: message, finish

    CHARACTER (*) :: name

!   CALL message (name, ' lookup table overflow')
!   lookupoverflow = .FALSE.
!    CALL finish (name, ' lookup table overflow')

  END SUBROUTINE lookuperror

END MODULE mo_convect_tables
