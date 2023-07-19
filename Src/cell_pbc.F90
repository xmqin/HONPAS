      MODULE  cell_pbc
      USE precision,   ONLY:  dp
      
      IMPLICIT NONE

      PRIVATE    

      PUBLIC :: pbc, real_to_scaled, scaled_to_real, trans_pbc
      
      INTERFACE pbc
        MODULE PROCEDURE pbc1,pbc2,pbc3
      END INTERFACE

CONTAINS

      FUNCTION pbc1(r,cell,rcell) RESULT(r_pbc)
        REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
        REAL(dp), DIMENSION(3,3), INTENT(IN)     :: cell, rcell
        REAL(KIND=dp), DIMENSION(3)              :: r_pbc
        REAL(KIND=dp), DIMENSION(3)              :: s

        s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
        s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
        s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)
        s(1) = s(1) - ANINT(s(1))
        s(2) = s(2) - ANINT(s(2))
        s(3) = s(3) - ANINT(s(3))
        r_pbc(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
        r_pbc(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
        r_pbc(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

      END FUNCTION pbc1

      FUNCTION pbc2(r,cell,rcell,nl) RESULT(r_pbc)
        REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
        REAL(dp), DIMENSION(3,3), INTENT(IN)     :: cell, rcell
        INTEGER, DIMENSION(3), INTENT(IN)        :: nl
        REAL(KIND=dp), DIMENSION(3)              :: r_pbc

        REAL(KIND=dp), DIMENSION(3)              :: s

      s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
      s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
      s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)
      s(1) = s(1) - REAL(NINT(s(1)) - nl(1),dp)
      s(2) = s(2) - REAL(NINT(s(2)) - nl(2),dp)
      s(3) = s(3) - REAL(NINT(s(3)) - nl(3),dp)
      r_pbc(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
      r_pbc(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
      r_pbc(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

      END FUNCTION pbc2

      FUNCTION pbc3(ra,rb,cell,rcell) RESULT(rab_pbc)
        REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb
        REAL(dp), DIMENSION(3,3), INTENT(IN)     :: cell, rcell
        REAL(KIND=dp), DIMENSION(3)              :: rab_pbc

        INTEGER                                  :: icell, jcell, kcell
        INTEGER, DIMENSION(3)                    :: periodic
        REAL(KIND=dp)                            :: rab2, rab2_pbc
        REAL(KIND=dp), DIMENSION(3)              :: r, ra_pbc, rab, rb_image, &
                                                    rb_pbc, s2r

!      CALL get_cell(cell=cell,periodic=periodic)

       ra_pbc(:) = pbc(ra(:),cell, rcell)
       rb_pbc(:) = pbc(rb(:),cell, rcell)

       rab2_pbc = HUGE(1.0_dp)

      DO icell=-periodic(1),periodic(1)
        DO jcell=-periodic(2),periodic(2)
          DO kcell=-periodic(3),periodic(3)
            r = REAL((/icell,jcell,kcell/),dp)
            CALL scaled_to_real(s2r,r,cell)
            rb_image(:) = rb_pbc(:) + s2r
            rab(:) = rb_image(:) - ra_pbc(:)
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            IF (rab2 < rab2_pbc) THEN
              rab2_pbc = rab2
              rab_pbc(:) = rab(:)
            END IF
          END DO
        END DO
      END DO
      END FUNCTION pbc3


      SUBROUTINE real_to_scaled(s,r,rcell)
       REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: s
       REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
       REAL(dp), DIMENSION(3,3), INTENT(IN)     :: rcell

      s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
      s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
      s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)

      END SUBROUTINE real_to_scaled


      SUBROUTINE scaled_to_real(r,s,cell)
       REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: r
       REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: s
       REAL(dp), DIMENSION(3,3), INTENT(IN)     :: cell

       r(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
       r(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
       r(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)
      END SUBROUTINE scaled_to_real





      subroutine trans_pbc(r,cell,cell_inv,r_pbc)
!  *****************************************************************************
!  Apply the periodic boundary conditions  to a position vector r.
!  author  shanghui
!  date    2010.01.18
! ****************************************************************************
!-------------input--------------------------!
!r(3)  : the vector
!scell(3,3) :  the supercell
!-------------output-------------------------!
!r_pbc(3): the vector in the pbc

       real(dp) r(3), cell(3,3), r_pbc(3), &
                cell_inv(3,3),s(3)

       real(dp), parameter  :: tiny=1.e-12_dp

!       call reclat(cell,cell_inv,0)
       s(1) = cell_inv(1,1)*r(1) + cell_inv(2,1)*r(2)  &
              + cell_inv(3,1)*r(3)
       s(2) = cell_inv(1,2)*r(1) + cell_inv(2,2)*r(2)  &
              + cell_inv(3,2)*r(3)
       s(3) = cell_inv(1,3)*r(1) + cell_inv(2,3)*r(2)  &
              + cell_inv(3,3)*r(3)
       if( (dabs(s(1))-0.5d0) .gt. tiny )  s(1) = s(1) - ANINT(s(1))
       if( (dabs(s(2))-0.5d0) .gt. tiny )  s(2) = s(2) - ANINT(s(2))
       if( (dabs(s(3))-0.5d0) .gt. tiny )  s(3) = s(3) - ANINT(s(3))

       r_pbc(1) =cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
       r_pbc(2) =cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
       r_pbc(3) =cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

      end subroutine trans_pbc

      END MODULE cell_pbc
