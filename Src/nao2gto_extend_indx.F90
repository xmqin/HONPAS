      SUBROUTINE extend_indx( indsc1,indsc2,nsc,nuotot,norb)
! *********************************************************************
! idea from Initializes listsc by J.M.Soler.
! little change  by shanghui,2008.11.20
! ************** INPUT ************************************************
! INTEGER NSC(3) : Number of cells in each supercell direction:
!                    SuperCell(ix,i) = CELL(ix,i) * NSC(i)
! INTEGER NUO    : Number of orbitals in unit cell
! *********************************************************************
!***************output*************************************************
! ind1(8*norb) : give the true supercell index from extend supercell
! ind2(norb)   : change the true supercell index to the extend big one,
!                in my vertion ,it is the index in the 8th extend
!                supercell,
!                and in Soler's vertion ,it is the index in the 1st
!                extend supercell.
!**************************************************************************************
      IMPLICIT NONE

      integer, intent(in)   :: nsc(3)
      integer, intent(in)   :: nuotot
      integer, intent(in)   :: norb
      integer, intent(inout)  :: indsc1(8*norb)
      integer, intent(inout)  :: indsc2(norb)
!      integer, dimension(:), allocatable, save :: IND1, IND2

      INTEGER I1, I2, I3, IC, J1, J2, J3, JC, &
             KUO, LASTIO, LASTJO, NCELLS !NO


      NCELLS = NSC(1) * NSC(2) * NSC(3)
!      NO=NCELLS*nuotot
!     (twice as large in each direction)
      DO J3 = 0,2*NSC(3)-1
      DO J2 = 0,2*NSC(2)-1
      DO J1 = 0,2*NSC(1)-1

!       I1,I2,I3 fold the extended supercell into the normal supercell
        I1 = MOD(J1,NSC(1))
        I2 = MOD(J2,NSC(2))
        I3 = MOD(J3,NSC(3))

!       IC,JC are the cell indexes in the normal and extended supercells
        IC = I1 +   NSC(1)*I2 +   NSC(1)*NSC(2)*I3
        JC = J1 + 2*NSC(1)*J2 + 4*NSC(1)*NSC(2)*J3

!       LASTIO,LASTJO are the indexes of last orbitals of
!       previous unit cells
        LASTIO = IC * NUOTOT
        LASTJO = JC * NUOTOT

!       IND1 folds the extended supercell onto the single supercell
        DO KUO = 1,NUOTOT
          INDSC1(LASTJO+KUO) = LASTIO + KUO
        ENDDO

!       IND2 inserts the single supercell into the extended supercell
        IF ( (I1+nsc(1)).EQ.J1 .AND. (I2+nsc(2)).EQ.J2  &
            .AND. (I3+nsc(3)).EQ.J3 ) THEN
          DO KUO = 1,NUOTOT
            INDSC2(LASTIO+KUO) = LASTJO + KUO
          ENDDO
        ENDIF

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE extend_indx
