! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      MODULE LISTSC_MODULE

      private 
      public :: LISTSC_INIT, LISTSC

      integer, pointer, save :: IND1(:), IND2(:)
      logical, save          :: nullified_pointers = .false.

      CONTAINS

      SUBROUTINE LISTSC_INIT( NSC, NUO )
      use alloc, only: re_alloc

C *********************************************************************
C Initializes listsc.
C Written by J.M.Soler, July 1999
C ************** INPUT ************************************************
C INTEGER NSC(3) : Number of cells in each supercell direction:
C                    SuperCell(ix,i) = CELL(ix,i) * NSC(i)
C INTEGER NUO    : Number of orbitals in unit cell
C *********************************************************************

      IMPLICIT NONE

      integer, intent(in)   :: nsc(3)
      integer, intent(in)   :: nuo

      INTEGER I1, I2, I3, IC, J1, J2, J3, JC,
     .        KUO, LASTIO, LASTJO, NCELLS, NO
      EXTERNAL MEMORY
      
      NCELLS = NSC(1) * NSC(2) * NSC(3)
      NO = NUO * NCELLS

      if (.not. nullified_pointers) then
         nullify(ind1, ind2)
         nullified_pointers = .true.
      endif

      call re_alloc(IND1,1,8*NO,name="ind1",routine="listsc_init")
      call re_alloc(IND2,1,NO,name="ind2",routine="listsc_init")

C     Loop on unit cells of an extended supercell 
C     (twice as large in each direction)
      DO J3 = 0,2*NSC(3)-1
      DO J2 = 0,2*NSC(2)-1
      DO J1 = 0,2*NSC(1)-1

C       I1,I2,I3 fold the extended supercell into the normal supercell
        I1 = MOD(J1,NSC(1))
        I2 = MOD(J2,NSC(2))
        I3 = MOD(J3,NSC(3))

C       IC,JC are the cell indexes in the normal and extended supercells
        IC = I1 +   NSC(1)*I2 +   NSC(1)*NSC(2)*I3
        JC = J1 + 2*NSC(1)*J2 + 4*NSC(1)*NSC(2)*J3

C       LASTIO,LASTJO are the indexes of last orbitals of
C       previous unit cells
        LASTIO = IC * NUO
        LASTJO = JC * NUO

C       IND1 folds the extended supercell onto the single supercell
        DO KUO = 1,NUO
          IND1(LASTJO+KUO) = LASTIO + KUO
        ENDDO

C       IND2 inserts the single supercell into the extended supercell
        IF (I1.EQ.J1 .AND. I2.EQ.J2 .AND. I3.EQ.J3) THEN
          DO KUO = 1,NUO
            IND2(LASTIO+KUO) = LASTJO + KUO
          ENDDO
        ENDIF

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE LISTSC_INIT



      FUNCTION LISTSC( IO, IUO, JUO )
C *********************************************************************
C Translates a listed neighbour from unit cell to supercell
C Written by J.M.Soler, July 1999
C ************** INPUT ************************************************
C INTEGER IO  : Orbital index
C INTEGER IUO : Orbital within unit cell equivalent to IO
C INTEGER JUO : A neighbor orbital of IUO (not necessarily in unit cell)
C ************** OUTPUT ***********************************************
C INTEGER LISTSC : Orbital, equivalent to JUO, which is related to IO
C                  exactly as JUO is related to IUO
C ************** BEHAVIOR *********************************************
C Subroutine LISTSC_INIT must be called once before any call to LISTSC,
C   and it must be called again if NSC or NUO change.
C Arguments must be positive and IO.LE.NO, IUO.LE.NUO, JUO.LE.NO. This 
C   is NOT checked, and a core dump is likely if not true.
C *********************************************************************

      IMPLICIT NONE

      integer listsc
      INTEGER, intent(in) ::  IO, IUO, JUO

      LISTSC = IND1( IND2(IO) + IND2(JUO) - IND2(IUO) )

      END FUNCTION LISTSC

      END MODULE LISTSC_MODULE
