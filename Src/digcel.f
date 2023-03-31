! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
C $Id: digcel.f,v 1.2 1999/01/31 10:53:52 emilio Exp $

      SUBROUTINE DIGCEL( CELL, MSCELL, UCELL, SCELL, NSC, ISDIAG )

C *********************************************************************
C Inputs a unit cell and a general supercell:
C    SuperCell(ix,i) = Sum_j( CELL(ix,j) * MSCELL(j,i) )
C and outputs an equivalent unit cell and an equivalent supercell 
C that are 'diagonal', ie:
C    SCELL(ix,i) = UCELL(ix,i) * NSC(i)
C Written by J.M.Soler, June 1998.
C Based on an algorithm of Juana Moreno
C Ref: J.Moreno and J.M.Soler, PRB 45, 13891 (1992)
C *********** Input ***************************************************
C Real*8  CELL(3,3)   : Unit cell vectors CELL(Ixyz,Ivector)
C Integer MSCELL(3,3) : Supercell vectors in units of CELL:
C                         SuperCell(ix,i) = Sum_j( CELL(ix,j) *
C                                                 MSCELL(j,i) )
C *********** Output **************************************************
C Real*8  UCELL(3,3)  : New equivalent unit cell vectors
C Real*8  SCELL(3,3)  : New supercell vectors
C Integer NSC(3)      : New diagonal elements of MSCELL, i.e.:
C                         SCELL(ix,i) = UCELL(ix,i) * NSC(i)
C Logical ISDIAG      : Is MSCELL already diagonal?
C *********************************************************************
      IMPLICIT          NONE
      LOGICAL           ISDIAG
      INTEGER           MSCELL(3,3), NSC(3)
      DOUBLE PRECISION  CELL(3,3), SCELL(3,3), UCELL(3,3)
      EXTERNAL          IDIAG

C Internal variables
      INTEGER           I, IX, J,
     .                  MAUX(3,3,2), MLEFT(3,3), MRIGHT(3,3), MSC(3,3)
      DOUBLE PRECISION  SCIN(3,3)

C Check if MSCELL is already diagonal
      ISDIAG = .TRUE.
      DO 20 J = 1,3
        DO 10 I = 1,3
          IF (I.NE.J .AND. MSCELL(I,J).NE.0) ISDIAG = .FALSE.
   10   CONTINUE
   20 CONTINUE

C If input MSCELL was already diagonal
      IF (ISDIAG) THEN

C       Copy input unit cell and supercell to output arrays
        DO 40 I = 1,3
          NSC(I) = MSCELL(I,I)
          DO 30 IX = 1,3
            UCELL(IX,I) = CELL(IX,I)
            SCELL(IX,I) = CELL(IX,I) * NSC(I)
   30     CONTINUE
   40   CONTINUE

      ELSE

C       Find initial supercell
        DO 60 I = 1,3
          DO 50 IX = 1,3
            SCIN(IX,I) = CELL(IX,1) * MSCELL(1,I) +
     .                   CELL(IX,2) * MSCELL(2,I) +
     .                   CELL(IX,3) * MSCELL(3,I)
   50     CONTINUE
   60  CONTINUE

C       Diagonalize MSCELL matrix
        CALL IDIAG( 3, MSCELL, MSC, MLEFT, MRIGHT, MAUX )

C       Find new supercell and unit cell
        DO 80 I = 1,3
          NSC(I) = MSC(I,I)
          DO 70 IX = 1,3
            SCELL(IX,I) = SCIN(IX,1) * MRIGHT(1,I) +
     .                    SCIN(IX,2) * MRIGHT(2,I) +
     .                    SCIN(IX,3) * MRIGHT(3,I)
            UCELL(IX,I) = SCELL(IX,I) / NSC(I)
   70     CONTINUE
   80   CONTINUE

      ENDIF
      END


                  
