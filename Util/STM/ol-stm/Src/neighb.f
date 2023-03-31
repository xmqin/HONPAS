! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE NEIGHB( CELL, RANGE, NA, XA, XPO, IA, ISC,
     .                   NNA, JAN, XIJ, R2IJ, FIRST )

C ********************************************************************
C Finds the neighbours of an atom in a cell with periodic boundary 
C conditions. This is an interface to routine ranger, which has
C extended functionalities.
C Written by J.M.Soler. March 1997.
C *********** INPUT **************************************************
C REAL*8  CELL(3,3) : Unit cell vectors CELL(IXYZ,IVECT)
C REAL*8  RANGE     : Maximum distance of neighbours required
C INTEGER NA        : Number of atoms
C REAL*8  XA(3,NA)  : Atomic positions in cartesian coordinates
C REAL*8  XPO(3)    : Position of a point in real space (cartesian coor)
C INTEGER IA        : Atom whose neighbours are needed.
C                     A routine initialization must be done by
C                     a first call with IA = 0
C INTEGER ISC       : Single-counting switch (0=No, 1=Yes). If ISC=1,
C                     only neighbours with JA.LE.IA are included in JAN
C INTEGER NNA       : Size of arrays JAN, XIJ and R2IJ
C *********** OUTPUT *************************************************
C INTEGER NNA        : Number of neighbour atoms within RANGE of IA
C INTEGER JAN(NNA)   : Atom-index of neighbours
C REAL*8  XIJ(3,NNA) : Vectors from atom IA to neighbours
C REAL*8  R2IJ(NNA)  : Squared distances to neighbours
C *********** UNITS **************************************************
C Units of CELL, RANGE and XA are arbitrary but must be equal
C *********** SUBROUTINES USED ***************************************
C CHKDIM, DISMIN, DOT, PRMEM, RECLAT
C *********** BEHAVIOUR **********************************************
C CPU time and memory scale linearly with the number of atoms, for
C   sufficiently large numbers.
C If internal dimension variables are too small, an include file named
C   neighb.h is printed with the required dimensions and the routine
C   stops, asking to be recompiled. Then, neighb.h is automatically
C   included in the new compilation, if the source file is in the same 
C   directory where the program has run. Initially, you can make all 
C   the parameters in neighb.h equal to 1.
C Different ranges can be used for different atoms, but for good 
C   performance, the largest range should be used in the initial
C   call (with IA=0).
C There are no limitations regarding cell shape or size. The range may
C   be larger than the cell size, in which case many 'images' of the
C   same atom will be included in the neighbour list, with different
C   interatomic vectors XIJ and distances R2IJ.
C The atom IA itself is included in the neighbour list, with zero
C   distance. You have to discard it if you want so.
C If the number of neighbour atoms found is larger than the size of
C   the arrays JAN, XIJ and R2IJ, i.e. if NNAout > NNAin, these arrays
C   are filled only up to their size NNAin. With dynamic memory
C   allocation, this allows to find first the required array sizes
C   and then find the neighbours. Notice however that no warning is
C   given, so that you should always check that NNAout.LE.NNAin.
C *********** USAGE **************************************************
C Sample usage for a molecular dynamics simulation:
C    DIMENSION JAN(MAXNNA), XIJ(3,MAXNNA)
C    Define CELL and initial positions XA
C    DO ITER = 1,NITER                (Molecular dynamics iteration)
C      NNA = MAXNNA
C      CALL NEIGHB( CELL, RANGE, NA, XA, 0, 1, NNA, JAN, XIJ, R2IJ )
C      IF (NNA .GT. MAXNNA) STOP 'MAXNNA too small'
C      Initialize to zero all atomic forces FA(IX,IA)
C      DO IA = 1,NA                   (Loop on atoms)
C        NNA = MAXNNA
C        CALL NEIGHB( CELL, RANGE, NA, XA, IA, 1, NNA, JAN, XIJ, R2IJ )
C        IF (NNA .GT. MAXNNA) STOP 'MAXNNA too small'
C        DO IN = 1,NNA                (Loop on neighbours of IA)
C          JA = JAN(IN)               (Atomic index of neighbour)
C          RIJ = SQRT(R2IJ(IN))       (Interatomic distance)
C          IF (RIJ .GT. 1.D-12) THEN  (Discard atom itself)
C            Find interatomic force FIJ( RIJ )
C            DO IX = 1,3
C              FA(IX,IA) = FA(IX,IA) - FIJ * XIJ(IX,IN) / RIJ
C              FA(IX,JA) = FA(IX,JA) + FIJ * XIJ(IX,IN) / RIJ
C            ENDDO
C          ENDIF
C        ENDDO
C      ENDDO
C      Move atomic positions XA       (Molecular dynamics step)
C    ENDDO
C ********************************************************************

C Next line is non-standard but may be supressed
      IMPLICIT NONE

C Define space dimension
      INTEGER NX
      PARAMETER ( NX = 3 )

C Argument types and dimensions
      INTEGER           IA, ISC, JAN(*), NA, NNA
      DOUBLE PRECISION  CELL(NX,NX), RANGE, R2IJ(*),
     .                  XA(NX,NA), XIJ(NX,*)
      DOUBLE PRECISION  XPO(3)

C Internal variables
      LOGICAL           FRSTME, SAMCEL
      INTEGER           IAMOVE, IX, JX
      DOUBLE PRECISION  CELAST(NX,NX), RGLAST, X0(NX)
      LOGICAL FIRST
      SAVE  
      DATA FRSTME / .TRUE. /
      DATA IAMOVE / 0 /
      DATA X0 / NX*0.D0 /

C Initialization section
      IF ( FRSTME .OR. IA.LE.0 .OR. RANGE.GT.RGLAST ) THEN

C       Find if cell or range have changed
        SAMCEL = .TRUE.
        DO 20 IX = 1,NX
          DO 10 JX = 1,NX
            IF (CELL(JX,IX) .NE. CELAST(JX,IX)) SAMCEL = .FALSE.
   10     CONTINUE
   20   CONTINUE
        IF (RANGE .NE. RGLAST) SAMCEL = .FALSE.

C       Cell initializations
        IF (.NOT.SAMCEL) THEN

C         Store cell and range for comparison in subsequent calls
          DO 40 IX = 1,NX
            DO 30 JX = 1,NX
              CELAST(JX,IX) = CELL(JX,IX)
   30       CONTINUE
   40     CONTINUE
          RGLAST = RANGE
          FRSTME = .FALSE.

C         Notify to RANGER that CELL has changed
          CALL RANGER( 'CELL', NX, CELL, RANGE, NA, XA,
     .                 NA, IAMOVE,
     .                 IA, ISC, X0,
     .                 NNA, JAN, XIJ, R2IJ )
        ENDIF

C       Notify to RANGER that atoms have moved
        CALL RANGER( 'MOVE', NX, CELL, RANGE, NA, XA,
     .               NA, IAMOVE,
     .               IA, ISC, X0,
     .               NNA, JAN, XIJ, R2IJ )

      ENDIF

C Find neighbours of atom IA
      IF ( .NOT. FIRST ) THEN

        DO IX = 1,NX
          X0(IX) = XPO(IX)
        ENDDO

        CALL RANGER( 'FIND', NX, CELL, RANGE, NA, XA,
     .               NA, IAMOVE,
     .               IA, ISC, X0,
     .               NNA, JAN, XIJ, R2IJ )

      ENDIF

      END



