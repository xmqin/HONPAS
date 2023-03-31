! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
C $Id: dismin.f,v 1.4 2000/01/11 10:13:50 jgale Exp $

      DOUBLE PRECISION FUNCTION DISMIN( CELL, X )

C FINDS THE MINIMUM DISTANCE FROM A PARALLELEPIPED CELL TO A POINT X.
C THIS IS A PROVISIONAL VERSION WICH RETURNS ONLY AN APPROXIMATE
C VALUE: THE MINIMUM DISTANCE TO A FINITE SET OF POINTS ON THE SURFACE,
C WHERE PARAMETER N BELOW DETERMINES HOW MANY POINTS ARE USED.
C WRITTEN BY J.M.SOLER. NOV'96.

      IMPLICIT NONE
      INTEGER I, I1, I2, I3, N
      DOUBLE PRECISION CELL(3,3), X(3), DX, D2, D2MIN
      PARAMETER ( N = 1 )

      D2MIN = 1.D30
      DO I1 = 0,N
        DO I2 = 0,N
          DO I3 = 0,N
            IF ( I1.EQ.0 .OR. I1.EQ.N .OR.
     .           I2.EQ.0 .OR. I2.EQ.N .OR.
     .           I3.EQ.0 .OR. I3.EQ.N      ) THEN
              D2 = 0.0D0
              DO I = 1,3
                DX = (CELL(I,1)*I1+CELL(I,2)*I2+CELL(I,3)*I3)/N
                DX = DX - X(I)
                D2 = D2 + DX * DX
              ENDDO
              D2MIN = MIN( D2, D2MIN )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DISMIN = SQRT( D2MIN )
      END

