! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE PLSURF( ORIGIN, CELL, NMESH, NSPAN,
     .                   F, FVALUE, NCOLOR, COLOR,
     .                   NT, COLMIN, COLAVG, COLMAX )

C *******************************************************************
C Plots the 3-D shape of a surface with constant function value.
C The surface is determined by the condition function=value, and
C it is ploted as a set of triangles in 3-D. The function must
C be given in a regular 3-D grid of points.
C Optionally, other function(s) determine the color(s) of the surface.
C Notice single precision in this version
C Written by J.M.Soler. October 1997.
C ************************* INPUT ***********************************
C REAL    ORIGIN(3)    : Coordinates of cell origin
C REAL    CELL(3,3)    : Unit cell vectors CELL(ixyz,ivector)
C INTEGER NMESH(3)     : Number of mesh divisions of each vector
C INTEGER NSPAN(3)     : Physical dimensions of arrays F and
C                        COLOR (or memory span between array elements)
C                        See usage section for more information
C REAL    F(*)         : Function such that F=FVALUE determines
C                        the shape of the solid surface.
C REAL    FVALUE       : Value such that F=FVALUE
C                        determines the shape of the solid surface.
C INTEGER NCOLOR       : Number of color functions (may be zero).
C REAL    COLOR(NCOLOR,*) : Function(s) which determines the
C                        color(s) of the solid surface.
C ************************* OUTPUT **********************************
C INTEGER NT           : Number of ploted triangles
C REAL    COLMIN(NCOLOR) : Minimun value(s) of color function(s)
C                          at the vertex triangles
C REAL    COLAVG(NCOLOR) : Average value(s) of color function(s)
C REAL    COLMAX(NCOLOR) : Maximun value(s) of color function(s)
C ************************ USAGE ************************************
C Typical calls for different array dimensions:
C     PARAMETER ( MAXP = 1000000 )
C     DIMENSION NMESH(3), F(MAXP), COLOR(MAXP)
C     Find CELL vectors
C     Find shape and color functions at N1*N2*N3 mesh points (less than 
C     MAXP) and place them consecutively in arrays F and COLOR
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL PLSURF( ORIGIN, CELL, NMESH, NMESH, F, FVALUE, 1, COLOR, NT)
C Or alternatively:
C     PARAMETER ( M1=100, M2=100, M3=100 )
C     DIMENSION NMESH(3), NSPAN(3), F(M1,M2,M3), COLOR(M1,M2,M3)
C     DATA NSPAN / M1, M2, M3 /
C     Find CELL vectors
C     Find F and COLOR at N1*N2*N3 mesh points
C     NMESH(1) = N1
C     NMESH(2) = N2
C     NMESH(3) = N3
C     CALL PLSURF( ORIGIN, CELL, NMESH, NSPAN, F, FVALUE, 1, COLOR, NT)
C ********* ROUTINES CALLED *****************************************
C CUBE, PLTR3D, RECLAT
C *******************************************************************
C    6  10        20        30        40        50        60        7072

C Next line is nonstandard but may be suppressed
      IMPLICIT NONE

C Argument types and dimensions
      INTEGER 
     .  NCOLOR, NMESH(3), NSPAN(3), NT
      REAL
     .  CELL(3,3), COLAVG(NCOLOR), COLMAX(NCOLOR),
     .  COLMIN(NCOLOR), COLOR(NCOLOR,*),
     .  F(*), FVALUE, ORIGIN(3)
      EXTERNAL
     .  CUBE, RECLAT

C Fix the order of the numerical derivatives
C NN is the number of points used in each coordinate and direction,
C i.e. a total of 6*NN neighbour points is used to find the gradients
      INTEGER NN
      PARAMETER ( NN = 1 )

C 	MAXCOL must be the same in routines solid, cube and pyramid
      INTEGER MAXCOL
      PARAMETER ( MAXCOL = 10 )

C Local variables and arrays
      LOGICAL
     .  ALLHGH, ALLLOW, HIGH
      INTEGER
     .  I1, I2, I3, IC, ICOLOR, IN, IP, IPAR, IT, IV, IX,
     .  J1, J2, J3, JN, JP(3,-NN:NN), JX,
     .  K1, K2, K3, L1, L2, L3, NP, NTC
      REAL
     .  CT(3,MAXCOL), CTC(MAXCOL,3,10),
     .  DGDM(-NN:NN), DGIDFJ(3,3,-NN:NN), DMDX(3,3), DXDM(3,3),
     .  F1, F2, FC(MAXCOL,0:1,0:1,0:1), FS(0:1,0:1,0:1),
     .  GF(3,0:1,0:1,0:1), GTC(3,3,10), XF(3,0:1,0:1,0:1), XTC(3,3,10)

C Avoid invading undeclared memory
      NCOLOR = MIN( NCOLOR, MAXCOL )

C Find weights of numerical derivation from Lagrange interp. formula
      DO IN = -NN,NN
        F1 = 1
        F2 = 1
        DO JN = -NN,NN
          IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
          IF (JN.NE.IN)               F2 = F2 * (IN - JN)
        ENDDO
        DGDM(IN) = F1 / F2
      ENDDO
      DGDM(0) = 0

C Find total number of mesh points
      NP = NMESH(1) * NMESH(2) * NMESH(3)

C Find Jacobian matrix dx/dmesh and its inverse
      DO IC = 1,3
        DO IX = 1,3
          DXDM(IX,IC) = CELL(IX,IC) / NMESH(IC)
        ENDDO
      ENDDO
      CALL RECLAT( DXDM, DMDX, 0 )

C Find the weights for the derivative d(gradF(i))/d(F(j)) of
C the gradient at point i with respect to the value at point j
      DO IN = -NN,NN
        DO IC = 1,3
          DO IX = 1,3
            DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
          ENDDO
        ENDDO
      ENDDO

C Initialize output
      NT = 0
      DO ICOLOR = 1,NCOLOR
        COLMIN(ICOLOR) =  1.E30
        COLMAX(ICOLOR) = -1.E30
        COLAVG(ICOLOR) =  0.
      ENDDO

C Loop on mesh points
      DO K3 = 0,NMESH(3)-1
      DO K2 = 0,NMESH(2)-1
      DO K1 = 0,NMESH(1)-1

C       Check if all cube vertices are above or below equi-surface.
        ALLHGH = .TRUE.
        ALLLOW = .TRUE.
        DO L3 = 0,1
        DO L2 = 0,1
        DO L1 = 0,1

C         Find mesh index of this point
          I1 = MOD( K1+L1, NMESH(1) )
          I2 = MOD( K2+L2, NMESH(2) )
          I3 = MOD( K3+L3, NMESH(3) )
          IP = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3

C         Find if this point is above FVALUE
          HIGH = (F(IP) .GT. FVALUE)
          ALLHGH = (ALLHGH .AND. HIGH)
          ALLLOW = (ALLLOW .AND. .NOT.HIGH)
        ENDDO
        ENDDO
        ENDDO

C       Skip the rest if the equi-surface does not cross this cube.
        IF (ALLHGH .OR. ALLLOW) GOTO 10
        
        DO L3 = 0,1
        DO L2 = 0,1
        DO L1 = 0,1
          I1 = K1 + L1
          I2 = K2 + L2
          I3 = K3 + L3

C         Find vertex coordinates
          DO IX = 1,3
            XF(IX,L1,L2,L3) = DXDM(IX,1) * I1 +
     .                        DXDM(IX,2) * I2 +
     .                        DXDM(IX,3) * I3 + ORIGIN(IX)
          ENDDO

C         Find mesh indexes of neighbour points
          I1 = MOD( I1, NMESH(1) )
          I2 = MOD( I2, NMESH(2) )
          I3 = MOD( I3, NMESH(3) )
          DO IN = -NN,NN
            J1 = MOD( I1+IN+100*NMESH(1), NMESH(1) )
            J2 = MOD( I2+IN+100*NMESH(2), NMESH(2) )
            J3 = MOD( I3+IN+100*NMESH(3), NMESH(3) )
            JP(1,IN) = 1 + J1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * I3
            JP(2,IN) = 1 + I1 + NSPAN(1) * J2 + NSPAN(1) * NSPAN(2) * I3
            JP(3,IN) = 1 + I1 + NSPAN(1) * I2 + NSPAN(1) * NSPAN(2) * J3
          ENDDO

C         Find mesh index of this point
          IP = JP(1,0)

C         Find COLOR and F and its gradient at this point
          DO ICOLOR = 1,NCOLOR
            FC(ICOLOR,L1,L2,L3) = COLOR(ICOLOR,IP)
          ENDDO
          FS(L1,L2,L3) = F(IP)
          DO IX = 1,3
            GF(IX,L1,L2,L3) = 0
            DO IN = -NN,NN
              GF(IX,L1,L2,L3) = GF(IX,L1,L2,L3) +
     .                 DGIDFJ(IX,1,IN) * F(JP(1,IN)) +
     .                 DGIDFJ(IX,2,IN) * F(JP(2,IN)) +
     .                 DGIDFJ(IX,3,IN) * F(JP(3,IN))
            ENDDO
          ENDDO

        ENDDO
        ENDDO
        ENDDO

        IPAR = MOD( K1+K2+K3, 2 )
        CALL CUBE( IPAR, FS, XF, GF, FVALUE, NCOLOR, FC,
     .             NTC, XTC, GTC, CTC )

        DO IT = 1,NTC
          NT = NT + 1
          DO IV = 1,3
            DO ICOLOR = 1,NCOLOR
              CT(IV,ICOLOR) = CTC(ICOLOR,IV,IT)
              COLMIN(ICOLOR) = MIN( COLMIN(ICOLOR), CT(IV,ICOLOR) )
              COLMAX(ICOLOR) = MAX( COLMAX(ICOLOR), CT(IV,ICOLOR) )
              COLAVG(ICOLOR) = COLAVG(ICOLOR) + CT(IV,ICOLOR)
            ENDDO
          ENDDO
          CALL PLTR3D( XTC(1,1,IT), GTC(1,1,IT), NCOLOR, CT )
        ENDDO

   10 ENDDO
      ENDDO
      ENDDO

      DO ICOLOR = 1,NCOLOR
        COLAVG(ICOLOR) = COLAVG(ICOLOR) / (3*NT)
      ENDDO
      END




      SUBROUTINE CUBE( IPAR, F, XF, GF, FV, NCOL, FCOL,
     .                 NT, XT, GT, TCOL )

      IMPLICIT NONE
C 	    MAXCOL must be the same in routines solid, cube and pyramid
      INTEGER     MAXCOL
      PARAMETER ( MAXCOL = 10 )
      INTEGER
     .  IPAR, NCOL, NT
      REAL
     .  TCOL(MAXCOL,3,*), FCOL(MAXCOL,8), F(8), FV,
     .  GF(3,8), GT(3,3,*), XF(3,8), XT(3,3,*)

      INTEGER
     .  ICOL, INDXV(4,5,0:1), IP, IVC, IVP, IX, NTP
      REAL
     .  FP(4), GP(3,4), PCOL(MAXCOL,4), XP(3,4)

      DATA (((INDXV(IVP,IP,IPAR),IVP=1,4),IPAR=0,1),IP=1,5)
     .  / 1,2,3,5,   1,2,4,6,
     .    2,3,4,8,   1,3,4,7,
     .    2,5,6,8,   1,5,6,7,
     .    3,5,7,8,   4,6,7,8,
     .    2,3,5,8,   1,4,6,7 /

      NT = 0

C     For each of the five pyramids inside the cube
      DO IP = 1,5

C       Copy coordinates, gradients and colors for this pyramid
        DO IVP = 1,4
          IVC = INDXV(IVP,IP,IPAR)
          FP(IVP) = F(IVC)
          DO IX = 1,3
            XP(IX,IVP) = XF(IX,IVC)
            GP(IX,IVP) = GF(IX,IVC)
          ENDDO
          DO ICOL = 1,NCOL
            PCOL(ICOL,IVP) = FCOL(ICOL,IVC)
          ENDDO
        ENDDO

C       Find the tringles covering the cut plane for this pyramid
        CALL PYRAMID( FP, XP, GP, FV, NCOL, PCOL,
     .                NTP, XT(1,1,NT+1), GT(1,1,NT+1),
     .                TCOL(1,1,NT+1) )
        NT = NT + NTP
      ENDDO
      END



      SUBROUTINE PYRAMID( F, X, G, VAL, NC, C, NT, XT, GT, CT )
      
      IMPLICIT NONE
C 	    MAXCOL must be the same in routines solid, cube and pyramid
      INTEGER     MAXCOL
      PARAMETER ( MAXCOL = 10 )
      INTEGER  
     .  NC, NT
      REAL
     .  C(MAXCOL,4), CT(MAXCOL,3,*), F(4), G(3,4), GT(3,3,*),
     .  VAL, X(3,4), XT(3,3,*)

      INTEGER
     .  IC, IP, IT, IV, IV1, IV2, IVT, IVTP(2,3,2,0:7),
     .  IX, K(4), KK(0:7), NTP(0:7)
      REAL
     .  W1, W2

      DATA ( KK(IP), NTP(IP),
     .      (((IVTP(IV,IVT,IT,IP),IV=1,2),IVT=1,3),IT=1,2), IP=0,7)
     .  / 0000,  0,  0,0, 0,0, 0,0,  0,0, 0,0, 0,0,
     .    1000,  1,  1,2, 1,3, 1,4,  0,0, 0,0, 0,0,
     .    0100,  1,  1,2, 2,3, 2,4,  0,0, 0,0, 0,0,
     .    1100,  2,  1,3, 1,4, 2,3,  1,4, 2,3, 2,4,
     .    0010,  1,  1,3, 2,3, 3,4,  0,0, 0,0, 0,0,
     .    1010,  2,  1,2, 1,4, 2,3,  1,4, 2,3, 3,4,
     .    0110,  2,  1,2, 1,3, 2,4,  1,3, 2,4, 3,4,
     .    1110,  1,  1,4, 2,4, 3,4,  0,0, 0,0, 0,0 /

C     Find index IP which determines how the pyramid is cut
      DO IV = 1,4
        IF (F(IV) .GT. VAL) THEN
          K(IV) = 1
        ELSE
          K(IV) = 0
        ENDIF
      ENDDO
      IP = K(1) + K(2)*2 + K(3)*4 + K(4)*8
      IF (IP .GE. 8) IP = 15 - IP

C     Find zero, one or two triangles covering the cut plane
      NT = NTP(IP)
      DO IT = 1,NT
        DO IVT = 1,3
          IV1 = IVTP(1,IVT,IT,IP)
          IV2 = IVTP(2,IVT,IT,IP)
          W1 = ABS(F(IV2)-VAL) / ABS(F(IV2)-F(IV1))
          W2 = 1.D0 - W1
          DO IX = 1,3
            XT(IX,IVT,IT) = X(IX,IV1)*W1 + X(IX,IV2)*W2
            GT(IX,IVT,IT) = G(IX,IV1)*W1 + G(IX,IV2)*W2
          ENDDO
          DO IC = 1,NC
            CT(IC,IVT,IT) = C(IC,IV1)*W1 + C(IC,IV2)*W2
          ENDDO
        ENDDO
      ENDDO
      END




      SUBROUTINE RECLAT (A,B,IOPT)

C  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
C  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      IMPLICIT NONE
      INTEGER I, IOPT
      REAL
     .  A(3,3), B(3,3), C, CI, PI
      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
      ENDDO
      END
