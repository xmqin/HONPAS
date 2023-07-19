c
      SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
C
C
C     .. Parameters ..
      INTEGER NTRY
      PARAMETER (NTRY=50)
      DOUBLE PRECISION FACTOR
      PARAMETER (FACTOR=1.6D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      LOGICAL SUCCES
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F1,F2
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (X1.EQ.X2) THEN
          WRITE (*,FMT=*) 'You have to guess an initial range'
          STOP 'range'

      END IF

      F1 = FUNC(X1)
      F2 = FUNC(X2)
      SUCCES = .TRUE.
      DO 11 J = 1,NTRY
c
c     AG: avoid overflow
c
          IF (sign(1.d0,F1)*sign(1.d0,F2) .LT. 0.D0) RETURN
          IF (ABS(F1).LT.ABS(F2)) THEN
              X1 = X1 + FACTOR* (X1-X2)
              F1 = FUNC(X1)

          ELSE
              X2 = X2 + FACTOR* (X2-X1)
              F2 = FUNC(X2)
          END IF

   11 CONTINUE
      SUCCES = .FALSE.
      RETURN

      END
      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      INTEGER N,NB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XB1(1),XB2(1)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FX
      EXTERNAL FX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DX,FC,FP,X
      INTEGER I,NBB
C     ..
      NBB = NB
      NB = 0
      X = X1
      DX = (X2-X1)/N
      FP = FX(X)
      DO 11 I = 1,N
          X = X + DX
          FC = FX(X)
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FC)*sign(1.d0,FP) .LT. 0.D0) THEN
              NB = NB + 1
              XB1(NB) = X - DX
              XB2(NB) = X
          END IF

          FP = FC
          IF (NBB.EQ.NB) RETURN
   11 CONTINUE
      RETURN

      END
      DOUBLE PRECISION FUNCTION ZBRENT(FUNC,X1,X2,TOL)
C
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-8)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL,X1,X2
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      INTEGER ITER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      A = X1
      B = X2
      FA = FUNC(A)
      FB = FUNC(B)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,FB)*sign(1.d0,FA) .GT. 0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for ZBRENT.'
          STOP 'range'
      END IF

      FC = FB
      DO 11 ITER = 1,ITMAX
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FB)*sign(1.d0,FC) .GT. 0.D0) THEN
CAG          IF (FB*FC.GT.0.D0) THEN
              C = A
              FC = FA
              D = B - A
              E = D
          END IF

          IF (ABS(FC).LT.ABS(FB)) THEN
              A = B
              B = C
              C = A
              FA = FB
              FB = FC
              FC = FA
          END IF

          TOL1 = 2.D0*EPS*ABS(B) + 0.5D0*TOL
          XM = .5D0* (C-B)
          IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0) THEN
              ZBRENT = B
              RETURN

          END IF

          IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S = FB/FA
              IF (A.EQ.C) THEN
                  P = 2.D0*XM*S
                  Q = 1.D0 - S

              ELSE
                  Q = FA/FC
                  R = FB/FC
                  P = S* (2.D0*XM*Q* (Q-R)- (B-A)* (R-1.D0))
                  Q = (Q-1.D0)* (R-1.D0)* (S-1.D0)
              END IF

              IF (P.GT.0.D0) Q = -Q
              P = ABS(P)
              IF (2.D0*P.LT.MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E = D
                  D = P/Q

              ELSE
                  D = XM
                  E = D
              END IF

          ELSE
              D = XM
              E = D
          END IF

          A = B
          FA = FB
          IF (ABS(D).GT.TOL1) THEN
              B = B + D

          ELSE
              B = B + SIGN(TOL1,XM)
          END IF

          FB = FUNC(B)
   11 CONTINUE
      WRITE (*,FMT=*) 'ZBRENT exceeding maximum iterations.'
      ZBRENT = B
      STOP 'ITER'

      END

      DOUBLE PRECISION FUNCTION RTBIS(FUNC,X1,X2,XACC)
C
C
C     .. Parameters ..
      INTEGER JMAX
      PARAMETER (JMAX=40)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2,XACC
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DX,F,FMID,XMID
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      FMID = FUNC(X2)
      F = FUNC(X1)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,F)*sign(1.d0,FMID) .GE. 0.D0) THEN
CAG      IF (F*FMID.GE.0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for bisection.'
          STOP 'RTBIS'

      END IF

      IF (F.LT.0.D0) THEN
          RTBIS = X1
          DX = X2 - X1

      ELSE
          RTBIS = X2
          DX = X1 - X2
      END IF

      DO 11 J = 1,JMAX
          DX = DX*.5D0
          XMID = RTBIS + DX
          FMID = FUNC(XMID)
          IF (FMID.LE.0.D0) RTBIS = XMID
          IF (ABS(DX).LT.XACC .OR. FMID.EQ.0.D0) RETURN
   11 CONTINUE
      WRITE (*,FMT=*) 'too many bisections'
      STOP 'ITER'
      END
      SUBROUTINE BRAC(FUNC,X1,X2,SUCCES)
C
c     This is ZBRAC frm Numerical Recipes.
c     Changed name to brac to avoid recursion...
c
C
C     .. Parameters ..
      INTEGER NTRY
      PARAMETER (NTRY=50)
      DOUBLE PRECISION FACTOR
      PARAMETER (FACTOR=1.6D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      LOGICAL SUCCES
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F1,F2
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (X1.EQ.X2) THEN
          WRITE (*,FMT=*) 'You have to guess an initial range'
          STOP 'range'

      END IF

      F1 = FUNC(X1)
      F2 = FUNC(X2)
      SUCCES = .TRUE.
      DO 11 J = 1,NTRY
c
c     AG: avoid overflow
c
          IF (sign(1.d0,F1)*sign(1.d0,F2) .LT. 0.D0) RETURN
CAG          IF (F1*F2.LT.0.D0) RETURN
          IF (ABS(F1).LT.ABS(F2)) THEN
              X1 = X1 + FACTOR* (X1-X2)
              F1 = FUNC(X1)

          ELSE
              X2 = X2 + FACTOR* (X2-X1)
              F2 = FUNC(X2)
          END IF

   11 CONTINUE
      SUCCES = .FALSE.
      RETURN

      END
C
      DOUBLE PRECISION FUNCTION BRENT(FUNC,X1,X2,TOL)
C
c     This is ZBRENT from Numerical Recipes.
c     Changed name to avoid recursion...
c
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-8)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL,X1,X2
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      INTEGER ITER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      A = X1
      B = X2
      FA = FUNC(A)
      FB = FUNC(B)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,FB)*sign(1.d0,FA) .GT. 0.D0) THEN
CAG      IF (FB*FA.GT.0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for BRENT.'
          STOP 'range'

      END IF

      FC = FB
      DO 11 ITER = 1,ITMAX
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FB)*sign(1.d0,FC) .GT. 0.D0) THEN
CAG          IF (FB*FC.GT.0.D0) THEN
              C = A
              FC = FA
              D = B - A
              E = D
          END IF

          IF (ABS(FC).LT.ABS(FB)) THEN
              A = B
              B = C
              C = A
              FA = FB
              FB = FC
              FC = FA
          END IF

          TOL1 = 2.D0*EPS*ABS(B) + 0.5D0*TOL
          XM = .5D0* (C-B)
          IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0) THEN
              BRENT = B
              RETURN

          END IF

          IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S = FB/FA
              IF (A.EQ.C) THEN
                  P = 2.D0*XM*S
                  Q = 1.D0 - S

              ELSE
                  Q = FA/FC
                  R = FB/FC
                  P = S* (2.D0*XM*Q* (Q-R)- (B-A)* (R-1.D0))
                  Q = (Q-1.D0)* (R-1.D0)* (S-1.D0)
              END IF

              IF (P.GT.0.D0) Q = -Q
              P = ABS(P)
              IF (2.D0*P.LT.MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E = D
                  D = P/Q

              ELSE
                  D = XM
                  E = D
              END IF

          ELSE
              D = XM
              E = D
          END IF

          A = B
          FA = FB
          IF (ABS(D).GT.TOL1) THEN
              B = B + D

          ELSE
              B = B + SIGN(TOL1,XM)
          END IF

          FB = FUNC(B)
   11 CONTINUE
      WRITE (*,FMT=*) 'BRENT exceeding maximum iterations.'
      BRENT = B
      RETURN

      END
