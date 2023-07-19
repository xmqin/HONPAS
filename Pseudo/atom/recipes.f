      SUBROUTINE HUNT(XX,N,X,JLO)
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER JLO,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XX(N)
C     ..
C     .. Local Scalars ..
      INTEGER INC,JHI,JM
      LOGICAL ASCND
C     ..
      ASCND = XX(N) .GT. XX(1)
      IF (JLO.LE.0 .OR. JLO.GT.N) THEN
          JLO = 0
          JHI = N + 1
          GO TO 3

      END IF

      INC = 1
      IF (X.GE.XX(JLO) .EQV. ASCND) THEN
    1     JHI = JLO + INC
          IF (JHI.GT.N) THEN
              JHI = N + 1

          ELSE IF (X.GE.XX(JHI) .EQV. ASCND) THEN
              JLO = JHI
              INC = INC + INC
              GO TO 1

          END IF

      ELSE
          JHI = JLO
    2     JLO = JHI - INC
          IF (JLO.LT.1) THEN
              JLO = 0

          ELSE IF (X.LT.XX(JLO) .EQV. ASCND) THEN
              JHI = JLO
              INC = INC + INC
              GO TO 2

          END IF

      END IF

    3 IF (JHI-JLO.EQ.1) RETURN
      JM = (JHI+JLO)/2
      IF (X.GT.XX(JM) .EQV. ASCND) THEN
          JLO = JM

      ELSE
          JHI = JM
      END IF

      GO TO 3

      END
      SUBROUTINE LOCATE(XX,N,X,J)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER J,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XX(N)
C     ..
C     .. Local Scalars ..
      INTEGER JL,JM,JU
C     ..
      JL = 0
      JU = N + 1
   10 IF (JU-JL.GT.1) THEN
          JM = (JU+JL)/2
          IF ((XX(N).GT.XX(1)) .EQV. (X.GT.XX(JM))) THEN
              JL = JM

          ELSE
              JU = JM
          END IF

          GO TO 10

      END IF

      J = JL
      RETURN

      END
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,
     +                  RKQC)
C
      include 'ode_path.h'
C     ..
C     .. Parameters ..
      INTEGER MAXSTP
      PARAMETER (MAXSTP=10000)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1.D-30)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS,H1,HMIN,X1,X2
      INTEGER NBAD,NOK,NVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YSTART(NVAR)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS,RKQC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H,HDID,HNEXT,X,XSAV
      INTEGER I,NSTP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYDX(NMAX),Y(NMAX),YSCAL(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
C     ..
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO 11 I = 1,NVAR
          Y(I) = YSTART(I)
   11 CONTINUE
      XSAV = X - DXSAV*TWO
      DO 16 NSTP = 1,MAXSTP
          CALL DERIVS(X,Y,DYDX)
          DO 12 I = 1,NVAR
              YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
   12     CONTINUE
          IF (KMAX.GT.0) THEN
              IF (ABS(X-XSAV).GT.ABS(DXSAV)) THEN
                  IF (KOUNT.LT.KMAX-1) THEN
                      KOUNT = KOUNT + 1
                      XP(KOUNT) = X
                      DO 13 I = 1,NVAR
                          YP(I,KOUNT) = Y(I)
   13                 CONTINUE
                      XSAV = X
                  END IF

              END IF

          END IF

          IF ((X+H-X2)* (X+H-X1).GT.ZERO) H = X2 - X
          CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
          IF (HDID.EQ.H) THEN
              NOK = NOK + 1

          ELSE
              NBAD = NBAD + 1
          END IF

          IF ((X-X2)* (X2-X1).GE.ZERO) THEN
              DO 14 I = 1,NVAR
                  YSTART(I) = Y(I)
   14         CONTINUE
              IF (KMAX.NE.0) THEN
                  KOUNT = KOUNT + 1
                  XP(KOUNT) = X
                  DO 15 I = 1,NVAR
                      YP(I,KOUNT) = Y(I)
   15             CONTINUE
              END IF

              RETURN

          END IF

          IF (ABS(HNEXT).LT.HMIN) THEN
              WRITE (*,FMT=*) 'ODEINT - Stepsize smaller than minimum.'
              RETURN

          END IF

          H = HNEXT
   16 CONTINUE
      WRITE (*,FMT=*) 'ODEINT - Too many steps.'
      RETURN

      END
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DY,X,Y
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XA(N),YA(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DEN,DIF,DIFT,HO,HP,W
      INTEGER I,M,NS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(NMAX),D(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      NS = 1
      DIF = ABS(X-XA(1))
      DO 11 I = 1,N
          DIFT = ABS(X-XA(I))
          IF (DIFT.LT.DIF) THEN
              NS = I
              DIF = DIFT
          END IF

          C(I) = YA(I)
          D(I) = YA(I)
   11 CONTINUE
      Y = YA(NS)
      NS = NS - 1
      DO 13 M = 1,N - 1
          DO 12 I = 1,N - M
              HO = XA(I) - X
              HP = XA(I+M) - X
              W = C(I+1) - D(I)
              DEN = HO - HP
              IF (DEN.EQ.0.D0) THEN
                  WRITE (*,FMT=*) 'POLINT - DEN=0.0'
                  RETURN

              END IF

              DEN = W/DEN
              D(I) = HP*DEN
              C(I) = HO*DEN
   12     CONTINUE
          IF (2*NS.LT.N-M) THEN
              DY = C(NS+1)

          ELSE
              DY = D(NS)
              NS = NS - 1
          END IF

          Y = Y + DY
   13 CONTINUE
      RETURN

      END
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
      DOUBLE PRECISION FCOR
      PARAMETER (FCOR=.0666666667D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0)
      DOUBLE PRECISION SAFETY
      PARAMETER (SAFETY=0.9D0)
      DOUBLE PRECISION ERRCON
      PARAMETER (ERRCON=6.D-4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYDX(N),Y(N),YSCAL(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERRMAX,H,HH,PGROW,PSHRNK,XSAV
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYSAV(NMAX),YSAV(NMAX),YTEMP(NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL RK4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      PGROW = -0.20D0
      PSHRNK = -0.25D0
      XSAV = X
      DO 11 I = 1,N
          YSAV(I) = Y(I)
          DYSAV(I) = DYDX(I)
   11 CONTINUE
      H = HTRY
    1 HH = 0.5D0*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF (X.EQ.XSAV) THEN
          WRITE (*,FMT=*) 'Stepsize not significant in RKQC.'
          RETURN

      END IF

      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO 12 I = 1,N
          YTEMP(I) = Y(I) - YTEMP(I)
          ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
   12 CONTINUE
      ERRMAX = ERRMAX/EPS
      IF (ERRMAX.GT.ONE) THEN
          H = SAFETY*H* (ERRMAX**PSHRNK)
          GO TO 1

      ELSE
          HDID = H
          IF (ERRMAX.GT.ERRCON) THEN
              HNEXT = SAFETY*H* (ERRMAX**PGROW)

          ELSE
              HNEXT = 4.D0*H
          END IF

      END IF

      DO 13 I = 1,N
          Y(I) = Y(I) + YTEMP(I)*FCOR
   13 CONTINUE
      RETURN

      END
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
C
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION H,X
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYDX(N),Y(N),YOUT(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H6,HH,XH
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYM(NMAX),DYT(NMAX),YT(NMAX)
C     ..
      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO 11 I = 1,N
          YT(I) = Y(I) + HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I = 1,N
          YT(I) = Y(I) + HH*DYT(I)
   12 CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I = 1,N
          YT(I) = Y(I) + H*DYM(I)
          DYM(I) = DYT(I) + DYM(I)
   13 CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I = 1,N
          YOUT(I) = Y(I) + H6* (DYDX(I)+DYT(I)+2.D0*DYM(I))
   14 CONTINUE
      RETURN

      END
