!
C ******************************************************************
C *  HERMDP.FOR  SOLVES EIGENVALUE PROBLEM                          *
C *                                       H X = E X                 *
C *  WHERE H IS A HERMITIAN MATRIX, X ARE THE EIGENVECTORS AND E    *
C *  THE EIGENVALUES.                                               *
C *                             OTTO F. SANKEY                      *
C *                             DEPT. OF PHYSICS                    *
C *                             ARIZONA STATE UNIVERSITY            *
C *                             602 965-4334                        *
C *                                                                 *
C *             DOUBLE PRECISION   5/15/85                          *
c *             updated            2/10/89 ofs.                     *
C *******************************************************************
C
C TO USE:
C               SUBROUTINE HERMDP(H,ZR,ZI,E,IDIM,M,W1,W2,icall)
C H     = INPUT HAMILTONIAN. 
c         it contains real and imaginary parts of h as follows:
c         h(i,j) = real (h(i,j))  , [ i.ge.j, put Real(lower triangle)
c                                             in lower triangle].
c         h(j,i) = imag (h(i,j))  , [ i.ge.j, put imag(lower triangle)
c                                             in lower triangle].
c        
C  IDIM = DIMENSION OF H EXACTLY AS IN CALLING ROUTINE.
C     M = DIMENSION OF MATRIX YOU WANT DIAGONALIZED (M.LE.IDIM)
C ZR,ZI = OUTPUT EIGENVECTORS (R MEANS REAL PART AND I IMAG. PART).
c         zr(i,j) is the real part of the i'th component of the j'th
c         eigenvector -- zr(component, eigenvector).
C     E = OUTPUT EIGENVALUES
C  icall= 3 for eigenvectors AND eigenvalues, icall.ne.3 for eigenvalues
c                                             only.
C W1, W2 ARE WORK VECTORS, W1(IDIM),W2(2,IDIM)
c aid is a work vector dimensioned 501. for bigger matrices, please
c                                       redimension.
c
c program computes trace & compares with sum of eigenvalues as a check.
C
C CALL SUBROUTINES: DPHTRI,DPIMTQ,DPHTBK
c ==============================================================
        SUBROUTINE hERMDP(H,ZR,ZI,E,IDIM,M,W1,W2,icall)
c ==============================================================
c        DOUBLE PRECISION H(IDIM,M),E(M),AID(501)
        integer :: idim, m, icall
        DOUBLE PRECISION H(IDIM,M),E(M),AID(5501),w3(5501)
        DOUBLE PRECISION ZR(IDIM,M),ZI(IDIM,M),W1(M),W2(2,M)
        DOUBLE PRECISION EIGSUM,TRACE
c ==============================================================
        COMMON /AIDIAG/AID
c     ==============================================================
        integer :: i, j, ier
        
c        IF(M.GT.501)WRITE(*,*)' ERROR IN HERMDP ******* M>501'
c        IF(M.GT.501)GO TO 999
c      write(*,*)' welcome to hermdp.f ......'
        IF(M.GT.5501)WRITE(*,*)' ERROR IN HERMDP ******* M>5501'
        IF(M.GT.5501)GO TO 999
C FIND TRACE OF MATRIX, to test sum of eigenvalues.
       TRACE=0.D0
       DO 10 I=1,M
10     TRACE=TRACE+H(I,I)
C INITIALIZE ZR,ZI. must set zr to identity.
        DO 1 I=1,M
        DO 2 J=1,M
2       ZR(I,J)=0.D0
1       ZR(I,I)=1.D0
c ==============================================================
c first call.
        CALL DPHTRI(IDIM,M,H,E,W1,w3,W2)
       IF(M.LT.300)GO TO 3
       WRITE(*,*)' FIRST CALL FINISHED'
3      continue 
c ==============================================================
c second call.
       CALL DPIMTQ(IDIM,M,E,W1,ZR,IER)
       IF(M.LT.300)GO TO 4
       WRITE(*,*)' SECOND CALL FINISHED'
4      continue 
c ==============================================================
c third call: for eigenvectors.
       if(icall.eq.3)then
           CALL DPHTBK(IDIM,M,H,W2,M,ZR,ZI)
           IF(M.GT.300)WRITE(*,*)' THIRD CALL FINISHED'
       else 
           write(*,*)' **caution** in hermdp, eigenvcs not desired'
       end if
       IF(M.GT.300.and.icall.eq.3)WRITE(*,*)' THIRD CALL FINISHED'
c ==============================================================
C SUM EIGENVALUES AND COMPARE WITH TRACE.
       EIGSUM=0.D0
       DO 11 I=1,M
11     EIGSUM=EIGSUM+E(I)
       IF(100.D0*DABS(EIGSUM-TRACE).LT.DABS(TRACE)*1.D-9)GOTO 12
       WRITE(*,*)' ERROR: hERMDP SUM OF EIGENVALUES NOT EQUAL TO TRACE'
       WRITE(*,*)' EIGSUM=',EIGSUM,' TRACE=',TRACE
       IF(DABS(EIGSUM-TRACE).GT.DABS(TRACE)*1.D-08)
     X  WRITE(*,*)' BAD ERROR IN hERMDP ******* !!!!!! ******'
12      IF(IER.NE.0) WRITE(*,*)' ERROR IN HERMDP IER IN HERMDP=',IER
c ==============================================================
        RETURN
999     STOP
        END
c ==============================================================
C
      SUBROUTINE DPHTRI(NM,N,AR,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION AR(NM,N),D(N),E(N),E2(N),TAU(2,N)
c      DOUBLE PRECISION F,FI,G,GI,H,HH,SI,SCALE,AID(501)
      DOUBLE PRECISION F,FI,G,GI,H,HH,SI,SCALE,AID(5501)
      DOUBLE PRECISION DSQRT,CDABS,DABS
      double precision ajunk
       COMMON /AIDIAG/AID
C      COMPLEX*16 DCMPLX
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX
C     TO A REAL SYMMETRIC TRIDIAGONAL MATRIX USING
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        AR AND AI CONTAIN THE REAL AND AIMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX HERMITIAN INPUT MATRIX.
C          ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT-
C
C       AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION IN THEIR FULL LOWER
C          TRIANGLES.  THEIR STRICT UPPER TRIANGLES AND THE
C          DIAGONAL OF AR ARE UNALTERED,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED,
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE USE OF THE SUBROUTINES
C     CABS AND CMPLX IN COMPUTING COMPLEX ABSOLUTE VALUES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0
c       DO 324 I=1,501
       DO 324 I=1,5501
324      AID(I)=0.D0
C
C
      DO 100 I = 1, N
  100 D(I) = AR(I,I)
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(AR(I,K)) + DABS(AR(K,I))
C
         IF (SCALE .NE. 0.0) GO TO 140
         TAU(1,L) = 1.0D0
         TAU(2,L) = 0.0D0
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AR(K,I) = AR(K,I) / SCALE
            H = H + AR(I,K) * AR(I,K) + AR(K,I) * AR(K,I)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = DSQRT(H)
         E(I) = SCALE * G
         F=DSQRT(AR(I,L)**2+AR(L,I)**2)
C         F = CDABS(DCMPLX(AR(I,L),AI(I,L)))
C     ********** FORM NEXT DIAGONAL ELEMENT OF MATRIX T **********
         IF (F .EQ. 0.0D0) GO TO 160
         TAU(1,L) = (AR(L,I) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AR(L,I) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0D0 + G / F
         AR(I,L) = G * AR(I,L)
         AR(L,I) = G * AR(L,I)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
            GI = 0.0D0
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, J
       AJUNK=AR(K,J)
       IF(K.EQ.J)AJUNK=AID(K)
               G = G + AR(J,K) * AR(I,K) + AJUNK * AR(K,I)
               GI = GI - AR(J,K) * AR(K,I) + AJUNK * AR(I,K)
  180       CONTINUE
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AR(J,K) * AR(K,I)
               GI = GI - AR(K,J) * AR(K,I)- AR(J,K) * AR(I,K)
  200       CONTINUE
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AR(J,I)
  240    CONTINUE
C
         HH = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AR(J,I)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
C
            DO 260 K = 1, J
       AJUNK=AR(K,J)
       IF(J.EQ.K)AJUNK=AID(K)
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     X                           + FI * TAU(2,K) + GI * AR(K,I)
       IF(K.EQ.J)GO TO 17
               AR(K,J) = AR(K,J) - F * TAU(2,K) - G * AR(K,I)
     X                           - FI * E(K) - GI * AR(I,K)
       GO TO 260
17               AID(K)=AID(K) - F * TAU(2,K) - G * AR(K,I)
     X                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AR(K,I) = SCALE * AR(K,I)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AID(I) = SCALE*DSQRT(H)
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF DPHTRI **********
      END
c ==============================================================
      SUBROUTINE DPIMTQ(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,F,G,P,R,S,MACHEP
      DOUBLE PRECISION DSQRT,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE DPIMTQ,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C     MACHEP=2**(-26) SP.    2**(-56) DP. (26 AND 56 ARE MAXIMUMS!)
C     WE USE 54 FOR SAFETY. (1.D0+MACHEP.GT.1.D0) DEFINES MINIMUM MACHEP.
C                **********
      MACHEP=2.D0**(-54)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (DABS(E(M)) .LE. MACHEP * (DABS(D(M)) + DABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = DSQRT(G*G+1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (DABS(F)  .LT. DABS(G))  GO TO 150
            C = G / F
            R =DSQRT(C*C+1.0D0)
            E(I+1) = F * R
            S = 1.0D0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = DSQRT(S*S+1.0D0)
            E(I+1) = G * R
            C = 1.0D0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     ********** FORM VECTOR **********
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
C
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C               EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF DPIMTQ **********
      END
c ==============================================================
      SUBROUTINE DPHTBK(NM,N,AR,TAU,M,ZR,ZI)
C
      INTEGER I,J,K,L,M,N,NM
      DOUBLE PRECISION AR(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
c      DOUBLE PRECISION H,S,SI,AID(501)
      DOUBLE PRECISION H,S,SI,AID(5501)
       COMMON /AIDIAG/AID
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  DPHTRI.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  DPHTRI  IN THEIR
C          FULL LOWER TRIANGLES EXCEPT FOR THE DIAGONAL OF AR.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        ZR AND ZI CONTAIN THE REAL AND AIMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
        IF(M .EQ. 0 ) GOTO 200
C     ********** TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. **********
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = - ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     ********** RECOVER AND APPLY THE HOUSEHOLDER MATRICES **********
      DO 140 I = 2, N
         L = I - 1
         H = AID(I)
         IF (H .EQ. 0.0D0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0D0
            SI = 0.0D0
C
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AR(K,I) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AR(K,I) * ZR(K,J)
  110       CONTINUE
C
        S=(S/H)/H
        SI=(SI/H)/H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AR(K,I)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AR(K,I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
C     ********** LAST CARD OF DPHTBK **********
      END
