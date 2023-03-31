!Copyright (c) 1992-2017 The University of Tennessee and The University
!                        of Tennessee Research Foundation.  All rights
!                        reserved.
!Copyright (c) 2000-2017 The University of California Berkeley. All
!                        rights reserved.
!Copyright (c) 2006-2017 The University of Colorado Denver.  All rights
!                        reserved.
!
!$COPYRIGHT$
!
!Additional copyrights may follow
!
!$HEADER$
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!- Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!- Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer listed
!  in this license in the documentation and/or other materials
!  provided with the distribution.
!
!- Neither the name of the copyright holders nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!The copyright holders provide no reassurances that the source code
!provided does not infringe any patent, copyright, or any other
!intellectual property rights of third parties.  The copyright holders
!disclaim any liability to any recipient for claims brought against
!recipient by any third party for infringement of that parties
!intellectual property rights.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/caxpy.f
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
      COMPLEX CA
      INTEGER INCX,INCY,N
      COMPLEX CX(*),CY(*)
      INTEGER I,IX,IY
      REAL SCABS1
      EXTERNAL SCABS1
      IF (N.LE.0) RETURN
      IF (SCABS1(CA).EQ.0.0E+0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CY(I) = CY(I) + CA*CX(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CY(IY) = CY(IY) + CA*CX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ccopy.f
      SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX CX(*),CY(*)
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CY(I) = CX(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CY(IY) = CX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cdotc.f
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX CX(*),CY(*)
      COMPLEX CTEMP
      INTEGER I,IX,IY
      INTRINSIC CONJG
      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(I))*CY(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTC = CTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cdotu.f
      COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX CX(*),CY(*)
      COMPLEX CTEMP
      INTEGER I,IX,IY
      CTEMP = (0.0,0.0)
      CDOTU = (0.0,0.0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CTEMP = CTEMP + CX(I)*CY(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CX(IX)*CY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      CDOTU = CTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cgemv.f
      SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      COMPLEX A(LDA,*),X(*),Y(*)
      COMPLEX ONE
      PARAMETER (ONE= (1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
      COMPLEX TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC CONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('CGEMV ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      NOCONJ = LSAME(TRANS,'T')
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + CONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + CONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cgeru.f
      SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      COMPLEX ALPHA
      INTEGER INCX,INCY,LDA,M,N
      COMPLEX A(LDA,*),X(*),Y(*)
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
      COMPLEX TEMP
      INTEGER I,INFO,IX,J,JY,KX
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('CGERU ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cscal.f
      SUBROUTINE CSCAL(N,CA,CX,INCX)
      COMPLEX CA
      INTEGER INCX,N
      COMPLEX CX(*)
      INTEGER I,NINCX
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            CX(I) = CA*CX(I)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            CX(I) = CA*CX(I)
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/csscal.f
      SUBROUTINE CSSCAL(N,SA,CX,INCX)
      REAL SA
      INTEGER INCX,N
      COMPLEX CX(*)
      INTEGER I,NINCX
      INTRINSIC AIMAG,CMPLX,REAL
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/cswap.f
      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX CX(*),CY(*)
      COMPLEX CTEMP
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            CTEMP = CX(I)
            CX(I) = CY(I)
            CY(I) = CTEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CX(IX)
            CX(IX) = CY(IY)
            CY(IY) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ctbsv.f
      SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX A(LDA,*),X(*)
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
      COMPLEX TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC CONJG,MAX,MIN
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('CTBSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      IF (NOCONJ) THEN
                          DO 90 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      ELSE
                          DO 100 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(KPLUS1,J))
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      IF (NOCONJ) THEN
                          DO 120 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      ELSE
                          DO 130 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(KPLUS1,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      IF (NOCONJ) THEN
                          DO 150 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      ELSE
                          DO 160 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      IF (NOCONJ) THEN
                          DO 180 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      ELSE
                          DO 190 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ctrsm.f
      SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      COMPLEX A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC CONJG,MAX
      COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
      COMPLEX ONE
      PARAMETER (ONE= (1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO= (0.0E+0,0.0E+0))
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('CTRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
  110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - CONJG(A(K,I))*B(K,J)
  120                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/CONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  130                 CONTINUE
  140             CONTINUE
              ELSE
                  DO 180 J = 1,N
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
  150                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - CONJG(A(K,I))*B(K,J)
  160                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/CONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  170                 CONTINUE
  180             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 230 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 190 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  190                     CONTINUE
                      END IF
                      DO 210 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 220 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  220                     CONTINUE
                      END IF
  230             CONTINUE
              ELSE
                  DO 280 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 240 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  240                     CONTINUE
                      END IF
                      DO 260 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 270 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 330 K = N,1,-1
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/CONJG(A(K,K))
                          END IF
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
                      DO 310 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = CONJG(A(J,K))
                              END IF
                              DO 300 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  300                         CONTINUE
                          END IF
  310                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 320 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  320                     CONTINUE
                      END IF
  330             CONTINUE
              ELSE
                  DO 380 K = 1,N
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/CONJG(A(K,K))
                          END IF
                          DO 340 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  340                     CONTINUE
                      END IF
                      DO 360 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = CONJG(A(J,K))
                              END IF
                              DO 350 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  350                         CONTINUE
                          END IF
  360                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 370 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  370                     CONTINUE
                      END IF
  380             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dasum.f
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
      INTRINSIC DABS,MOD
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) +
     $              DABS(DX(I+2)) + DABS(DX(I+3)) +
     $              DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM = DTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/daxpy.f
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dcabs1.f
      DOUBLE PRECISION FUNCTION DCABS1(Z)
      COMPLEX*16 Z
      INTRINSIC ABS,DBLE,DIMAG
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dcopy.f
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ddot.f
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dgemm.f
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (NOTB) THEN
          IF (NOTA) THEN
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dgemv.f
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dger.f
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGER  ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dnrm2.f
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
      INTEGER INCX,N
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
      INTRINSIC ABS,SQRT
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
      DNRM2 = NORM
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/drot.f
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dscal.f
      SUBROUTINE DSCAL(N,DA,DX,INCX)
      DOUBLE PRECISION DA
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      INTEGER I,M,MP1,NINCX
      INTRINSIC MOD
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsdot.f
      DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      INTEGER I,KX,KY,NS
      INTRINSIC DBLE
      DSDOT = 0.0D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dspmv.f
      SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      DOUBLE PRECISION AP(*),X(*),Y(*)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 6
      ELSE IF (INCY.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSPMV ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  K = KK
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  KK = KK + J
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 K = KK,KK + J - 2
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   80         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*AP(KK)
                  K = KK + 1
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(I)
                      K = K + 1
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
                  KK = KK + (N-J+1)
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*AP(KK)
                  IX = JX
                  IY = JY
                  DO 110 K = KK + 1,KK + N - J
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + AP(K)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + (N-J+1)
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dspr.f
      SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,N
      CHARACTER UPLO
      DOUBLE PRECISION AP(*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSPR  ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      K = KK
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                  END IF
                  KK = KK + J
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + J
   40         CONTINUE
          END IF
      ELSE
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      K = KK
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  END IF
                  KK = KK + N - J + 1
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dspr2.f
      SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      DOUBLE PRECISION AP(*),X(*),Y(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSPR2 ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      DO 10 I = 1,J
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   10                 CONTINUE
                  END IF
                  KK = KK + J
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      DO 30 K = KK,KK + J - 1
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   40         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      K = KK
                      DO 50 I = J,N
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   50                 CONTINUE
                  END IF
                  KK = KK + N - J + 1
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      DO 70 K = KK,KK + N - J
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dswap.f
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsymm.f
      SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,K,NROWA
      LOGICAL UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      IF (LSAME(SIDE,'L')) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (LSAME(SIDE,'L')) THEN
          IF (UPPER) THEN
              DO 70 J = 1,N
                  DO 60 I = 1,M
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   50                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) +
     +                             ALPHA*TEMP2
                      END IF
   60             CONTINUE
   70         CONTINUE
          ELSE
              DO 100 J = 1,N
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   80                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) +
     +                             ALPHA*TEMP2
                      END IF
   90             CONTINUE
  100         CONTINUE
          END IF
      ELSE
          DO 170 J = 1,N
              TEMP1 = ALPHA*A(J,J)
              IF (BETA.EQ.ZERO) THEN
                  DO 110 I = 1,M
                      C(I,J) = TEMP1*B(I,J)
  110             CONTINUE
              ELSE
                  DO 120 I = 1,M
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  120             CONTINUE
              END IF
              DO 140 K = 1,J - 1
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(K,J)
                  ELSE
                      TEMP1 = ALPHA*A(J,K)
                  END IF
                  DO 130 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  130             CONTINUE
  140         CONTINUE
              DO 160 K = J + 1,N
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(J,K)
                  ELSE
                      TEMP1 = ALPHA*A(K,J)
                  END IF
                  DO 150 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsymv.f
      SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYMV ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*A(J,J)
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*A(J,J)
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + A(I,J)*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsyr.f
      SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
      DOUBLE PRECISION A(LDA,*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR  ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = KX
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*X(JX)
                      IX = JX
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsyr2.f
      SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2 ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 10 I = 1,J
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                 CONTINUE
                  END IF
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = KX
                      IY = KY
                      DO 30 I = 1,J
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   40         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(J)
                      TEMP2 = ALPHA*X(J)
                      DO 50 I = J,N
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                 CONTINUE
                  END IF
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*Y(JY)
                      TEMP2 = ALPHA*X(JX)
                      IX = JX
                      IY = JY
                      DO 70 I = J,N
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   70                 CONTINUE
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsyr2k.f
      SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYR2K',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dsyrk.f
      SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
      DOUBLE PRECISION A(LDA,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DSYRK ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP = ZERO
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtbsv.f
      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      DOUBLE PRECISION A(LDA,*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX,MIN
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTBSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      DO 90 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      DO 110 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      DO 130 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      DO 150 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtpmv.f
      SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      INTEGER INCX,N
      CHARACTER DIAG,TRANS,UPLO
      DOUBLE PRECISION AP(*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTPMV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          K = KK
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*AP(K)
                              K = K + 1
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*AP(KK+J-1)
                      END IF
                      KK = KK + J
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 K = KK,KK + J - 2
                              X(IX) = X(IX) + TEMP*AP(K)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*AP(KK+J-1)
                      END IF
                      JX = JX + INCX
                      KK = KK + J
   40             CONTINUE
              END IF
          ELSE
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          K = KK
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*AP(K)
                              K = K - 1
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*AP(KK-N+J)
                      END IF
                      KK = KK - (N-J+1)
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 K = KK,KK - (N- (J+1)),-1
                              X(IX) = X(IX) + TEMP*AP(K)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*AP(KK-N+J)
                      END IF
                      JX = JX - INCX
                      KK = KK - (N-J+1)
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*AP(KK)
                      K = KK - 1
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + AP(K)*X(I)
                          K = K - 1
   90                 CONTINUE
                      X(J) = TEMP
                      KK = KK - J
  100             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*AP(KK)
                      DO 110 K = KK - 1,KK - J + 1,-1
                          IX = IX - INCX
                          TEMP = TEMP + AP(K)*X(IX)
  110                 CONTINUE
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - J
  120             CONTINUE
              END IF
          ELSE
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 140 J = 1,N
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*AP(KK)
                      K = KK + 1
                      DO 130 I = J + 1,N
                          TEMP = TEMP + AP(K)*X(I)
                          K = K + 1
  130                 CONTINUE
                      X(J) = TEMP
                      KK = KK + (N-J+1)
  140             CONTINUE
              ELSE
                  JX = KX
                  DO 160 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*AP(KK)
                      DO 150 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          TEMP = TEMP + AP(K)*X(IX)
  150                 CONTINUE
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + (N-J+1)
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtpsv.f
      SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      INTEGER INCX,N
      CHARACTER DIAG,TRANS,UPLO
      DOUBLE PRECISION AP(*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTPSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK - 1
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*AP(K)
                              K = K - 1
   10                     CONTINUE
                      END IF
                      KK = KK - J
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*AP(K)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
                      KK = KK - J
   40             CONTINUE
              END IF
          ELSE
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK + 1
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*AP(K)
                              K = K + 1
   50                     CONTINUE
                      END IF
                      KK = KK + (N-J+1)
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 K = KK + 1,KK + N - J
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*AP(K)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
                      KK = KK + (N-J+1)
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      K = KK
                      DO 90 I = 1,J - 1
                          TEMP = TEMP - AP(K)*X(I)
                          K = K + 1
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      X(J) = TEMP
                      KK = KK + J
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      DO 110 K = KK,KK + J - 2
                          TEMP = TEMP - AP(K)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + J
  120             CONTINUE
              END IF
          ELSE
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      K = KK
                      DO 130 I = N,J + 1,-1
                          TEMP = TEMP - AP(K)*X(I)
                          K = K - 1
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      X(J) = TEMP
                      KK = KK - (N-J+1)
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      DO 150 K = KK,KK - (N- (J+1)),-1
                          TEMP = TEMP - AP(K)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - (N-J+1)
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtrmm.f
      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 110 J = 1,N
                      DO 100 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 90 K = 1,I - 1
                              TEMP = TEMP + A(K,I)*B(K,J)
   90                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  100                 CONTINUE
  110             CONTINUE
              ELSE
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = B(I,J)
                          IF (NOUNIT) TEMP = TEMP*A(I,I)
                          DO 120 K = I + 1,M
                              TEMP = TEMP + A(K,I)*B(K,J)
  120                     CONTINUE
                          B(I,J) = ALPHA*TEMP
  130                 CONTINUE
  140             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 180 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  150                 CONTINUE
                      DO 170 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 160 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  160                         CONTINUE
                          END IF
  170                 CONTINUE
  180             CONTINUE
              ELSE
                  DO 220 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 190 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  190                 CONTINUE
                      DO 210 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
  220             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 260 K = 1,N
                      DO 240 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 250 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              ELSE
                  DO 300 K = N,1,-1
                      DO 280 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = ALPHA*A(J,K)
                              DO 270 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  270                         CONTINUE
                          END IF
  280                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      IF (TEMP.NE.ONE) THEN
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
  300             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtrmv.f
      SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      DOUBLE PRECISION A(LDA,*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRMV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = N,1,-1
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 90 I = J - 1,1,-1
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 110 I = J - 1,1,-1
                          IX = IX - INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  110                 CONTINUE
                      X(JX) = TEMP
                      JX = JX - INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = 1,N
                      TEMP = X(J)
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 130 I = J + 1,N
                          TEMP = TEMP + A(I,J)*X(I)
  130                 CONTINUE
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  JX = KX
                  DO 160 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 150 I = J + 1,N
                          IX = IX + INCX
                          TEMP = TEMP + A(I,J)*X(IX)
  150                 CONTINUE
                      X(JX) = TEMP
                      JX = JX + INCX
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtrsm.f
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dtrsv.f
      SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      DOUBLE PRECISION A(LDA,*),X(*)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTRSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      DO 90 I = 1,J - 1
                          TEMP = TEMP - A(I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      DO 110 I = 1,J - 1
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(JX) = TEMP
                      JX = JX + INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      DO 130 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      DO 150 I = N,J + 1,-1
                          TEMP = TEMP - A(I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(J,J)
                      X(JX) = TEMP
                      JX = JX - INCX
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dzasum.f
      DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)
      INTEGER INCX,N
      COMPLEX*16 ZX(*)
      DOUBLE PRECISION STEMP
      INTEGER I,NINCX
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
      DZASUM = 0.0d0
      STEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      END IF
      DZASUM = STEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/dznrm2.f
      DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
      INTEGER INCX,N
      COMPLEX*16 X(*)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
      INTEGER IX
      INTRINSIC ABS,DBLE,DIMAG,SQRT
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (DBLE(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DBLE(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (DIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
      DZNRM2 = NORM
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/icamax.f
      INTEGER FUNCTION ICAMAX(N,CX,INCX)
      INTEGER INCX,N
      COMPLEX CX(*)
      REAL SMAX
      INTEGER I,IX
      REAL SCABS1
      EXTERNAL SCABS1
      ICAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
         SMAX = SCABS1(CX(1))
         DO I = 2,N
            IF (SCABS1(CX(I)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(I))
            END IF
         END DO
      ELSE
         IX = 1
         SMAX = SCABS1(CX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (SCABS1(CX(IX)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/idamax.f
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION DMAX
      INTEGER I,IX
      INTRINSIC DABS
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/isamax.f
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
      INTEGER INCX,N
      REAL SX(*)
      REAL SMAX
      INTEGER I,IX
      INTRINSIC ABS
      ISAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ISAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
         SMAX = ABS(SX(1))
         DO I = 2,N
            IF (ABS(SX(I)).GT.SMAX) THEN
               ISAMAX = I
               SMAX = ABS(SX(I))
            END IF
         END DO
      ELSE
         IX = 1
         SMAX = ABS(SX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (ABS(SX(IX)).GT.SMAX) THEN
               ISAMAX = I
               SMAX = ABS(SX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/izamax.f
      INTEGER FUNCTION IZAMAX(N,ZX,INCX)
      INTEGER INCX,N
      COMPLEX*16 ZX(*)
      DOUBLE PRECISION DMAX
      INTEGER I,IX
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
      IZAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
         DMAX = DCABS1(ZX(1))
         DO I = 2,N
            IF (DCABS1(ZX(I)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(I))
            END IF
         END DO
      ELSE
         IX = 1
         DMAX = DCABS1(ZX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DCABS1(ZX(IX)).GT.DMAX) THEN
               IZAMAX = I
               DMAX = DCABS1(ZX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/lsame.f
      LOGICAL FUNCTION LSAME(CA,CB)
      CHARACTER CA,CB
      INTRINSIC ICHAR
      INTEGER INTA,INTB,ZCODE
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
      ZCODE = ICHAR('Z')
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sasum.f
      REAL FUNCTION SASUM(N,SX,INCX)
      INTEGER INCX,N
      REAL SX(*)
      REAL STEMP
      INTEGER I,M,MP1,NINCX
      INTRINSIC ABS,MOD
      SASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = STEMP + ABS(SX(I))
            END DO
            IF (N.LT.6) THEN
               SASUM = STEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            STEMP = STEMP + ABS(SX(I)) + ABS(SX(I+1)) +
     $              ABS(SX(I+2)) + ABS(SX(I+3)) +
     $              ABS(SX(I+4)) + ABS(SX(I+5))
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(SX(I))
         END DO
      END IF
      SASUM = STEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/scabs1.f
      REAL FUNCTION SCABS1(Z)
      COMPLEX Z
      INTRINSIC ABS,AIMAG,REAL
      SCABS1 = ABS(REAL(Z)) + ABS(AIMAG(Z))
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/scasum.f
      REAL FUNCTION SCASUM(N,CX,INCX)
      INTEGER INCX,N
      COMPLEX CX(*)
      REAL STEMP
      INTEGER I,NINCX
      INTRINSIC ABS,AIMAG,REAL
      SCASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      END IF
      SCASUM = STEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/scnrm2.f
      REAL FUNCTION SCNRM2(N,X,INCX)
      INTEGER INCX,N
      COMPLEX X(*)
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL NORM,SCALE,SSQ,TEMP
      INTEGER IX
      INTRINSIC ABS,AIMAG,REAL,SQRT
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (REAL(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(REAL(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (AIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(AIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
      SCNRM2 = NORM
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/scopy.f
      SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               SY(I) = SX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            SY(I) = SX(I)
            SY(I+1) = SX(I+1)
            SY(I+2) = SX(I+2)
            SY(I+3) = SX(I+3)
            SY(I+4) = SX(I+4)
            SY(I+5) = SX(I+5)
            SY(I+6) = SX(I+6)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            SY(IY) = SX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sdot.f
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      STEMP = 0.0e0
      SDOT = 0.0e0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = STEMP + SX(I)*SY(I)
            END DO
            IF (N.LT.5) THEN
               SDOT=STEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) +
     $            SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            STEMP = STEMP + SX(IX)*SY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      SDOT = STEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sdsdot.f
      REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
      REAL SB
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      DOUBLE PRECISION DSDOT
      INTEGER I,KX,KY,NS
      INTRINSIC DBLE
      DSDOT = SB
      IF (N.LE.0) THEN
         SDSDOT = DSDOT
         RETURN
      END IF
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
      SDSDOT = DSDOT
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sgemv.f
      SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      REAL ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      REAL A(LDA,*),X(*),Y(*)
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('SGEMV ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sger.f
      SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      REAL ALPHA
      INTEGER INCX,INCY,LDA,M,N
      REAL A(LDA,*),X(*),Y(*)
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,J,JY,KX
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('SGER  ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/snrm2.f
      REAL FUNCTION SNRM2(N,X,INCX)
      INTEGER INCX,N
      REAL X(*)
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL ABSXI,NORM,SCALE,SSQ
      INTEGER IX
      INTRINSIC ABS,SQRT
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
      SNRM2 = NORM
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sscal.f
      SUBROUTINE SSCAL(N,SA,SX,INCX)
      REAL SA
      INTEGER INCX,N
      REAL SX(*)
      INTEGER I,M,MP1,NINCX
      INTRINSIC MOD
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               SX(I) = SA*SX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            SX(I) = SA*SX(I)
            SX(I+1) = SA*SX(I+1)
            SX(I+2) = SA*SX(I+2)
            SX(I+3) = SA*SX(I+3)
            SX(I+4) = SA*SX(I+4)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            SX(I) = SA*SX(I)
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/sswap.f
      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
      INTEGER INCX,INCY,N
      REAL SX(*),SY(*)
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
      INTRINSIC MOD
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = SX(I)
               SX(I) = SY(I)
               SY(I) = STEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            STEMP = SX(I)
            SX(I) = SY(I)
            SY(I) = STEMP
            STEMP = SX(I+1)
            SX(I+1) = SY(I+1)
            SY(I+1) = STEMP
            STEMP = SX(I+2)
            SX(I+2) = SY(I+2)
            SY(I+2) = STEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            STEMP = SX(IX)
            SX(IX) = SY(IY)
            SY(IY) = STEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/stbsv.f
      SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      REAL A(LDA,*),X(*)
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      REAL TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX,MIN
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('STBSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      DO 90 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(J) = TEMP
  100             CONTINUE
              ELSE
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      DO 110 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      DO 130 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(J) = TEMP
  140             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      DO 150 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/strsm.f
      SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      REAL ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      REAL A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      REAL TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('STRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/xerbla.f
      SUBROUTINE XERBLA( SRNAME, INFO )
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
      INTRINSIC          LEN_TRIM
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
      STOP
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zaxpy.f
      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
      COMPLEX*16 ZX(*),ZY(*)
      INTEGER I,IX,IY
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
      IF (N.LE.0) RETURN
      IF (DCABS1(ZA).EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            ZY(I) = ZY(I) + ZA*ZX(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZY(IY) = ZY(IY) + ZA*ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zcopy.f
      SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX*16 ZX(*),ZY(*)
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
          ZY(I) = ZX(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZY(IY) = ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zdotc.f
      COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX*16 ZX(*),ZY(*)
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
      INTRINSIC DCONJG
      ZTEMP = (0.0d0,0.0d0)
      ZDOTC = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTC = ZTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zdotu.f
      COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX*16 ZX(*),ZY(*)
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
      ZTEMP = (0.0d0,0.0d0)
      ZDOTU = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(I)*ZY(I)
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      ZDOTU = ZTEMP
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zdrot.f
      SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C, S
      COMPLEX*16         CX( * ), CY( * )
      INTEGER            I, IX, IY
      COMPLEX*16         CTEMP
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
         DO I = 1, N
            CTEMP = C*CX( I ) + S*CY( I )
            CY( I ) = C*CY( I ) - S*CX( I )
            CX( I ) = CTEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF( INCX.LT.0 )
     $      IX = ( -N+1 )*INCX + 1
         IF( INCY.LT.0 )
     $      IY = ( -N+1 )*INCY + 1
         DO I = 1, N
            CTEMP = C*CX( IX ) + S*CY( IY )
            CY( IY ) = C*CY( IY ) - S*CX( IX )
            CX( IX ) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zdscal.f
      SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
      DOUBLE PRECISION DA
      INTEGER INCX,N
      COMPLEX*16 ZX(*)
      INTEGER I,NINCX
      INTRINSIC DCMPLX
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zgemm.f
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
          NROWA = M
          NCOLA = K
      ELSE
          NROWA = K
          NCOLA = M
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (NOTB) THEN
          IF (NOTA) THEN
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE IF (CONJA) THEN
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
              DO 150 J = 1,N
                  DO 140 I = 1,M
                      TEMP = ZERO
                      DO 130 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  130                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (NOTA) THEN
          IF (CONJB) THEN
              DO 200 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 160 I = 1,M
                          C(I,J) = ZERO
  160                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 170 I = 1,M
                          C(I,J) = BETA*C(I,J)
  170                 CONTINUE
                  END IF
                  DO 190 L = 1,K
                      TEMP = ALPHA*DCONJG(B(J,L))
                      DO 180 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          ELSE
              DO 250 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 210 I = 1,M
                          C(I,J) = ZERO
  210                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 220 I = 1,M
                          C(I,J) = BETA*C(I,J)
  220                 CONTINUE
                  END IF
                  DO 240 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 230 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  230                 CONTINUE
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (CONJA) THEN
          IF (CONJB) THEN
              DO 280 J = 1,N
                  DO 270 I = 1,M
                      TEMP = ZERO
                      DO 260 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  260                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (CONJB) THEN
              DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
              DO 370 J = 1,N
                  DO 360 I = 1,M
                      TEMP = ZERO
                      DO 350 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  350                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zgemv.f
      SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMV ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      NOCONJ = LSAME(TRANS,'T')
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 110 J = 1,N
                  TEMP = ZERO
                  IF (NOCONJ) THEN
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  ELSE
                      DO 100 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          ELSE
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  IF (NOCONJ) THEN
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  ELSE
                      DO 130 I = 1,M
                          TEMP = TEMP + DCONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  END IF
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zgerc.f
      SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
      COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERC ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*DCONJG(Y(JY))
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*DCONJG(Y(JY))
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zgeru.f
      SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
      COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
      EXTERNAL XERBLA
      INTRINSIC MAX
      INFO = 0
      IF (M.LT.0) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGERU ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF (INCY.GT.0) THEN
          JY = 1
      ELSE
          JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zhemm.f
      SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      COMPLEX*16 ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG,MAX
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,J,K,NROWA
      LOGICAL UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      IF (LSAME(SIDE,'L')) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHEMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (LSAME(SIDE,'L')) THEN
          IF (UPPER) THEN
              DO 70 J = 1,N
                  DO 60 I = 1,M
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*DCONJG(A(K,I))
   50                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*DBLE(A(I,I)) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*DBLE(A(I,I)) +
     +                             ALPHA*TEMP2
                      END IF
   60             CONTINUE
   70         CONTINUE
          ELSE
              DO 100 J = 1,N
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*DCONJG(A(K,I))
   80                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*DBLE(A(I,I)) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*DBLE(A(I,I)) +
     +                             ALPHA*TEMP2
                      END IF
   90             CONTINUE
  100         CONTINUE
          END IF
      ELSE
          DO 170 J = 1,N
              TEMP1 = ALPHA*DBLE(A(J,J))
              IF (BETA.EQ.ZERO) THEN
                  DO 110 I = 1,M
                      C(I,J) = TEMP1*B(I,J)
  110             CONTINUE
              ELSE
                  DO 120 I = 1,M
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  120             CONTINUE
              END IF
              DO 140 K = 1,J - 1
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(K,J)
                  ELSE
                      TEMP1 = ALPHA*DCONJG(A(J,K))
                  END IF
                  DO 130 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  130             CONTINUE
  140         CONTINUE
              DO 160 K = J + 1,N
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*DCONJG(A(J,K))
                  ELSE
                      TEMP1 = ALPHA*A(K,J)
                  END IF
                  DO 150 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zhemv.f
      SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 5
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      ELSE IF (INCY.EQ.0) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHEMV ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + DCONJG(A(I,J))*X(I)
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*DBLE(A(J,J)) + ALPHA*TEMP2
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 I = 1,J - 1
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + DCONJG(A(I,J))*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*DBLE(A(J,J)) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*DBLE(A(J,J))
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + DCONJG(A(I,J))*X(I)
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*DBLE(A(J,J))
                  IX = JX
                  IY = JY
                  DO 110 I = J + 1,N
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*A(I,J)
                      TEMP2 = TEMP2 + DCONJG(A(I,J))*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zher.f
      SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,LDA,N
      CHARACTER UPLO
      COMPLEX*16 A(LDA,*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHER  ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.DBLE(ZERO))) RETURN
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      DO 10 I = 1,J - 1
                          A(I,J) = A(I,J) + X(I)*TEMP
   10                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(J)*TEMP)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 I = 1,J - 1
                          A(I,J) = A(I,J) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) + DBLE(X(JX)*TEMP)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
   40         CONTINUE
          END IF
      ELSE
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(J))
                      DO 50 I = J + 1,N
                          A(I,J) = A(I,J) + X(I)*TEMP
   50                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      A(J,J) = DBLE(A(J,J)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 I = J + 1,N
                          IX = IX + INCX
                          A(I,J) = A(I,J) + X(IX)*TEMP
   70                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zher2.f
      SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      COMPLEX*16 A(LDA,*),X(*),Y(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHER2 ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      DO 10 I = 1,J - 1
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) +
     +                         DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      IX = KX
                      IY = KY
                      DO 30 I = 1,J - 1
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                      A(J,J) = DBLE(A(J,J)) +
     +                         DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   40         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      A(J,J) = DBLE(A(J,J)) +
     +                         DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                      DO 50 I = J + 1,N
                          A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      A(J,J) = DBLE(A(J,J)) +
     +                         DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                      IX = JX
                      IY = JY
                      DO 70 I = J + 1,N
                          IX = IX + INCX
                          IY = IY + INCY
                          A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
   70                 CONTINUE
                  ELSE
                      A(J,J) = DBLE(A(J,J))
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zher2k.f
      SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      COMPLEX*16 ALPHA
      DOUBLE PRECISION BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG,MAX
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHER2K',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.DBLE(ZERO)) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.DBLE(ZERO)) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*DCONJG(B(J,L))
                          TEMP2 = DCONJG(ALPHA*A(J,L))
                          DO 110 I = 1,J - 1
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                          C(J,J) = DBLE(C(J,J)) +
     +                             DBLE(A(J,L)*TEMP1+B(J,L)*TEMP2)
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.DBLE(ZERO)) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*DCONJG(B(J,L))
                          TEMP2 = DCONJG(ALPHA*A(J,L))
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                          C(J,J) = DBLE(C(J,J)) +
     +                             DBLE(A(J,L)*TEMP1+B(J,L)*TEMP2)
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                          TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  190                 CONTINUE
                      IF (I.EQ.J) THEN
                          IF (BETA.EQ.DBLE(ZERO)) THEN
                              C(J,J) = DBLE(ALPHA*TEMP1+
     +                                 DCONJG(ALPHA)*TEMP2)
                          ELSE
                              C(J,J) = BETA*DBLE(C(J,J)) +
     +                                 DBLE(ALPHA*TEMP1+
     +                                 DCONJG(ALPHA)*TEMP2)
                          END IF
                      ELSE
                          IF (BETA.EQ.DBLE(ZERO)) THEN
                              C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                          ELSE
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                                 DCONJG(ALPHA)*TEMP2
                          END IF
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + DCONJG(A(L,I))*B(L,J)
                          TEMP2 = TEMP2 + DCONJG(B(L,I))*A(L,J)
  220                 CONTINUE
                      IF (I.EQ.J) THEN
                          IF (BETA.EQ.DBLE(ZERO)) THEN
                              C(J,J) = DBLE(ALPHA*TEMP1+
     +                                 DCONJG(ALPHA)*TEMP2)
                          ELSE
                              C(J,J) = BETA*DBLE(C(J,J)) +
     +                                 DBLE(ALPHA*TEMP1+
     +                                 DCONJG(ALPHA)*TEMP2)
                          END IF
                      ELSE
                          IF (BETA.EQ.DBLE(ZERO)) THEN
                              C(I,J) = ALPHA*TEMP1 + DCONJG(ALPHA)*TEMP2
                          ELSE
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                                 DCONJG(ALPHA)*TEMP2
                          END IF
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zherk.f
      SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
      COMPLEX*16 A(LDA,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCMPLX,DCONJG,MAX
      COMPLEX*16 TEMP
      DOUBLE PRECISION RTEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'C'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHERK ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                      C(J,J) = BETA*DBLE(C(J,J))
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.DCMPLX(ZERO)) THEN
                          TEMP = ALPHA*DCONJG(A(J,L))
                          DO 110 I = 1,J - 1
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(I,L))
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      C(J,J) = BETA*DBLE(C(J,J))
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  ELSE
                      C(J,J) = DBLE(C(J,J))
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.DCMPLX(ZERO)) THEN
                          TEMP = ALPHA*DCONJG(A(J,L))
                          C(J,J) = DBLE(C(J,J)) + DBLE(TEMP*A(J,L))
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 220 J = 1,N
                  DO 200 I = 1,J - 1
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
                  RTEMP = ZERO
                  DO 210 L = 1,K
                      RTEMP = RTEMP + DCONJG(A(L,J))*A(L,J)
  210             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                      C(J,J) = ALPHA*RTEMP
                  ELSE
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  END IF
  220         CONTINUE
          ELSE
              DO 260 J = 1,N
                  RTEMP = ZERO
                  DO 230 L = 1,K
                      RTEMP = RTEMP + DCONJG(A(L,J))*A(L,J)
  230             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                      C(J,J) = ALPHA*RTEMP
                  ELSE
                      C(J,J) = ALPHA*RTEMP + BETA*DBLE(C(J,J))
                  END IF
                  DO 250 I = J + 1,N
                      TEMP = ZERO
                      DO 240 L = 1,K
                          TEMP = TEMP + DCONJG(A(L,I))*A(L,J)
  240                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  250             CONTINUE
  260         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zhpmv.f
      SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      COMPLEX*16 AP(*),X(*),Y(*)
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 6
      ELSE IF (INCY.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHPMV ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (N-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (N-1)*INCY
      END IF
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,N
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,N
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,N
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,N
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  K = KK
                  DO 50 I = 1,J - 1
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + DCONJG(AP(K))*X(I)
                      K = K + 1
   50             CONTINUE
                  Y(J) = Y(J) + TEMP1*DBLE(AP(KK+J-1)) + ALPHA*TEMP2
                  KK = KK + J
   60         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 80 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  IX = KX
                  IY = KY
                  DO 70 K = KK,KK + J - 2
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + DCONJG(AP(K))*X(IX)
                      IX = IX + INCX
                      IY = IY + INCY
   70             CONTINUE
                  Y(JY) = Y(JY) + TEMP1*DBLE(AP(KK+J-1)) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   80         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 100 J = 1,N
                  TEMP1 = ALPHA*X(J)
                  TEMP2 = ZERO
                  Y(J) = Y(J) + TEMP1*DBLE(AP(KK))
                  K = KK + 1
                  DO 90 I = J + 1,N
                      Y(I) = Y(I) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + DCONJG(AP(K))*X(I)
                      K = K + 1
   90             CONTINUE
                  Y(J) = Y(J) + ALPHA*TEMP2
                  KK = KK + (N-J+1)
  100         CONTINUE
          ELSE
              JX = KX
              JY = KY
              DO 120 J = 1,N
                  TEMP1 = ALPHA*X(JX)
                  TEMP2 = ZERO
                  Y(JY) = Y(JY) + TEMP1*DBLE(AP(KK))
                  IX = JX
                  IY = JY
                  DO 110 K = KK + 1,KK + N - J
                      IX = IX + INCX
                      IY = IY + INCY
                      Y(IY) = Y(IY) + TEMP1*AP(K)
                      TEMP2 = TEMP2 + DCONJG(AP(K))*X(IX)
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP2
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + (N-J+1)
  120         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zhpr.f
      SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
      DOUBLE PRECISION ALPHA
      INTEGER INCX,N
      CHARACTER UPLO
      COMPLEX*16 AP(*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHPR  ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.DBLE(ZERO))) RETURN
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF (INCX.EQ.1) THEN
              DO 20 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      K = KK
                      DO 10 I = 1,J - 1
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   10                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(J)*TEMP)
                  ELSE
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  KK = KK + J
   20         CONTINUE
          ELSE
              JX = KX
              DO 40 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      IX = KX
                      DO 30 K = KK,KK + J - 2
                          AP(K) = AP(K) + X(IX)*TEMP
                          IX = IX + INCX
   30                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) + DBLE(X(JX)*TEMP)
                  ELSE
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  JX = JX + INCX
                  KK = KK + J
   40         CONTINUE
          END IF
      ELSE
          IF (INCX.EQ.1) THEN
              DO 60 J = 1,N
                  IF (X(J).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(J))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(J))
                      K = KK + 1
                      DO 50 I = J + 1,N
                          AP(K) = AP(K) + X(I)*TEMP
                          K = K + 1
   50                 CONTINUE
                  ELSE
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  KK = KK + N - J + 1
   60         CONTINUE
          ELSE
              JX = KX
              DO 80 J = 1,N
                  IF (X(JX).NE.ZERO) THEN
                      TEMP = ALPHA*DCONJG(X(JX))
                      AP(KK) = DBLE(AP(KK)) + DBLE(TEMP*X(JX))
                      IX = JX
                      DO 70 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          AP(K) = AP(K) + X(IX)*TEMP
   70                 CONTINUE
                  ELSE
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  JX = JX + INCX
                  KK = KK + N - J + 1
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zhpr2.f
      SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      COMPLEX*16 AP(*),X(*),Y(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,IX,IY,J,JX,JY,K,KK,KX,KY
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DBLE,DCONJG
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (N.LT.0) THEN
          INFO = 2
      ELSE IF (INCX.EQ.0) THEN
          INFO = 5
      ELSE IF (INCY.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZHPR2 ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
          IF (INCY.GT.0) THEN
              KY = 1
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
          JX = KX
          JY = KY
      END IF
      KK = 1
      IF (LSAME(UPLO,'U')) THEN
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 20 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      K = KK
                      DO 10 I = 1,J - 1
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   10                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) +
     +                             DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                  ELSE
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  KK = KK + J
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      IX = KX
                      IY = KY
                      DO 30 K = KK,KK + J - 2
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                          IX = IX + INCX
                          IY = IY + INCY
   30                 CONTINUE
                      AP(KK+J-1) = DBLE(AP(KK+J-1)) +
     +                             DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                  ELSE
                      AP(KK+J-1) = DBLE(AP(KK+J-1))
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + J
   40         CONTINUE
          END IF
      ELSE
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
              DO 60 J = 1,N
                  IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(J))
                      TEMP2 = DCONJG(ALPHA*X(J))
                      AP(KK) = DBLE(AP(KK)) +
     +                         DBLE(X(J)*TEMP1+Y(J)*TEMP2)
                      K = KK + 1
                      DO 50 I = J + 1,N
                          AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                          K = K + 1
   50                 CONTINUE
                  ELSE
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  KK = KK + N - J + 1
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                      TEMP1 = ALPHA*DCONJG(Y(JY))
                      TEMP2 = DCONJG(ALPHA*X(JX))
                      AP(KK) = DBLE(AP(KK)) +
     +                         DBLE(X(JX)*TEMP1+Y(JY)*TEMP2)
                      IX = JX
                      IY = JY
                      DO 70 K = KK + 1,KK + N - J
                          IX = IX + INCX
                          IY = IY + INCY
                          AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
   70                 CONTINUE
                  ELSE
                      AP(KK) = DBLE(AP(KK))
                  END IF
                  JX = JX + INCX
                  JY = JY + INCY
                  KK = KK + N - J + 1
   80         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zscal.f
      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
      COMPLEX*16 ZA
      INTEGER INCX,N
      COMPLEX*16 ZX(*)
      INTEGER I,NINCX
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         DO I = 1,N
            ZX(I) = ZA*ZX(I)
         END DO
      ELSE
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = ZA*ZX(I)
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zswap.f
      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX*16 ZX(*),ZY(*)
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
         DO I = 1,N
            ZTEMP = ZX(I)
            ZX(I) = ZY(I)
            ZY(I) = ZTEMP
         END DO
      ELSE
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZX(IX)
            ZX(IX) = ZY(IY)
            ZY(IY) = ZTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zsymm.f
      SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      COMPLEX*16 ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,J,K,NROWA
      LOGICAL UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      IF (LSAME(SIDE,'L')) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSAME(SIDE,'L')) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF (M.LT.0) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,M)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSYMM ',INFO)
          RETURN
      END IF
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
      IF (LSAME(SIDE,'L')) THEN
          IF (UPPER) THEN
              DO 70 J = 1,N
                  DO 60 I = 1,M
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 50 K = 1,I - 1
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   50                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) +
     +                             ALPHA*TEMP2
                      END IF
   60             CONTINUE
   70         CONTINUE
          ELSE
              DO 100 J = 1,N
                  DO 90 I = M,1,-1
                      TEMP1 = ALPHA*B(I,J)
                      TEMP2 = ZERO
                      DO 80 K = I + 1,M
                          C(K,J) = C(K,J) + TEMP1*A(K,I)
                          TEMP2 = TEMP2 + B(K,J)*A(K,I)
   80                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) +
     +                             ALPHA*TEMP2
                      END IF
   90             CONTINUE
  100         CONTINUE
          END IF
      ELSE
          DO 170 J = 1,N
              TEMP1 = ALPHA*A(J,J)
              IF (BETA.EQ.ZERO) THEN
                  DO 110 I = 1,M
                      C(I,J) = TEMP1*B(I,J)
  110             CONTINUE
              ELSE
                  DO 120 I = 1,M
                      C(I,J) = BETA*C(I,J) + TEMP1*B(I,J)
  120             CONTINUE
              END IF
              DO 140 K = 1,J - 1
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(K,J)
                  ELSE
                      TEMP1 = ALPHA*A(J,K)
                  END IF
                  DO 130 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  130             CONTINUE
  140         CONTINUE
              DO 160 K = J + 1,N
                  IF (UPPER) THEN
                      TEMP1 = ALPHA*A(J,K)
                  ELSE
                      TEMP1 = ALPHA*A(K,J)
                  END IF
                  DO 150 I = 1,M
                      C(I,J) = C(I,J) + TEMP1*B(I,K)
  150             CONTINUE
  160         CONTINUE
  170     CONTINUE
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zsyr2k.f
      SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      COMPLEX*16 TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSYR2K',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 +
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 +
     +                             ALPHA*TEMP2
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/zsyrk.f
      SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO
      COMPLEX*16 A(LDA,*),C(LDC,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC MAX
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 10
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSYRK ',INFO)
          RETURN
      END IF
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (UPPER) THEN
              DO 130 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 90 I = 1,J
                          C(I,J) = ZERO
   90                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 100 I = 1,J
                          C(I,J) = BETA*C(I,J)
  100                 CONTINUE
                  END IF
                  DO 120 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 110 I = 1,J
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                  END IF
                  DO 170 L = 1,K
                      IF (A(J,L).NE.ZERO) THEN
                          TEMP = ALPHA*A(J,L)
                          DO 160 I = J,N
                              C(I,J) = C(I,J) + TEMP*A(I,L)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP = ZERO
                      DO 190 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  190                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP = ZERO
                      DO 220 L = 1,K
                          TEMP = TEMP + A(L,I)*A(L,J)
  220                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztbsv.f
      SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX*16 A(LDA,*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX,MIN
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTBSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      IF (NOCONJ) THEN
                          DO 90 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      ELSE
                          DO 100 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - DCONJG(A(L+I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(KPLUS1,J))
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      IF (NOCONJ) THEN
                          DO 120 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      ELSE
                          DO 130 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - DCONJG(A(L+I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(KPLUS1,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      IF (NOCONJ) THEN
                          DO 150 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      ELSE
                          DO 160 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - DCONJG(A(L+I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(1,J))
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      IF (NOCONJ) THEN
                          DO 180 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      ELSE
                          DO 190 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - DCONJG(A(L+I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(1,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztpmv.f
      SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      INTEGER INCX,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX*16 AP(*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTPMV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          K = KK
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*AP(K)
                              K = K + 1
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*AP(KK+J-1)
                      END IF
                      KK = KK + J
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 K = KK,KK + J - 2
                              X(IX) = X(IX) + TEMP*AP(K)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*AP(KK+J-1)
                      END IF
                      JX = JX + INCX
                      KK = KK + J
   40             CONTINUE
              END IF
          ELSE
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          K = KK
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*AP(K)
                              K = K - 1
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*AP(KK-N+J)
                      END IF
                      KK = KK - (N-J+1)
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 K = KK,KK - (N- (J+1)),-1
                              X(IX) = X(IX) + TEMP*AP(K)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*AP(KK-N+J)
                      END IF
                      JX = JX - INCX
                      KK = KK - (N-J+1)
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      K = KK - 1
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*AP(KK)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + AP(K)*X(I)
                              K = K - 1
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(AP(KK))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(AP(K))*X(I)
                              K = K - 1
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
                      KK = KK - J
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*AP(KK)
                          DO 120 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              TEMP = TEMP + AP(K)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(AP(KK))
                          DO 130 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(AP(K))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - J
  140             CONTINUE
              END IF
          ELSE
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      K = KK + 1
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*AP(KK)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + AP(K)*X(I)
                              K = K + 1
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(AP(KK))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(AP(K))*X(I)
                              K = K + 1
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
                      KK = KK + (N-J+1)
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*AP(KK)
                          DO 180 K = KK + 1,KK + N - J
                              IX = IX + INCX
                              TEMP = TEMP + AP(K)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(AP(KK))
                          DO 190 K = KK + 1,KK + N - J
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(AP(K))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + (N-J+1)
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztpsv.f
      SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      INTEGER INCX,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX*16 AP(*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,K,KK,KX
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (INCX.EQ.0) THEN
          INFO = 7
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTPSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK - 1
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*AP(K)
                              K = K - 1
   10                     CONTINUE
                      END IF
                      KK = KK - J
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*AP(K)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
                      KK = KK - J
   40             CONTINUE
              END IF
          ELSE
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK + 1
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*AP(K)
                              K = K + 1
   50                     CONTINUE
                      END IF
                      KK = KK + (N-J+1)
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 K = KK + 1,KK + N - J
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*AP(K)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
                      KK = KK + (N-J+1)
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              KK = 1
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      K = KK
                      IF (NOCONJ) THEN
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - AP(K)*X(I)
                              K = K + 1
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      ELSE
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - DCONJG(AP(K))*X(I)
                              K = K + 1
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(AP(KK+J-1))
                      END IF
                      X(J) = TEMP
                      KK = KK + J
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      IF (NOCONJ) THEN
                          DO 120 K = KK,KK + J - 2
                              TEMP = TEMP - AP(K)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      ELSE
                          DO 130 K = KK,KK + J - 2
                              TEMP = TEMP - DCONJG(AP(K))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(AP(KK+J-1))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + J
  140             CONTINUE
              END IF
          ELSE
              KK = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      K = KK
                      IF (NOCONJ) THEN
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - AP(K)*X(I)
                              K = K - 1
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      ELSE
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - DCONJG(AP(K))*X(I)
                              K = K - 1
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(AP(KK-N+J))
                      END IF
                      X(J) = TEMP
                      KK = KK - (N-J+1)
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      IF (NOCONJ) THEN
                          DO 180 K = KK,KK - (N- (J+1)),-1
                              TEMP = TEMP - AP(K)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      ELSE
                          DO 190 K = KK,KK - (N- (J+1)),-1
                              TEMP = TEMP - DCONJG(AP(K))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(AP(KK-N+J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - (N-J+1)
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztrmm.f
      SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 50 J = 1,N
                      DO 40 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              DO 30 I = 1,K - 1
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   30                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP*A(K,K)
                              B(K,J) = TEMP
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      DO 70 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*B(K,J)
                              B(K,J) = TEMP
                              IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                              DO 60 I = K + 1,M
                                  B(I,J) = B(I,J) + TEMP*A(I,K)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 120 J = 1,N
                      DO 110 I = M,1,-1
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
   90                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 100 K = 1,I - 1
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  100                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = 1,M
                          TEMP = B(I,J)
                          IF (NOCONJ) THEN
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 130 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
  130                         CONTINUE
                          ELSE
                              IF (NOUNIT) TEMP = TEMP*DCONJG(A(I,I))
                              DO 140 K = I + 1,M
                                  TEMP = TEMP + DCONJG(A(K,I))*B(K,J)
  140                         CONTINUE
                          END IF
                          B(I,J) = ALPHA*TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 200 J = N,1,-1
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 170 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  170                 CONTINUE
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 J = 1,N
                      TEMP = ALPHA
                      IF (NOUNIT) TEMP = TEMP*A(J,J)
                      DO 210 I = 1,M
                          B(I,J) = TEMP*B(I,J)
  210                 CONTINUE
                      DO 230 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              TEMP = ALPHA*A(K,J)
                              DO 220 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 280 K = 1,N
                      DO 260 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 K = N,1,-1
                      DO 300 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = ALPHA*A(J,K)
                              ELSE
                                  TEMP = ALPHA*DCONJG(A(J,K))
                              END IF
                              DO 290 I = 1,M
                                  B(I,J) = B(I,J) + TEMP*B(I,K)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      TEMP = ALPHA
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = TEMP*A(K,K)
                          ELSE
                              TEMP = TEMP*DCONJG(A(K,K))
                          END IF
                      END IF
                      IF (TEMP.NE.ONE) THEN
                          DO 310 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztrmv.f
      SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX*16 A(LDA,*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 10 I = 1,J - 1
                              X(I) = X(I) + TEMP*A(I,J)
   10                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX
                  DO 40 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 30 I = 1,J - 1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX + INCX
   30                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX + INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          TEMP = X(J)
                          DO 50 I = N,J + 1,-1
                              X(I) = X(I) + TEMP*A(I,J)
   50                     CONTINUE
                          IF (NOUNIT) X(J) = X(J)*A(J,J)
                      END IF
   60             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          TEMP = X(JX)
                          IX = KX
                          DO 70 I = N,J + 1,-1
                              X(IX) = X(IX) + TEMP*A(I,J)
                              IX = IX - INCX
   70                     CONTINUE
                          IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                      END IF
                      JX = JX - INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 90 I = J - 1,1,-1
                              TEMP = TEMP + A(I,J)*X(I)
   90                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 100 I = J - 1,1,-1
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 140 J = N,1,-1
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 120 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  120                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 130 I = J - 1,1,-1
                              IX = IX - INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  130                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = J + 1,N
                              TEMP = TEMP + A(I,J)*X(I)
  150                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 160 I = J + 1,N
                              TEMP = TEMP + DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  JX = KX
                  DO 200 J = 1,N
                      TEMP = X(JX)
                      IX = JX
                      IF (NOCONJ) THEN
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 180 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + A(I,J)*X(IX)
  180                     CONTINUE
                      ELSE
                          IF (NOUNIT) TEMP = TEMP*DCONJG(A(J,J))
                          DO 190 I = J + 1,N
                              IX = IX + INCX
                              TEMP = TEMP + DCONJG(A(I,J))*X(IX)
  190                     CONTINUE
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztrsm.f
      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
      COMPLEX*16 A(LDA,*),B(LDB,*)
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
      COMPLEX*16 ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
          NROWA = M
      ELSE
          NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
          INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
          INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
          INFO = 4
      ELSE IF (M.LT.0) THEN
          INFO = 5
      ELSE IF (N.LT.0) THEN
          INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRSM ',INFO)
          RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      IF (ALPHA.EQ.ZERO) THEN
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 140 J = 1,N
                      DO 130 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 110 K = 1,I - 1
                                  TEMP = TEMP - A(K,I)*B(K,J)
  110                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 120 K = 1,I - 1
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  120                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  130                 CONTINUE
  140             CONTINUE
              ELSE
                  DO 180 J = 1,N
                      DO 170 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          IF (NOCONJ) THEN
                              DO 150 K = I + 1,M
                                  TEMP = TEMP - A(K,I)*B(K,J)
  150                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/A(I,I)
                          ELSE
                              DO 160 K = I + 1,M
                                  TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  160                         CONTINUE
                              IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                          END IF
                          B(I,J) = TEMP
  170                 CONTINUE
  180             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
              IF (UPPER) THEN
                  DO 230 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 190 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  190                     CONTINUE
                      END IF
                      DO 210 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 200 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 220 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  220                     CONTINUE
                      END IF
  230             CONTINUE
              ELSE
                  DO 280 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 240 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  240                     CONTINUE
                      END IF
                      DO 260 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 250 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 270 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              END IF
          ELSE
              IF (UPPER) THEN
                  DO 330 K = N,1,-1
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 290 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  290                     CONTINUE
                      END IF
                      DO 310 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 300 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  300                         CONTINUE
                          END IF
  310                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 320 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  320                     CONTINUE
                      END IF
  330             CONTINUE
              ELSE
                  DO 380 K = 1,N
                      IF (NOUNIT) THEN
                          IF (NOCONJ) THEN
                              TEMP = ONE/A(K,K)
                          ELSE
                              TEMP = ONE/DCONJG(A(K,K))
                          END IF
                          DO 340 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  340                     CONTINUE
                      END IF
                      DO 360 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              IF (NOCONJ) THEN
                                  TEMP = A(J,K)
                              ELSE
                                  TEMP = DCONJG(A(J,K))
                              END IF
                              DO 350 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  350                         CONTINUE
                          END IF
  360                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 370 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  370                     CONTINUE
                      END IF
  380             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
! SOURCE-FILE = /home/nicpa/LA/lapack/BLAS/SRC/ztrsv.f
      SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
      COMPLEX*16 A(LDA,*),X(*)
      COMPLEX*16 ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA
      INTRINSIC DCONJG,MAX
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRSV ',INFO)
          RETURN
      END IF
      IF (N.EQ.0) RETURN
      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      IF (LSAME(TRANS,'N')) THEN
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      IF (X(JX).NE.ZERO) THEN
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 110 J = 1,N
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - DCONJG(A(I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  110             CONTINUE
              ELSE
                  JX = KX
                  DO 140 J = 1,N
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 120 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 130 I = 1,J - 1
                              TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX + INCX
  140             CONTINUE
              END IF
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      IF (NOCONJ) THEN
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - DCONJG(A(I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(J) = TEMP
  170             CONTINUE
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      IX = KX
                      TEMP = X(JX)
                      IF (NOCONJ) THEN
                          DO 180 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      ELSE
                          DO 190 I = N,J + 1,-1
                              TEMP = TEMP - DCONJG(A(I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/DCONJG(A(J,J))
                      END IF
                      X(JX) = TEMP
                      JX = JX - INCX
  200             CONTINUE
              END IF
          END IF
      END IF
      RETURN
      END
