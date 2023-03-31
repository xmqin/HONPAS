! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE IDIAG (NN,MOLD,MNEW,MLEFT,MRIGHT,MAUX)

      use sys, only: die
      use parallel, only: ionode

C GIVEN A SQUARE INTEGER MATRIX MOLD, FINDS A DIAGONAL MATRIX MNEW,
C AND TWO MATRICES MLEFT AND MRIGHT OF DETERMINANT ONE, SUCH THAT
C MNEW = MLEFT * MOLD * MRIGHT
C Written by J.Moreno and J.M.Soler

      integer :: niter, ibig, i, j, nn, n, iter, mmin, imin
      integer :: jmin

      INTEGER MOLD(NN,NN),MNEW(NN,NN),MLEFT(NN,NN),MRIGHT(NN,NN),
     .        MAUX(NN,NN,2),ISUM
      PARAMETER (NITER=50,IBIG=9999999)
      ISUM = 0
      DO 20 J=1,NN
        DO 10 I=1,NN
          ISUM = ISUM + ABS(MOLD(I,J))
          MNEW(I,J)=MOLD(I,J)
          MLEFT(I,J)=0
          MRIGHT(I,J)=0
          MAUX(I,J,1)=0
   10   CONTINUE
        MLEFT(J,J)=1
        MRIGHT(J,J)=1
        MAUX(J,J,1)=1
   20 CONTINUE
      IF (ISUM.EQ.0) RETURN
      DO 60 N=NN,2,-1
        DO 50 ITER=1,NITER
          MMIN=IBIG
          DO J=1,N
          DO 30 I=1,N
            IF ((I.NE.N.AND.J.NE.N) .OR. (I.EQ.N.AND.J.EQ.N)) GOTO 30
            IF (MNEW(I,J).EQ.0 .OR. ABS(MNEW(I,J)).GE.MMIN) GOTO 30
               IMIN=I
               JMIN=J
               MMIN=ABS(MNEW(I,J))
   30     CONTINUE
          ENDDO
          IF (MMIN.EQ.IBIG) GOTO 60
          I=MIN(IMIN,JMIN)
          MAUX(I,I,1)=0
          MAUX(N,N,1)=0
          MAUX(I,N,1)=SIGN(1,MNEW(IMIN,JMIN))
          MAUX(N,I,1)=SIGN(1,MNEW(IMIN,JMIN))
          IF (IMIN.LT.JMIN) THEN
            CALL IMXM (MAUX,MNEW,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MAUX,MLEFT,NN,NN,NN,MLEFT,MAUX(1,1,2))
          ELSE
            CALL IMXM (MNEW,MAUX,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MRIGHT,MAUX,NN,NN,NN,MRIGHT,MAUX(1,1,2))
          ENDIF
          MAUX(I,N,1)=0
          MAUX(N,I,1)=0
          MAUX(I,I,1)=1
          MAUX(N,N,1)=1
          DO 40 I=1,N-1
            MAUX(I,N,1)=-(MNEW(I,N)/MNEW(N,N))
            CALL IMXM (MAUX,MNEW,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MAUX,MLEFT,NN,NN,NN,MLEFT,MAUX(1,1,2))
            MAUX(I,N,1)=0
            MAUX(N,I,1)=-(MNEW(N,I)/MNEW(N,N))
            CALL IMXM (MNEW,MAUX,NN,NN,NN,MNEW,MAUX(1,1,2))
            CALL IMXM (MRIGHT,MAUX,NN,NN,NN,MRIGHT,MAUX(1,1,2))
            MAUX(N,I,1)=0
   40     CONTINUE
   50   CONTINUE
          if (IOnode)
     $       WRITE(6,*)'IDIAG: ERROR. ITERATION HAS NOT CONVERGED. N=',N
          call die('error in idiag')
   60 CONTINUE
      END



      SUBROUTINE IMXM (MA,MB,N1,N2,N3,MC,MAUX)

C MULTIPLIES TWO INTEGER MATRICES. ARGUMENTS:
C   MA(N1,N2),MB(N2,N3) : INPUT MATRICES OF DIMENSIONS AS INDICATED
C   MC(N1,N3) : OUTPUT PRODUCT MATRIX. MC MAY BE THE SAME AS MA OR MB
C              (OR BOTH) FOR 'IN PLACE' MULTIPLICATION
C   MAUX : AUXILIARY ARRAY OF MINIMUM SIZE N1*N3. IF MC IS DIFFERENT
C          FROM BOTH MA AND MB, YOU CAN MAKE MAUX=MC
C WRITTEN BY JOSE SOLER. 22/5/90

      integer :: n1, n2, n3, i1, i2, i3, i, j
      INTEGER MA(N1,N2),MB(N2,N3),MC(N1,N3),MAUX(N1,N3)
      DO I3=1,N3
        DO I1=1,N1
          MAUX(I1,I3)=0
          DO I2=1,N2
            MAUX(I1,I3)=MAUX(I1,I3)+MA(I1,I2)*MB(I2,I3)
          ENDDO
        ENDDO
      ENDDO
      DO I=N3,1,-1
        DO J=N1,1,-1
          MC(J,I)=MAUX(J,I)
        ENDDO
      ENDDO
      END
