! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE ROTATE (N,X,ALPHA,BETA,GAMMA)

C ROTATES N VECTORS X(3,N) ACCORDING TO THE EULER ANGLES:
C   ALPHA: FIRST ROTATION AROUND AXIS Z
C   BETA : SECOND ROTATION AROUND Y
C   GAMMA: THIRD ROTATION AROUND Z AGAIN
C WRITTEN BY J.SOLER.  AUG/91

*     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3,N),RMAT(3,3),XNEW(3)
      ZERO=0.D0
      PI=ACOS(-1.D0)
      CONS=PI/180.D0
      SA=SIN(-ALPHA*CONS)
      SB=SIN(-BETA *CONS)
      SG=SIN(-GAMMA*CONS)
      CA=COS(-ALPHA*CONS)
      CB=COS(-BETA *CONS)
      CG=COS(-GAMMA*CONS)
      RMAT(1,1)= CB*CA*CG-SA*SG
      RMAT(1,2)= CB*SA*CG+CA*SG
      RMAT(1,3)=-SB*CG
      RMAT(2,1)=-CB*CA*SG-SA*CG
      RMAT(2,2)=-CB*SA*SG+CA*CG
      RMAT(2,3)= SB*SG
      RMAT(3,1)= SB*CA
      RMAT(3,2)= SB*SA
      RMAT(3,3)= CB
      DO 40 I=1,N
        DO 20 J=1,3
          XNEW(J)=ZERO
          DO 10 K=1,3
            XNEW(J)=XNEW(J)+RMAT(J,K)*X(K,I)
   10     CONTINUE
   20   CONTINUE
        DO 30 J=1,3
          X(J,I)=XNEW(J)
   30   CONTINUE
   40 CONTINUE
      END

