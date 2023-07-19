! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      DOUBLE PRECISION FUNCTION DOT(A,B,N)
C
C     RETURNS REAL*8 SCALAR PRODUCT OF TWO REAL*8 VECTORS
C
      integer :: n, i
      DOUBLE PRECISION A(N),B(N),SUM
      SUM=0.D0
      DO I=1,N
         SUM = SUM + A(I) * B(I)
      ENDDO
      DOT=SUM
      RETURN
      END
