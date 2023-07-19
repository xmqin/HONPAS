C
c $Id: splift.f,v 1.2 1997/05/22 17:32:32 wdpgaara Exp $
c
c $Log: splift.f,v $
c Revision 1.2  1997/05/22 17:32:32  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine splift(x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)
C
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
C
C     .. Parameters ..
      double precision four
      parameter (four=4.D0)
C     ..
C     .. Scalar Arguments ..
      double precision a1, an, b1, bn
      integer ierr, isx, n
C     ..
C     .. Array Arguments ..
      double precision w(n,3), x(n), y(n), yp(n), ypp(n)
C     ..
C     .. Local Scalars ..
      double precision dnew, dold
      integer i, j, nm1, nm2
C     ..
      if (n .lt. 4) then
         ierr = 2
c
         return
c
      end if
      nm1 = n - 1
      nm2 = n - 2
      if (isx .gt. 0) go to 40
      do 10 i = 2, n
         if (x(i)-x(i-1) .le. 0) then
            ierr = 3
c
            return
c
         end if
   10 continue
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
      w(1,3) = x(2) - x(1)
      do 20 i = 2, nm1
         w(i,2) = w(i-1,3)
         w(i,3) = x(i+1) - x(i)
         w(i,1) = 2*(w(i,2)+w(i,3))
   20 continue
      w(1,1) = four
      w(1,3) = -4*a1
      w(n,1) = four
      w(n,2) = -4*an
C
C     L U DECOMPOSITION
C
      do 30 i = 2, n
         w(i-1,3) = w(i-1,3)/w(i-1,1)
         w(i,1) = w(i,1) - w(i,2)*w(i-1,3)
   30 continue
C
C     DEFINE *CONSTANT* VECTOR
C
   40 continue
      ypp(1) = 4*b1
      dold = (y(2)-y(1))/w(2,2)
      do 50 i = 2, nm2
         dnew = (y(i+1)-y(i))/w(i+1,2)
         ypp(i) = 6*(dnew-dold)
         yp(i) = dold
         dold = dnew
   50 continue
      dnew = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(nm1) = 6*(dnew-dold)
      ypp(n) = 4*bn
      yp(nm1) = dold
      yp(n) = dnew
C
C     FORWARD SUBSTITUTION
C
      ypp(1) = ypp(1)/w(1,1)
      do 60 i = 2, n
         ypp(i) = (ypp(i)-w(i,2)*ypp(i-1))/w(i,1)
   60 continue
C
C     BACKWARD SUBSTITUTION
C
      do 70 j = 1, nm1
         i = n - j
         ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
   70 continue
C
C     COMPUTE FIRST DERIVATIVES
C
      yp(1) = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2*ypp(1)+ypp(2))/6
      do 80 i = 2, nm1
         yp(i) = yp(i) + w(i,2)*(ypp(i-1)+2*ypp(i))/6
   80 continue
      yp(n) = yp(n) + (x(n)-x(nm1))*(ypp(nm1)+2*ypp(n))/6
C
      ierr = 1
c
      return
c
      end
