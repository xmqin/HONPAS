c
      double precision function divdif(f,x,n,z,m)
c
      implicit double precision (a-h,o-z)
c
C     CERN LIBRARY PROGRAM NO E-105.
C     REVISED VERSION JULY 1973.
C     PURPOSE = TO INTERPOLATE IN TABLE OF GIVEN FUNCTION VALUES WHICH
C               ARE STORED AFTER INCREASING OR DECREASING VALUES OF THE
C               ARGUMENTS.NEWTONS GENERAL INTERPOLATION FORMULA IS USED.
C     PARAMETERS ( IN LIST ).
C     F       = THE ARRAY OF THE GIVEN FUNCTION VALUES.F(K)=F(X(,)).
C     X       = THE ARRAY OF GIVEN ARGUMENTS.
C     N       = DIMENSION OF THE ARRAYS F AND X,I.D.THE NUMBER OF POINTS
C               TABLE.
C     Z       = ARGUMENT FOR WHICH THE INTERPOLATION IS WANTED.
C     M       = ORDER OF INTERPOLATION.
C
C     PARAMETERS ( IN COMMON BLOCK / DIVCOF / ).
C     MM      = THE NUMBER OF ELEMENTS STORED IN THE FOLLOWING ARRAYS
C               (MM=M+1).
C     ARG     = AN ARRAY USED FOR STORING THE ARGUMENTS USED IN THE IN-
C               TERPOLATION.
C     VAL     = AN ARRAY USED FOR STORING THE FUNCTION VALUES USED IN
C               THE INTERPOLATION.
C     COF     = AN ARRAY USED FOR STORING THE COEFFICIENTS IN NEWTONS
C               INTERPOLATION FORMULA.
C
      integer n, m
      integer i, ipoint, il, iu, jl, ju, i1, i2, ii, j
      integer index, jndex
      double precision z, zero
      double precision f(n), x(n)
      integer mm, mmax
c
      common /divcof/ arg(11), val(11), cof(11)
      common /divint/ mm
      data zero, mmax/0.0D0, 10/
c
C
C     INTERNAL PARAMETER.
C     MMAX    = THE MAXIMUM ORDER OF INTERPOLATION PERMITTED.THE DIMEN-
C               SIONS OF THE ARRAYS ARG , VAL AND COF IN THE COMMON
C               BLOCK / DIVCOF / SHOULD BE MMAX+1.
C
   10 continue
      if ((z-x(1))*(x(n)-z) .ge. zero) go to 30
C     Z-VALUE OUTSIDE RANGE,PRINT ERROR MESSAGE.
   20 continue
      write(*,9000) z
      divdif = zero
c
      return
C
   30 continue
      if ((m.le.(n-1)) .and. (m.le.mmax)) go to 40
      mm = m
      if (m .gt. (n-1)) m = n - 1
      if (m .gt. mmax) m = mmax
C     REQUIRED ORDER OF INTERPOLATION TOO HIGH.PRINT ERROR MESSAGE AND
C     REDUCE ORDER.
      write(*,9010) mm, m
C
C     START ACTUAL CALCULATION.
C     COMPUTE POINTER,IPOINT,FOR THE LEFT BOUNDARY OF THE INTERVAL IN
C     WHICH WE HAVE Z I.D. Z IN THE INTERVAL X(IPOINT),X(IPOINT+1).
   40 continue
      cof1 = z - x(1)
      do 50 i = 2, n
         ipoint = i - 1
         cof2 = z - x(i)
         if (cof1*cof2 .le. zero) go to 60
         cof1 = cof2
   50 continue
C     CONSTRUCT TABLE TO BE USED IN THE INTERPOLATION.
   60 continue
      il = ipoint
      iu = il + 1
      jl = 1
      ju = 1
      mm = m + 1
      do 80 i = 1, mm
         i1 = 1
         i2 = 1
         if ((jl.eq.0) .or. (ju.eq.0)) go to 70
         cof1 = dabs(z-x(il))
         cof2 = dabs(x(iu)-z)
         if (cof1 .gt. cof2) i1 = 0
         if (i1 .eq. 1) i2 = 0
   70    continue
         if ((jl.eq.0) .or. (i1.eq.0)) ii = iu
         if ((ju.eq.0) .or. (i2.eq.0)) ii = il
         arg(i) = x(ii)
         cof(i) = f(ii)
         val(i) = f(ii)
         if ((jl.eq.1) .and. (i1.eq.1)) il = il - 1
         if ((ju.eq.1) .and. (i2.eq.1)) iu = iu + 1
         if (il .lt. 1) jl = 0
         if (iu .gt. n) ju = 0
   80 continue
C
      do 100 i = 1, m
         do 90 j = i, m
            index = m + 1 + i - j
            jndex = index - i
            cof(index) = (cof(index)-cof(index-1))/
     1                   (arg(index)-arg(jndex))
   90    continue
  100 continue
C
      sum = cof(m+1)
      do 110 i = 1, m
         index = m + 1 - i
         sum = (z-arg(index))*sum + cof(index)
  110 continue
C
      divdif = sum
c
      return
C
 9000 format(//5x,'*** ERROR MESSAGE FUNCTION DIVDIF *** , ARGUMENT Z ='
     1      ,e21.14,' OUTSIDE RANGE.',//)
 9010 format(//5x,
     1'*** ERROR MESSAGE FUNCTION DIVDIF *** , ORDER OF INTERPOLATION M
     2=',i3,' TOO HIGH,M IS REDUCED TO ',i2,//)
C     OPPO
      end
