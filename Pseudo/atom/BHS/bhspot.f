c     
      program bhspot
c     
c     Generates VPS-compatible file from BHS parameters
c     
      implicit none
      
      integer nrp
      parameter (nrp=1000)
      
      real*8 r(nrp), rab(nrp)
      real*8 v(nrp,0:2), vcore(nrp)
c     
      character*2 nameat
      character irel*3, nicore*4, icorr*2
      character ray(6)*10, title(7)*10
      real*8 zv, ztot, znuc, zion, ar, b, aa, bb, rmax
      
      real*8 alfcore(2), ccore(2)
      real*8 alf(3), c(6), a(6)
      
      integer i,j,k,l
      integer nr, npotd, npotu
      
      real*8 erf
      external erf
      
      include 'symbols.h'
      
c---------------------------------
      open(5,file='BHSINP',form='formatted')
      open(1,file='VPS',form='unformatted')
c     
      read(5,'(a2)') nameat
      read(5,*) ztot, zv
      
      znuc = ztot - zv
      
      rmax = set_value('RMAX','80.')
      aa = set_value('AA_R','6.')
      bb = set_value('BB_R','40.')
      
      ar = exp (-aa) / znuc
      b = 1/bb
      do 30 i = 1, nrp
         r(i) = ar*(exp(b*(i-1))-1)
         rab(i) = (r(i)+ar)*b
         if (r(i) .gt. rmax) go to 40
 30      continue
c     
 40      continue
         nr = i - 1
c
         npotd = 3
         npotu = 0
         irel = 'bhs'
         nicore = 'nc'
         icorr = 'ca'
c
         write(ray,*) 'BHS fit'
         write(title,*) 'Random'
c
      write(1) nameat, icorr, irel, nicore, (ray(i),i=1,6),
     &  (title(i),i=1,7), npotd, npotu, nr - 1, ar, b, zv
      write(1) (r(i),i=2,nr)
c
         read(5,*) (alfcore(i),i=1,2)
         read(5,*) (ccore(i),i=1,2)
c     
c Note factor of 2 to convert to rydberg units...
c
         do i=2,nr
            vcore(i) = (ccore(1)*erf(sqrt(alfcore(1))*r(i))) +
     $           (ccore(2)*erf(sqrt(alfcore(2))*r(i))) 
            vcore(i) = -2.d0*zv*vcore(i) / r(i)
         enddo
         
         do i=0,2
            read(5,*) l, (alf(j),j=1,3), (c(j),j=1,6)
            if (l.ne.i) stop 'L'
            call bhsfit(alf,c,a)
            do j=1,nr
               v(j,i) = vcore(j)
               do k = 1, 3
                  v(j,i) = v(j,i) + 
     $                 2.d0 * (a(k) + r(j)*r(j)*a(k+3))
     $                 * exp(-alf(k)*r(j)*r(j))
               enddo
               if (j.ne.1) write(8,*) r(j), v(j,i)
            enddo
            write(1) i, (r(j)*v(j,i),j=2,nr)
            call potran(i+1,v(1,i),r,nr,zv)
         enddo
c         
c        Charges... Zero for now.
c
         write(1) (0.d0,i=2,nr)
         write(1) (0.d0,i=2,nr)
c
         end
c
      subroutine bhsfit(alf,c,a)
c
c   Transformation of BHS fitting parameters as given by Warren Pickett
c   Input: Alf and C. Output: A
c
      implicit none

      real*8 alf(3), a(6)
      real*8 c(6)
c
      real*8 q(6,6), s(6,6)
c
      integer i,j,k, ir, ii, jj
      real*8 alfinv
      real*8 pi

      data pi / 3.1415926535898d0 /
c
      do i = 1, 6
         do j = 1, 6
            q(i,j) = 0.d0
            ii = i
            jj = j
            if (ii .gt. 3) ii = i-3
            if (jj .gt. 3) jj = j-3
            alfinv = 1.d0 / (alf(ii) + alf(jj))
            s(i,j) = 0.25d0*alfinv * sqrt(pi*alfinv)

            if ((i .gt. 3 .and. j.le.3) .or. (i.le.3 .and. j.gt.3) )
     +         s(i,j) = s(i,j) * 1.5d0 * alfinv

            if (i.gt.3 .and. j.gt.3)
     $           s(i,j) = s(i,j) * 3.75d0*alfinv*alfinv
         enddo
      enddo
c
      q(1,1) = sqrt(s(1,1)) 
      do j=2,6
         q(1,j) = s(1,j) / q(1,1)
      enddo

      q(2,2) = sqrt( s(2,2) - q(1,2)**2 ) 
      q(2,3) = ( s(2,3) - q(1,2)*q(1,3) ) / q(2,2)
      q(3,3) = sqrt( s(3,3) - q(1,3)**2 - q(2,3)**2 ) 
      q(2,4) = ( s(2,4) - q(1,2)*q(1,4) ) / q(2,2)
      q(3,4) = ( s(3,4) - q(1,3)*q(1,4) - q(2,3)*q(2,4) ) / q(3,3)
      q(4,4) = sqrt( s(4,4) - q(1,4)**2 - q(2,4)**2 - q(3,4)**2 ) 
      q(2,5) = ( s(2,5) - q(1,2)*q(1,5) ) / q(2,2)
      q(3,5) = ( s(3,5) - q(1,3)*q(1,5) - q(2,3)*q(2,5) ) / q(3,3)
      q(4,5) = ( s(4,5) - q(1,4)*q(1,5) - q(2,4)*q(2,5)
     $                  - q(3,4)*q(3,5) ) / q(4,4)
      q(5,5) = sqrt( s(5,5) - q(1,5)**2 - q(2,5)**2 - q(3,5)**2
     $                      - q(4,5)**2 ) 
      q(2,6) = ( s(2,6) - q(1,2)*q(1,6) ) / q(2,2)
      q(3,6) = ( s(3,6) - q(1,3)*q(1,6) - q(2,3)*q(2,6) ) / q(3,3)
      q(4,6) = ( s(4,6) - q(1,4)*q(1,6) - q(2,4)*q(2,6)
     $                  - q(3,4)*q(3,6) ) / q(4,4)
      q(5,6) = ( s(5,6) - q(1,5)*q(1,6) - q(2,5)*q(2,6)
     $                  - q(3,5)*q(3,6) - q(4,5)*q(4,6) ) / q(5,5)
      q(6,6) = sqrt( s(6,6) - q(1,6)**2 - q(2,6)**2 -q(3,6)**2
     $                      - q(4,6)**2 - q(5,6)**2 ) 
c
c     Backsubstitute
c
      a(6) = - c(6) / q(6,6)
      a(5) = - ( c(5) + q(5,6)*a(6) ) / q(5,5)
      a(4) = - ( c(4) + q(4,5)*a(5) + q(4,6)*a(6) ) / q(4,4)
      a(3) = - ( c(3) + q(3,4)*a(4) + q(3,5)*a(5)
     $                + q(3,6)*a(6) ) / q(3,3)
      a(2) = - ( c(2) + q(2,3)*a(3) + q(2,4)*a(4)
     $                + q(2,5)*a(5) + q(2,6)*a(6) ) / q(2,2)
      a(1) = - ( c(1) + q(1,2)*a(2) + q(1,3)*a(3) + q(1,4)*a(4)
     $                + q(1,5)*a(5) + q(1,6)*a(6) ) / q(1,1)
c
      return

      end
c
c
      double precision function erfc(xx)
c
c
c     sandia mathematical program library
c     applied mathematics division 2613
c     sandia laboratories
c     albuquerque, new mexico  87185
c     control data 6600/7600  version 7.2  may 1978
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * the primary document for the library of which this routine is     *
c  * part is sand77-1441.                                              *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     written by j.e. vogel from approximations derived by w.j. cody .
c
c     abstract
c
c         erfc(x) computes 2.0/sqrt(pi) times the integral from x to
c         infinity of exp(-x**2). this is done using rational approx-
c         imations.  eleven correct significant figures are provided.
c
c     description of parameters
c
c         x may be any real value
c
c     erfc is documented completely in sc-m-70-275.
c
C     .. Scalar Arguments ..
      double precision xx
C     ..
C     .. Local Scalars ..
      double precision a, r, sqpi, x, x2, xi2
      integer i
C     ..
C     .. Local Arrays ..
      double precision p1(4), p2(6), p3(4), q1(4), q2(6), q3(4)
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp
C     ..
C     .. Data statements ..
      data (p1(i),i=1,4)/242.6679552305318D0, 21.97926161829415D0,
     1     6.996383488619136D0, -3.560984370181539D-02/,
     2     (q1(i),i=1,4)/215.0588758698612D0, 91.16490540451490D0,
     3     15.08279763040779D0, 1.0D0/, (p2(i),i=1,6)/22.898992851659D0,
     4     26.094746956075D0, 14.571898596926D0, 4.2677201070898D0,
     5     .56437160686381D0, -6.0858151959688D-06/,
     6     (q2(i),i=1,6)/22.898985749891D0, 51.933570687552D0,
     7     50.273202863803D0, 26.288795758761D0, 7.5688482293618D0,
     8     1.0D0/, (p3(i),i=1,4)/-1.21308276389978D-2,
     9     -.1199039552681460D0, -.243911029488626D0,
     A     -3.24319519277746D-2/, (q3(i),i=1,4)/4.30026643452770D-02,
     B     .489552441961437D0, 1.43771227937118D0, 1.0D0/
      data sqpi/.564189583547756D0/
C     ..
c
      x = abs(xx)
      x2 = x*x
      if (xx .lt. -6.0D0) then
         erfc = 2.0D0
      else if (x .gt. 25.8D0) then
         erfc = 0.0D0
      else if (x .gt. 4.0D0) then
         xi2 = 1.D0/x2
         r = xi2*(p3(1)+xi2*(p3(2)+xi2*(p3(3)+xi2*p3(4))))/
     1       (q3(1)+xi2*(q3(2)+xi2*(q3(3)+xi2*q3(4))))
         a = exp(-x2)*(sqpi+r)/x
         if (xx .lt. 0.0D0) then
            erfc = 2.0D0 - a
         else
            erfc = a
         end if
      else if (x .gt. .46875D0) then
         a = exp(-x2)*(p2(1)+x*(p2(2)+x*(p2(3)+x*(p2(4)+x*(p2(5)+
     1       x*p2(6))))))
         a = a/(q2(1)+x*(q2(2)+x*(q2(3)+x*(q2(4)+x*(q2(5)+x*q2(6))))))
         if (xx .le. 0.0D0) then
            erfc = 2.0D0 - a
         else
            erfc = a
         end if
      else
         a = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))
         a = a/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         if (xx .lt. 0.D0) a = -a
         erfc = 1.0D0 - a
      end if
c
      return
c
      end
c
c
      real*8 function erf(x)

      real*8 x
      real*8 erfc
c
      erf = 1 - erfc(x)
     
      return
      end
