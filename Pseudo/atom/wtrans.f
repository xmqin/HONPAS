C
c $Id: wtrans.f,v 1.3 1997/05/22 18:05:43 wdpgaara Exp $
c
      subroutine wtrans(vd,r,nr,i,ist)
c
      implicit none
c
      include 'plot.h'
c
c **********************************************************
c *  The wave function is fitted
c *  with a second degree polynomial which is muliplied
c *  with the appropriate functions and then integrated
c *  by parts in taking the fourier transform. 
c **********************************************************
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Scalar Arguments ..
      integer i, ist, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision d1, d2, d3, q, q2, r0, rm, rp, v0, vm, vp
      integer j, k
      character*7 filename
C     ..
C     .. Local Arrays ..
      double precision vql(100)
C     ..
C     .. Intrinsic Functions ..
      intrinsic cos, sin
C     ..
      rm = zero
      vm = zero
      do 10 k = 2, nr-1, 2
         r0 = r(k)
         v0 = vd(k)
         rp = r(k+1)
         vp = vd(k+1)
         d1 = 1/((rp-rm)*(r0-rm))
         d2 = 1/((rp-r0)*(rm-r0))
         d3 = 1/((r0-rp)*(rm-rp))
         a(k) = vm*d1 + v0*d2 + vp*d3
         b(k) = -vm*(r0+rp)*d1 - v0*(rm+rp)*d2 - vp*(rm+r0)*d3
         c(k) = vm*r0*rp*d1 + v0*rm*rp*d2 + vp*rm*r0*d3
         rm = rp
         vm = vp
   10 continue
c
c  Find the fourier transform-vql.
c
      do 30 j = 1, 54
         q = one/4*j
         q2 = q*q
         vql(j) = zero
         rm = zero
         do 20 k = 2, nr - 1, 2
            rp = r(k+1)
            vql(j) = vql(j) + (2*a(k)*rp+b(k))/q*sin(q*rp) -
     &               ((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp) -
     &               (2*a(k)*rm+b(k))/q*sin(q*rm) +
     &               ((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
            rm = rp
   20    continue
         vql(j) = vql(j)/q2
   30 continue
c
c  Print out the transform vql(q) to the current plot.dat
c  file (unit=3) for latter plotting.
c
      write(filename,9900) i
 9900 format('PSWFNQ',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
c
      do 40 j = 1, 48
         write(3,9000) one/4*j, ist*vql(j)
   40 continue
 9000 format(1x,f7.4,3x,f10.6)
c
      close(unit=3)
cag
c      write(3,9010) i
c 9010 format(1x,'marker fw',i1)
c
      return
c
      end
