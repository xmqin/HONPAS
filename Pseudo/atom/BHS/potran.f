C
c $Id: potran.f,v 1.1 1997/05/22 18:12:33 wdpgaara Exp $
c
c $Log: potran.f,v $
c Revision 1.1  1997/05/22 18:12:33  wdpgaara
c *** empty log message ***
c
c Revision 1.1.1.1  1997/01/07 08:38:56  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine potran(i,vd,r,nr,zion)
c
      implicit none
c
      include 'plot.h'
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should adjust   *
c *  for their own needs.  The potential is fitted with a   *
c *  second degree polynomial, which is muliplied with the  *
c *  appropriate functions and then integrated by parts     *
c *  to find the fourier transform.  The result is then     *
c *  printed to the current plot.dat file (unit=3) for      *
c *  later plotting.  A marker(marker fn#) is placed at     *
c *  the end of each set of data.                           *
c *                                                         *
c ***********************************************************
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Scalar Arguments ..
      double precision zion
      integer i, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision d1, d2, d3, q, q2, r0, rm, rp, v0, vline, vm, vp
      integer j, k
      character*9 filename
C     ..
C     .. Local Arrays ..
      double precision vql(100)
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, cos, sin
C     ..
      rm = zero
      vm = 2*zion
      do 10 k = 2, nr-1, 2
         r0 = r(k)
         v0 = r0*vd(k) + 2*zion
         rp = r(k+1)
         vp = rp*vd(k+1) + 2*zion
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
c  Find the fourier transform q^2/4pi/zion*vql. Everything is
c  rescaled  by zion.  Integration is from q=1 to q=25.
c
      do 30 j = 1, 100
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
         vql(j) = vql(j)/2/zion - one
   30 continue
c
c  Print out the transforms( really q^2/(4pi*zion)*v(q) ) to
c  the current plot.dat file (unit=3) for latter plotting.
c
      write(filename,9900) i-1
 9900 format('PSPOTQ',i1)
      open(unit=3,file=filename,form='formatted')
c
      do 40 j = 1, 100
         write(3,9000) one/4*j, vql(j), 0.d0
   40 continue
 9000 format(1x,f7.4,3x,2f10.6)
c
      close(3)
c
cag      write(3,9010) i
c 9010 format(1x,'marker fn',i1)
c
c     Compute the absolute area 
c
      vline = 7*one + 32*abs(vql(1)) + 12*abs(vql(2)) + 32*abs(vql(3)) +
     &        7*abs(vql(4))
      do 50 j = 4, 88, 4
         vline = vline + 7*abs(vql(j)) + 32*abs(vql(j+1)) +
     &           12*abs(vql(j+2)) + 32*abs(vql(j+3)) + 7*abs(vql(j+4))
   50 continue
      write(6,9020) i - 1, vline/90
 9020 format(1x,'The Fourier(q^2/(4pi*zion)*V(q)) absolute',
     &      ' area for l=',i1,' is ',f10.6)
c
      return
c
      end
