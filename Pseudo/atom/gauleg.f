c
      subroutine gauleg(x1,x2,x,w,n)
C
C $Id: gauleg.f,v 1.1 2000/02/10 17:56:11 wdpgaara Exp $
C
C $Log: gauleg.f,v $
C Revision 1.1  2000/02/10 17:56:11  wdpgaara
C Implement Fourier transform of core charge.
C
C Revision 1.1.1.1  1997/01/07 08:37:19  wdpgaara
C PS fourier transform package
C
c Revision 1.3  1991/12/13  23:32:01  alberto
c More twiddling
c
c Revision 1.2  1991/12/13  23:12:50  alberto
c Cosmetic changes only
c
c Revision 1.1  1991/12/13  22:56:11  alberto
c Initial revision
c
C     Taken from Numerical Recipes (W.H. Press et al., 1986, p.125)
C
c     Given the lower and upper limits of integration x1 and x2, 
c     and given n, this routine returns arrays x and w of length
c     n, containing the abscissas and weights of the Gauss-Legendre
c     quadrature formula.
c
C     .. Parameters ..
      double precision eps
      parameter (eps=3.D-14)
C     ..
C     .. Scalar Arguments ..
      double precision x1, x2
      integer n
C     ..
C     .. Array Arguments ..
      double precision w(n), x(n)
C     ..
C     .. Local Scalars ..
      double precision p1, p2, p3, pp, xl, xm, z, z1
      integer i, j, m
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, cos
C     ..
c
c     The roots are symmetric in the interval, so we only have to
c     find half of them.
c
      m = (n+1)/2
      xm = 0.5D0*(x2+x1)
      xl = 0.5D0*(x2-x1)
c
c     Loop over the desired roots
c
      do 30 i = 1, m
         z = cos(3.141592654D0*(i-.25D0)/(n+.5D0))
c
c        Starting with the above approximation to the Ith root, we 
c        enter the main loop of refinement by Newton's method
c
   10    continue
         p1 = 1.0D0
         p2 = 0.0D0
c
c        Look up the recurrence relation to get the Legendre 
c        polynomial evaluated at z
c
         do 20 j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
   20    continue
c
c        P1 is now the desired Legendre polynomial. We next compute
c        pp, its derivative, by a standard relation involving also p2,
c        the polynomial of one lower order.
c
         pp = n*(z*p1-p2)/(z*z-1.0D0)
         z1 = z
c
c        Newton's method
c
         z = z1 - p1/pp
         if (abs(z-z1) .gt. eps) go to 10
c
c        Scale the root to the desired interval and put in its
c        symmetric counterpart
c
         x(i) = xm - xl*z
         x(n+1-i) = xm + xl*z
c
c       Compute the weight and its symmetric counterpart
c
         w(i) = 2.0D0*xl/((1.0D0-z*z)*pp*pp)
         w(n+1-i) = w(i)
c
   30 continue
c
      return
c
      end
