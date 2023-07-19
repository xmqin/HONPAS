c
c $Id: genrhs.f,v 1.2 1997/05/22 17:32:13 wdpgaara Exp $
c
c $Log: genrhs.f,v $
c Revision 1.2  1997/05/22 17:32:13  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine genrhs(gamma,delta,y)
c
c     Generates the right-hand sides for set of five linear equations
c     that determines alpha(i) given gamma and delta.
c
      implicit none
c
      include 'nonlinear.h'
c
      double precision gamma, delta
      double precision y(5)
c
      double precision a2, a3
c
      y(1) = log(arc/rc1**lp) - gamma * rc2 - delta
      y(2) = brc - lp/rc1 - 2 * gamma * rc1
      a2 = y(2) + 2*gamma*rc1
      y(3) = vrc - eigv - 2*lp/rc1 * a2 - a2**2 - 2*gamma
      a3 = y(3) + 2*gamma
      y(4) = vap + 2*lp/rc2 * a2 - 2*lp/rc1 * a3 - 2*a2*a3
      y(5) = vapp - 4*lp/rc3 * a2 + 4*lp/rc2 * a3 - 2*lp/rc1 * y(4) -
     &       2 * a3**2 - 2 * a2 * y(4)
c
      return
c
      end
