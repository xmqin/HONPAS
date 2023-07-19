c $Id: radial.h,v 1.3 1997/05/22 18:07:41 wdpgaara Exp $
c
c     nrmax: Maximum number of points in the radial mesh
c
      integer nrmax
      parameter (nrmax=1500)
c
c     The radial coordinate is reparametrized to allow a more
c     precise sampling of the area near the origin:
c
c     r (x) = a { exp(beta*(x-1)) - 1 }
c
c     r'(x) =  (r(x) + a) * beta    (rab in the program)
c
      integer nr
      double precision a, b
      double precision r(nrmax), rab(nrmax)
c
      common /radial/ a, b, r, rab
      common /rad_int/ nr
      save /radial/, /rad_int/
c------





