c------
c $Id: energy.h,v 1.2 1997/05/22 17:32:10 wdpgaara Exp $
c
c  ev: Eigenvalues
c  ep: Potential energy
c  ek: Kinetic energy
c
c  etot(10): total energy array
c
      double precision ev(norbmx), ep(norbmx), ek(norbmx)
      double precision etot(10)
c
      common /energy/ ev, ep, ek, etot
      save /energy/
c------
