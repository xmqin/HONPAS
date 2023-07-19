c------
c $Id: orbital.h,v 1.3 1999/02/26 14:26:44 wdpgaara Exp $
c
      integer norbmx
      parameter (norbmx=40)
c
c     norb:     Total number of orbitals.
c     ncp :     First valence orbital.
c
      integer norb, ncp
      character il(5)*1
c
      integer no(norbmx), lo(norbmx)
      double precision so(norbmx), zo(norbmx)
      logical down(norbmx)
c
      common /orbital/ so, zo
      common /orb_int/ norb, ncp, no, lo
      common /orb_char/ il
      common /orb_log/  down
      save /orbital/, /orb_int/, /orb_char/, /orb_log/
c------
