c------
c $Id: plot.h,v 1.3 1997/09/17 17:33:39 wdpgaara Exp $
c
c $Log: plot.h,v $
c Revision 1.3  1997/09/17 17:33:39  wdpgaara
c nrmax updated to 1500... shame!
c
c Revision 1.2  1997/05/22  17:32:24  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
c  Scratch arrays
c
      integer  nrmax
      parameter (nrmax=1500)
c
      double precision a(nrmax), b(nrmax), c(nrmax)
c 
      common /plot/ a, b, c
      save /plot/
c------
