c------
c $Id: ion.h,v 1.2 1997/05/22 17:32:17 wdpgaara Exp $
c
c $Log: ion.h,v $
c Revision 1.2  1997/05/22 17:32:17  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      integer lmax
      parameter (lmax=4)
c
c     Viod/Viou is the ionic potential times r
c
      double precision viod(lmax,nrmax), viou(lmax,nrmax)
c
      common /ion/ viod, viou
      save /ion/
c------
