c------
c $Id: plot.h,v 1.1 1997/05/22 18:12:33 wdpgaara Exp $
c
c $Log: plot.h,v $
c Revision 1.1  1997/05/22 18:12:33  wdpgaara
c *** empty log message ***
c
c Revision 1.1.1.1  1997/01/07 08:38:56  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
c  Scratch arrays
c
      integer  nrmax
      parameter (nrmax=1000)
c
      double precision a(nrmax), b(nrmax), c(nrmax)
c 
      common /plot/ a, b, c
      save /plot/
c------
