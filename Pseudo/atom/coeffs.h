c------
c $Id: coeffs.h,v 1.2 1997/05/22 17:32:05 wdpgaara Exp $
c
c $Log: coeffs.h,v $
c Revision 1.2  1997/05/22 17:32:05  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      double precision alpha, alpha1, alpha2, alpha3, alpha4, 
     &                 delta, gamma
      integer ang_moment
c
      common /coeffs/ alpha, alpha1, alpha2, alpha3, alpha4,
     &                delta, gamma,
     &                ang_moment
      save /coeffs/
c------
