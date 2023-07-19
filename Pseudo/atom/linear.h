c-------
c $Id: linear.h,v 1.3 1997/05/22 17:32:18 wdpgaara Exp $
c
c $Log: linear.h,v $
c Revision 1.3  1997/05/22 17:32:18  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1992/03/17  00:24:07  alberto
c Aligned parity.
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      double precision parity
      double precision alin(5,5)
      integer indx(5)
c
      common /linear/ alin, parity, indx
      save /linear/
c-------
