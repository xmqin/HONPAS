c------
c $Id: nonlinear.h,v 1.2 1997/05/22 17:32:19 wdpgaara Exp $
c
c $Log: nonlinear.h,v $
c Revision 1.2  1997/05/22 17:32:19  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      double precision arc, brc, eigv, cdrc, vrc, vap, vapp, 
     &                 rc1, rc2, rc3
      integer lp, jrc
c
      common /nonlinear/ arc, brc, eigv, cdrc, vrc, vap, vapp,  
     &                   rc1, rc2, rc3, lp, jrc
      save /nonlinear/
c------
