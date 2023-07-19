c------
c $Id: tm2_blk.h,v 1.2 1997/05/22 17:32:34 wdpgaara Exp $
c
c $Log: tm2_blk.h,v $
c Revision 1.2  1997/05/22 17:32:34  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      double precision alpha, alpha1, alpha2, alpha3, alpha4, arc, brc,
     &                 cdrc, delta, eigv, rc1, rc2, rc3, rc4, rc5,
     &                 rc6, rc7, rc8, vap, vapp, vrc
      integer jrc, lp
      double precision ar(nrmax)
c
      common /tm2_blk/ alpha, alpha1, alpha2, alpha3, alpha4,
     &                 arc, brc, cdrc, delta, eigv,
     &                 rc1, rc2, rc3, rc4, rc5, rc6, rc7, rc8,
     &                 vap,  vapp, vrc, 
     &                 jrc, lp,
     &                 ar
      save /tm2_blk/
c------
