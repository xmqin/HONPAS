c------
c $Id: charge.h,v 1.2 1997/05/22 17:32:03 wdpgaara Exp $
c
c $Log: charge.h,v $
c Revision 1.2  1997/05/22 17:32:03  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
c
c     Cdd/Cdu : Down/Up valence charge  (4pi r^2 rho(r))
c     Cdc     :         core charge           "
c
      double precision cdd(nrmax), cdu(nrmax), cdc(nrmax)
c
      common /charge/ cdd, cdu, cdc
      save /charge/
c------
