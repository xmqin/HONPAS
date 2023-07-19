c-----
c $Id: elecpot.h,v 1.2 1997/05/22 17:32:10 wdpgaara Exp $
c
c $Log: elecpot.h,v $
c Revision 1.2  1997/05/22 17:32:10  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  02:06:34  alberto
c Initial revision
c
c
c     Vid/Viu : Down/Up "input" potentials
c     Vod/Vou : Down/Up "output" potentials
c
      double precision vid(nrmax), viu(nrmax),
     &                 vod(nrmax), vou(nrmax)
c
      common /elecpot/ vid, viu, vod, vou
      save /elecpot/
c------
