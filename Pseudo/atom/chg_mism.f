c
c $Id: chg_mism.f,v 1.2 1997/05/22 17:32:04 wdpgaara Exp $
c
c $Log: chg_mism.f,v $
c Revision 1.2  1997/05/22 17:32:04  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      double precision function chg_mism(xdelta)
c
c     chg_mism gives the discrepancy (one half the log of their
c     ratio) between the pseudo and all-electron charge densities
c     inside rc when the first polynomial coefficient is xdelta.
c
      implicit none
c
      include 'radial.h'
      include 'nonlinear.h'
      include 'coeffs.h'
      include 'linear.h'
c
      double precision xdelta
c
      double precision ar(nrmax), bj(5)
c
      double precision polyr, rp, r2, cdps
      integer k, ll
c
c     Find the {alpha} set for this xdelta and the current gamma.
c
      call genrhs(gamma,xdelta,bj)
      call sgesl(alin,5,5,indx,bj,0)
c
         alpha = bj(1)
         alpha1 = bj(2)
         alpha2 = bj(3)
         alpha3 = bj(4)
         alpha4 = bj(5)
c
c   Generate the pseudo wavefunction
c   (note that exp(xdelta) is put at the end )
c
         do 30 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &              alpha)*r2+gamma)
            ar(k) = rp**lp*exp(polyr)
   30    continue
c
c   Integrate pseudo charge density from r = 0 to rc
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 40 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   50       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
         chg_mism = log(cdrc/cdps) - 2*xdelta
c
         return
c  
         end
