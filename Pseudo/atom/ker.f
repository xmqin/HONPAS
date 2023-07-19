c
c $Id: ker.f,v 1.2 1997/05/22 17:32:17 wdpgaara Exp $
c
c $Log: ker.f,v $
c Revision 1.2  1997/05/22 17:32:17  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine ker(i,ar,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      double precision zero, one, pfive, smtol
      parameter (zero=0.D0,one=1.D0,pfive=0.5D0,smtol=1.D-12)
      double precision small, ai
      parameter (small=1.D-12,ai=2*137.0360411D0)
c
      integer i
      double precision ar(nrmax), br(nrmax)
c
      integer ist, iswtch, j, jrc, lp, ka, ll, k
      logical spin_down
      double precision polyr, eigv, cdrc, rc2, rc3, rc4, arc, arp, 
     &                 brc, vrc, alpha, beta, gamma, delta, cdps,
     &                 fdnew, ddelta, fdold, expd, xlamda, vj
c
      integer isrchfgt
c
      lp = lo(i) + 1
      ka = lo(i) + 1
      eigv = ev(i)
      if (down(i) .and. lo(i) .ne. 0) ka = -lo(i)
c
c     Reset rc to grid point r(j) such that r(j) <= rc < r(j+1)
c
      jrc = isrchfgt(nr,r,1,rc_input(lp)) - 1
      rc(lp) = r(jrc)
c
c  Find the integrated charge inside rc.
c
      ll = 2
      if (relativistic) then
         cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 40 k = jrc, 1, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   50       continue
            cdrc = cdrc - (ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc + 9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     &             3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     &             3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     &             (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
         end if
         cdrc = cdrc/3
      else
         cdrc = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 60 k = jrc, 1, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   60       continue
         else
            do 70 k = jrc, 4, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   70       continue
            cdrc = cdrc - ar(4)*ar(4)*rab(4)
            cdrc = cdrc + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdrc = cdrc/3
      end if
c
c   The initial values for alpha, beta, gamma and delta.
c
      rc2 = r(jrc)*r(jrc)
      rc3 = r(jrc)*rc2
      rc4 = r(jrc)*rc3
      iswtch = 1
      if (ar(jrc) .lt. zero) iswtch = -1
      arc = iswtch*ar(jrc)
      arp = br(jrc)
c
      if (relativistic) then
         if (down(i)) then
            arp = ka*ar(jrc)/r(jrc) + (eigv-viod(lp,jrc)/r(jrc)-
     &            vid(jrc)+ai*ai)*br(jrc)/ai
         else
            arp = ka*ar(jrc)/r(jrc) + (eigv-viou(lp,jrc)/r(jrc)-
     &            viu(jrc)+ai*ai)*br(jrc)/ai
         end if
      end if
c
      brc = arp/ar(jrc)
      if (down(i)) then
         vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
      else
         vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
      end if
      alpha = (3*log(arc/r(jrc)**lp)-2*(r(jrc)*brc-lp)+
     &        (rc2*vrc+lp*lp-rc2*(eigv+brc*brc))/2)/rc4
      beta = (-8*log(arc/r(jrc)**lp)+5*(r(jrc)*brc-lp)-
     &       (rc2*vrc+lp*lp-rc2*(eigv+brc*brc)))/rc3
      gamma = (6*log(arc/r(jrc)**lp)-3*(r(jrc)*brc-lp)+
     &        (rc2*vrc+lp*lp-rc2*(eigv+brc*brc))/2)/rc2
      delta = zero
c
c  Start the iteration loop to find delta.
c
      do 110 j = 1, 50
c
c  Generate the pseudo-wavefunction (note missing factor exp(delta)).
c
         do 80 k = 1, jrc
            polyr = r(k)*r(k)*((alpha*r(k)+beta)*r(k)+gamma)
            ar(k) = iswtch*r(k)**lp*exp(polyr)
   80    continue
c
c  Integrate  the pseudo charge density from r = 0 to rc.
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 90 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   90       continue
         else
            do 100 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  100       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
c  Find the new delta.
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. smtol) go to 120
         if (j .eq. 1) then
            ddelta = pfive
         else
            ddelta = -fdnew*ddelta/(fdnew-fdold)
         end if
         alpha = alpha - 3*ddelta/rc4
         beta = beta + 8*ddelta/rc3
         gamma = gamma - 6*ddelta/rc2
         delta = delta + ddelta
         fdold = fdnew
  110 continue
c
c  End the iteration loop for delta.
c
      call ext(820+lp)
c
c    Augment the charge density and invert schroedinger equation
c  to find new potential.
c
  120 continue
      expd = exp(delta)
      if (down(i)) then
         do 130 j = 1, jrc
            ar(j) = expd*ar(j)
            xlamda = (4*alpha*r(j)+3*beta)*r(j) + 2*gamma
            vj = eigv + xlamda*(2*lp+xlamda*r(j)**2) +
     &           (12*alpha*r(j)+6*beta)*r(j) + 2*gamma
            viod(lp,j) = (vj-vid(j))*r(j)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
  130    continue
         do 140 j = jrc + 1, nr
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
  140    continue
      else
         do 150 j = 1, jrc
            ar(j) = expd*ar(j)
            xlamda = (4*alpha*r(j)+3*beta)*r(j) + 2*gamma
            vj = eigv + xlamda*(2*lp+xlamda*r(j)**2) +
     &           (12*alpha*r(j)+6*beta)*r(j) + 2*gamma
            viou(lp,j) = (vj-viu(j))*r(j)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
  150    continue
         do 160 j = jrc + 1, nr
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
  160    continue
      end if
      write(6,9000) lo(i)+1, il(lp), so(i), eigv, rc(lp), cdrc, delta
 9000 format(1x,i1,a1,f6.1,5f12.6)
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = 1
      if (ar(nr-85) .lt. zero) ist = -1
      call potrw(ar,r,nr-85,lo(i),0,ist,rc(lp))
      call wtrans(ar,r,nr,lo(i),ist)
c
c     Convention mandates that the pseudo for the s channel is
c     "down" in a relativistic calculation. To make the wavefunctions
c     adhere to the same convention (encoded in the indd and indu arrays)
c     we do the same here.
c
      spin_down = down(i)
      if (relativistic .and. (lo(i).eq.0)) spin_down = .true.
      call pswf_store(ar,lo(i),spin_down)
c
      return
c
      end
