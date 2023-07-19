C
c $Id: difrel.f,v 1.3 1997/05/22 17:32:07 wdpgaara Exp $
c
c $Log: difrel.f,v $
c Revision 1.3  1997/05/22 17:32:07  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1994/02/18  01:26:18  garcia
c *** empty log message ***
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine difrel(iter,iorb,v,ar,br,n,l,spin,eigv)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c     
c  difrel integrates the relativistic Dirac equation
c  it finds the eigenvalue ev, the major and minor component
c  of the wavefunction, ar and br.  It uses an intial guess
c  for the eigenvalues from dsolv1
c
c  njtj  ***  modifications  ***
c    This routine has major modifications.
c    1)The data needed inside the loops has been calculated
c    outside the main loop(increases speed for non-opt
c    compilers, i.e. dumb compilers).
c    2)The predict/correct values are placed in an array.
c    Output is unchanged
c  njtj  ***  modifications  ***
c
c  twb   ***  modifications  ***
c    Underflow trap fixed
c  twb   ***  modifications  ***
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
c  Tolerance
c
c
c  Integration coefficients
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
C     .. Parameters ..
      double precision zero, pnine, one, ai
      parameter (zero=0.D0,pnine=0.9D0,one=1.D0,ai=2*137.0360411D0)
      double precision etol
      parameter (etol=-1.D-7)
      double precision tol
      parameter (tol=1.D-10)
      double precision abc1, abc2, abc3, abc4, abc5, amc0, amc1, amc2,
     &                 amc3, amc4
      parameter (abc1=190.1D0/72,abc2=-138.7D0/36,abc3=10.9D0/3,
     &          abc4=-63.7D0/36,abc5=25.1D0/72,amc0=25.1D0/72,
     &          amc1=32.3D0/36,amc2=-1.1D0/3,amc3=5.3D0/36,
     &          amc4=-1.9D0/72)
C     ..
C     .. Scalar Arguments ..
      double precision eigv, spin
      integer iorb, iter, n, l
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*), v(*)
C     ..
C     .. Local Scalars ..
      double precision a1, a2, ai2, alf, arc, arin, arout, arp, arpin,
     &                 arpout, az, b0, b1, b2, brc, brp, dev, emax,
     &                 emin, evold, evv, evvai2, expzer, factor, faj,
     &                 fbj, s, temp, vzero
      integer icount, istop, itmax, j, juflow, ka, ll, nctp, ninf, nodes
C     ..
C     .. Local Arrays ..
      double precision rs(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log, sign, sqrt
C     ..
C     .. Arrays in Common ..
      double precision fa(nrmax), fb(nrmax), rabai(nrmax), rabkar(nrmax)
C     ..
C     .. Common blocks ..
      common  rabkar, rabai, fa, fb
C     ..
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
      itmax = 100
      ai2 = ai*ai
      az = znuc/(2*ai)
      ka = l + 1
      if (spin .lt. 0.1D0 .and. l .ne. 0) ka = -l
c
c  determine effective charge and vzero for startup of
c  outward integration
c  ar = r**s * (1  + a1 r + a2 r**2 + ... )
c  br = r**s * (b0 + b1 r + b2 r**2 + ... )
c  s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c  an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c        / (n ai (2 s + n))
c  bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) then
         b0 = -az/(s+ka)
      else
         b0 = (s-ka)/az
      end if
      if (spin .lt. 0.1D0) then
         vzero = vid(2)
      else
         vzero = viu(2)
      end if
c
c  njtj  ***  start major modification  ***
c    Loop data calculated only once.
c    Set ar() and br() to zero.
c
      do 10 j = 1, nr
         ar(j) = zero
         br(j) = zero
   10 continue
      do 20 j = 2, nr
         rabkar(j) = rab(j)*ka/r(j)
   20 continue
      do 30 j = 2, nr
         rabai(j) = rab(j)/ai
   30 continue
      do 40 j = 2, 5
         rs(j) = r(j)**s
   40 continue
c
c  set the underflow trap.
c  twb *** begin modification ***
c
      juflow = 1
      do 50 j = 2, nr
         if (s*abs(log(r(j))) .lt. sqrt(expzer)) go to 60
         juflow = j
   50 continue
   60 continue
c
c  twb *** end modification ***
c  njtj *** end major modification  ***
c
      emax = zero
      emin = -one*100000
      if (eigv .gt. emax) eigv = emax
   70 continue
      if (itmax .lt. 2) write(6,9000) iorb, iter, eigv, nodes
 9000 format(' iorb =',i3,' iter =',i3,' ev =',d18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (eigv .gt. zero) then
         write(6,9010) iorb
         call ext(620+iorb)
      end if
 9010 format(//' error in difrel - ev(',i2,') greater then v(infinty)')
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
c
      icount = 0
   80 continue
      icount = icount + 1
      do 90 j = nr, 2, -1
         temp = v(j) - eigv
         if (temp .lt. zero) temp = zero
         if (r(j)*sqrt(temp) .lt. expzer) go to 100
   90 continue
  100 continue
      ninf = j
      nctp = ninf - 5
      do 110 j = 2, ninf - 5
         if (v(j) .lt. eigv) nctp = j
  110 continue
      if (eigv .ge. etol*100) nctp = ninf - 5
      if (eigv .ge. etol) eigv = zero
      if (nctp .le. 6) then
         eigv = pnine*eigv
         if (icount .gt. 100) then
            write(6,9020) iorb
            call ext(650+iorb)
         end if
c
         go to 80
c
      end if
 9020 format(//'error in difrel - cannot find classical',
     &      /'turning point in orbital ',i2)
c
c  Outward integration from 1 to nctp, startup.
c
      a1 = (az*(vzero-eigv)-(s+1+ka)*(vzero-eigv-ai2)*b0)/
     &     (ai*(2*s+1))
      b1 = ((vzero-eigv)-2*znuc*a1)/(ai*(s+1+ka))
      a2 = (az*(vzero-eigv)*a1-(s+2+ka)*(vzero-eigv-ai2)*b1)/
     &     (2*ai*(2*s+2))
      b2 = ((vzero-eigv)*a1-2*znuc*a2)/(ai*(s+2+ka))
      do 120 j = 2, 5
         ar(j) = rs(j)*(1+(a1+a2*r(j))*r(j))
         br(j) = rs(j)*(b0+(b1+b2*r(j))*r(j))
  120 continue
      fa(1) = zero
      fb(1) = zero
      fa(2) = rabkar(2)*ar(2) + (eigv-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2) - (eigv-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3) + (eigv-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3) - (eigv-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4) + (eigv-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4) - (eigv-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5) + (eigv-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5) - (eigv-v(5))*ar(5)*rabai(5)
c
c  Intergration loop.
c
      nodes = 0
      do 130 j = 6, nctp
c
c  Predictor (Adams-Bashforth).
c
         evvai2 = eigv - v(j) + ai2
         evv = eigv - v(j)
         arp = ar(j-1) + abc1*fa(j-1) + abc2*fa(j-2) + abc3*fa(j-3) +
     &         abc4*fa(j-4) + abc5*fa(j-5)
         brp = br(j-1) + abc1*fb(j-1) + abc2*fb(j-2) + abc3*fb(j-3) +
     &         abc4*fb(j-4) + abc5*fb(j-5)
         fa(j) = rabkar(j)*arp + evvai2*brp*rabai(j)
         fb(j) = -rabkar(j)*brp - evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
         arc = ar(j-1) + amc0*fa(j) + amc1*fa(j-1) + amc2*fa(j-2) +
     &         amc3*fa(j-3) + amc4*fa(j-4)
         brc = br(j-1) + amc0*fb(j) + amc1*fb(j-1) + amc2*fb(j-2) +
     &         amc3*fb(j-3) + amc4*fb(j-4)
         faj = rabkar(j)*arc + evvai2*brc*rabai(j)
         fbj = -rabkar(j)*brc - evv*arc*rabai(j)
c
c  Error reduction step.
c
         ar(j) = arc + amc0*(faj-fa(j))
         br(j) = brc + amc0*(fbj-fb(j))
         fa(j) = rabkar(j)*ar(j) + evvai2*br(j)*rabai(j)
         fb(j) = -rabkar(j)*br(j) - evv*ar(j)*rabai(j)
c
c  Count nodes - if no underflow.
c
         if (j .gt. juflow .and. ar(j)*ar(j-1) .lt.
     &       zero) nodes = nodes + 1
  130 continue
      arout = ar(nctp)
      arpout = fa(nctp)
c
c  End outward integration.
c  If number of nodes correct, start inward integration
c  else modify energy stepwise and try again.
c
      if (nodes .ne. n-l-1) then
c
c  too many nodes decrease ev
c
         if (nodes .gt. n-l-1) then
            if (eigv .lt. emax) emax = eigv
            eigv = eigv + eigv/10
c
c  too few nodes increase ev
c
         else
            if (eigv .gt. emin) emin = eigv
            eigv = eigv - eigv/10
         end if
         itmax = itmax - 1
c
         go to 70
c
      end if
c
c  Inward integration from ninf to nctp startup.
c
      do 140 j = ninf, ninf - 4, -1
         alf = v(j) - eigv
         if (alf .lt. zero) alf = zero
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-eigv-ai2)
  140 continue
      fa(ninf) = rabkar(ninf)*ar(ninf) +
     &           (eigv-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf) -
     &           (eigv-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1) +
     &             (eigv-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1) -
     &             (eigv-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2) +
     &             (eigv-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2) -
     &             (eigv-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3) +
     &             (eigv-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3) -
     &             (eigv-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4) +
     &             (eigv-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4) -
     &             (eigv-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)
c
c  Integration loop.
c
      istop = ninf - nctp
      if (istop .lt. 5) go to 160
      do 150 j = ninf - 5, nctp, -1
c
c  Predictor (Adams-Bashforth).
c
         evvai2 = eigv - v(j) + ai2
         evv = eigv - v(j)
         arp = ar(j+1) - (abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)+
     &         abc4*fa(j+4)+abc5*fa(j+5))
         brp = br(j+1) - (abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)+
     &         abc4*fb(j+4)+abc5*fb(j+5))
         fa(j) = rabkar(j)*arp + evvai2*brp*rabai(j)
         fb(j) = -rabkar(j)*brp - evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
         arc = ar(j+1) - (amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)+
     &         amc3*fa(j+3)+amc4*fa(j+4))
         brc = br(j+1) - (amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)+
     &         amc3*fb(j+3)+amc4*fb(j+4))
         faj = rabkar(j)*arc + evvai2*brc*rabai(j)
         fbj = -rabkar(j)*brc - evv*arc*rabai(j)
c
c  Error reduction step.
c
         ar(j) = arc + amc0*(faj-fa(j))
         br(j) = brc + amc0*(fbj-fb(j))
         fa(j) = rabkar(j)*ar(j) + evvai2*br(j)*rabai(j)
         fb(j) = -rabkar(j)*br(j) - evv*ar(j)*rabai(j)
  150 continue
  160 continue
      arin = ar(nctp)
      arpin = fa(nctp)
c
c  End inward integration
c  Rescale ar and br outside nctp to match ar(nctp) from
c  outward integration.
c
      factor = arout/arin
      do 170 j = nctp, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  170 continue
      arpin = factor*arpin
c
c  Find the normalizing factor.
c
      factor = zero
      ll = 4
      do 180 j = 2, ninf
         factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
         ll = 6 - ll
  180 continue
      factor = factor/3
c
c  Modify the eigenvalue ev.
c
      dev = arout*(arpout-arpin)/(factor*rab(nctp))
      if (5*abs(dev) .gt. -eigv) dev = sign(eigv,dev)/5
      itmax = itmax - 1
      evold = eigv
      eigv = eigv + dev
      if (eigv .gt. emax) then
         eigv = (evold+emax)/2
      else if (eigv .lt. emin) then
         eigv = (evold+emin)/2
      end if
      if (abs(dev) .gt. tol*(1-eigv)) go to 70
c
c  Normalize the wavefunction.
c
      factor = 1/sqrt(factor)
      do 190 j = 1, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  190 continue
c
      return
c
      end
