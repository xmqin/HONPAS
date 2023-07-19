c
c $Id: difnrl.f,v 1.3 1997/05/22 17:32:06 wdpgaara Exp $
c
c $Log: difnrl.f,v $
c Revision 1.3  1997/05/22 17:32:06  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1994/02/18  01:26:11  garcia
c *** empty log message ***
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine difnrl(iter,iorb,v,ar,br,n,l,spin,eigv,iflag)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue eigv, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
c    iorb:       orbital number
c    n, l, spin: orbital's quantum numbers
c    eigv:       eigenvalue
c
c  njtj  ***  modifications  ***
c    This routine has had major modifications.  Some
c    of the data used inside the main loop has been
c    calculated outside the main loop to reduce the number
c    of operations(uses extra array space to gain speed)
c    The predictor-corrector functions have been put
c    into a array.
c    The iflag variable was added to indicate nonconvergence
c    for other programs.  It has no use in the atom program
c    and can be removed by the user.
c    All output from the routine is compatible to
c    the Berkeley/Sverre Froyen version.
c  njtj  ***  modifications  ***
c
c  twb   ***  modifications  ***
c    Underflow trap fixed
c  twb   ***  modifications  ***
c
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
c  Integration coefficients
c
c
c  njtj  *** start modification  ***
c    Arrays added to gain speed.
c
c
c  njtj  ***  end modification  ***
c
C     .. Parameters ..
      double precision zero, pnine, two, etol
      parameter (zero=0.D0,pnine=0.9D0,two=2.D0,etol=-1.D-7)
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
      integer iflag, iorb, iter, n, l
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*), v(*)
C     ..
C     .. Local Scalars ..
      double precision aa, alf, arc, arctp, arp, bb, brc, brctp, brp,
     &                 dev, emax, emin, evold, expzer, factor, fb0, fb1,
     &                 temp, var0, vev, vzero, zeff
      integer icount, istop, itmax, j, j1, j2, j3, j4, j5, juflow, ll,
     &        lp, nctp, ninf, ninf1, ninf2, ninf3, ninf4, nodes
C     ..
C     .. Local Arrays ..
      double precision rabrlo(5), rlp(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log, sign, sqrt
C     ..
C     .. Arrays in Common ..
      double precision fa(nrmax), fb(nrmax), rab2(nrmax)
C     ..
C     .. Common blocks ..
      common  rab2, fa, fb
C     ..
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
c      AG: Let TINY be the smallest normalized number:
c
c      Then: expzer = log(1.d0/sqrt(TINY))
c
c      We can obtain TINY by a call to smach (LINPACK/SCILIB):
c
c      TINY = smach(2)
c
c      Beware of machines with IEEE arithmetic (Sun, etc). 
c
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer =  2.8E3
c
c  njtj  *** major modification start  ***
c
c    Loop data calculated outside loop to gain speed.
c
      itmax = 100
      iflag = 0
      lp = l + 1
      ar(1) = zero
      if (l .eq. 0) then
         br(1) = b*a
      else
         br(1) = zero
      end if
      do 10 j = 2, nr
         ar(j) = zero
         br(j) = zero
   10 continue
c
c     Startup for predictor-corrector (rR goes as r^l near r=0)
c
      do 30 j = 2, 5
         rlp(j) = r(j)**lp
         rabrlo(j) = rab(j)*r(j)**l
   30 continue
c
      do 50 j = 1, nr
         rab2(j) = rab(j)*rab(j)
   50 continue
c
c   set underflow trap
c   twb *** begin modification ***
c
      juflow = 1
      do 60 j = 2, nr
         if (lp*abs(log(r(j))) .lt. sqrt(expzer)) go to 70
         juflow = j
   60 continue
   70 continue
c
c  twb *** end modification ***
c  njtj  *** end major modification  ***
c
c   determine effective charge and vzero for startup of
c   outward integration
c
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = zero
      if (spin .lt. 0.1D0 .and. viod(lp,2) .lt. -0.1D0) zeff = znuc
      if (spin .gt. 0.1D0 .and. viou(lp,2) .lt. -0.1D0) zeff = znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. zero) then
         if (spin .lt. 0.1D0) then
            vzero = vzero + viod(lp,2)/r(2)
         else
            vzero = vzero + viou(lp,2)/r(2)
         end if
      end if
      if (spin .lt. 0.1D0) then
         vzero = vzero + vid(2)
      else
         vzero = vzero + viu(2)
      end if
      var0 = zero
      if (l .eq. 0) var0 = -2*zeff
      if (l .eq. 1) var0 = two
c
      emax = zero
      emin = -two*100000
      if (eigv .gt. emax) eigv = emax
   80 continue
      if (itmax .lt. 2) write(6,9000) iorb, iter, eigv, nodes
 9000 format(' iorb =',i3,' iter =',i3,' ev =',d18.10,' nodes =',i2)
      if (itmax .eq. 0) then
         iflag = 1
c
         return
c
      end if
      if (eigv .gt. zero) then
         write(6,9010) iorb
         call ext(620+iorb)
      end if
 9010 format(//' error in difnrl - ev(',i2,') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount = 0
   90 continue
      icount = icount + 1
      do 100 j = nr, 2, -1
         temp = v(j) - eigv
         if (temp .lt. zero) temp = zero
         if (r(j)*sqrt(temp) .lt. expzer) go to 110
  100 continue
  110 continue
      ninf = j
      nctp = ninf - 5
      do 120 j = 2, ninf - 5
         if (v(j) .lt. eigv) nctp = j
  120 continue
      if (eigv .ge. etol*10) nctp = ninf - 5
      if (eigv .ge. etol) eigv = zero
      if (nctp .le. 6) then
         eigv = pnine*eigv
         if (icount .gt. 100) then
            write(6,9020) iorb
            call ext(650+iorb)
         end if
c
         go to 90
c
      end if
 9020 format(//'error in difnrl - cannot find the classical ',
     &      /' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-eigv)/(4*lp+2)
      do 130 j = 2, 5
         ar(j) = rlp(j)*(1+(aa+bb*r(j))*r(j))
         br(j) = rabrlo(j)*(lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
  130 continue
c
c  njtj  ***  start major modification  ***
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-eigv)*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-eigv)*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-eigv)*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-eigv)*ar(5)
c
c   integration loop
c
      nodes = 0
      do 140 j = 6, nctp
c
c   predictor (Adams-Bashforth)
c
         j1 = j - 1
         j2 = j - 2
         j3 = j - 3
         j4 = j - 4
         j5 = j - 5
         vev = v(j) - eigv
         arp = ar(j1) + abc1*fa(j1) + abc2*fa(j2) + abc3*fa(j3) +
     &         abc4*fa(j4) + abc5*fa(j5)
         brp = br(j1) + abc1*fb(j1) + abc2*fb(j2) + abc3*fb(j3) +
     &         abc4*fb(j4) + abc5*fb(j5)
         fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
         arc = ar(j1) + amc0*brp + amc1*fa(j1) + amc2*fa(j2) +
     &         amc3*fa(j3) + amc4*fa(j4)
         brc = br(j1) + amc0*fb1 + amc1*fb(j1) + amc2*fb(j2) +
     &         amc3*fb(j3) + amc4*fb(j4)
         fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
         ar(j) = arc + amc0*(brc-brp)
         br(j) = brc + amc0*(fb0-fb1)
         fa(j) = br(j)
         fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
         if (j .gt. juflow .and. ar(j)*ar(j-1) .lt.
     &       zero) nodes = nodes + 1
  140 continue
c
c  njtj  ***  end major modification  ***
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (nodes .ne. n-l-1) then
         if (nodes .lt. n-l-1) then
c
c  too few nodes; increase ev
c
            if (eigv .gt. emin) emin = eigv
            eigv = eigv - eigv/10
         else
c
c  too many nodes; decrease ev
c
            if (eigv .lt. emax) emax = eigv
            eigv = eigv + eigv/10
         end if
         itmax = itmax - 1
c
         go to 80
c
      end if
c
c   inward integration from ninf to nctp
c   startup
c
      do 150 j = ninf, ninf - 4, -1
         alf = v(j) - eigv
         if (alf .lt. zero) alf = zero
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = -rab(j)*alf*ar(j)
  150 continue
c
c  njtj  ***  start major modification  ***
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*(v(ninf)-eigv)*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*(v(ninf1)-eigv)*
     &            ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*(v(ninf2)-eigv)*
     &            ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*(v(ninf3)-eigv)*
     &            ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*(v(ninf4)-eigv)*
     &            ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) go to 170
      do 160 j = ninf - 5, nctp, -1
c
c   predictor (Adams-Bashforth)
c
         j1 = j + 1
         j2 = j + 2
         j3 = j + 3
         j4 = j + 4
         j5 = j + 5
         vev = v(j) - eigv
         arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     &         abc4*fa(j4)+abc5*fa(j5))
         brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     &         abc4*fb(j4)+abc5*fb(j5))
         fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
         arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+amc3*fa(j3)+
     &         amc4*fa(j4))
         brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+amc3*fb(j3)+
     &         amc4*fb(j4))
c
         fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
         ar(j) = arc - amc0*(brc-brp)
         br(j) = brc - amc0*(fb1-fb0)
         fa(j) = br(j)
         fb(j) = b*br(j) + rab2(j)*vev*ar(j)
  160 continue
c
c   end inward integration
c
c  njtj  *** end major modification  ***
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  170 continue
      factor = arctp/ar(nctp)
      do 180 j = nctp, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  180 continue
c
c   find normalizing factor
c
      factor = zero
      ll = 4
      do 190 j = 2, ninf
         factor = factor + ll*ar(j)*ar(j)*rab(j)
         ll = 6 - ll
  190 continue
      factor = factor/3
c
c   modify eigenvalue ev
c
      dev = arctp*(brctp-br(nctp))/(factor*rab(nctp))
      if (5*abs(dev) .gt. -eigv) dev = sign(eigv,dev)/5
      itmax = itmax - 1
      evold = eigv
      eigv = eigv + dev
      if (eigv .gt. emax) eigv = (evold+emax)/2
      if (eigv .lt. emin) eigv = (evold+emin)/2
      if (abs(dev) .gt. tol*(1-eigv)) go to 80
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1/sqrt(factor)
      do 200 j = 1, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)/rab(j)
  200 continue
c
      return
c
      end
