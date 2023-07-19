c
c $Id: hsc.f,v 1.2 1997/05/22 17:32:15 wdpgaara Exp $
c
c $Log: hsc.f,v $
c Revision 1.2  1997/05/22 17:32:15  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine hsc(i,ar,br)
c
c     Implements the HSC method.
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
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
      double precision small, small2, small3, pzfive
      parameter (small=1.D-13,small2=1.D-10,small3=1.D-18,pzfive=.05D0)
      double precision pfive, small4
      parameter (pfive=0.5D0,small4=1.D-6)
c
      integer i
      double precision ar(*), br(*)
c
      double precision aa, ag, cl, dcl, delta, dev, devold, eviae, fjm1,
     &                 gamma, gg, gpp, rra, rrc, rrp
c
      integer iflag, ist, j, j3rc, k, ll, llp, lp
      logical spin_down
      character id*1
      integer nops(norbmx)
c
      double precision f(nrmax), g(nrmax), v(nrmax), arps(nrmax)
c
      external difnrl, dsolv2, potrw, wtrans
      intrinsic nint, sign
c
c     Do everything non-relativistically...
c
      if (polarized) then
         id = 's'
      else
         id = ' '
      endif
c
c  Reset the n quantum number to give the proper number of
c  nodes (0) for the pseudowavefunction. Zero out the rest.
c
      do 10 j = 1, norb
         nops(j) = 0
   10 continue
      lp = lo(i) + 1
      llp = lo(i) * lp
      rc(lp) = rc_input(lp)
      nops(i) = lp
c
c  njtj  ***  modification start  ***
c  Set up the functions f(r/rc) and g(r/rc) and
c  modify the ionic potential.
c
      aa = 4*one
      dcl = -6*one*lp
      cl = dcl
c
      do 20 j = 1, nr
         rrc = r(j)/rc(lp)
         rra = rrc**aa
         f(j) = zero
         if (rra .lt. 88*one) f(j) = exp(-rra)
         g(j) = rrc**lp*f(j)
         fjm1 = one - f(j)
         if (fjm1 .lt. small4) fjm1 = (one-pfive*rra)*rra
         if (down(i)) then
            viod(lp,j) = fjm1*viod(lp,j) - f(j)*r(j)*vid(j) +
     &                   dcl*r(j)*f(j)
         else
            viou(lp,j) = fjm1*viou(lp,j) - f(j)*r(j)*viu(j) +
     &                   dcl*r(j)*f(j)
         end if
         if (rrc .lt. 3*one) j3rc = j
   20 continue
      dcl = dcl/2
c
c   Start the iteration loop to find cl.
c
      eviae = ev(i)
      devold = zero
      do 60 j = 1, 100
c
         call dsolv2(j,2,id,1,norb,ncore,nops)
         dev = eviae - ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c    division by zero errors in the calculation of
c    dcl = -dev*dcl / (dev-devold).
c
         if (((abs(dev).lt.small2).or.(abs(dev-devold).lt.small3)) .and.
     &       (j.ne.1)) then
c
            go to 70
c
         else
            if (j .gt. 20 .or. abs(dev) .lt. 0.001D0) then
c
c                 Use newton-raphson iteration to change cl.
c
               dcl = -dev*dcl/(dev-devold)
            else
               if (dev*dcl .lt. zero) dcl = -dcl/3
            end if
         end if
c
c  njtj  ***  modification end  ***
c
c  Find the new potential.
c
   30    continue
         if (down(i)) then
            do 40 k = 2, nr
               viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
   40       continue
         else
            do 50 k = 2, nr
               viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
   50       continue
         end if
c
c        Update...
c
         cl = cl + dcl
         devold = dev
c
   60 continue
c
c  End the iteration loop for cl.
c
      call ext(820+lp)
c
c   Find the pseudo-wavefunction.
c
   70 continue
      if (down(i)) then
         do 80 j = 2, nr
            v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
   80    continue
      else
         do 90 j = 2, nr
            v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
   90    continue
      end if
c
      call difnrl(0,i,v,arps,br,nops(i),lo(i),so(i),ev(i),iflag)
c
c  Compute delta and gamma.
c
      gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
      ag = zero
      gg = zero
      ll = 4
      do 100 j = 2, nr
         ag = ag + ll*arps(j)*g(j)*rab(j)
         gg = gg + ll*g(j)*g(j)*rab(j)
         ll = 6 - ll
  100 continue
      ag = ag/3
      gg = gg/3
      delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c     Modify the pseudo-wavefunction and pseudo-potential and
c     add to charge density.
c
      if (down(i)) then
         do 110 j = 2, nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vod(j) = vod(j) + zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j) = small
            rrp = r(j)/rc(lp)
            gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &            r(j)**2
            viod(lp,j) = viod(lp,j) + gamma*delta*
     &                   ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  110    continue
      else
         do 120 j = 2, nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vou(j) = vou(j) + zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j) = small
            rrp = r(j)/rc(lp)
            gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &            r(j)**2
            viou(lp,j) = viou(lp,j) + gamma*delta*
     &                   ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  120    continue
      end if
c
c  wtrans is called to fourier transform the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = nint(sign(1.d0,arps(nr-85)))
      call potrw(arps,r,nr-85,lo(i),0,ist,rc(lp))
      call wtrans(arps,r,nr,lo(i),ist)
c
      write(6,9000) nops(i), il(lp), so(i), ev(i), rc(lp), cl, gamma,
     &  delta
 9000 format(1x,i1,a1,f6.1,5f12.6)
c
c     Convention mandates that the pseudo for the s channel is
c     "down" in a relativistic calculation. To make the wavefunctions
c     adhere to the same convention (encoded in the indd and indu arrays)
c     we do the same here.
c
      spin_down = down(i)
      if (relativistic .and. (lo(i).eq.0)) spin_down = .true.
      call pswf_store(arps,lo(i),spin_down)
c
      return
c
      end
