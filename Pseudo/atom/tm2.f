c
c $Id: tm2.f,v 1.5 2002/07/04 18:29:34 wdpgaara Exp $
c
      subroutine tm2(i,wfr,br)
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
      include 'tm2_blk.h'
c
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
      double precision small, ai
      parameter (small=1.D-12,ai=2*137.0360411D0)
      double precision accuracy
      parameter (accuracy=1.d-10)
c
      integer i
      double precision wfr(nrmax), br(nrmax)
c
      double precision work(5)
c
C     .. Local Scalars ..
      double precision arp, bj1, bj2, bj3, bj4, bj5,
     &                 cdps, ddelta, expd, fdnew, fdold, gamma,
     &                 poly, polyr, r2, rc10, rc9, rp, vj, x1, x2,
     &                 xlamda, rcond
      integer ierr, ist, j, k, ka, ll
      logical spin_down, bracketed
C     ..
C     .. Local Arrays ..
      double precision aj(5,5), bj(5)
      double precision aa(nrmax), aap(nrmax), aapp(nrmax), 
     &                 wspl(3*nrmax)
c
      integer indx(5)
C     ..
      integer isrchfgt
      double precision zbrent, v0pp
      external zbrent, v0pp
c
c AG
      fdold = 0.d0
c
      lp = lo(i) + 1
      ka = lo(i) + 1
      eigv = ev(i)
      spin_down = down(i)
      if (spin_down .and. lo(i) .ne. 0) ka = -lo(i)
c
c     Need this to put ar in a common block...
c
      call scopy(nr,wfr,1,ar,1)
c
c     Reset rc to grid point r(j) such that r(j) <= rc < r(j+1)
c
      jrc = isrchfgt(nr,r,1,rc_input(lp)) - 1
      rc(lp) = r(jrc)
c
c  Find the integrated charge inside rc (1-charge outside).
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
c  Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
c  d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
c
      rc1 = r(jrc)
      rc2 = rc1*rc1
      rc3 = rc2*rc1
      rc4 = rc2*rc2
      rc5 = rc4*rc1
      rc6 = rc4*rc2
      rc7 = rc4*rc3
      rc8 = rc4*rc4
      rc9 = rc4*rc5
      rc10 = rc4*rc6
c
      arc = ar(jrc)
      arp = br(jrc)
c
      if (relativistic) then
         if (spin_down) then
            arp = ka*ar(jrc)/r(jrc) + (eigv-viod(lp,jrc)/r(jrc)-
     &            vid(jrc)+ai*ai)*br(jrc)/ai
         else
            arp = ka*ar(jrc)/r(jrc) + (eigv-viou(lp,jrc)/r(jrc)-
     &            viu(jrc)+ai*ai)*br(jrc)/ai
         end if
      end if
cag?         arp = arp
      brc = arp/arc
c
c
      do 500 j = 2, nr
         if (spin_down) then
            aa(j) = viod(lp,j)/r(j) + vid(j)
         else
            aa(j) = viou(lp,j)/r(j) + viu(j)
         endif
  500 continue   
c
c      Use splines to compute V' and V'' 
c   
      aa(1) = aa(2) - (aa(3)-aa(2))*r(2)/(r(3)-r(2))
      call splift(r,aa,aap,aapp,nr,wspl,ierr,0,zero,zero,zero,zero)
c       
      vrc = aa(jrc)
      vap = aap(jrc)
      vapp = aapp(jrc)
c
c
c   Set up matrix without the d2(potential(0)/dr2=0 condition
c   to find an initial guess for gamma.
c   Note that the equations in the paper have been modified to
c   account for the use of rydberg units.
c
      delta = zero
c
      bj(1) = log(arc/rc1**lp)
      bj1 = bj(1)
      bj(2) = brc - lp/rc1
      bj2 = bj(2)
      bj(3) = vrc - eigv - 2*lp/rc1*bj2 - bj2**2
      bj3 = bj(3)
      bj(4) = vap + 2*lp/rc2*bj2 - 2*lp/rc1*bj3 - 2*bj2*bj3
      bj4 = bj(4)
      bj(5) = vapp - 4*lp/rc3*bj2 + 4*lp/rc2*bj3 - 2*lp/rc1*bj4 -
     &        2*bj3**2 - 2*bj2*bj4
      bj5 = bj(5)
c
      aj(1,1) = rc2
      aj(1,2) = rc4
      aj(1,3) = rc6
      aj(1,4) = rc8
      aj(1,5) = rc10
      aj(2,1) = 2*rc1
      aj(2,2) = 4*rc3
      aj(2,3) = 6*rc5
      aj(2,4) = 8*rc7
      aj(2,5) = 10*rc9
      aj(3,1) = 2*one
      aj(3,2) = 12*rc2
      aj(3,3) = 30*rc4
      aj(3,4) = 56*rc6
      aj(3,5) = 90*rc8
      aj(4,1) = zero
      aj(4,2) = 24*rc1
      aj(4,3) = 120*rc3
      aj(4,4) = 336*rc5
      aj(4,5) = 720*rc7
      aj(5,1) = zero
      aj(5,2) = 24*one
      aj(5,3) = 360*rc2
      aj(5,4) = 1680*rc4
      aj(5,5) = 5040*rc6
c
c     Use LU decomposition to solve the linear system of
c     equations. We can re-use the decomposition inside the
c     loop. This is not the case if gaussian elimination is
c     used ! (Alberto Garcia, April 1991) See Numerical Recipes, p. 31.
c
      call sgeco(aj,5,5,indx,rcond,work)
      if (rcond .lt. 1.d-7) write(6,*) ' rcond too small:' ,rcond
c
      call sgesl(aj,5,5,indx,bj,0)
c
      gamma = bj(1)
      alpha = bj(2)
      alpha1 = bj(3)
      alpha2 = bj(4)
      alpha3 = bj(5)
c
c  Start iteration loop to find delta, uses false position.
c
      do 150 j = 1, 50
c
c  Generate pseudo wavefunction-note missing factor exp(delta).
c
         do 100 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*((((alpha3*r2+alpha2)*r2+alpha1)*r2+alpha)*r2+
     &              gamma)
            ar(k) = rp**lp*exp(polyr)
  100    continue
c
c           Integrate pseudo charge density from r = 0 to rc.
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 110 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  110       continue
         else
            do 120 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  120       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
c           Calculate new delta
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. small) go to 160
         if (j .eq. 1) then
            ddelta = -one/2
         else
            ddelta = -fdnew*ddelta/(fdnew-fdold)
         endif
         delta = delta + ddelta
c
         bj(1) = bj1 - delta
         bj(2) = bj2
         bj(3) = bj3
         bj(4) = bj4
         bj(5) = bj5
c
         call sgesl(aj,5,5,indx,bj,0)
c
         gamma = bj(1)
         alpha = bj(2)
         alpha1 = bj(3)
         alpha2 = bj(4)
         alpha3 = bj(5)
c
         fdold = fdnew
c
  150 continue
c
c  End iteration loop for delta.
c
      write(6,9000) lp - 1
      call ext(820+lp)
 9000 format(//'error in pseud2 - nonconvergence in finding',
     &      /' starting delta for angular momentum ',i1)
c
c  Bracket the correct gamma, use gamma and -gamma
c  from above as initial brackets, expands brackets
c  until a root is found..
c
  160 continue
      alpha4 = zero
      x1 = gamma
      x2 = -gamma
c
      call zbrac(v0pp,x1,x2,bracketed)
c
      if ( .not. bracketed) then
         write(6,9010) lp
         call ext(830+lp)
 9010    format(//'Error in zbractk - cannot bracket orbital ',i2)
      end if
c
c  Iteration loop to find correct gamma, uses
c  bisection to find gamma.
c
      gamma = zbrent(v0pp,x1,x2,accuracy)
c
c  Augment charge density and invert schroedinger equation
c  to find new potential.
c
  170 continue
      expd = exp(delta)
         do 180 j = 1, jrc
            r2 = r(j)*r(j)
            poly = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &             alpha)*r2+gamma)
            ar(j) = r(j)**lp*expd*exp(poly)
c
c           For those of us with inferior minds, xlamda = p'(r)/r
c           and the thing goes like this (eq. 23 in TM2, in rydberg):
c
c                                               2
c           V = eigv + 2 p'* l(l+1)/r + p'' + p'   ==>
c
c           V = eigv + (p'/r) * [ l(l+1) + rp' ] + p''
c                        |                  |
c                     "xlamda"         "xlamda*r2"
c
            xlamda = ((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     &               6*alpha1)*r2+4*alpha)*r2 + 2*gamma
c
            vj = eigv + xlamda*(2*lp+xlamda*r2) +
     &           ((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+
     &           30*alpha1)*r2+12*alpha)*r2 + 2*gamma
c
          if (spin_down) then
            viod(lp,j) = (vj-vid(j))*r(j)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
          else
            viou(lp,j) = (vj-viu(j))*r(j)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
          endif
c
  180    continue
c
         do 190 j = jrc + 1, nr
          if (spin_down) then
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
          else
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
          endif
  190    continue
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a useful number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = 1
      call potrw(ar,r,nr-85,lo(i),0,ist,rc(lp))
      call wtrans(ar,r,nr,lo(i),ist)
c
      write(6,9020) lp, il(lp), so(i), eigv, rc(lp), cdrc, delta
 9020 format(1x,i1,a1,f6.1,5f12.6)
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
