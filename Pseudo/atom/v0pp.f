c
      double precision function v0pp(gamma)
c
      implicit none
c
      include 'radial.h'
      include 'tm2_blk.h'
c
      double precision gamma
c
C     .. Parameters ..
      double precision zero, pfive, one, errmin
      parameter (zero=0.D0,pfive=0.5D0,one=1.D0,errmin=1.D-12)
C     ..
C     .. Local Scalars ..
      double precision bj1, bj2, bj2a, bj3, bj3a, bj4, bj5, cdps,
     &                 ddelta, fdnew, fdold, polyr, r2, rc10, rc11,
     &                 rc12, rc9, rp, rcond
      integer j, k, ll
C     ..
C     .. Local Arrays ..
      double precision aj(5,5), bj(5), work(5)
      integer indx(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log
C     ..

      fdold = 0.d0

      rc9 = rc8*rc1
      rc10 = rc8*rc2
      rc11 = rc8*rc3
      rc12 = rc8*rc4
c
      delta = zero
c
      bj(1) = log(arc/rc1**lp) - gamma*rc2
      bj1 = bj(1)
      bj(2) = brc - lp/rc1 - 2*gamma*rc1
      bj2a = bj(2) + 2*gamma*rc1
      bj2 = bj(2)
      bj(3) = vrc - eigv - 2*lp/rc1*bj2a - bj2a**2 - 2*gamma
      bj3 = bj(3)
      bj3a = bj(3) + 2*gamma
      bj(4) = vap + 2*lp/rc2*bj2a - 2*lp/rc1*bj3a - 2*bj2a*bj3a
      bj4 = bj(4)
      bj(5) = vapp - 4*lp/rc3*bj2a + 4*lp/rc2*bj3a - 2*lp/rc1*bj4 -
     &        2*bj3a**2 - 2*bj2a*bj4
      bj5 = bj(5)
c
      aj(1,1) = rc4
      aj(1,2) = rc6
      aj(1,3) = rc8
      aj(1,4) = rc10
      aj(1,5) = rc12
      aj(2,1) = 4*rc3
      aj(2,2) = 6*rc5
      aj(2,3) = 8*rc7
      aj(2,4) = 10*rc9
      aj(2,5) = 12*rc11
      aj(3,1) = 12*rc2
      aj(3,2) = 30*rc4
      aj(3,3) = 56*rc6
      aj(3,4) = 90*rc8
      aj(3,5) = 132*rc10
      aj(4,1) = 24*rc1
      aj(4,2) = 120*rc3
      aj(4,3) = 336*rc5
      aj(4,4) = 720*rc7
      aj(4,5) = 1320*rc9
      aj(5,1) = 24*one
      aj(5,2) = 360*rc2
      aj(5,3) = 1680*rc4
      aj(5,4) = 5040*rc6
      aj(5,5) = 11880*rc8
c
c     Use LU decomposition (AG, April 1991) See Numerical Recipes.
c
      call sgeco(aj,5,5,indx,rcond,work)
      if (rcond .lt. 1.d-7) write(6,*) ' rcond too small:' ,rcond
c
      call sgesl(aj,5,5,indx,bj,0)
c
      alpha = bj(1)
      alpha1 = bj(2)
      alpha2 = bj(3)
      alpha3 = bj(4)
      alpha4 = bj(5)
c
c   start iteration loop to find delta (with gamma fixed)
c
      do 80 j = 1, 100
c
c   generate pseudo wavefunction-note missing factor exp(delta)
c
         do 30 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &              alpha)*r2+gamma)
            ar(k) = rp**lp*exp(polyr)
   30    continue
c
c   integrate pseudo charge density from r = 0 to rc
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (jrc .ne. 2*(jrc/2)) then
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
c        Calculate new delta (with gamma fixed), uses false position
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. errmin) then
            v0pp = 8*((2*one*(lp-one)+5*one)*alpha+gamma**2)
c
            return
c
         end if
c
         if (j .eq. 1) then
            ddelta = -pfive
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
         alpha = bj(1)
         alpha1 = bj(2)
         alpha2 = bj(3)
         alpha3 = bj(4)
         alpha4 = bj(5)
c
         fdold = fdnew
c
   80 continue
      v0pp = 1.d-60
      write(6,9000)
 9000 format(//'error in gamfind (aka v0pp) - delta not found')
      call ext(860+lp)
c
      end

