(Message inbox:52)
Received: from savanna.Berkeley.EDU by jungle.Berkeley.EDU (AIX 3.1/UCB 5.61/4.03(?))
          id AA32236; Sat, 11 Jul 92 17:06:41 -0700
Received: from jungle.Berkeley.EDU by savanna.Berkeley.EDU (AIX 3.1/UCB 5.61/4.03(?))
          id AA14657; Sat, 11 Jul 92 17:06:42 -0700
Received: by jungle.Berkeley.EDU (AIX 3.1/UCB 5.61/4.03(?))
          id AA34793; Sat, 11 Jul 92 17:06:37 -0700
Date: Sat, 11 Jul 92 17:06:37 -0700
From: elsae@jungle.Berkeley.EDU (C. Elsasser)
Message-Id: <9207120006.AA34793@jungle.Berkeley.EDU>
To: alberto@savanna.Berkeley.EDU

      subroutine parcor(icore,r,cdc,ac,gc,bc)
c
c     calculation of partial core density:
c     - outside: cdc
c     - inside : r**2 * ( ac*j_0(gc*r) + bc ) 
c     with continuity in the density and the two first derivatives
c
      implicit real*8 (a-h,o-z)
c
      dimension r(*), cdc(*)
c
      pi = 4.d0*atan(1.d0)
c
c     calculate derivatives of density by parabola fit:
c     f(x) = aa + bb * x + cc * x**2
c
      x1 = r(icore-1)
      x2 = r(icore)
      x3 = r(icore+1)
c
      f1 = cdc(icore-1)
      f2 = cdc(icore)
      f3 = cdc(icore+1)
c
      print*, ' data points : '
      print*, ' icore-1 : ', x1, f1
      print*, ' icore   : ', x2, f2
      print*, ' icore+1 : ', x3, f3
c
      cc = (f3-f1)/(x3-x1)/(x3-x2) + (f2-f1)/(x2-x1)/(x2-x3)
      bbnum = (f3-f1)/(x3*x3-x1*x1) - (f2-f1)/(x2*x2-x1*x1)
      bbden = 1.d0/(x3+x1) - 1.d0/(x2+x1)
      bb = bbnum/bbden
      aa = f1 - bb*x1 - cc*x1*x1
c
      den = f2
      denp = bb + 2.d0*cc*x2
      denpp = 2.d0*cc
      rc = x2
c
      print*, ' parabola fit: '
      print*, ' icore-1 : ', x1, aa + bb*x1 + cc*x1*x1
      print*, ' icore   : ', x2, aa + bb*x2 + cc*x2*x2
      print*, ' icore+1 : ', x3, aa + bb*x3 + cc*x3*x3
c
      ttnum = 0.5*rc*denp - den
      ttden = rc*denp - 0.5*rc*rc*denpp - den
      tt = ttnum/ttden
c 
      xold = 2.5d0
      do 10 i=1,50
        xnew = pi + atan(xold/(tt*xold*xold+1.d0))
        print*, i, xnew
        if (abs(xnew-xold) .lt. 1.d-4) go to 20
        xold = xnew
   10 continue
      call ext(830)
   20 continue
      gc = xnew/rc
      xx = sin(xnew)/gc - rc*cos(xnew) + 0.5d0*rc*xnew*sin(xnew)
      ac = (den-0.5d0*rc*rc*denpp)/rc/xx
      bc = 0.5d0*(denpp - 2.d0*ac*cos(xnew) + ac*xnew*sin(xnew))
c
      do 30 i=1,icore
        cdc(i) = ac*r(i)*sin(gc*r(i))/gc + bc*r(i)*r(i)
   30 continue
c
      write(6,9050) rc, ac, bc, gc
 9050 format(//'core correction used',
     1        /'pseudo core inside r =',f6.3,
     2        /' ac ='f8.3,' bc =',f8.3,' gc =',f8.3,/)
c
      return
c
      end
