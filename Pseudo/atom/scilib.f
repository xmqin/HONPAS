c
      subroutine saxpy(n,sa,sx,incx,sy,incy)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      double precision sa
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (sa .eq. 0.0D0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         sy(iy) = sy(iy) + sa*sx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,4)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sy(i) = sy(i) + sa*sx(i)
   30 continue
      if (n .lt. 4) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 4
         sy(i) = sy(i) + sa*sx(i)
         sy(i+1) = sy(i+1) + sa*sx(i+1)
         sy(i+2) = sy(i+2) + sa*sx(i+2)
         sy(i+3) = sy(i+3) + sa*sx(i+3)
   50 continue
c
      return
c
      end
c
      double precision function sdot(n,sx,incx,sy,incy)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      double precision stemp
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      stemp = 0.0D0
      sdot = 0.0D0
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         stemp = stemp + sx(ix)*sy(iy)
         ix = ix + incx
         iy = iy + incy
   10 continue
      sdot = stemp
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         stemp = stemp + sx(i)*sy(i)
   30 continue
      if (n .lt. 5) go to 60
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +
     &           sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 continue
   60 continue
      sdot = stemp
c
      return
c
      end
c
      subroutine scopy(n,sx,incx,sy,incy)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         sy(iy) = sx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,7)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sy(i) = sx(i)
   30 continue
      if (n .lt. 7) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 7
         sy(i) = sx(i)
         sy(i+1) = sx(i+1)
         sy(i+2) = sx(i+2)
         sy(i+3) = sx(i+3)
         sy(i+4) = sx(i+4)
         sy(i+5) = sx(i+5)
         sy(i+6) = sx(i+6)
   50 continue
c
      return
c
      end
c
      subroutine sscal(n,sa,sx,incx)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      double precision sa
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      integer i, m, mp1, nincx
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (incx .eq. 1) go to 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      nincx = n*incx
      do 10 i = 1, nincx, incx
         sx(i) = sa*sx(i)
   10 continue
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sx(i) = sa*sx(i)
   30 continue
      if (n .lt. 5) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         sx(i) = sa*sx(i)
         sx(i+1) = sa*sx(i+1)
         sx(i+2) = sa*sx(i+2)
         sx(i+3) = sa*sx(i+3)
         sx(i+4) = sa*sx(i+4)
   50 continue
c
      return
c
      end
      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)
c
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      double precision d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      double precision u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,epslon
      integer ind(m)
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        m11 specifies the lower boundary index for the desired
c          eigenvalues.
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        lb and ub define an interval containing exactly the desired
c          eigenvalues.
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection impossible.
c
c        rv4 and rv5 are temporary storage arrays.
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0d0
c     .......... look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0d0
         if (i .ne. n) u = dabs(e(i+1))
         xu = dmin1(d(i)-(x1+u),xu)
         x0 = dmax1(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         tst1 = dabs(d(i)) + dabs(d(i-1))
         tst2 = tst1 + dabs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c
      x1 = n
      x1 = x1 * epslon(dmax1(dabs(xu),dabs(x0)))
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     .......... determine an interval containing exactly
c                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5d0
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = dabs(e(q+1))
         v = e2(q+1)
  110    xu = dmin1(d(q)-(x1+u),xu)
         x0 = dmax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = epslon(dmax1(dabs(xu),dabs(x0)))
      if (eps1 .le. 0.0d0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = dmax1(t1,xu-x1)
      ub = dmin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5d0
         if ((x0 - xu) .le. dabs(eps1)) go to 420
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = dabs(e(i)) / epslon(1.0d0)
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- interval cannot be found containing
c                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      double precision d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,
     x       pythag
      integer ind(m)
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here.
c
c        m is the number of specified eigenvalues.
c
c        w contains the m eigenvalues in ascending or descending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output
c
c        all input arrays are unaltered.
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0d0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0d0) go to 140
  120 continue
c     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     .......... check for isolated root ..........
         xu = 1.0d0
         if (p .ne. q) go to 490
         rv6(p) = 1.0d0
         go to 870
  490    norm = dabs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0d-3 * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / dsqrt(uk)
         s = p
  505    group = 0
         go to 520
c     .......... look for close or coincident roots ..........
  510    if (dabs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (dabs(e(i)) .lt. dabs(u)) go to 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
      end
C  LINPACKS.FOR  02 December 1989
C  The single precision version of LINPACK.
C
C***********************************************************************
C
      subroutine sgeco(a,lda,n,ipvt,rcond,z)
C
C     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SGEFA
C     BLAS SAXPY,SDOT,SSCAL,SASUM
C     FORTRAN ABS,AMAX1,SIGN
C
C     INTERNAL VARIABLES
C
C
C
C     COMPUTE 1-NORM OF A
C
C     .. Scalar Arguments ..
      double precision rcond
      integer lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*), z(*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision anorm, ek, s, sm, t, wk, wkm, ynorm
      integer info, j, k, kb, kp1, l
C     ..
C     .. External Functions ..
      double precision sasum, sdot
      external sasum, sdot
C     ..
C     .. External Subroutines ..
      external saxpy, sgefa, sscal
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, dmax1, sign
C     ..
      anorm = 0.0D0
      do 10 j = 1, n
         anorm = dmax1(anorm,sasum(n,a(1,j),1))
   10 continue
C
C     FACTOR
C
      call sgefa(a,lda,n,ipvt,info)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      ek = 1.0D0
      do 20 j = 1, n
         z(j) = 0.0D0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0D0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30
         s = abs(a(k,k))/abs(ek-z(k))
         call sscal(n,s,z,1)
         ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (a(k,k) .eq. 0.0D0) go to 40
         wk = wk/a(k,k)
         wkm = wkm/a(k,k)
c
         go to 50
c
   40    continue
         wk = 1.0D0
         wkm = 1.0D0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
         do 60 j = kp1, n
            sm = sm + abs(z(j)+wkm*a(k,j))
            z(j) = z(j) + wk*a(k,j)
            s = s + abs(z(j))
   60    continue
         if (s .ge. sm) go to 80
         t = wkm - wk
         wk = wkm
         do 70 j = kp1, n
            z(j) = z(j) + t*a(k,j)
   70    continue
   80    continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0D0) go to 110
         s = 1.0D0/abs(z(k))
         call sscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
C
      ynorm = 1.0D0
C
C     SOLVE L*V = Y
C
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0D0) go to 130
         s = 1.0D0/abs(z(k))
         call sscal(n,s,z,1)
         ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
C
C     SOLVE  U*Z = V
C
      do 160 kb = 1, n
         k = n + 1 - kb
         if (abs(z(k)) .le. abs(a(k,k))) go to 150
         s = abs(a(k,k))/abs(z(k))
         call sscal(n,s,z,1)
         ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0D0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0D0) z(k) = 1.0D0
         t = -z(k)
         call saxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
C     MAKE ZNORM = 1.0
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
C
      if (anorm .ne. 0.0D0) rcond = ynorm/anorm
      if (anorm .eq. 0.0D0) rcond = 0.0D0
c
      return
c
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
      double precision function sasum(n,sx,incx)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      double precision stemp
      integer i, m, mp1, nincx
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, mod
C     ..
      sasum = 0.0D0
      stemp = 0.0D0
      if (n .le. 0) return
      if (incx .eq. 1) go to 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      nincx = n*incx
      do 10 i = 1, nincx, incx
         stemp = stemp + abs(sx(i))
   10 continue
      sasum = stemp
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,6)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         stemp = stemp + abs(sx(i))
   30 continue
      if (n .lt. 6) go to 60
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 6
         stemp = stemp + abs(sx(i)) + abs(sx(i+1)) + abs(sx(i+2)) +
     &           abs(sx(i+3)) + abs(sx(i+4)) + abs(sx(i+5))
   50 continue
   60 continue
      sasum = stemp
c
      return
c
      end
c
      subroutine sgefa(a,lda,n,ipvt,info)
C
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,ISAMAX
C
C     INTERNAL VARIABLES
C
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C     .. Scalar Arguments ..
      integer info, lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision t
      integer j, k, kp1, l, nm1
C     ..
C     .. External Functions ..
      integer isamax
      external isamax
C     ..
C     .. External Subroutines ..
      external saxpy, sscal
C     ..
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
C
C        FIND L = PIVOT INDEX
C
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         if (a(l,k) .eq. 0.0D0) go to 40
C
C           INTERCHANGE IF NECESSARY
C
         if (l .eq. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
   10    continue
C
C           COMPUTE MULTIPLIERS
C
         t = -1.0D0/a(k,k)
         call sscal(n-k,t,a(k+1,k),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
         do 30 j = kp1, n
            t = a(l,j)
            if (l .eq. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
   20       continue
            call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30    continue
c
         go to 50
c
   40    continue
         info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0D0) info = n
c
      return
c
      end
c
      integer function isamax(n,sx,incx)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      double precision smax
      integer i, ix
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      isamax = 0
      if (n .lt. 1) return
      isamax = 1
      if (n .eq. 1) return
      if (incx .eq. 1) go to 30
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 20 i = 2, n
         if (abs(sx(ix)) .le. smax) go to 10
         isamax = i
         smax = abs(sx(ix))
   10    continue
         ix = ix + incx
   20 continue
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   30 continue
      smax = abs(sx(1))
      do 40 i = 2, n
         if (abs(sx(i)) .le. smax) go to 40
         isamax = i
         smax = abs(sx(i))
   40 continue
c
      return
c
      end
c
      subroutine sgesl(a,lda,n,ipvt,b,job)
C
C     SGESL SOLVES THE REAL SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C
C     INTERNAL VARIABLES
C
C
C     .. Scalar Arguments ..
      integer job, lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*), b(*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision t
      integer k, kb, l, nm1
C     ..
C     .. External Functions ..
      double precision sdot
      external sdot
C     ..
C     .. External Subroutines ..
      external saxpy
C     ..
      nm1 = n - 1
      if (job .ne. 0) go to 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
      if (nm1 .lt. 1) go to 30
      do 20 k = 1, nm1
         l = ipvt(k)
         t = b(l)
         if (l .eq. k) go to 10
         b(l) = b(k)
         b(k) = t
   10    continue
         call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20 continue
   30 continue
C
C        NOW SOLVE  U*X = Y
C
      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call saxpy(k-1,t,a(1,k),1,b(1),1)
   40 continue
c
      go to 100
c
   50 continue
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
      do 60 k = 1, n
         t = sdot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k)-t)/a(k,k)
   60 continue
C
C        NOW SOLVE TRANS(L)*X = Y
C
      if (nm1 .lt. 1) go to 90
      do 80 kb = 1, nm1
         k = n - kb
         b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l .eq. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
   70    continue
   80 continue
   90 continue
  100 continue
c
      return
c
      end
c
      integer function isrchfgt(n,array,inc,target)
C
C     Returns the location of the first element in a real ARRAY
C          that is greater than a real TARGET.  Returns N+1 if
C          TARGET is not found and 0 if N < 1.
C     Martin J. McBride.  7/17/85.
C     General Electric CRD, Information System Operation.
C
c
C     .. Scalar Arguments ..
      double precision target
      integer inc, n
C     ..
C     .. Array Arguments ..
      double precision array(*)
C     ..
C     .. Local Scalars ..
      integer i, ix
C     ..
      isrchfgt = 0
      if (n .lt. 1) return
      ix = 1
      if (inc .lt. 0) ix = (-inc)*(n-1) + 1
      do 10 i = 1, n
         if (array(ix) .gt. target) go to 20
         ix = ix + inc
   10 continue
   20 continue
      isrchfgt = i
c
      return
c
      end
