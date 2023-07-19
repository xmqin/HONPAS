c
c THIS CODE COMES FROM THE BERKELEY PLANE-WAVE PROGRAM
c (Maintained by Alberto Garcia)
c
c It is the responsibility of the user to make sure that the
c only prime factors of the fft lenghts are 2, 3, 5.
c
      subroutine c1d_fft(n1,f,mode)
c
c      Computes complex 1d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1               grid size. 
c
c I/O  f(i)                i=1,n1
c                          array to be fourier transformed
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
      implicit none
c
C     .. Scalar Arguments ..
      integer mode, n1
C     ..
C     .. Array Arguments ..
      complex f(*)
c
C     ..
      integer nmax
      parameter (nmax=20000)
c
c I/O  wrk(i)              i=1,4*max(n1)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex wrk(4*nmax+15)
c
      real factor
      integer i, j, k, np

      external cfftb, cfftf, cffti

      np = n1 + 1

      if (n1 .gt. nmax) STOP 'FFT'
c
c     Transform the first dimension
c
      call cffti(n1,wrk(np))
            do 90 k = 1, n1
               wrk(k) = f(k)
   90       continue
            if (mode .gt. 0) call cfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k) = wrk(k)
  100       continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.E0/n1
      do 150 i = 1, n1
               f(i) = factor*f(i)
  150 continue
c
      return
c
      end
c
      subroutine c2d_fft(n1,n2,f,mode,np1)
c
c      Computes complex 2d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1,n2            grid size
c Inp  np1              leading physical dimensions of array f
c
c I/O  f(i,j)              i=1,n1,j=1,n2
c                          array to be fourier transformed
c
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
      implicit none

C     .. Scalar Arguments ..
      integer mode, n1, n2, np1
C     ..
C     .. Array Arguments ..
      complex f(np1,*)
c
C     ..
      integer nmax
      parameter (nmax=1024)
c
c I/O  wrk(i)              i=1,4*max(n1,n2)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex wrk(4*nmax+15)
c
      real factor
      integer i, j, k, np

      external cfftb, cfftf, cffti

      intrinsic max

      np = max(n1,n2) + 1
c
      if (np-1 .gt. nmax) STOP 'FFT'
c
c     Transform the second dimension
c
      call cffti(n2,wrk(np))

         do 70 j = 1, n1
            do 50 k = 1, n2
               wrk(k) = f(j,k)
   50       continue
            if (mode .gt. 0) call cfftf(n2,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n2,wrk,wrk(np))
            do 60 k = 1, n2
               f(j,k) = wrk(k)
   60       continue
   70    continue
c
c     Transform the first dimension
c
      if (n1 .ne. n2) call cffti(n1,wrk(np))
      do 120 i = 1, n2
            do 90 k = 1, n1
               wrk(k) = f(k,i)
   90       continue
            if (mode .gt. 0) call cfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k,i) = wrk(k)
  100       continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.E0/(n1*n2)
      do 150 i = 1, n1
         do 140 j = 1, n2
               f(i,j) = factor*f(i,j)
  140    continue
  150 continue
c
      return
c
      end
c
      subroutine c3d_fft(n1,n2,n3,f,mode,np1,np2)
c
c      Computes complex 3d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1,n2,n3            grid size
c Inp  np1, np2            leading physical dimensions of array f
c
c I/O  f(i,j,k)            i=1,n1,j=1,n2,k=1,n3.
c                          array to be fourier transformed
c
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
C     .. Scalar Arguments ..
      integer mode, n1, n2, n3, np1, np2
C     ..
C     .. Array Arguments ..
      complex f(np1,np2,*)
c
C     ..
      integer nmax
      parameter (nmax=1024)
c
c I/O  wrk(i)              i=1,4*max(n1,n2,n3)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex wrk(4*nmax+15)
c
      real factor
      integer i, j, k, np

      external cfftb, cfftf, cffti

      intrinsic max

      np = max(n1,n2,n3) + 1
c
      if (np-1 .gt. nmax) STOP 'FFT'
c
c     First call to initialize the work arrays
c
      call cffti(n3,wrk(np))
c
c     Transform the third dimension
c
      do 40 i = 1, n1
         do 30 j = 1, n2
            do 10 k = 1, n3
               wrk(k) = f(i,j,k)
   10       continue
            if (mode .gt. 0) call cfftf(n3,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n3,wrk,wrk(np))
            do 20 k = 1, n3
               f(i,j,k) = wrk(k)
   20       continue
   30    continue
   40 continue
c
c     Transform the second dimension
c
      if (n2 .ne. n3) call cffti(n2,wrk(np))
      do 80 i = 1, n3
         do 70 j = 1, n1
            do 50 k = 1, n2
               wrk(k) = f(j,k,i)
   50       continue
            if (mode .gt. 0) call cfftf(n2,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n2,wrk,wrk(np))
            do 60 k = 1, n2
               f(j,k,i) = wrk(k)
   60       continue
   70    continue
   80 continue
c
c     Transform the first dimension
c
      if (n1 .ne. n2) call cffti(n1,wrk(np))
      do 120 i = 1, n2
         do 110 j = 1, n3
            do 90 k = 1, n1
               wrk(k) = f(k,i,j)
   90       continue
            if (mode .gt. 0) call cfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call cfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k,i,j) = wrk(k)
  100       continue
  110    continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.E0/(n1*n2*n3)
      do 150 i = 1, n1
         do 140 j = 1, n2
            do 130 k = 1, n3
               f(i,j,k) = factor*f(i,j,k)
  130       continue
  140    continue
  150 continue
c
      return
c
      end
c fftpack    from ulib                                     06/05/81
c
cAG: Names begin with z or d to emphasize that this is a double precision
c    implementation of fftpack.
c
c purpose  this package consists of programs which perform fast fourier
c          transforms for both complex and real periodic sequences and
c          certain other symmetric sequences that are listed below.
c
c
c            cffti     initialize cfftf and cfftb
c            cfftf     forward transform of a complex periodic sequence
c            cfftb     unnormalized inverse of cfftf
c
c
c     subroutine cffti(n,wsave)
c
c     ******************************************************************
c
c     subroutine cffti initializes the array wsave which is used in
c     both cfftf and cfftb. the prime factorization of n together with
c     a tabulation of the trigonometric functions are computed and
c     stored in wsave.
c
c     input parameter
c
c     n       the length of the sequence to be transformed
c
c     output parameter
c
c     wsave   a work array which must be dimensioned at least 4*n+15
c             the same work array can be used for both cfftf and cfftb
c             as long as n remains unchanged. different wsave arrays
c             are required for different values of n. the contents of
c             wsave must not be changed between calls of cfftf or cfftb.
c
c     ******************************************************************
c
c     subroutine cfftf(n,c,wsave)
c
c     ******************************************************************
c
c     subroutine cfftf computes the forward complex discrete fourier
c     transform (the fourier analysis). equivalently , cfftf computes
c     the fourier coefficients of a complex periodic sequence.
c     the transform is defined below at output parameter c.
c
c     the transform is not normalized. to obtain a normalized transform
c     the output must be divided by n. otherwise a call of cfftf
c     followed by a call of cfftb will multiply the sequence by n.
c
c     the array wsave which is used by subroutine cfftf must be
c     initialized by calling subroutine cffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes. n
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls cfftf. the wsave array must be
c             initialized by calling subroutine cffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by cfftf and cfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(-i*j*k*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine cfftf or cfftb
c
c     ******************************************************************
c
c     subroutine cfftb(n,c,wsave)
c
c     ******************************************************************
c
c     subroutine cfftb computes the backward complex discrete fourier
c     transform (the fourier synthesis). equivalently , cfftb computes
c     a complex periodic sequence from its fourier coefficients.
c     the transform is defined below at output parameter c.
c
c     a call of cfftf followed by a call of cfftb will multiply the
c     sequence by n.
c
c     the array wsave which is used by subroutine cfftb must be
c     initialized by calling subroutine cffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes.
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls cfftb. the wsave array must be
c             initialized by calling subroutine cffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by cfftf and cfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(i*j*k*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine cfftf or cfftb
c___________________________________________________________________
c
      subroutine cffti(n,wsave)
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      real wsave(*)
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external cffti1
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cffti1(n,wsave(iw1),wsave(iw2))
c
      return
c
      end
c____________________________________________________________________
c
      subroutine cfftf(n,c,wsave)
c
C     .. Array Arguments ..
      real c(*), wsave(*)
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external cfftf1
C     ..
C     .. Scalar Arguments ..
      integer n
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cfftf1(n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
c
      end
c____________________________________________________________________
c
      subroutine cfftb(n,c,wsave)
c
C     .. Array Arguments ..
C
c
      real c(*), wsave(*)
C     ..
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external cfftb1
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call cfftb1(n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
c
      end
c_____________________________________________________________
c
      subroutine cffti1(n,wa,ifac)
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      real wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      real arg, argh, argld, fi, tpi
      integer i, i1, ib, ido, idot, ii, ip, ipm, j, k1, l1, l2, ld, nf,
     1        nl, nq, nr, ntry
C     ..
C     .. Local Arrays ..
      integer ntryh(4)
C     ..
C     .. Intrinsic Functions ..
      intrinsic cos, float, sin
C     ..
C     .. Data statements ..
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4)/3, 4, 2, 5/
      data tpi/6.28318530717958647692528676655900577e0/
C     ..
c
      nl = n
      nf = 0
      j = 0
c
   10 continue
      j = j + 1
      if (j .le. 4) ntry = ntryh(j)
      if (j .gt. 4) ntry = ntry + 2
   20 continue
      nq = nl/ntry
      nr = nl - ntry*nq
      if (nr .ne. 0) go to 10
c
      nf = nf + 1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 40
      if (nf .eq. 1) go to 40
      do 30 i = 2, nf
         ib = nf - i + 2
         ifac(ib+2) = ifac(ib+1)
   30 continue
      ifac(3) = 2
c
   40 continue
      if (nl .ne. 1) go to 20
c
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/float(n)
      i = 2
      l1 = 1
      do 70 k1 = 1, nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido + ido + 2
         ipm = ip - 1
c
         do 60 j = 1, ipm
            i1 = i
            wa(i-1) = 1.e0
            wa(i) = 0.e0
            ld = ld + l1
            fi = 0.e0
            argld = float(ld)*argh
            do 50 ii = 4, idot, 2
               i = i + 2
               fi = fi + 1.e0
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
   50       continue
            if (ip .le. 5) go to 60
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
   60    continue
c
         l1 = l2
   70 continue
c
      return
c
      end
c________________________________________________________________
c
      subroutine cfftf1(n,c,ch,wa,ifac)
c Dummy dimension compatible with semantic analyzer.
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      real c(*), ch(*), wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      integer i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2, n2,
     1        na, nac, nf
C     ..
C     .. External Subroutines ..
      external spssf, spssf2, spssf3, spssf4, spssf5
C     ..
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 160 k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido + ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 30
         ix2 = iw + idot
         ix3 = ix2 + idot
         if (na .ne. 0) go to 10
         call spssf4(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
c
         go to 20
c
   10    continue
         call spssf4(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
   20    continue
         na = 1 - na
c
         go to 150
c
   30    continue
         if (ip .ne. 2) go to 60
         if (na .ne. 0) go to 40
         call spssf2(idot,l1,c,ch,wa(iw))
c
         go to 50
c
   40    continue
         call spssf2(idot,l1,ch,c,wa(iw))
   50    continue
         na = 1 - na
c
         go to 150
c
   60    continue
         if (ip .ne. 3) go to 90
         ix2 = iw + idot
         if (na .ne. 0) go to 70
         call spssf3(idot,l1,c,ch,wa(iw),wa(ix2))
c
         go to 80
c
   70    continue
         call spssf3(idot,l1,ch,c,wa(iw),wa(ix2))
   80    continue
         na = 1 - na
c
         go to 150
c
   90    continue
         if (ip .ne. 5) go to 120
         ix2 = iw + idot
         ix3 = ix2 + idot
         ix4 = ix3 + idot
         if (na .ne. 0) go to 100
         call spssf5(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
c
         go to 110
c
  100    continue
         call spssf5(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    continue
         na = 1 - na
c
         go to 150
c
  120    continue
         if (na .ne. 0) go to 130
         call spssf(nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
c
         go to 140
c
  130    continue
         call spssf(nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  140    continue
         if (nac .ne. 0) na = 1 - na
c
  150    continue
         l1 = l2
         iw = iw + (ip-1)*idot
  160 continue
      if (na .eq. 0) return
c
      n2 = n + n
      do 170 i = 1, n2
         c(i) = ch(i)
  170 continue
c
      return
c
      end
c___________________________________________________________________
c
      subroutine cfftb1(n,c,ch,wa,ifac)
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      real c(*), ch(*), wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      integer i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2, n2,
     1        na, nac, nf
C     ..
C     .. External Subroutines ..
      external spssb, spssb2, spssb3, spssb4, spssb5
C     ..
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 160 k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido + ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 30
         ix2 = iw + idot
         ix3 = ix2 + idot
         if (na .ne. 0) go to 10
         call spssb4(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
c
         go to 20
c
   10    continue
         call spssb4(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
   20    continue
         na = 1 - na
c
         go to 150
c
   30    continue
         if (ip .ne. 2) go to 60
         if (na .ne. 0) go to 40
         call spssb2(idot,l1,c,ch,wa(iw))
c
         go to 50
c
   40    continue
         call spssb2(idot,l1,ch,c,wa(iw))
   50    continue
         na = 1 - na
c
         go to 150
c
   60    continue
         if (ip .ne. 3) go to 90
         ix2 = iw + idot
         if (na .ne. 0) go to 70
         call spssb3(idot,l1,c,ch,wa(iw),wa(ix2))
c
         go to 80
c
   70    continue
         call spssb3(idot,l1,ch,c,wa(iw),wa(ix2))
   80    continue
         na = 1 - na
c
         go to 150
c
   90    continue
         if (ip .ne. 5) go to 120
         ix2 = iw + idot
         ix3 = ix2 + idot
         ix4 = ix3 + idot
         if (na .ne. 0) go to 100
         call spssb5(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
c
         go to 110
c
  100    continue
         call spssb5(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    continue
         na = 1 - na
c
         go to 150
c
  120    continue
         if (na .ne. 0) go to 130
         call spssb(nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
c
         go to 140
c
  130    continue
         call spssb(nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  140    continue
         if (nac .ne. 0) na = 1 - na
c
  150    continue
         l1 = l2
         iw = iw + (ip-1)*idot
  160 continue
      if (na .eq. 0) return
c
      n2 = n + n
      do 170 i = 1, n2
         c(i) = ch(i)
  170 continue
c
      return
c
      end
c_____________________________________________________________________
c
      subroutine spssb(nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
C     .. Scalar Arguments ..
      integer idl1, ido, ip, l1, nac
C     ..
C     .. Array Arguments ..
      real c1(ido,l1,ip), c2(idl1,ip), cc(ido,ip,l1),
     1                 ch(ido,l1,ip), ch2(idl1,ip), wa(*)
C     ..
C     .. Local Scalars ..
      real wai, war
      integer i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph,
     1        j, jc, k, l, lc, nt
C     ..
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip + 2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 60
      do 30 j = 2, ipph
         jc = ipp2 - j
         do 20 k = 1, l1
            do 10 i = 1, ido
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   10       continue
   20    continue
   30 continue
c
      do 50 k = 1, l1
         do 40 i = 1, ido
            ch(i,k,1) = cc(i,1,k)
   40    continue
   50 continue
c
      go to 120
c
   60 continue
      do 90 j = 2, ipph
         jc = ipp2 - j
         do 80 i = 1, ido
            do 70 k = 1, l1
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   70       continue
   80    continue
   90 continue
c
      do 110 i = 1, ido
         do 100 k = 1, l1
            ch(i,k,1) = cc(i,1,k)
  100    continue
  110 continue
c
  120 continue
      idl = 2 - ido
      inc = 0
      do 160 l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do 130 ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
  130    continue
         idlj = idl
         inc = inc + ido
         do 150 j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp) idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 140 ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) + wai*ch2(ik,jc)
  140       continue
  150    continue
  160 continue
c
      do 180 j = 2, ipph
         do 170 ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
  170    continue
  180 continue
c
      do 200 j = 2, ipph
         jc = ipp2 - j
         do 190 ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
  190    continue
  200 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 210 ik = 1, idl1
         c2(ik,1) = ch2(ik,1)
  210 continue
c
      do 230 j = 2, ip
         do 220 k = 1, l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  220    continue
  230 continue
c
      if (idot .gt. l1) go to 270
      idij = 0
      do 260 j = 2, ip
         idij = idij + 2
         do 250 i = 4, ido, 2
            idij = idij + 2
            do 240 k = 1, l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) + wa(idij)*ch(i-1,k,j)
  240       continue
  250    continue
  260 continue
c
      return
c
  270 continue
      idj = 2 - ido
      do 300 j = 2, ip
         idj = idj + ido
         do 290 k = 1, l1
            idij = idj
            do 280 i = 4, ido, 2
               idij = idij + 2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) + wa(idij)*ch(i-1,k,j)
  280       continue
  290    continue
  300 continue
c
      return
c
      end
c
      subroutine spssb2(ido,l1,cc,ch,wa1)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,2,l1), ch(ido,l1,2), wa1(*)
C     ..
C     .. Local Scalars ..
      real ti2, tr2
      integer i, k
C     ..
      if (ido .gt. 2) go to 20
      do 10 k = 1, l1
         ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
         ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
         ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
            tr2 = cc(i-1,1,k) - cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
            ti2 = cc(i,1,k) - cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2 + wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2 - wa1(i)*ti2
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssb3(ido,l1,cc,ch,wa1,wa2)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui,
     1                 taur, ti2, tr2
      integer i, k
C     ..
C     .. Data statements ..
      data taur/-0.5E0/
      data taui/0.86602540378443864676372317075293618E0/
C     ..
c
c     one half sqrt(3) = .866025.....  .
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         tr2 = cc(1,2,k) + cc(1,3,k)
         cr2 = cc(1,1,k) + taur*tr2
         ch(1,k,1) = cc(1,1,k) + tr2
         ti2 = cc(2,2,k) + cc(2,3,k)
         ci2 = cc(2,1,k) + taur*ti2
         ch(2,k,1) = cc(2,1,k) + ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2 - ci3
         ch(1,k,3) = cr2 + ci3
         ch(2,k,2) = ci2 + cr3
         ch(2,k,3) = ci2 - cr3
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            tr2 = cc(i-1,2,k) + cc(i-1,3,k)
            cr2 = cc(i-1,1,k) + taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k) + tr2
            ti2 = cc(i,2,k) + cc(i,3,k)
            ci2 = cc(i,1,k) + taur*ti2
            ch(i,k,1) = cc(i,1,k) + ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssb4(ido,l1,cc,ch,wa1,wa2,wa3)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1                 wa3(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4,
     1                 tr1, tr2, tr3, tr4
      integer i, k
C     ..
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti1 = cc(2,1,k) - cc(2,3,k)
         ti2 = cc(2,1,k) + cc(2,3,k)
         tr4 = cc(2,4,k) - cc(2,2,k)
         ti3 = cc(2,2,k) + cc(2,4,k)
         tr1 = cc(1,1,k) - cc(1,3,k)
         tr2 = cc(1,1,k) + cc(1,3,k)
         ti4 = cc(1,2,k) - cc(1,4,k)
         tr3 = cc(1,2,k) + cc(1,4,k)
         ch(1,k,1) = tr2 + tr3
         ch(1,k,3) = tr2 - tr3
         ch(2,k,1) = ti2 + ti3
         ch(2,k,3) = ti2 - ti3
         ch(1,k,2) = tr1 + tr4
         ch(1,k,4) = tr1 - tr4
         ch(2,k,2) = ti1 + ti4
         ch(2,k,4) = ti1 - ti4
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti1 = cc(i,1,k) - cc(i,3,k)
            ti2 = cc(i,1,k) + cc(i,3,k)
            ti3 = cc(i,2,k) + cc(i,4,k)
            tr4 = cc(i,4,k) - cc(i,2,k)
            tr1 = cc(i-1,1,k) - cc(i-1,3,k)
            tr2 = cc(i-1,1,k) + cc(i-1,3,k)
            ti4 = cc(i-1,2,k) - cc(i-1,4,k)
            tr3 = cc(i-1,2,k) + cc(i-1,4,k)
            ch(i-1,k,1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i,k,1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1,k,2) = wa1(i-1)*cr2 - wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2 + wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3 - wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3 + wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4 - wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4 + wa3(i)*cr4
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssb5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1                 wa3(*), wa4(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3,
     1                 di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2,
     2                 ti3, ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      integer i, k
C     ..
C     .. Data statements ..
      data tr11/0.30901699437494742410229341718281906E0/
      data ti11/0.95105651629515357211643933337938214E0/
      data tr12/-0.80901699437494742410229341718281906E0/
      data ti12/0.58778525229247312916870595463907277E0/
C     ..
c
c     sin(pi/10) = .30901699....    .
c     cos(pi/10) = .95105651....    .
c     sin(pi/5 ) = .58778525....    .
c     cos(pi/5 ) = .80901699....    .
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti5 = cc(2,2,k) - cc(2,5,k)
         ti2 = cc(2,2,k) + cc(2,5,k)
         ti4 = cc(2,3,k) - cc(2,4,k)
         ti3 = cc(2,3,k) + cc(2,4,k)
         tr5 = cc(1,2,k) - cc(1,5,k)
         tr2 = cc(1,2,k) + cc(1,5,k)
         tr4 = cc(1,3,k) - cc(1,4,k)
         tr3 = cc(1,3,k) + cc(1,4,k)
         ch(1,k,1) = cc(1,1,k) + tr2 + tr3
         ch(2,k,1) = cc(2,1,k) + ti2 + ti3
         cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
         ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
         cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
         ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         ch(1,k,2) = cr2 - ci5
         ch(1,k,5) = cr2 + ci5
         ch(2,k,2) = ci2 + cr5
         ch(2,k,3) = ci3 + cr4
         ch(1,k,3) = cr3 - ci4
         ch(1,k,4) = cr3 + ci4
         ch(2,k,4) = ci3 - cr4
         ch(2,k,5) = ci2 - cr5
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti5 = cc(i,2,k) - cc(i,5,k)
            ti2 = cc(i,2,k) + cc(i,5,k)
            ti4 = cc(i,3,k) - cc(i,4,k)
            ti3 = cc(i,3,k) + cc(i,4,k)
            tr5 = cc(i-1,2,k) - cc(i-1,5,k)
            tr2 = cc(i-1,2,k) + cc(i-1,5,k)
            tr4 = cc(i-1,3,k) - cc(i-1,4,k)
            tr3 = cc(i-1,3,k) + cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
            ch(i,k,1) = cc(i,1,k) + ti2 + ti3
            cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4 - wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4 + wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5 - wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5 + wa4(i)*dr5
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssf(nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
C     .. Scalar Arguments ..
      integer idl1, ido, ip, l1, nac
C     ..
C     .. Array Arguments ..
      real c1(ido,l1,ip), c2(idl1,ip), cc(ido,ip,l1),
     1                 ch(ido,l1,ip), ch2(idl1,ip), wa(*)
C     ..
C     .. Local Scalars ..
      real wai, war
      integer i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph,
     1        j, jc, k, l, lc, nt
C     ..
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip + 2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 60
      do 30 j = 2, ipph
         jc = ipp2 - j
         do 20 k = 1, l1
            do 10 i = 1, ido
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   10       continue
   20    continue
   30 continue
c
      do 50 k = 1, l1
         do 40 i = 1, ido
            ch(i,k,1) = cc(i,1,k)
   40    continue
   50 continue
c
      go to 120
c
   60 continue
      do 90 j = 2, ipph
         jc = ipp2 - j
         do 80 i = 1, ido
            do 70 k = 1, l1
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   70       continue
   80    continue
   90 continue
c
      do 110 i = 1, ido
         do 100 k = 1, l1
            ch(i,k,1) = cc(i,1,k)
  100    continue
  110 continue
c
  120 continue
      idl = 2 - ido
      inc = 0
      do 160 l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do 130 ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  130    continue
         idlj = idl
         inc = inc + ido
         do 150 j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp) idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 140 ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) - wai*ch2(ik,jc)
  140       continue
  150    continue
  160 continue
c
      do 180 j = 2, ipph
         do 170 ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
  170    continue
  180 continue
c
      do 200 j = 2, ipph
         jc = ipp2 - j
         do 190 ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
  190    continue
  200 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 210 ik = 1, idl1
         c2(ik,1) = ch2(ik,1)
  210 continue
c
      do 230 j = 2, ip
         do 220 k = 1, l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  220    continue
  230 continue
c
      if (idot .gt. l1) go to 270
      idij = 0
      do 260 j = 2, ip
         idij = idij + 2
         do 250 i = 4, ido, 2
            idij = idij + 2
            do 240 k = 1, l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) + wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
  240       continue
  250    continue
  260 continue
c
      return
c
  270 continue
      idj = 2 - ido
      do 300 j = 2, ip
         idj = idj + ido
         do 290 k = 1, l1
            idij = idj
            do 280 i = 4, ido, 2
               idij = idij + 2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) + wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
  280       continue
  290    continue
  300 continue
c
      return
c
      end
c
      subroutine spssf2(ido,l1,cc,ch,wa1)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,2,l1), ch(ido,l1,2), wa1(*)
C     ..
C     .. Local Scalars ..
      real ti2, tr2
      integer i, k
C     ..
      if (ido .gt. 2) go to 20
      do 10 k = 1, l1
         ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
         ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
         ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
            tr2 = cc(i-1,1,k) - cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
            ti2 = cc(i,1,k) - cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2 - wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2 + wa1(i)*ti2
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssf3(ido,l1,cc,ch,wa1,wa2)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui,
     1                 taur, ti2, tr2
      integer i, k
C     ..
C     .. Data statements ..
      data taur/-0.5E0/
      data taui/-0.86602540378443864676372317075293618E0/
C     ..
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         tr2 = cc(1,2,k) + cc(1,3,k)
         cr2 = cc(1,1,k) + taur*tr2
         ch(1,k,1) = cc(1,1,k) + tr2
         ti2 = cc(2,2,k) + cc(2,3,k)
         ci2 = cc(2,1,k) + taur*ti2
         ch(2,k,1) = cc(2,1,k) + ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2 - ci3
         ch(1,k,3) = cr2 + ci3
         ch(2,k,2) = ci2 + cr3
         ch(2,k,3) = ci2 - cr3
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            tr2 = cc(i-1,2,k) + cc(i-1,3,k)
            cr2 = cc(i-1,1,k) + taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k) + tr2
            ti2 = cc(i,2,k) + cc(i,3,k)
            ci2 = cc(i,1,k) + taur*ti2
            ch(i,k,1) = cc(i,1,k) + ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssf4(ido,l1,cc,ch,wa1,wa2,wa3)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1                 wa3(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4,
     1                 tr1, tr2, tr3, tr4
      integer i, k
C     ..
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti1 = cc(2,1,k) - cc(2,3,k)
         ti2 = cc(2,1,k) + cc(2,3,k)
         tr4 = cc(2,2,k) - cc(2,4,k)
         ti3 = cc(2,2,k) + cc(2,4,k)
         tr1 = cc(1,1,k) - cc(1,3,k)
         tr2 = cc(1,1,k) + cc(1,3,k)
         ti4 = cc(1,4,k) - cc(1,2,k)
         tr3 = cc(1,2,k) + cc(1,4,k)
         ch(1,k,1) = tr2 + tr3
         ch(1,k,3) = tr2 - tr3
         ch(2,k,1) = ti2 + ti3
         ch(2,k,3) = ti2 - ti3
         ch(1,k,2) = tr1 + tr4
         ch(1,k,4) = tr1 - tr4
         ch(2,k,2) = ti1 + ti4
         ch(2,k,4) = ti1 - ti4
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti1 = cc(i,1,k) - cc(i,3,k)
            ti2 = cc(i,1,k) + cc(i,3,k)
            ti3 = cc(i,2,k) + cc(i,4,k)
            tr4 = cc(i,2,k) - cc(i,4,k)
            tr1 = cc(i-1,1,k) - cc(i-1,3,k)
            tr2 = cc(i-1,1,k) + cc(i-1,3,k)
            ti4 = cc(i-1,4,k) - cc(i-1,2,k)
            tr3 = cc(i-1,2,k) + cc(i-1,4,k)
            ch(i-1,k,1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i,k,1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1,k,2) = wa1(i-1)*cr2 + wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2 - wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3 + wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3 - wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4 + wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4 - wa3(i)*cr4
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine spssf5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      real cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1                 wa3(*), wa4(*)
C     ..
C     .. Local Scalars ..
      real ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3,
     1                 di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2,
     2                 ti3, ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      integer i, k
C     ..
C     .. Data statements ..
      data tr11/0.30901699437494742410229341718281906E0/
      data ti11/-0.95105651629515357211643933337938214E0/
      data tr12/-0.80901699437494742410229341718281906E0/
      data ti12/-0.58778525229247312916870595463907277E0/
C     ..
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti5 = cc(2,2,k) - cc(2,5,k)
         ti2 = cc(2,2,k) + cc(2,5,k)
         ti4 = cc(2,3,k) - cc(2,4,k)
         ti3 = cc(2,3,k) + cc(2,4,k)
         tr5 = cc(1,2,k) - cc(1,5,k)
         tr2 = cc(1,2,k) + cc(1,5,k)
         tr4 = cc(1,3,k) - cc(1,4,k)
         tr3 = cc(1,3,k) + cc(1,4,k)
         ch(1,k,1) = cc(1,1,k) + tr2 + tr3
         ch(2,k,1) = cc(2,1,k) + ti2 + ti3
         cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
         ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
         cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
         ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         ch(1,k,2) = cr2 - ci5
         ch(1,k,5) = cr2 + ci5
         ch(2,k,2) = ci2 + cr5
         ch(2,k,3) = ci3 + cr4
         ch(1,k,3) = cr3 - ci4
         ch(1,k,4) = cr3 + ci4
         ch(2,k,4) = ci3 - cr4
         ch(2,k,5) = ci2 - cr5
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti5 = cc(i,2,k) - cc(i,5,k)
            ti2 = cc(i,2,k) + cc(i,5,k)
            ti4 = cc(i,3,k) - cc(i,4,k)
            ti3 = cc(i,3,k) + cc(i,4,k)
            tr5 = cc(i-1,2,k) - cc(i-1,5,k)
            tr2 = cc(i-1,2,k) + cc(i-1,5,k)
            tr4 = cc(i-1,3,k) - cc(i-1,4,k)
            tr3 = cc(i-1,3,k) + cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
            ch(i,k,1) = cc(i,1,k) + ti2 + ti3
            cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4 + wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4 - wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5 + wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5 - wa4(i)*dr5
   30    continue
   40 continue
c
      return
c
      end
c
      integer function nfft(n)
      
      integer n

c     Given a natural number n, this function returns 
c     another natural number nfft such that nfft>=n and nfft 
c     satisfies:
c
c     (for IBM ESSL routines):
c     The only prime factors of nfft are 2,3,5, 7, 11
c     2 appears at least once, 3 appears at most twice, 
c     and 5, 7, and 11 appear at most once.
c
c     (for CONVEX and DEC's DXML)
c     set235 contains the possible nfft's when
c     only factors of 2, 3, and 5 (no other restrictions) are allowed
c
c     (for completeness, set2 contains the powers of 2)
c
      integer nessl, n235, npower2
      parameter (nessl=52,n235=51,npower2=8)

      integer i, max_ind

      integer setessl(nessl)
      integer set235(n235)
      integer set2(npower2)

      data setessl /2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24,
     &              28, 30, 32, 36, 40, 42, 44, 48, 56, 60, 64, 
     &              66, 70, 72, 80, 84, 88, 90, 96, 110, 112, 120,
     &              126, 128, 132, 140, 144, 154, 160, 168, 176,
     &              180, 192, 198, 210, 220, 224, 240, 252, 256/

      data set235 /2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20,
     &            24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60,
     &            64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125,
     &            128, 135, 144, 150, 160, 162, 180, 192, 200,
     &            216, 225, 240, 243, 250, 256/

      data set2 /2, 4, 8, 16, 32, 64, 128, 256/

C     ..
c
c     Choose appropriate settings below!!!
c
c      max_ind = nessl
c      max_ind = npower2
      max_ind = n235
c
      do i = 1, max_ind
         if ( set235(i) .ge. n) then
            nfft = set235(i)
            go to 10
         endif
      enddo
      
      write (6,'(/,a,i4,/)') '*** n>256 in nfft!  n=', n
      STOP

 10   continue

      return

      end
c
c It is the responsibility of the user to make sure that the
c only prime factors of the fft lenghts are 2, 3, 5.
c
c
      subroutine z1d_fft(n1,f,mode)
c
c      Computes (double)complex 1d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1               grid size.
c
c I/O  f(i)                i=1,n1
c                          array to be fourier transformed
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
      implicit none
c
C     .. Scalar Arguments ..
      integer mode, n1
C     ..
C     .. Array Arguments ..
      complex*16 f(*)
c
C     ..
      integer nmax
      parameter (nmax=20000)
c
c I/O  wrk(i)              i=1,4*max(n1)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex*16 wrk(4*nmax+15)
c
      double precision factor
      integer i, j, k, np

      external zfftb, zfftf, zffti

      np = n1 + 1

      if (n1 .gt. nmax) STOP 'FFT'
c
c     Transform the first dimension
c
      call zffti(n1,wrk(np))
            do 90 k = 1, n1
               wrk(k) = f(k)
   90       continue
            if (mode .gt. 0) call zfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k) = wrk(k)
  100       continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.D0/n1
      do 150 i = 1, n1
               f(i) = factor*f(i)
  150 continue
c
      return
c
      end
c
      subroutine z2d_fft(n1,n2,f,mode,np1)
c
c      Computes (double)complex 2d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1,n2            grid size
c Inp  np1              leading physical dimensions of array f
c
c I/O  f(i,j)              i=1,n1,j=1,n2
c                          array to be fourier transformed
c
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
      implicit none

C     .. Scalar Arguments ..
      integer mode, n1, n2, np1
C     ..
C     .. Array Arguments ..
      complex*16 f(np1,*)
c
C     ..
      integer nmax
      parameter (nmax=1024)
c
c I/O  wrk(i)              i=1,4*max(n1,n2)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex*16 wrk(4*nmax+15)
c
      double precision factor
      integer i, j, k, np

      external zfftb, zfftf, zffti

      intrinsic max

      np = max(n1,n2) + 1
c
      if (np-1 .gt. nmax) STOP 'FFT'
c
c     Transform the second dimension
c
      call zffti(n2,wrk(np))

         do 70 j = 1, n1
            do 50 k = 1, n2
               wrk(k) = f(j,k)
   50       continue
            if (mode .gt. 0) call zfftf(n2,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n2,wrk,wrk(np))
            do 60 k = 1, n2
               f(j,k) = wrk(k)
   60       continue
   70    continue
c
c     Transform the first dimension
c
      if (n1 .ne. n2) call zffti(n1,wrk(np))
      do 120 i = 1, n2
            do 90 k = 1, n1
               wrk(k) = f(k,i)
   90       continue
            if (mode .gt. 0) call zfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k,i) = wrk(k)
  100       continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.D0/(n1*n2)
      do 150 i = 1, n1
         do 140 j = 1, n2
               f(i,j) = factor*f(i,j)
  140    continue
  150 continue
c
      return
c
      end
c
      subroutine z3d_fft(n1,n2,n3,f,mode,np1,np2)
c
c      Computes (double)complex 3d fft of f
c      NOTE: A COMPLEX array is passed
c
c Inp  n1,n2,n3            grid size
c Inp  np1, np2            leading physical dimensions of array f
c
c I/O  f(i,j,k)            i=1,n1,j=1,n2,k=1,n3.
c                          array to be fourier transformed
c
c
c Inp  mode                determines transform mode
c                          mode > 0   forward
c                          mode < 0   backward
c
C     .. Scalar Arguments ..
      integer mode, n1, n2, n3, np1, np2
C     ..
C     .. Array Arguments ..
      complex*16 f(np1,np2,*)
c
C     ..
      integer nmax
      parameter (nmax=1024)
c
c I/O  wrk(i)              i=1,4*max(n1,n2,n3)+15
c                          work array. Treated as complex here for
c                          convenience.
c 
      complex*16 wrk(4*nmax+15)
c
      double precision factor
      integer i, j, k, np

      external zfftb, zfftf, zffti

      intrinsic max

      np = max(n1,n2,n3) + 1
c
      if (np-1 .gt. nmax) STOP 'FFT'
c
c     First call to initialize the work arrays
c
      call zffti(n3,wrk(np))
c
c     Transform the third dimension
c
      do 40 i = 1, n1
         do 30 j = 1, n2
            do 10 k = 1, n3
               wrk(k) = f(i,j,k)
   10       continue
            if (mode .gt. 0) call zfftf(n3,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n3,wrk,wrk(np))
            do 20 k = 1, n3
               f(i,j,k) = wrk(k)
   20       continue
   30    continue
   40 continue
c
c     Transform the second dimension
c
      if (n2 .ne. n3) call zffti(n2,wrk(np))
      do 80 i = 1, n3
         do 70 j = 1, n1
            do 50 k = 1, n2
               wrk(k) = f(j,k,i)
   50       continue
            if (mode .gt. 0) call zfftf(n2,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n2,wrk,wrk(np))
            do 60 k = 1, n2
               f(j,k,i) = wrk(k)
   60       continue
   70    continue
   80 continue
c
c     Transform the first dimension
c
      if (n1 .ne. n2) call zffti(n1,wrk(np))
      do 120 i = 1, n2
         do 110 j = 1, n3
            do 90 k = 1, n1
               wrk(k) = f(k,i,j)
   90       continue
            if (mode .gt. 0) call zfftf(n1,wrk,wrk(np))
            if (mode .lt. 0) call zfftb(n1,wrk,wrk(np))
            do 100 k = 1, n1
               f(k,i,j) = wrk(k)
  100       continue
  110    continue
  120 continue
c
c      Normalize
c
      if (mode .lt. 0) return
      factor = 1.D0/(n1*n2*n3)
      do 150 i = 1, n1
         do 140 j = 1, n2
            do 130 k = 1, n3
               f(i,j,k) = factor*f(i,j,k)
  130       continue
  140    continue
  150 continue
c
      return
c
      end
c fftpack    from ulib                                     06/05/81
c
cAG: Names begin with z or d to emphasize that this is a double precision
c    implementation of fftpack.
c
c purpose  this package consists of programs which perform fast fourier
c          transforms for both complex and real periodic sequences and
c          certain other symmetric sequences that are listed below.
c
c
c            zffti     initialize zfftf and zfftb
c            zfftf     forward transform of a complex periodic sequence
c            zfftb     unnormalized inverse of zfftf
c
c
c     subroutine zffti(n,wsave)
c
c     ******************************************************************
c
c     subroutine zffti initializes the array wsave which is used in
c     both zfftf and zfftb. the prime factorization of n together with
c     a tabulation of the trigonometric functions are computed and
c     stored in wsave.
c
c     input parameter
c
c     n       the length of the sequence to be transformed
c
c     output parameter
c
c     wsave   a work array which must be dimensioned at least 4*n+15
c             the same work array can be used for both zfftf and zfftb
c             as long as n remains unchanged. different wsave arrays
c             are required for different values of n. the contents of
c             wsave must not be changed between calls of zfftf or zfftb.
c
c     ******************************************************************
c
c     subroutine zfftf(n,c,wsave)
c
c     ******************************************************************
c
c     subroutine zfftf computes the forward complex discrete fourier
c     transform (the fourier analysis). equivalently , zfftf computes
c     the fourier coefficients of a complex periodic sequence.
c     the transform is defined below at output parameter c.
c
c     the transform is not normalized. to obtain a normalized transform
c     the output must be divided by n. otherwise a call of zfftf
c     followed by a call of zfftb will multiply the sequence by n.
c
c     the array wsave which is used by subroutine zfftf must be
c     initialized by calling subroutine zffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes. n
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls zfftf. the wsave array must be
c             initialized by calling subroutine zffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by zfftf and zfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(-i*j*k*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine zfftf or zfftb
c
c     ******************************************************************
c
c     subroutine zfftb(n,c,wsave)
c
c     ******************************************************************
c
c     subroutine zfftb computes the backward complex discrete fourier
c     transform (the fourier synthesis). equivalently , zfftb computes
c     a complex periodic sequence from its fourier coefficients.
c     the transform is defined below at output parameter c.
c
c     a call of zfftf followed by a call of zfftb will multiply the
c     sequence by n.
c
c     the array wsave which is used by subroutine zfftb must be
c     initialized by calling subroutine zffti(n,wsave).
c
c     input parameters
c
c
c     n      the length of the complex sequence c. the method is
c            more efficient when n is the product of small primes.
c
c     c      a complex array of length n which contains the sequence
c
c     wsave   a real work array which must be dimensioned at least 4n+15
c             in the program that calls zfftb. the wsave array must be
c             initialized by calling subroutine zffti(n,wsave) and a
c             different wsave array must be used for each different
c             value of n. this initialization does not have to be
c             repeated so long as n remains unchanged thus subsequent
c             transforms can be obtained faster than the first.
c             the same wsave array can be used by zfftf and zfftb.
c
c     output parameters
c
c     c      for j=1,...,n
c
c                c(j)=the sum from k=1,...,n of
c
c                      c(k)*exp(i*j*k*2*pi/n)
c
c                            where i=sqrt(-1)
c
c     wsave   contains initialization calculations which must not be
c             destroyed between calls of subroutine zfftf or zfftb
c___________________________________________________________________
c
      subroutine zffti(n,wsave)
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision wsave(*)
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external zffti1
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call zffti1(n,wsave(iw1),wsave(iw2))
c
      return
c
      end
c____________________________________________________________________
c
      subroutine zfftf(n,c,wsave)
c
C     .. Array Arguments ..
      double precision c(*), wsave(*)
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external zfftf1
C     ..
C     .. Scalar Arguments ..
      integer n
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call zfftf1(n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
c
      end
c____________________________________________________________________
c
      subroutine zfftb(n,c,wsave)
c
C     .. Array Arguments ..
C
c
      double precision c(*), wsave(*)
C     ..
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Local Scalars ..
      integer iw1, iw2
C     ..
C     .. External Subroutines ..
      external zfftb1
C     ..
      if (n .eq. 1) return
c
      iw1 = n + n + 1
      iw2 = iw1 + n + n
      call zfftb1(n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
c
      end
c_____________________________________________________________
c
      subroutine zffti1(n,wa,ifac)
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      double precision arg, argh, argld, fi, tpi
      integer i, i1, ib, ido, idot, ii, ip, ipm, j, k1, l1, l2, ld, nf,
     1        nl, nq, nr, ntry
C     ..
C     .. Local Arrays ..
      integer ntryh(4)
C     ..
C     .. Intrinsic Functions ..
      intrinsic cos, dfloat, sin
C     ..
C     .. Data statements ..
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4)/3, 4, 2, 5/
      data tpi/6.28318530717958647692528676655900577d0/
C     ..
c
      nl = n
      nf = 0
      j = 0
c
   10 continue
      j = j + 1
      if (j .le. 4) ntry = ntryh(j)
      if (j .gt. 4) ntry = ntry + 2
   20 continue
      nq = nl/ntry
      nr = nl - ntry*nq
      if (nr .ne. 0) go to 10
c
      nf = nf + 1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 40
      if (nf .eq. 1) go to 40
      do 30 i = 2, nf
         ib = nf - i + 2
         ifac(ib+2) = ifac(ib+1)
   30 continue
      ifac(3) = 2
c
   40 continue
      if (nl .ne. 1) go to 20
c
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/dfloat(n)
      i = 2
      l1 = 1
      do 70 k1 = 1, nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido + ido + 2
         ipm = ip - 1
c
         do 60 j = 1, ipm
            i1 = i
            wa(i-1) = 1.d0
            wa(i) = 0.d0
            ld = ld + l1
            fi = 0.d0
            argld = dfloat(ld)*argh
            do 50 ii = 4, idot, 2
               i = i + 2
               fi = fi + 1.d0
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
   50       continue
            if (ip .le. 5) go to 60
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
   60    continue
c
         l1 = l2
   70 continue
c
      return
c
      end
c________________________________________________________________
c
      subroutine zfftf1(n,c,ch,wa,ifac)
c Dummy dimension compatible with semantic analyzer.
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision c(*), ch(*), wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      integer i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2, n2,
     1        na, nac, nf
C     ..
C     .. External Subroutines ..
      external dpssf, dpssf2, dpssf3, dpssf4, dpssf5
C     ..
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 160 k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido + ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 30
         ix2 = iw + idot
         ix3 = ix2 + idot
         if (na .ne. 0) go to 10
         call dpssf4(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
c
         go to 20
c
   10    continue
         call dpssf4(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
   20    continue
         na = 1 - na
c
         go to 150
c
   30    continue
         if (ip .ne. 2) go to 60
         if (na .ne. 0) go to 40
         call dpssf2(idot,l1,c,ch,wa(iw))
c
         go to 50
c
   40    continue
         call dpssf2(idot,l1,ch,c,wa(iw))
   50    continue
         na = 1 - na
c
         go to 150
c
   60    continue
         if (ip .ne. 3) go to 90
         ix2 = iw + idot
         if (na .ne. 0) go to 70
         call dpssf3(idot,l1,c,ch,wa(iw),wa(ix2))
c
         go to 80
c
   70    continue
         call dpssf3(idot,l1,ch,c,wa(iw),wa(ix2))
   80    continue
         na = 1 - na
c
         go to 150
c
   90    continue
         if (ip .ne. 5) go to 120
         ix2 = iw + idot
         ix3 = ix2 + idot
         ix4 = ix3 + idot
         if (na .ne. 0) go to 100
         call dpssf5(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
c
         go to 110
c
  100    continue
         call dpssf5(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    continue
         na = 1 - na
c
         go to 150
c
  120    continue
         if (na .ne. 0) go to 130
         call dpssf(nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
c
         go to 140
c
  130    continue
         call dpssf(nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  140    continue
         if (nac .ne. 0) na = 1 - na
c
  150    continue
         l1 = l2
         iw = iw + (ip-1)*idot
  160 continue
      if (na .eq. 0) return
c
      n2 = n + n
      do 170 i = 1, n2
         c(i) = ch(i)
  170 continue
c
      return
c
      end
c___________________________________________________________________
c
      subroutine zfftb1(n,c,ch,wa,ifac)
c
C     .. Scalar Arguments ..
      integer n
C     ..
C     .. Array Arguments ..
      double precision c(*), ch(*), wa(*)
      integer ifac(*)
C     ..
C     .. Local Scalars ..
      integer i, idl1, ido, idot, ip, iw, ix2, ix3, ix4, k1, l1, l2, n2,
     1        na, nac, nf
C     ..
C     .. External Subroutines ..
      external dpssb, dpssb2, dpssb3, dpssb4, dpssb5
C     ..
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 160 k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido + ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 30
         ix2 = iw + idot
         ix3 = ix2 + idot
         if (na .ne. 0) go to 10
         call dpssb4(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
c
         go to 20
c
   10    continue
         call dpssb4(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
   20    continue
         na = 1 - na
c
         go to 150
c
   30    continue
         if (ip .ne. 2) go to 60
         if (na .ne. 0) go to 40
         call dpssb2(idot,l1,c,ch,wa(iw))
c
         go to 50
c
   40    continue
         call dpssb2(idot,l1,ch,c,wa(iw))
   50    continue
         na = 1 - na
c
         go to 150
c
   60    continue
         if (ip .ne. 3) go to 90
         ix2 = iw + idot
         if (na .ne. 0) go to 70
         call dpssb3(idot,l1,c,ch,wa(iw),wa(ix2))
c
         go to 80
c
   70    continue
         call dpssb3(idot,l1,ch,c,wa(iw),wa(ix2))
   80    continue
         na = 1 - na
c
         go to 150
c
   90    continue
         if (ip .ne. 5) go to 120
         ix2 = iw + idot
         ix3 = ix2 + idot
         ix4 = ix3 + idot
         if (na .ne. 0) go to 100
         call dpssb5(idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
c
         go to 110
c
  100    continue
         call dpssb5(idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    continue
         na = 1 - na
c
         go to 150
c
  120    continue
         if (na .ne. 0) go to 130
         call dpssb(nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
c
         go to 140
c
  130    continue
         call dpssb(nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  140    continue
         if (nac .ne. 0) na = 1 - na
c
  150    continue
         l1 = l2
         iw = iw + (ip-1)*idot
  160 continue
      if (na .eq. 0) return
c
      n2 = n + n
      do 170 i = 1, n2
         c(i) = ch(i)
  170 continue
c
      return
c
      end
c_____________________________________________________________________
c
      subroutine dpssb(nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
C     .. Scalar Arguments ..
      integer idl1, ido, ip, l1, nac
C     ..
C     .. Array Arguments ..
      double precision c1(ido,l1,ip), c2(idl1,ip), cc(ido,ip,l1),
     1                 ch(ido,l1,ip), ch2(idl1,ip), wa(*)
C     ..
C     .. Local Scalars ..
      double precision wai, war
      integer i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph,
     1        j, jc, k, l, lc, nt
C     ..
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip + 2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 60
      do 30 j = 2, ipph
         jc = ipp2 - j
         do 20 k = 1, l1
            do 10 i = 1, ido
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   10       continue
   20    continue
   30 continue
c
      do 50 k = 1, l1
         do 40 i = 1, ido
            ch(i,k,1) = cc(i,1,k)
   40    continue
   50 continue
c
      go to 120
c
   60 continue
      do 90 j = 2, ipph
         jc = ipp2 - j
         do 80 i = 1, ido
            do 70 k = 1, l1
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   70       continue
   80    continue
   90 continue
c
      do 110 i = 1, ido
         do 100 k = 1, l1
            ch(i,k,1) = cc(i,1,k)
  100    continue
  110 continue
c
  120 continue
      idl = 2 - ido
      inc = 0
      do 160 l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do 130 ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
  130    continue
         idlj = idl
         inc = inc + ido
         do 150 j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp) idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 140 ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) + wai*ch2(ik,jc)
  140       continue
  150    continue
  160 continue
c
      do 180 j = 2, ipph
         do 170 ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
  170    continue
  180 continue
c
      do 200 j = 2, ipph
         jc = ipp2 - j
         do 190 ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
  190    continue
  200 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 210 ik = 1, idl1
         c2(ik,1) = ch2(ik,1)
  210 continue
c
      do 230 j = 2, ip
         do 220 k = 1, l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  220    continue
  230 continue
c
      if (idot .gt. l1) go to 270
      idij = 0
      do 260 j = 2, ip
         idij = idij + 2
         do 250 i = 4, ido, 2
            idij = idij + 2
            do 240 k = 1, l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) + wa(idij)*ch(i-1,k,j)
  240       continue
  250    continue
  260 continue
c
      return
c
  270 continue
      idj = 2 - ido
      do 300 j = 2, ip
         idj = idj + ido
         do 290 k = 1, l1
            idij = idj
            do 280 i = 4, ido, 2
               idij = idij + 2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) - wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) + wa(idij)*ch(i-1,k,j)
  280       continue
  290    continue
  300 continue
c
      return
c
      end
c
      subroutine dpssb2(ido,l1,cc,ch,wa1)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(*)
C     ..
C     .. Local Scalars ..
      double precision ti2, tr2
      integer i, k
C     ..
      if (ido .gt. 2) go to 20
      do 10 k = 1, l1
         ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
         ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
         ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
            tr2 = cc(i-1,1,k) - cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
            ti2 = cc(i,1,k) - cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2 + wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2 - wa1(i)*ti2
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssb3(ido,l1,cc,ch,wa1,wa2)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui,
     1                 taur, ti2, tr2
      integer i, k
C     ..
C     .. Data statements ..
      data taur/-0.5D0/
      data taui/0.86602540378443864676372317075293618D0/
C     ..
c
c     one half sqrt(3) = .866025.....  .
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         tr2 = cc(1,2,k) + cc(1,3,k)
         cr2 = cc(1,1,k) + taur*tr2
         ch(1,k,1) = cc(1,1,k) + tr2
         ti2 = cc(2,2,k) + cc(2,3,k)
         ci2 = cc(2,1,k) + taur*ti2
         ch(2,k,1) = cc(2,1,k) + ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2 - ci3
         ch(1,k,3) = cr2 + ci3
         ch(2,k,2) = ci2 + cr3
         ch(2,k,3) = ci2 - cr3
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            tr2 = cc(i-1,2,k) + cc(i-1,3,k)
            cr2 = cc(i-1,1,k) + taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k) + tr2
            ti2 = cc(i,2,k) + cc(i,3,k)
            ci2 = cc(i,1,k) + taur*ti2
            ch(i,k,1) = cc(i,1,k) + ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssb4(ido,l1,cc,ch,wa1,wa2,wa3)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1                 wa3(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4,
     1                 tr1, tr2, tr3, tr4
      integer i, k
C     ..
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti1 = cc(2,1,k) - cc(2,3,k)
         ti2 = cc(2,1,k) + cc(2,3,k)
         tr4 = cc(2,4,k) - cc(2,2,k)
         ti3 = cc(2,2,k) + cc(2,4,k)
         tr1 = cc(1,1,k) - cc(1,3,k)
         tr2 = cc(1,1,k) + cc(1,3,k)
         ti4 = cc(1,2,k) - cc(1,4,k)
         tr3 = cc(1,2,k) + cc(1,4,k)
         ch(1,k,1) = tr2 + tr3
         ch(1,k,3) = tr2 - tr3
         ch(2,k,1) = ti2 + ti3
         ch(2,k,3) = ti2 - ti3
         ch(1,k,2) = tr1 + tr4
         ch(1,k,4) = tr1 - tr4
         ch(2,k,2) = ti1 + ti4
         ch(2,k,4) = ti1 - ti4
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti1 = cc(i,1,k) - cc(i,3,k)
            ti2 = cc(i,1,k) + cc(i,3,k)
            ti3 = cc(i,2,k) + cc(i,4,k)
            tr4 = cc(i,4,k) - cc(i,2,k)
            tr1 = cc(i-1,1,k) - cc(i-1,3,k)
            tr2 = cc(i-1,1,k) + cc(i-1,3,k)
            ti4 = cc(i-1,2,k) - cc(i-1,4,k)
            tr3 = cc(i-1,2,k) + cc(i-1,4,k)
            ch(i-1,k,1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i,k,1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1,k,2) = wa1(i-1)*cr2 - wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2 + wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3 - wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3 + wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4 - wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4 + wa3(i)*cr4
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssb5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1                 wa3(*), wa4(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3,
     1                 di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2,
     2                 ti3, ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      integer i, k
C     ..
C     .. Data statements ..
      data tr11/0.30901699437494742410229341718281906D0/
      data ti11/0.95105651629515357211643933337938214D0/
      data tr12/-0.80901699437494742410229341718281906D0/
      data ti12/0.58778525229247312916870595463907277D0/
C     ..
c
c     sin(pi/10) = .30901699....    .
c     cos(pi/10) = .95105651....    .
c     sin(pi/5 ) = .58778525....    .
c     cos(pi/5 ) = .80901699....    .
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti5 = cc(2,2,k) - cc(2,5,k)
         ti2 = cc(2,2,k) + cc(2,5,k)
         ti4 = cc(2,3,k) - cc(2,4,k)
         ti3 = cc(2,3,k) + cc(2,4,k)
         tr5 = cc(1,2,k) - cc(1,5,k)
         tr2 = cc(1,2,k) + cc(1,5,k)
         tr4 = cc(1,3,k) - cc(1,4,k)
         tr3 = cc(1,3,k) + cc(1,4,k)
         ch(1,k,1) = cc(1,1,k) + tr2 + tr3
         ch(2,k,1) = cc(2,1,k) + ti2 + ti3
         cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
         ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
         cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
         ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         ch(1,k,2) = cr2 - ci5
         ch(1,k,5) = cr2 + ci5
         ch(2,k,2) = ci2 + cr5
         ch(2,k,3) = ci3 + cr4
         ch(1,k,3) = cr3 - ci4
         ch(1,k,4) = cr3 + ci4
         ch(2,k,4) = ci3 - cr4
         ch(2,k,5) = ci2 - cr5
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti5 = cc(i,2,k) - cc(i,5,k)
            ti2 = cc(i,2,k) + cc(i,5,k)
            ti4 = cc(i,3,k) - cc(i,4,k)
            ti3 = cc(i,3,k) + cc(i,4,k)
            tr5 = cc(i-1,2,k) - cc(i-1,5,k)
            tr2 = cc(i-1,2,k) + cc(i-1,5,k)
            tr4 = cc(i-1,3,k) - cc(i-1,4,k)
            tr3 = cc(i-1,3,k) + cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
            ch(i,k,1) = cc(i,1,k) + ti2 + ti3
            cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4 - wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4 + wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5 - wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5 + wa4(i)*dr5
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssf(nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
c
C     .. Scalar Arguments ..
      integer idl1, ido, ip, l1, nac
C     ..
C     .. Array Arguments ..
      double precision c1(ido,l1,ip), c2(idl1,ip), cc(ido,ip,l1),
     1                 ch(ido,l1,ip), ch2(idl1,ip), wa(*)
C     ..
C     .. Local Scalars ..
      double precision wai, war
      integer i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph,
     1        j, jc, k, l, lc, nt
C     ..
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip + 2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 60
      do 30 j = 2, ipph
         jc = ipp2 - j
         do 20 k = 1, l1
            do 10 i = 1, ido
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   10       continue
   20    continue
   30 continue
c
      do 50 k = 1, l1
         do 40 i = 1, ido
            ch(i,k,1) = cc(i,1,k)
   40    continue
   50 continue
c
      go to 120
c
   60 continue
      do 90 j = 2, ipph
         jc = ipp2 - j
         do 80 i = 1, ido
            do 70 k = 1, l1
               ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
   70       continue
   80    continue
   90 continue
c
      do 110 i = 1, ido
         do 100 k = 1, l1
            ch(i,k,1) = cc(i,1,k)
  100    continue
  110 continue
c
  120 continue
      idl = 2 - ido
      inc = 0
      do 160 l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do 130 ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  130    continue
         idlj = idl
         inc = inc + ido
         do 150 j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp) idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 140 ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) - wai*ch2(ik,jc)
  140       continue
  150    continue
  160 continue
c
      do 180 j = 2, ipph
         do 170 ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
  170    continue
  180 continue
c
      do 200 j = 2, ipph
         jc = ipp2 - j
         do 190 ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
  190    continue
  200 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 210 ik = 1, idl1
         c2(ik,1) = ch2(ik,1)
  210 continue
c
      do 230 j = 2, ip
         do 220 k = 1, l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  220    continue
  230 continue
c
      if (idot .gt. l1) go to 270
      idij = 0
      do 260 j = 2, ip
         idij = idij + 2
         do 250 i = 4, ido, 2
            idij = idij + 2
            do 240 k = 1, l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) + wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
  240       continue
  250    continue
  260 continue
c
      return
c
  270 continue
      idj = 2 - ido
      do 300 j = 2, ip
         idj = idj + ido
         do 290 k = 1, l1
            idij = idj
            do 280 i = 4, ido, 2
               idij = idij + 2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j) + wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j) - wa(idij)*ch(i-1,k,j)
  280       continue
  290    continue
  300 continue
c
      return
c
      end
c
      subroutine dpssf2(ido,l1,cc,ch,wa1)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(*)
C     ..
C     .. Local Scalars ..
      double precision ti2, tr2
      integer i, k
C     ..
      if (ido .gt. 2) go to 20
      do 10 k = 1, l1
         ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
         ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
         ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
         ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
            tr2 = cc(i-1,1,k) - cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
            ti2 = cc(i,1,k) - cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2 - wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2 + wa1(i)*ti2
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssf3(ido,l1,cc,ch,wa1,wa2)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui,
     1                 taur, ti2, tr2
      integer i, k
C     ..
C     .. Data statements ..
      data taur/-0.5D0/
      data taui/-0.86602540378443864676372317075293618D0/
C     ..
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         tr2 = cc(1,2,k) + cc(1,3,k)
         cr2 = cc(1,1,k) + taur*tr2
         ch(1,k,1) = cc(1,1,k) + tr2
         ti2 = cc(2,2,k) + cc(2,3,k)
         ci2 = cc(2,1,k) + taur*ti2
         ch(2,k,1) = cc(2,1,k) + ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2 - ci3
         ch(1,k,3) = cr2 + ci3
         ch(2,k,2) = ci2 + cr3
         ch(2,k,3) = ci2 - cr3
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            tr2 = cc(i-1,2,k) + cc(i-1,3,k)
            cr2 = cc(i-1,1,k) + taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k) + tr2
            ti2 = cc(i,2,k) + cc(i,3,k)
            ci2 = cc(i,1,k) + taur*ti2
            ch(i,k,1) = cc(i,1,k) + ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2 - ci3
            dr3 = cr2 + ci3
            di2 = ci2 + cr3
            di3 = ci2 - cr3
            ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssf4(ido,l1,cc,ch,wa1,wa2,wa3)
c
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1                 wa3(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4,
     1                 tr1, tr2, tr3, tr4
      integer i, k
C     ..
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti1 = cc(2,1,k) - cc(2,3,k)
         ti2 = cc(2,1,k) + cc(2,3,k)
         tr4 = cc(2,2,k) - cc(2,4,k)
         ti3 = cc(2,2,k) + cc(2,4,k)
         tr1 = cc(1,1,k) - cc(1,3,k)
         tr2 = cc(1,1,k) + cc(1,3,k)
         ti4 = cc(1,4,k) - cc(1,2,k)
         tr3 = cc(1,2,k) + cc(1,4,k)
         ch(1,k,1) = tr2 + tr3
         ch(1,k,3) = tr2 - tr3
         ch(2,k,1) = ti2 + ti3
         ch(2,k,3) = ti2 - ti3
         ch(1,k,2) = tr1 + tr4
         ch(1,k,4) = tr1 - tr4
         ch(2,k,2) = ti1 + ti4
         ch(2,k,4) = ti1 - ti4
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti1 = cc(i,1,k) - cc(i,3,k)
            ti2 = cc(i,1,k) + cc(i,3,k)
            ti3 = cc(i,2,k) + cc(i,4,k)
            tr4 = cc(i,2,k) - cc(i,4,k)
            tr1 = cc(i-1,1,k) - cc(i-1,3,k)
            tr2 = cc(i-1,1,k) + cc(i-1,3,k)
            ti4 = cc(i-1,4,k) - cc(i-1,2,k)
            tr3 = cc(i-1,2,k) + cc(i-1,4,k)
            ch(i-1,k,1) = tr2 + tr3
            cr3 = tr2 - tr3
            ch(i,k,1) = ti2 + ti3
            ci3 = ti2 - ti3
            cr2 = tr1 + tr4
            cr4 = tr1 - tr4
            ci2 = ti1 + ti4
            ci4 = ti1 - ti4
            ch(i-1,k,2) = wa1(i-1)*cr2 + wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2 - wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3 + wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3 - wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4 + wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4 - wa3(i)*cr4
   30    continue
   40 continue
c
      return
c
      end
c
      subroutine dpssf5(ido,l1,cc,ch,wa1,wa2,wa3,wa4)
C     .. Scalar Arguments ..
      integer ido, l1
C     ..
C     .. Array Arguments ..
      double precision cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1                 wa3(*), wa4(*)
C     ..
C     .. Local Scalars ..
      double precision ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3,
     1                 di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2,
     2                 ti3, ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      integer i, k
C     ..
C     .. Data statements ..
      data tr11/0.30901699437494742410229341718281906D0/
      data ti11/-0.95105651629515357211643933337938214D0/
      data tr12/-0.80901699437494742410229341718281906D0/
      data ti12/-0.58778525229247312916870595463907277D0/
C     ..
c
      if (ido .ne. 2) go to 20
      do 10 k = 1, l1
         ti5 = cc(2,2,k) - cc(2,5,k)
         ti2 = cc(2,2,k) + cc(2,5,k)
         ti4 = cc(2,3,k) - cc(2,4,k)
         ti3 = cc(2,3,k) + cc(2,4,k)
         tr5 = cc(1,2,k) - cc(1,5,k)
         tr2 = cc(1,2,k) + cc(1,5,k)
         tr4 = cc(1,3,k) - cc(1,4,k)
         tr3 = cc(1,3,k) + cc(1,4,k)
         ch(1,k,1) = cc(1,1,k) + tr2 + tr3
         ch(2,k,1) = cc(2,1,k) + ti2 + ti3
         cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
         ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
         cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
         ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
         cr5 = ti11*tr5 + ti12*tr4
         ci5 = ti11*ti5 + ti12*ti4
         cr4 = ti12*tr5 - ti11*tr4
         ci4 = ti12*ti5 - ti11*ti4
         ch(1,k,2) = cr2 - ci5
         ch(1,k,5) = cr2 + ci5
         ch(2,k,2) = ci2 + cr5
         ch(2,k,3) = ci3 + cr4
         ch(1,k,3) = cr3 - ci4
         ch(1,k,4) = cr3 + ci4
         ch(2,k,4) = ci3 - cr4
         ch(2,k,5) = ci2 - cr5
   10 continue
c
      return
c
   20 continue
      do 40 k = 1, l1
         do 30 i = 2, ido, 2
            ti5 = cc(i,2,k) - cc(i,5,k)
            ti2 = cc(i,2,k) + cc(i,5,k)
            ti4 = cc(i,3,k) - cc(i,4,k)
            ti3 = cc(i,3,k) + cc(i,4,k)
            tr5 = cc(i-1,2,k) - cc(i-1,5,k)
            tr2 = cc(i-1,2,k) + cc(i-1,5,k)
            tr4 = cc(i-1,3,k) - cc(i-1,4,k)
            tr3 = cc(i-1,3,k) + cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
            ch(i,k,1) = cc(i,1,k) + ti2 + ti3
            cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            dr3 = cr3 - ci4
            dr4 = cr3 + ci4
            di3 = ci3 + cr4
            di4 = ci3 - cr4
            dr5 = cr2 + ci5
            dr2 = cr2 - ci5
            di5 = ci2 - cr5
            di2 = ci2 + cr5
            ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4 + wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4 - wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5 + wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5 - wa4(i)*dr5
   30    continue
   40 continue
c
      return
c
      end
