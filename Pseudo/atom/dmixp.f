C
c $Id: dmixp.f,v 1.4 1999/02/25 13:48:42 wdpgaara Exp $
c
      subroutine dmixp(a,b,beta,icy,id,nmsh,c,d,vn1,vn12,vn2,vn22)
c
      implicit none
c
C*    ADAPTED FROM K.C.PANDEY
C*    USING ANDERSON'S EXTRAPOLATION SCHEME
C*    EQS 4.1-4.9,4.15-4.18 OF
C*    D.G.ANDERSON J.ASSOC.COMPUTING MACHINERY,12,547(1965)
c
C*    COMPUTES A NEW VECTOR IN A ITERATIVE SCHEME
c
C*    INPUT A=NEWPOT B=OLDPOT
C*    OUTPUT A=A-B B=NEWPOT
C*    BETA=MIXING,IN=ITER. NUMBER
C*    ID=1,2 OR 3 DIFF CONV METH.
C*    ICY CYCLE NUMBER ,ICY=1 ON FIRST/ZEROTH CALL
C*    C,D WORK ARRAYS OF SIZE NMSH
C*    VN1,VN12,VN2,VN22 STORAGE ARRAYS OF SIZE NMSH
C
c
c     Modified by Alberto Garcia to make use of Level 1 BLAS
c
C     .. Parameters ..
      double precision uze, um, detol
      parameter (uze=0.0D0,um=1.0D0,detol=1.D-9)
C     ..
C     .. Scalar Arguments ..
      double precision beta
      integer icy, id, nmsh
C     ..
C     .. Array Arguments ..
      double precision a(*), b(*), c(*), d(*)
C     ..
C     .. Local Scalars ..
      double precision a2, bt1, bt2, d11, d12, d22, det, dett, r2, rd1m,
     &                 rd2m, t1, t2, x
      integer in
C     ..
C     .. External Subroutines ..
c     BLAS level 1 (SCILIB)
      double precision sdot
      external sdot, saxpy, sscal, scopy
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      double precision vn1(*), vn12(*), vn2(*), vn22(*)
C     ..
      in = icy - 1
      if (in .eq. 0) then
         call saxpy(nmsh,um,a,1,b,1)
c
         return
c
      end if
      call saxpy(nmsh,-um,b,1,a,1)
      r2 = sdot(nmsh,a,1,a,1)
      if (id .eq. 1) then
         call saxpy(nmsh,beta,a,1,b,1)
c
         return
c
      end if
      if (in .eq. 1) then
         call scopy(nmsh,a,1,vn1,1)
         call scopy(nmsh,b,1,vn2,1)
         call saxpy(nmsh,beta,a,1,b,1)
c
         return
c
      end if
      call scopy(nmsh,vn1,1,c,1)
      if (id .eq. 3 .and. in .gt. 2) then
         call scopy(nmsh,vn12,1,d,1)
      end if
      call scopy(nmsh,a,1,vn1,1)
      if (id .gt. 2 .and. in .gt. 1) then
         call scopy(nmsh,c,1,vn12,1)
      end if
      call saxpy(nmsh,-um,a,1,c,1)
      d11 = sdot(nmsh,c,1,c,1)
      rd1m = sdot(nmsh,a,1,c,1)
      if (in .le. 2 .or. id .le. 2) then
         t1 = -rd1m/d11
         x = um - t1
         bt1 = beta*t1
         call sscal(nmsh,beta,a,1)
         call saxpy(nmsh,bt1,c,1,a,1)
         call scopy(nmsh,vn2,1,d,1)
         call saxpy(nmsh,t1,d,1,a,1)
         call scopy(nmsh,b,1,vn2,1)
         if (id .gt. 2 .and. in .eq. 2) then
            call scopy(nmsh,d,1,vn22,1)
         end if
         call sscal(nmsh,x,b,1)
         call saxpy(nmsh,1.d0,a,1,b,1)
c
         return
c
      end if
      call saxpy(nmsh,-um,a,1,d,1)
      d22 = sdot(nmsh,d,1,d,1)
      d12 = sdot(nmsh,c,1,d,1)
      rd2m = sdot(nmsh,a,1,d,1)
      a2 = d11*d22
      det = a2 - d12*d12
      dett = det/a2
      if (abs(dett) .ge. detol) then
         t1 = (-rd1m*d22+rd2m*d12)/det
         t2 = (rd1m*d12-rd2m*d11)/det
      else
         t1 = -rd1m/d11
         t2 = uze
      end if
      x = um - t1 - t2
      bt1 = beta*t1
      bt2 = beta*t2
      call sscal(nmsh,beta,a,1)
      call saxpy(nmsh,bt1,c,1,a,1)
      call saxpy(nmsh,bt2,d,1,a,1)
      call saxpy(nmsh,t1,vn2,1,a,1)
      call saxpy(nmsh,t2,vn22,1,a,1)
      call scopy(nmsh,vn2,1,vn22,1)
      call scopy(nmsh,b,1,vn2,1)
      call sscal(nmsh,x,b,1)
      call saxpy(nmsh,1.d0,a,1,b,1)
c
      return
c
      end







