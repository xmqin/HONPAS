      module m_linpack
      public :: zgedi, zgefa
      
      integer, parameter, private :: dp = kind(1.d0)  ! to make it consistent

c     Updated (somewhat) by A. Garcia, Sep 2016
c
      CONTAINS

      subroutine zgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(n),job
      complex(dp) a(lda,n),det(2),work(n)
c
c     zgedi computes the determinant and inverse of a matrix
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex(dp)(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        work    complex(dp)(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex(dp)(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,zswap
c     fortran abs,cmplx,mod
c
c     internal variables
c
      complex(dp) t
      real(dp) ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex(dp) zdum
      real(dp) cabs1
      real(dp) real,aimag
      complex(dp) zdumr,zdumi
      real(dp) :: local_real, local_aimag
      local_real(zdumr) = zdumr
      local_aimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = abs(local_real(zdum)) + abs(local_aimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
#ifdef OLD_CRAY
            call cscal(k-1,t,a(1,k),1)
#else
            call zscal(k-1,t,a(1,k),1)
#endif
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
#ifdef OLD_CRAY
               call caxpy(k,t,a(1,k),1,a(1,j),1)
#else
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
#endif
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
#ifdef OLD_CRAY
               call caxpy(n,t,a(1,j),1,a(1,k),1)
#else
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
#endif
  120       continue
            l = ipvt(k)
#ifdef OLD_CRAY
            if (l .ne. k) call cswap(n,a(1,k),1,a(1,l),1)
#else
            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
#endif
  130    continue
  140    continue
  150 continue
      return
      end subroutine zgedi

      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      complex(dp) a(lda,n)
c
c     zgefa factors a complex(dp) matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex(dp)(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax_l
c     fortran abs
c
c     internal variables
c
      complex(dp) t
      integer j,k,kp1,l,nm1
c
c!    Statement functions...
c!
      complex(dp) zdum
      real(dp) cabs1
      real(dp) real,aimag
      complex(dp) zdumr,zdumi
      real(zdumr) = zdumr
      aimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax_l(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
#ifdef OLD_CRAY
            call cscal(n-k,t,a(k+1,k),1)
#else
            call zscal(n-k,t,a(k+1,k),1)
#endif
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
#ifdef OLD_CRAY
               call caxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
#else
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
#endif
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end subroutine zgefa

      integer function izamax_l(n,zx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 3/11/78.
c
      integer n,incx
      integer :: ix, i
      complex(dp) zx(*)
      real(dp) smax
c
      izamax_l = 1
      if(n.le.1)return
      if(incx.eq.1)go to 20
c
c        code for increments not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      smax = sdcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(sdcabs1(zx(ix)).le.smax) go to 5
         izamax_l = ix
         smax = sdcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increments equal to 1
c
   20 smax = sdcabs1(zx(1))
      do 30 i = 2,n
         if(sdcabs1(zx(i)).le.smax) go to 30
         izamax_l = i
         smax = sdcabs1(zx(i))
   30 continue
      return
      end function izamax_l
C
C  The code below is a modified form of the routine sdcabs1
C  changed to avoid problems on the Cray when compiled
C  with optimisation turned on. JDG, November 1999
C
      real(dp) function sdcabs1(z)
      complex(dp) z
      real(dp) t(2)
      t(1) = real(z)
      t(2) = aimag(z)
      sdcabs1 = abs(t(1)) + abs(t(2))
      return
      end function sdcabs1

      end module m_linpack

