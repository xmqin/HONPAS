c
      subroutine coreq
c
c  Compute and print the fourier transform of the pseudocore charge
c
c  Note: It is assumed that the function is reasonably smooth... otherwise
c        the Gauss-Legendre quadrature will not work!
c
      include 'radial.h'
      include 'charge.h'
c
      double precision pi
      parameter (pi=3.141592653589d0)
      double precision tol
      parameter ( tol = 1.d-6 )
c
c     Number of abscissae
c
      integer ngauss
      parameter (ngauss=100)
      double precision cored(ngauss), t(ngauss), w(ngauss)
c
c     How many points at what spacing? Should go up to q=20, say
c
      double precision delql
      parameter (delql = 0.05d0)
      integer nql
      parameter (nql = 400)

      integer j, nrpnew
      integer norder, iu, k
      double precision rmin, rmax, zc, q, dcq

      double precision divdif
      external divdif
c
c     Fourier transform the core charge
c
c     find smallest r above which 4pi*r^2*core charge is negligible
c
      do j = nr, 1, -1
         if (abs(cdc(j)) .gt. tol) go to 130
      enddo
 130  continue
      nrpnew = j + 1
      if (nrpnew .gt. nr) nrpnew = nr
c
      rmin = r(1)
      rmax = r(nrpnew)
c
c     Generate abscissas and weights for Gauss-Legendre integration
c
      call gauleg(rmin,rmax,t,w,ngauss)
c
c     Total core charge    (q=0)
c
c     We evaluate the function to be integrated at the calculated
c     abscissas (and multiply by the weight for convenience)
c     The Cern library function divdif ( interpolating to third
c     order ) is used.
c
      norder = 3

      zc = 0.d0
      do  j = 1, ngauss
         cored(j) = w(j)*divdif(cdc,r,nr,t(j),norder)
         zc = zc + cored(j)
      enddo
      write(6,'(a,f8.4)') 'Total pseudocore charge: ', zc

      call get_unit(iu)
      open(iu,file='COREQ',form='formatted',status='unknown')
      rewind(iu)

      write(iu,'(f8.2,4x,f12.4)') 0.d0, zc
c
c      Rest of the Fourier transform
c
      do  k = 1, nql
         q = delql*k
         dcq = 0.d0
         do j = 1, ngauss
            dcq = dcq + cored(j)*sin(q*t(j))/t(j)
         enddo
         dcq = dcq/q
         write(iu,'(f8.2,4x,f12.4)') q, dcq
      enddo
      close(iu)

      end


