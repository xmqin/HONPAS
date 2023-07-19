C
c $Id: dsolv1.f,v 1.3 1999/02/25 13:49:13 wdpgaara Exp $
c
      subroutine dsolv1(nfirst,nlast,nn)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
c   dsolv1 finds the (non)-relativistic wave function
c   using finite differences and matrix diagonalization.
c   An initial guess for the eigenvalues need not be supplied.
c
C     .. Parameters ..
      double precision zero, one, pone, opf
      parameter (zero=0.D0,one=1.D0,pone=0.1D0,opf=1.5D0)
C     ..
C     .. Scalar Arguments ..
      integer nfirst, nlast
C     ..
c     nn will generally be 'no'
c
      integer nn(*)
C     ..
C     .. Local Scalars ..
      double precision bl, bu, c1, c2, denr, eps
      integer i, ierr, j, k, ki, kn, l, llp, nrm
C     ..
C     .. Local Arrays ..
      double precision e(10)
      integer ind(10), nmax(2,5)
C     ..
C     .. External Subroutines ..
      external tinvit, tridib
C     ..
      double precision d(nrmax), dk(nrmax), rv1(nrmax), rv2(nrmax),
     &                 rv3(nrmax), rv4(nrmax), rv5(nrmax), sd(nrmax),
     &                 sd2(nrmax), z(6*nrmax)
C     ..
c
c   Initialize the charge density arrays.
c
      do 10 i = 1, nr
         cdd(i) = zero
         cdu(i) = zero
   10 continue
c
c   Find the max n given l and s.
c   Zero spin is treated as down.
c
      do 40 i = 1, 2
         do 30 j = 1, lmax
            nmax(i,j) = 0
            do 20 k = nfirst, nlast
               if (nn(k) .le. 0) go to 20
               if (lo(k) .ne. j-1) go to 20
               if ((so(k)-pone)*(i-opf) .lt. zero) go to 20
               nmax(i,j) = nn(k)
   20       continue
   30    continue
   40 continue
c
c   Set up hamiltonian matrix for kinetic energy.
c   Only the diagonal depends on the potential.
c
      c2 = -one/b**2
      c1 = -2*one*c2 + one/4
      dk(1) = c1/(r(2)+a)**2
      sd(1) = zero
      sd2(1) = zero
      do 50 i = 3, nr
         dk(i-1) = c1/(r(i)+a)**2
         sd(i-1) = c2/((r(i)+a)*(r(i-1)+a))
         sd2(i-1) = sd(i-1)**2
   50 continue
c
c   Start loop over spin down=1 and spin up=2.
c
      nrm = nr - 1
      do 100 i = 1, 2
c
c   Start loop over s p d... states.
c
         do 90 j = 1, lmax
            if (nmax(i,j) .eq. 0) go to 90
            llp = j*(j-1)
            do 60 k = 2, nr
               if (i .eq. 1) then
                  d(k-1) = dk(k-1) + (viod(j,k)+llp/r(k))/r(k) + vid(k)
               else
                  d(k-1) = dk(k-1) + (viou(j,k)+llp/r(k))/r(k) + viu(k)
               end if
   60       continue
c
c   Diagonalize the matrix.
c
            eps = -one
            call tridib(nrm,eps,d,sd,sd2,bl,bu,1,nmax(i,j),e,ind,ierr,
     &                  rv4,rv5)
            if (ierr .ne. 0) write(6,9000) ierr
 9000       format(/' error in tridib ****** ierr =',i3,/)
            call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,rv1,rv2,
     &                  rv3,rv4,rv5)
            if (ierr .ne. 0) write(6,9010) ierr
 9010       format(/' error in tinvit ****** ierr =',i3,/)
c
c   Save the energy levels and add to charge density.
c
            ki = 1
            kn = 0
            do 80 k = nfirst, nlast
               if (nn(k) .le. 0) go to 80
               if (lo(k) .ne. j-1) go to 80
               if ((so(k)-pone)*(i-opf) .lt. zero) go to 80
               ev(k) = e(ki)
               do 70 l = 2, nr
                  denr = zo(k)*z(kn+l-1)**2/rab(l)
                  if (i .eq. 1) then
                     cdd(l) = cdd(l) + denr
                  else
                     cdu(l) = cdu(l) + denr
                  end if
   70          continue
               ki = ki + 1
               kn = kn + nrm
   80       continue
   90    continue
  100 continue
c
c   End loop over s p and d states.
c
      return
c
      end
