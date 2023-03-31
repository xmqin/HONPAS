! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine cart2frac(na,xc,yc,zc,rv,xyzfrac)
C
C  Converts Cartesian coordinates to fractional coordinates
C
C  On entry : 
C
C  na            = number of atoms
C  xc(na)        = Cartesian X coordinate
C  yc(na)        = Cartesian Y coordinate
C  zc(na)        = Cartesian Z coordinate
C  rv(3,3)       = cell vectors
C
C  On exit : 
C
C  xyzfrac(3,na) = fractional X coordinate
C
C  Julian Gale, NRI, Curtin University, May 2004
C
      use precision

      implicit none
C
C  Passed variables
C
      integer, intent(in)     :: na
      real(dp),  intent(inout)  :: xc(na)
      real(dp),  intent(inout)  :: yc(na)
      real(dp),  intent(inout)  :: zc(na)
      real(dp),  intent(out)    :: xyzfrac(3,na)
      real(dp),  intent(in)     :: rv(3,3)
C
C  Local variables
C
      integer                 :: i
      real(dp)                  :: rmat(3,3)
C
C  Copy lattice vectors to scratch array
C
      do i = 1,3
        rmat(1,i) = rv(1,i)
        rmat(2,i) = rv(2,i)
        rmat(3,i) = rv(3,i)
      enddo
C
C  Copy Cartesian coordinates to fractional array
C
      do i = 1,na
        xyzfrac(1,i) = xc(i)
        xyzfrac(2,i) = yc(i)
        xyzfrac(3,i) = zc(i)
      enddo
C
C  Convert cartesian coordinates to fractional for 3D case
C
      call gaussxyz(na,rmat,xyzfrac)
C
C  Place fractional coordinates in the range 0 to 1
C
      do i = 1,na
        xyzfrac(1,i) = dmod(xyzfrac(1,i)+10.0d0,1.0d0)
        xyzfrac(2,i) = dmod(xyzfrac(2,i)+10.0d0,1.0d0)
        xyzfrac(3,i) = dmod(xyzfrac(3,i)+10.0d0,1.0d0)
      enddo
C
C  Convert fractional back to Cartesian to place within cell
C
      do i = 1,na
        xc(i) = xyzfrac(1,i)*rv(1,1) + 
     .          xyzfrac(2,i)*rv(1,2) + 
     .          xyzfrac(3,i)*rv(1,3)
        yc(i) = xyzfrac(1,i)*rv(2,1) + 
     .          xyzfrac(2,i)*rv(2,2) + 
     .          xyzfrac(3,i)*rv(2,3)
        zc(i) = xyzfrac(1,i)*rv(3,1) + 
     .          xyzfrac(2,i)*rv(3,2) + 
     .          xyzfrac(3,i)*rv(3,3)
      enddo
C
      return
      end
C
      subroutine gaussxyz(m,a,xf)
C
C  Invert matrix A by Gaussian elimination and apply to
C  multiple vectors x
C
      use precision
      use sys, only : die
      implicit none
C
C  Passed variables
C
      integer, intent(in)     :: m
      real(dp),  intent(inout)  :: a(3,3)
      real(dp),  intent(inout)  :: xf(3,m)
C
C  Local variables
C
      integer                 :: i
      integer                 :: ie
      integer                 :: in
      integer                 :: ix
      integer                 :: j
      integer                 :: k
      integer                 :: kk
      real(dp)                  :: delt
      real(dp)                  :: u
      real(dp)                  :: x
C
      delt = 1.0d-10
      do k = 1,2
        u = abs(a(k,k))
        kk = k + 1
        in = k
C
C  Search for index in of maximum pivot value
C      
        do i = kk,3
          if (abs(a(i,k)).gt.u) then
            u = abs(a(i,k))
            in = i
          endif
        enddo
        if (k.ne.in) then
C
C  Interchange rows k and index in
C
          do j = k,3
            x = a(k,j)
            a(k,j) = a(in,j)
            a(in,j) = x
          enddo
          do j = 1,m
            x = xf(k,j)
            xf(k,j) = xf(in,j)
            xf(in,j) = x
          enddo
        endif
C
C  Check to see if pivot value is too small
C
        if (u.lt.delt) then
          call die('Cell matrix is singular')
        endif
C
C  Forward elimination step
C
        do j = kk,3
          do i = kk,3
            a(i,j) = a(i,j) - a(i,k)*a(k,j)/a(k,k)
          enddo
        enddo
        do j = 1,m
          do i = kk,3
            xf(i,j) = xf(i,j) - a(i,k)*xf(k,j)/a(k,k)
          enddo
        enddo
      enddo
      if (abs(a(3,3)).lt.delt) then
        call die('Cell matrix is singular')
      endif
C
C  Back substitution
C
      do k = 1,m
        xf(3,k) = xf(3,k)/a(3,3)
      	do ie = 1,2
      	  i = 3 - ie
      	  ix = i + 1
      	  do j = ix,3
            xf(i,k) = xf(i,k) - xf(j,k)*a(i,j)
          enddo
      	  xf(i,k) = xf(i,k)/a(i,i)
        enddo
      enddo
C
      return
      end
