! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine uncell(rv,a,b,c,alpha,beta,gamma,degtorad)
C
C  Convert cell vectors to parameters
C
C  Julian Gale, April 1992
C
      use precision
      implicit none
C
C  Passed variables
C
      real(dp)         :: a
      real(dp)         :: b
      real(dp)         :: c
      real(dp)         :: alpha
      real(dp)         :: beta
      real(dp)         :: gamma
      real(dp)         :: degtorad
      real(dp)         :: rv(3,3)
C
C  Local variables
C
      integer          :: i
      integer          :: j
      real(dp)         :: temp(6)
C
      do i = 1,3
        temp(i) = 0.0_dp
        do j = 1,3
          temp(i) = temp(i) + rv(j,i)**2
        enddo
        temp(i) = sqrt(temp(i))
      enddo
      a = abs(temp(1))
      b = abs(temp(2))
      c = abs(temp(3))
      do i = 1,3
        temp(3+i) = 0.0_dp
      enddo
      do j = 1,3
        temp(4) = temp(4) + rv(j,2)*rv(j,3)
        temp(5) = temp(5) + rv(j,1)*rv(j,3)
        temp(6) = temp(6) + rv(j,1)*rv(j,2)
      enddo
      temp(4) = temp(4)/(temp(2)*temp(3))
      temp(5) = temp(5)/(temp(1)*temp(3))
      temp(6) = temp(6)/(temp(1)*temp(2))
      alpha = acos(temp(4))/degtorad
      beta  = acos(temp(5))/degtorad
      gamma = acos(temp(6))/degtorad
C
C  Avoid round off errors for 90.0 and 120.0 degrees
C
      if (abs(alpha-90.0).lt.0.00001) alpha = 90.0_dp
      if (abs(alpha-120.0).lt.0.00001) alpha = 120.0_dp
      if (abs(beta-90.0).lt.0.00001) beta = 90.0_dp
      if (abs(beta-120.0).lt.0.00001) beta = 120.0_dp
      if (abs(gamma-90.0).lt.0.00001) gamma = 90.0_dp
      if (abs(gamma-120.0).lt.0.00001) gamma = 120.0_dp
C
      return
      end
C
      subroutine cellimagelist(ucell,xvec1cell,yvec1cell,zvec1cell,
     .  ixvec1cell,iyvec1cell,izvec1cell)
C
C  Store linear array of lattice vectors for ii/jj/kk=-1,1 
C  => 27 lattice vectors. This tidies up and speeds up the
C  loops over these indices.
C
C  Julian Gale, NRI, Curtin University, March 2004
C
      use precision
      implicit none
C
C  Passed variables
C
      integer            :: ixvec1cell(27)
      integer            :: iyvec1cell(27)
      integer            :: izvec1cell(27)
      real(dp)           :: ucell(3,3)
      real(dp)           :: xvec1cell(27)
      real(dp)           :: yvec1cell(27)
      real(dp)           :: zvec1cell(27)
C
C  Local variables
C
      integer            :: ii
      integer            :: iimax
      integer            :: jj
      integer            :: kk
      real(dp)           :: xcdi
      real(dp)           :: ycdi
      real(dp)           :: zcdi
      real(dp)           :: xcdj
      real(dp)           :: ycdj
      real(dp)           :: zcdj
      real(dp)           :: xcrd
      real(dp)           :: ycrd
      real(dp)           :: zcrd
C
!! Why this -2 ?
!! (and the rest below: Simply to modify 
!!  the indexes)
!! The final ordering is:
!!   1 --> (-1,-1,-1)
!!   2 --> (-1,-1, 0)
!!   3 --> (-1,-1, 1)
!!   4 --> (-1, 0,-1)
!    ...
!    So it might be clearer to do this:
c$$$      n = 0
c$$$      a1 => ucell(:,1) ... etc
c$$$      do kk = -1, 1
c$$$         do jj = -1, 1
c$$$            do ii = -1, 1
c$$$               n = n + 1
c$$$               vlin(1:3) = ii*a1 + jj*a2 + kk*a3
c$$$               xvec1(n)  = vlin(1) ... etc
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$            
!!
!! vector cdi = -2*a1
      xcdi = -2.0d0*ucell(1,1)
      ycdi = -2.0d0*ucell(2,1)
      zcdi = -2.0d0*ucell(3,1)
      iimax = 0
C
C  Loop over unit cells
C
      do ii = -1,1
        ! cdi = cdi + a1
        xcdi = xcdi + ucell(1,1)
        ycdi = ycdi + ucell(2,1)
        zcdi = zcdi + ucell(3,1)
        ! cdj = cdi - 2*a2
        xcdj = xcdi - 2.0d0*ucell(1,2)
        ycdj = ycdi - 2.0d0*ucell(2,2)
        zcdj = zcdi - 2.0d0*ucell(3,2)
        do jj = -1,1
          ! cdj = cdj + a2
          xcdj = xcdj + ucell(1,2)
          ycdj = ycdj + ucell(2,2)
          zcdj = zcdj + ucell(3,2)
          ! crd = cdj - 2*a3
          xcrd = xcdj - 2.0d0*ucell(1,3)
          ycrd = ycdj - 2.0d0*ucell(2,3)
          zcrd = zcdj - 2.0d0*ucell(3,3)
          do kk = -1,1
            iimax = iimax + 1
            ! crd = crd + a3
            xcrd = xcrd + ucell(1,3)
            ycrd = ycrd + ucell(2,3)
            zcrd = zcrd + ucell(3,3)
            ixvec1cell(iimax) = ii
            iyvec1cell(iimax) = jj
            izvec1cell(iimax) = kk
            ! Store successive points in linear arrays
            ! vec1(iimax) = crd
            xvec1cell(iimax) = xcrd
            yvec1cell(iimax) = ycrd
            zvec1cell(iimax) = zcrd
          enddo
        enddo
      enddo
C
      return
      end
