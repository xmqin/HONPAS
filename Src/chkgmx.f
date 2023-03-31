! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE

      MODULE m_chkgmx

! Used module procedures
      use sys,       only: die     ! Termination routine
      use m_minvec,  only: minvec  ! Finds a lattice basis of minimal length
      use cellsubs,  only: reclat  ! Finds reciprocal unit cell

! Used module parameters and variables
      use precision, only: dp      ! Double precision real kind
!      use parallel,  only: Node    ! My processor index

      implicit none

! Public procedures provided:
      PUBLIC::
     .  chkgmx,  ! Checks that a given cutoff is consistent with a mesh
     .  meshKcut ! Returns the planewave cutoff of a periodic mesh

! Public parameters, variables, and arrays:
!     none

      PRIVATE ! Nothing is declared public beyond this point

      CONTAINS
!---------------------------------------------------------------------------

      SUBROUTINE CHKGMX (K,BG,MESH,G2MAX)

      implicit none
      real(dp),intent(in) :: K(3)    ! Vector in reciprocal unit cell (BZ)
      real(dp),intent(in) :: BG(3,3) ! Reciprocal lattice vectors BG(ixyz,ivec)
      integer, intent(in) :: MESH(3) ! Number of mesh point divisions of each 
                                     ! real-space lattice vector.
      real(dp),intent(inout):: G2MAX ! Planewave cutoff, which is reduced to
                                     ! the maximum value supported by the mesh
                                     ! (without aliasing), if this is smaller.

      real(dp), parameter :: ZERO=0.0_dp,HALF=0.5_dp,
     .                       TINY=1.0e-8_dp,BIG=1.0e20_dp

      integer :: i, j, i1, i2, i3
      real(dp):: bm(3,3), baux(3,3), ctransf(3,3),
     .           g(3), gmod, gmax, g2msh, r

      DO I=1,3
        DO J=1,3
          BM(J,I)=BG(J,I)*MESH(I)
          Baux(J,I)=BG(J,I)*MESH(I)
        ENDDO
      ENDDO
!
!     Use Baux to avoid aliasing of arguments, as one is intent(in)
!     and the other intent(out)...
!
      CALL MINVEC (Baux,BM, ctransf)
      GMAX=BIG
      DO I3=-1,1
        DO I2=-1,1
          DO I1=-1,1
            IF (I1.EQ.0.AND.I2.EQ.0.AND.I3.EQ.0) GO TO 20
              G(1)=BM(1,1)*I1+BM(1,2)*I2+BM(1,3)*I3
              G(2)=BM(2,1)*I1+BM(2,2)*I2+BM(2,3)*I3
              G(3)=BM(3,1)*I1+BM(3,2)*I2+BM(3,3)*I3
              GMOD=SQRT(SUM(G**2))
              R=HALF*GMOD-SUM(K*G)/GMOD
              GMAX=MIN(GMAX,R)
  20        CONTINUE
          ENDDO
        ENDDO
      ENDDO
      IF (GMAX.LT.ZERO) call die('CHKGMX: K NOT IN FIRST BZ')
      G2MSH=GMAX*GMAX-TINY
      IF (G2MSH.LT.G2MAX) THEN
!       if (Node.eq.0) then
!         WRITE(6,*) 'CHKGMX: MESH TOO SPARSE. GMAX HAS BEEN REDUCED'
!         WRITE(6,*) 'CHKGMX: OLD G2MAX =',G2MAX
!         WRITE(6,*) 'CHKGMX: NEW G2MAX =',G2MSH
!       endif
        G2MAX=G2MSH
      ENDIF
      END SUBROUTINE CHKGMX

!----------------------------------------------------------------------

      FUNCTION meshKcut( cell, nMesh )

      ! Finds mesh planewave cutoff

      use precision, only: dp

      implicit none
      real(dp),intent(in):: cell(3,3) ! Unit cell vectors cell(iXYZ,iVector)
      integer, intent(in):: nMesh(3)  ! Mesh divisions of each vector
      real(dp)           :: meshKcut  ! Mesh wavevector cutoff

      real(dp):: g2max, k0(3), kcell(3,3)

      ! Find reciprocal cell vectors
      call reclat( cell, kcell, 1 )

      ! Initialize arguments to call chkgmx
      k0(1:3) = 0._dp
      g2max = huge(g2max)

      ! chkgmx will reduce g2max to square of mesh cutoff
      call chkgmx( k0, kcell, nMesh, g2max )

      ! Return wavevector cutoff
      meshKcut = sqrt(g2max)

      END FUNCTION meshKcut

      END MODULE m_chkgmx

