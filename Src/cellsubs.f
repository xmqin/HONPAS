! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE

      MODULE cellSubs

      USE precision, only: dp   ! Double precision real kind

      PUBLIC::  
     .  reclat, ! Reciprocal unit cell vectors
     .  volcel  ! Volume of the unit cell

      PRIVATE

      CONTAINS

!--------------------------------------------------------------------

      SUBROUTINE RECLAT (A,B,IOPT)

C  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
C  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      integer :: iopt, i
      real(dp):: A(3,3),B(3,3), pi, c , ci
      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO 20 I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
  20  CONTINUE
      END SUBROUTINE RECLAT

!--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION VOLCEL( C )

C  CALCULATES THE VOLUME OF THE UNIT CELL

      real(dp) C(3,3)
      VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) +
     .         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) +
     .         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL = ABS( VOLCEL )
      END FUNCTION VOLCEL

      END MODULE cellSubs

