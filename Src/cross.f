! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
C $Id: cross.f,v 1.2 1999/01/31 10:53:49 emilio Exp $

      SUBROUTINE CROSS(A,B,AXB)
C Finds the cross product AxB of vectors A and B
C Written by J.M.Soler
      IMPLICIT NONE
      DOUBLE PRECISION A(3),B(3),AXB(3)
      AXB(1)=A(2)*B(3)-A(3)*B(2)
      AXB(2)=A(3)*B(1)-A(1)*B(3)
      AXB(3)=A(1)*B(2)-A(2)*B(1)
      END
