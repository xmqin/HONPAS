! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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
