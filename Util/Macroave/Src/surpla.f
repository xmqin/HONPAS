! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
C $Id: surpla.f,v 1.1 2003/06/20 15:57:22 javier Exp $

      DOUBLE PRECISION FUNCTION SURPLA( C )

C  CALCULATES THE SRFACE OF THE UNIT CELL NORMAL TO THE INTERFACE

      DOUBLE PRECISION C(3,3)
      SURPLA = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) **2 +
     .         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) **2 +
     .         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) **2
      SURPLA = SQRT( ABS( SURPLA ) )
      END
