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
C $Id: surpla.f,v 1.1 2003/06/20 15:57:22 javier Exp $

      DOUBLE PRECISION FUNCTION SURPLA( C )

C  CALCULATES THE SRFACE OF THE UNIT CELL NORMAL TO THE INTERFACE

      DOUBLE PRECISION C(3,3)
      SURPLA = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) **2 +
     .         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) **2 +
     .         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) **2
      SURPLA = SQRT( ABS( SURPLA ) )
      END
