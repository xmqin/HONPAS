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
      SUBROUTINE CPUTIM (TIME)

C  RETURNS CPU TIME IN SECONDS SINCE PROGRAM START
C  WRITEN BY J.SOLER (JSOLER AT EMDUAM11)

      DOUBLE PRECISION TIME
      REAL TIMES(2)

C  NEXT LINE FOR CRAY
      TIME = SECOND()

      END
