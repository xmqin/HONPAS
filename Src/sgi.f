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

C  NEXT LINES FOR IBM SYSTEMS-370 (ASSEMBLE ROUTINE TASKTM REQUIRED)
*     CALL TASKTM (ITIME)
*     TIME = 1.D-4 * ITIME

C  NEXT LINES FOR IBM-3090
*     CALL CPUTIME (TIME,RCODE)
*     TIME = 1.D-6 * TIME

C  NEXT LINE FOR 1
*     TIME = SECOND()

C  NEXT LINE FOR SGI, SUN OR DEC WORKSTATIONS
      TIME = ETIME(TIMES)

C  NEXT LINE FOR IBM RS/6000 WORKSTATIONS
*     TIME = MCLOCK()*0.01D0

C  NEXT LINE FOR ANYTHING ELSE
*     TIME = 0.D0
      RETURN

      END
