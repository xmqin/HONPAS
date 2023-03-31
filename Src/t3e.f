! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE CPUTIM (TIME)

C  RETURNS CPU TIME IN SECONDS SINCE PROGRAM START
C  WRITEN BY J.SOLER (JSOLER AT EMDUAM11)

      DOUBLE PRECISION TIME
      REAL TIMES(2)

C  NEXT LINE FOR CRAY
      TIME = SECOND()

      END
