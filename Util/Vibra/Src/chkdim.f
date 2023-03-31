! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE CHKDIM (SUB,VAR,ND,N,IOPT)
C **********************************************************************
C Checks a dimension parameter, warning and stopping if it is too small
C Written by J.M.Soler. 1987.
C ************ INPUT ***************************************************
C CHARACTER SUB*(*) : Name of calling subroutine
C CHARACTER VAR*(*) : Name of dimension parameter
C INTEGER   ND      : Actual value of dimension parameter
C INTEGER   N       : Required value of dimension parameter
C INTEGER   IOPT    : Option switch: IOPT=0 => Require that ND.EQ.N
C                                    IOPT=1 => Require that ND.GE.N
C **********************************************************************
      integer :: nd,n,iopt
      CHARACTER SUB*(*),VAR*(*)
      IF ( IOPT.EQ.0 ) THEN
        IF ( ND.EQ.N ) RETURN
          WRITE (6,'(/5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .      ', dimension ',VAR,' =',ND,'. It must be exactly ',N
          STOP
      ELSE
        IF ( ND.GE.N ) RETURN
          WRITE (6,'(/5A,I8,A,I8)') 'chkdim: ERROR: In ', SUB,
     .      ', dimension ',VAR,' =',ND,'. It must be at least ',N
          STOP
      ENDIF
      END



      SUBROUTINE CHKDIME (ND,N,OVERFLOW,NM)
      integer :: nd, n, nm
      LOGICAL OVERFLOW
      NM = MAX(N,NM)
      IF ( ND.GE.N ) RETURN
      OVERFLOW = .TRUE.
      RETURN
      END

