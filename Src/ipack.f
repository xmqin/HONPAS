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


      SUBROUTINE IPACK( TASK, ND, N, I, IND )

C *********************************************************************
C Packs/unpacks several integer indexes into/out of one
C Written by J.M.Soler. July 1998
C ******** INPUT ******************************************************
C INTEGER TASK : Task switch: TASK=+1 => Pack I into IND
C                             TASK=-1 => Unpack I out of IND
C INTEGER ND   : Number of indexes to pack/unpack (dimension of I)
C INTEGER N(ND): Upper range of I. Must be N(j)>0 for all j.
C ******** INPUT or OUTPUT (Depending of TASK) ************************
C INTEGER I(ND): Indexes to pack/unpack
C INTEGER IND  : Combined (packed) index such that
C                IND = 1 + I(1) + N(1)*I(2) + N(1)*N(2)*I(3) + ...
C                where 0 <= I(j) < N(j)
C ******** BEHAVIOR ***************************************************
C Notice that the range of I(j) begins at 0, and that of IND at 1
C A modulus operation is done to bring I(j) to the range (0:N(j)-1)
C Does not check that N(j) > 0
C If TASK=-1 and IND is not in the range (1:N(1)*N(2)*...*N(ND)),
C   the program stops with an error message
C *********************************************************************
      
      use sys, only: die

      IMPLICIT NONE
      INTEGER  IND, ND, TASK
      INTEGER  I(ND), N(ND)

      INTEGER  IJ, J

      IF (TASK .GT. 0) THEN
        IND = 0
        DO J = ND,1,-1
          IJ = MOD( I(J)+1000*N(J), N(J) )
          IND = IJ + N(J) * IND
        ENDDO
        IND = 1 + IND
      ELSE
        IF (IND.LT.1) call die('IPACK: ERROR: IND < 1')
        IJ = IND - 1
        DO J = 1,ND
          I(J) = MOD( IJ, N(J) )
          IJ = (IJ-I(J)) / N(J)
        ENDDO
        IF (IJ.GT.0) call die('IPACK: ERROR: IND out of range')
      ENDIF

      END

