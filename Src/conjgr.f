! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_conjgr
      use precision, only: dp

      implicit none

      public :: conjgr
      private
      
      CONTAINS

      subroutine conjgr(N,X,G,dxmax,GTOL,cntrol,H)

C Directs a conjugate-gradient minimization of a function which
C is evaluated by the calling program.
C  N     : INPUT SPACE DIMENSIONALITY
C  X     : INPUT POSITION AT WHICH GRADIENT HAS BEEN EVALUATED AND
C          OUTPUT NEW POSITION AT WHICH GRADIENT MUST BE EVALUATED NEXT
C  G     : INPUT GRADIENT (WITH A MINUS SIGN, OR ACTIVATE A LINE BELOW)
C  dxmax : INPUT MAXIMUM ALLOWED DISPLACEMENT IN EACH COORDINATE
C  GTOL  : INPUT MAXIMUM FINAL VALUE OF EACH GRADIENT COMPONENT
C  cntrol: CONTROL ARRAY. FIRST ELEMENT MUST BE MADE ZERO BEFORE
C          FIRST CALL. IF IT IS ZERO ON OUTPUT, MINIMIZATION IS
C          CONVERGED. OTHERWISE, CALCULATE GRADIENT AT THE NEW
C          POSITION AND CALL AGAIN THIS ROUTINE. DO NOT MODIFY ANY
C          ARGUMENT OTHER THAN G BETWEEN CALLS.
C  H     : AUXILIARY ARRAY WHICH MUST NOT BE MODIFIED BETWEEN CALLS
C  IOPT  : PARAMETER BELOW WHICH DETERMINES METHOD USED AND AUXILIARY
C          STORAGE REQUIRED: IOP=1 => FLETCHER-REEVES. IOPT=2 =>
C          POLAK-RIBIERE. DETAILS IN SECT. 10.6 OF 'NUMERICAL RECIPES'
C WRITTEN BY J.SOLER. JAN/91. BASED ON ROUTINES IN 'NUMERICAL RECIPES'
C Modified by Rainer Hoft Mar/05 to allow for an array of gradient
C    tolerances and maximum displacements (a tolerance and displacement
C    for each variable)

C Local parameters
      integer,  parameter     :: iopt = 2

C Passed variables
      integer,  intent(in)    :: n

      real(dp), intent(in)    :: dxmax(n), gtol(n), G(N)

      real(dp), intent(inout) :: X(N),H(N,IOPT),cntrol(0:19)

C Local variables
      real(dp)                :: gg, gamma
      integer                 :: i, j

      real(dp)                :: ddot
      external  ddot

C If gradient is smaller than tolerence, return
      do j = 1,n
        if (abs(g(j)).gt.gtol(j)) then
          goto 70
        endif
      enddo
      cntrol(0) = 0
      goto 60
  70  continue

C First-call initializations
      if (nint(cntrol(0)).eq.0) then
        do I = 1,IOPT
          do J = 1,N
            H(J,I) = G(J)
          enddo
        enddo
        cntrol(0) = 1
        cntrol(1) = 1
        cntrol(2) = ddot(n,G,1,G,1)
        cntrol(10) = 0

        cntrol(18) = sqrt(ddot(n,dxmax,1,dxmax,1))

      endif

C Line minimization is always called
   40 call linmin1(N,X,H,G,dxmax,cntrol(10:19))

C If line minimization is finished, find new line direction
      if (nint(cntrol(10)).eq.0) then
        GG = ddot(n,G,1,G,1)
        if (IOPT.EQ.2) GG = GG - ddot(n,G,1,H(1,2),1)
        GAMMA = GG/cntrol(2)
        do J = 1,N
          H(J,1) = G(J) + GAMMA*H(J,1)
          IF (IOPT.EQ.2) H(J,2) = G(J)
        enddo
        cntrol(1) = cntrol(1) + 1
        cntrol(2)  = ddot(n,G,1,G,1)
*       WRITE(6,'(A,I4,F15.6)')
*    .     ' CONJGR: NEW LINE DIRECTION. N,DX=',N,cntrol(18)
        goto 40
      endif

   60 continue
      end subroutine conjgr
!-----------------------------------------------------------------
      subroutine linmin1(n,xvec,hvec,gvec,dxmax,cntrol)

      integer, intent(in)       ::  n
      real(dp), intent(inout)   ::  xvec(N),cntrol(0:9)
      real(dp), intent(in)      ::  hvec(N),gvec(N)
      real(dp), intent(in)      ::  dxmax(n)

      real(dp), parameter :: FACTOR=1.6D0

      integer  :: i, icntrl
      real(dp) :: x1, x2, y1, y2, x0, hmod, dx, x, y
      
      real(dp) :: dxhmin   ! = min_i(dxmax(i)/|h(i)|)

      real(dp) :: ddot
      external  ddot

C Translate control parameters
      icntrl = nint(cntrol(0))
      X1 = cntrol(1)
      X2 = cntrol(2)
      Y1 = cntrol(3)
      Y2 = cntrol(4)
      X0 = cntrol(5)
      hmod = cntrol(6)

      dxhmin = cntrol(7)

      DX = cntrol(8)
      X = X0
      Y = ddot(n,gvec,1,hvec,1)
*     WRITE(6,'(A,I4,2F12.6)') ' LINMIN: icntrl,X,Y=',icntrl,X,Y

      if (icntrl.eq.0) then
C Initialize x1,y1 on first call
        X1 = 0.0_dp
        Y1 = Y
C Prepare second point
        icntrl = 1
        X0 = 0.0_dp

        if (dx.eq.0.d0) then
          dx = sqrt(ddot(n,dxmax,1,dxmax,1))
        endif

        hmod = sqrt(ddot(n,hvec,1,hvec,1))

        dxhmin = abs(dxmax(1)/hvec(1))
        do 10 i = 2,n
          if (abs(hvec(i)).gt.1d-12) then
            dxhmin = min(dxhmin,abs(dxmax(i)/hvec(i)))
          endif
   10   continue
        if (hmod.gt.1d-12) then
          x = min(dx/hmod,dxhmin)
        else
          x = 0.0_dp
        endif

        goto 20
      elseif (icntrl.eq.1) then
C Initialize x2,y2 on second call
        X2=X
        Y2=Y
        icntrl=2
      elseif (icntrl.eq.2) then
C Shift interval using new point
        X1=X2
        Y1=Y2
        X2=X
        Y2=Y
      elseif (icntrl.eq.3) then
C If root was found in last call, all is done now
        icntrl=0
        goto 20
      endif

      if (Y2.gt.0.0_dp) then
C Root not bracketed yet. try new right bracket

        x = x2+min(factor*(x2-x1),dxhmin)

      else
C Interpolate for root and return to calculate last gradient
        X = (X1*Y2-X2*Y1)/(Y2-Y1)
        icntrl = 3
      endif

C Store control parameters and set new point
   20 cntrol(0) = icntrl
      cntrol(1) = X1
      cntrol(2) = X2
      cntrol(3) = Y1
      cntrol(4) = Y2
      cntrol(5) = X
      cntrol(6) = hmod

      cntrol(7) = dxhmin

      cntrol(8) = abs(X)*hmod
      do I = 1,N
        xvec(I) = xvec(I) + hvec(I)*(X-X0)
      enddo

      end subroutine linmin1

      end module m_conjgr
