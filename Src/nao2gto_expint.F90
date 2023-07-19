MODULE m_expint_cp2k
  USE kinds,                           ONLY: dp
  implicit none
  
  PRIVATE        
  PUBLIC expint_cp2k

  CONTAINS


  !USE mathlib,                         ONLY: expint_cp2k
! *****************************************************************************
!> \brief computes the exponential integral
!>      En(x) = Int(exp(-x*t)/t^n,t=1..infinity)  x>0, n=0,1,..
!>      Note: Ei(-x) = -E1(x)
!> \par History
!>      05.2007 Created
!> \author Manuel Guidon (adapted from Numerical recipies)
! *****************************************************************************
  FUNCTION expint_cp2k(n,x)
    INTEGER                                  :: n
    REAL(dp)                                 :: x, expint_cp2k

    INTEGER, PARAMETER                       :: maxit = 100
    REAL(dp), PARAMETER :: eps = 6.e-14_dp, &
      euler = 0.5772156649015328606065120_dp, fpmin = TINY(0.0_dp)

    INTEGER                                  :: i, ii, nm1
    REAL(dp)                                 :: a, b, c, d, del, fact, h, psi

    nm1=n-1

    IF(n.lt.0.OR.x.lt.0.0_dp.OR.(x.eq.0.0_dp.AND.(n.EQ.0.or.n.EQ.1))) THEN
      write(6,*) 'Invalid argument'
    ELSE IF(n.EQ.0) THEN       !Special case.
      expint_cp2k=EXP(-x)/x
    ELSE IF(x.EQ.0.0_dp) THEN  !Another special case.
      expint_cp2k=1.0_dp/nm1
    ELSE IF(x.GT.1.0_dp) THEN  !Lentz’s algorithm (§5.2).
      b=x+n
      c=1.0_dp/FPMIN
      d=1.0_dp/b
      h=d
      DO i = 1,MAXIT
        a=-i*(nm1+i)
        b=b+2.0_dp
        d=1.0_dp/(a*d+b)
        c=b+a/c
        del=c*d
        h=h*del
        IF(ABS(del-1.0_dp).LT.EPS) THEN
          expint_cp2k=h*EXP(-x)
          RETURN
        END IF
      END DO
        write(6,*) 'continued fraction failed in expint_cp2k'
    ELSE !Evaluate series.
      IF(nm1.NE.0)THEN  !Set first term.
        expint_cp2k=1.0_dp/nm1
      ELSE
        expint_cp2k=-LOG(x)-euler
      END IF
      fact=1.0_dp
      DO i=1,MAXIT
        fact=-fact*x/i
        IF(i.NE.nm1) THEN
          del=-fact/(i-nm1)
        ELSE
          psi=-euler !Compute ψ(n).
          DO ii=1,nm1
            psi=psi+1.0_dp/ii
          END DO
          del=fact*(-LOG(x)+psi)
        END IF
        expint_cp2k=expint_cp2k+del
        IF(ABS(del).LT.ABS(expint_cp2k)*EPS) RETURN
      END DO
        write(6,*) 'series failed in expint_cp2k'
    END IF
    RETURN
  END FUNCTION expint_cp2k

end module m_expint_cp2k
