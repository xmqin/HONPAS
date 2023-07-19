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
      logical function propor( N, A, B, TOL, AOVERB )
C **********************************************************************
C Checks if two vectors are proportional within a tolerance.
C Written by J.M.Soler. August 1996.
C **********************************************************************

      use precision
      use sys

      implicit none

      integer           i, imax, n
      real(dp)          a(n), aoverb, b(n), bmax, tol

      bmax = 0.0d0
      imax = 0
      do i = 1,n
        if ( abs(b(i)) .gt. bmax ) then
          imax = i
          bmax = abs(b(i))
        endif
      enddo
      if (imax .eq. 0) call die("propor: ERROR:  IMAX = 0")

      propor = .true.
      if (bmax .eq. 0.0d0) then
        aoverb = 0.0d0
        do i = 1,n
          if ( abs(a(i)) .gt. tol ) then
            propor = .false.
            return
          endif
        enddo
      else
        aoverb = a(imax) / b(imax)
        do i = 1,n
          if ( abs(a(i)-b(i)*aoverb) .gt. tol ) then
            propor = .false.
            return
          endif
        enddo
      endif
      end



