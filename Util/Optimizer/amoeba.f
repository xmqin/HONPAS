! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!------------------------------------------------------------------------------
! module m_amoeba
!   Wrapper for subroutine amoeba
! Used modules
!   use precision,   only: dp
!   use vars_module, only: constrained
! Public procedures provided by this module:
!   amoeba : Nelder-Mead minimization algorithm
! Public derived types:
!   none
! Public variables and parameters:
!   none
! Written by A.Garcia and J.M.Soler, May 2015
!------------------------------------------------------------------------------
! SUBROUTINE amoeba(x,y,ndim,ftol,funk,neval,maxeval)
!   Nelder-Mead minimization algorithm
! Input:
!   integer ndim           ! number of coordinates
!   interface
!     function funk(x) result(value)        ! function to be minimized
!       use precision, only: dp
!       real(dp),intent(in) :: x(:)         ! point coordinates
!       real(dp)            :: value        ! function value at x
!     end function funk
!   end interface
!   real(dp) ftol           ! function tolerance at minimum
! Optional input:
!   integer  maxeval        ! max. function evaluations
! Input/output:
!   real(dp) x(ndim+1,ndim) ! simplex values x(ipoint,icoord)
!   real(dp) y(ndim+1)      ! function values at x
! Output:
!   integer neval           ! number of function evaluations
! Reference:
!   W.H.Press et al, Numerical Recipes, Cambridge U.P.
!------------------------------------------------------------------------------

      MODULE m_amoeba

      use precision, only: dp
      use vars_module, only: constrained

      implicit none

      public :: amoeba
      private

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE amoeba(x,y,ndim,ftol,funk,neval,maxeval)

      ! Arguments
      implicit none
      real(dp),intent(inout):: x(ndim+1,ndim) ! simplex values x(ipoint,icoord)
      real(dp),intent(inout):: y(ndim+1)      ! function values at x
      integer, intent(in)   :: ndim           ! number of coordinates
      interface
        function funk(x) result(value)        ! function to be minimized
          use precision, only: dp
          real(dp),intent(in) :: x(:)         ! point coordinates
          real(dp)            :: value        ! function value at x
        end function funk
      end interface
      real(dp),intent(in)   :: ftol           ! function tolerance at minimum
      integer, intent(out)  :: neval          ! number of function evaluations
      integer, intent(in),optional:: maxeval  ! max. function evaluations

      ! Internal parameters and variables
      integer, parameter :: maxeval_default = 2000
      real(dp),parameter :: tiny = 1.0e-12_dp
      INTEGER :: i,ihi,ilo,inhi,max_eval,npoint
      REAL(dp):: df,x0(ndim),xtry(ndim),yhi,ytry

      ! Set parameters
      npoint=ndim+1                           ! number of simplex points
      if (present(maxeval)) then
         max_eval = maxeval                   ! max. allowed func. evaluations
      else
         max_eval = maxeval_default
      endif

      ! Minimization iteration
      neval=0
      do while (neval<=max_eval)
        ilo=minloc(y,dim=1)                   ! lowest point
        ihi=maxloc(y,dim=1)                   ! highest point
        inhi=maxloc(y,dim=1,mask=(y/=y(ihi))) ! next highest point

        ! If converged, place lowest vertex first and return
        df=2*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+tiny)
        print *,"Fractional dispersion: ", df
        if (df<ftol) then
          x((/1,ilo/),:)=x((/ilo,1/),:)       ! swap 1 and ilo
          y((/1,ilo/))=y((/ilo,1/))
          exit ! do while loop
        endif

        ! First try a reflection
        ytry=trial(-1.0_dp)
        if (ytry<=y(ilo)) then ! lower than best => double the reflected point
          ytry=trial(2.0_dp)
        elseif (ytry>=y(inhi)) then ! worse than 2nd highest => contract
          yhi=y(ihi)
          ytry=trial(0.5_dp)
          if (ytry>=yhi) then  ! contract all points, relative to the lowest
            do i=1,npoint
              if(i/=ilo)then
                x(i,:)=0.5_dp*(x(i,:)+x(ilo,:)) ! factor-2 contraction
                y(i)=funk(x(i,:))
                neval=neval+1
              endif
            enddo ! i
          endif ! (ytry>=yhi)
        endif ! (ytry<=y(ilo))

      enddo ! while (neval<=max_eval). This is the return point

      CONTAINS ! trial is a contained function of amoeba

!------------------------------------------------------------------------------

      FUNCTION trial(fac) result(ytry)

      ! Tries a new position for highest point and substitutes it if better

      REAL(dp),intent(in):: fac  ! expansion factor of highest point,
                                 ! relative to opposite face
      REAL(dp)           :: ytry ! function value at new position

      INTEGER :: j

      ! Find center of face opposite to highest point ihi
      x0=(sum(x,dim=1)-x(ihi,:))/(npoint-1)

      ! Find unconstrained new position of ihi
      xtry=x0+fac*(x(ihi,:)-x0)

      ! Find constrained new position of ihi
      do j=1,ndim
        xtry(j)=constrained(xtry(j),j)
      enddo

      ! Find function at new position
      ytry=funk(xtry)
      neval=neval+1
      print *, "New point: ", xtry, " --- ",  ytry

      ! If succesful, substitute highest point by trial point
      if (ytry<y(ihi)) then
        x(ihi,:)=xtry(:)
        y(ihi)=ytry
      endif

      END function trial

!------------------------------------------------------------------------------

      END subroutine amoeba

      END MODULE m_amoeba

