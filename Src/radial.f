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
      module radial

      use precision
      use xml

      implicit none

      private

      public :: rad_alloc, rad_get, rad_setup_d2, rad_zero
      public :: radial_read_ascii, radial_dump_ascii
      public :: radial_dump_xml

      type, public :: rad_func
         integer          n
         double precision cutoff         
         double precision delta
         double precision, dimension(:), pointer :: f   ! Actual data
         double precision, dimension(:), pointer :: d2  ! Second derivative
      end type rad_func

      private :: splint, spline

      CONTAINS

      subroutine rad_alloc(func,n)
!
!     Sets the 'size' n of the arrays and allocates f and d2.
!
      type(rad_func), intent(inout)    :: func
      integer, intent(in)        :: n
      func%n = n
      allocate(func%f(n),func%d2(n))
      end subroutine rad_alloc

      subroutine rad_get(func,r,fr,dfdr)
      type(rad_func), intent(in) :: func
      real(dp), intent(in)         :: r
      real(dp), intent(out)        :: fr
      real(dp), intent(out)        :: dfdr

      if (func%n .eq. 0) then
          fr = 0._dp
          dfdr = 0._dp
       else
          call splint(func%delta,func%f,func%d2,func%n,r,fr,dfdr)
       endif
      
      end subroutine rad_get
!
!     Set up second derivative in a radial function
!
      subroutine rad_setup_d2(func,yp1,ypn)
      type(rad_func), intent(inout) :: func
      real(dp), intent(in)          :: yp1, ypn

      if (func%n .eq. 0) return
      call spline(func%delta,func%f,func%n,yp1,ypn,func%d2)
      
      end subroutine rad_setup_d2

      subroutine rad_zero(func)
      type(rad_func), intent(inout) :: func
      func%n      = 0
      end subroutine rad_zero
!
!     Do not use yet... interface in need of fuller specification
!
      function rad_rvals(func) result (r)
      real(dp), dimension(:), pointer :: r
      type(rad_func), intent(in) :: func

      integer i

      nullify(r)
      if (func%n .eq. 0) return
      allocate(r(func%n))
      do i=1,func%n
         r(i) = func%delta *(i-1)
      enddo
      end function rad_rvals


      SUBROUTINE SPLINT(DELT,YA,Y2A,N,X,Y,DYDX) 
C Cubic Spline Interpolation.
C Adapted from Numerical Recipes for a uniform grid.

      implicit none

      real(dp), intent(in)  :: delt
      real(dp), intent(in)  :: ya(:), y2a(:)
      integer, intent(in)   :: n
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: y
      real(dp), intent(out) :: dydx
      
      integer nlo, nhi
      real(dp) a, b

      NLO=INT(X/DELT)+1
      NHI=NLO+1

      A=NHI-X/DELT-1
      B=1.0_DP-A

      Y=A*YA(NLO)+B*YA(NHI)+
     *      ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DELT**2)/6._DP

      DYDX=(YA(NHI)-YA(NLO))/DELT +
     * (-((3*(A**2)-1._DP)*Y2A(NLO))+
     * (3*(B**2)-1._DP)*Y2A(NHI))*DELT/6._DP

      end subroutine splint

      SUBROUTINE SPLINE(DELT,Y,N,YP1,YPN,Y2) 

C Cubic Spline Interpolation.
C Adapted from Numerical Recipes routines for a uniform grid
C D. Sanchez-Portal, Oct. 1996.
! Alberto Garcia, June 2000

      implicit none

      integer, intent(in)    :: n
      real(dp), intent(in)   :: delt
      real(dp), intent(in)   :: yp1, ypn
      real(dp), intent(in)   :: Y(:)
      real(dp), intent(out)  :: Y2(:)
 
      integer i, k
      real(dp) sig, p, qn, un

      real(dp)  :: u(n)            ! Automatic array

      IF (YP1.eq. huge(1._dp)) THEN
        Y2(1)=0.0_dp
        U(1)=0.0_dp
      ELSE
        Y2(1)=-0.5_dp
        U(1)=(3.0_DP/DELT)*((Y(2)-Y(1))/DELT-YP1)
      ENDIF

      DO I=2,N-1
        SIG=0.5_DP
        P=SIG*Y2(I-1)+2.0_DP
        Y2(I)=(SIG-1.0_DP)/P
        U(I)=(3.0_DP*( Y(I+1)+Y(I-1)-2.0_DP*Y(I) )/(DELT*DELT)
     $      -SIG*U(I-1))/P
      ENDDO

      IF (YPN.eq.huge(1._dp)) THEN
        QN=0.0_DP
        UN=0.0_DP
      ELSE
        QN=0.5_DP
        UN=(3.0_DP/DELT)*(YPN-(Y(N)-Y(N-1))/DELT)
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1._DP)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      enddo

      END subroutine spline

      subroutine radial_read_ascii(op,lun,yp1,ypn)
      type(rad_func)    :: op 
      real(dp), intent(in)          :: yp1, ypn

      integer lun
      integer j, npts
      real(dp) dummy

      read(lun,*) npts, op%delta, op%cutoff
      call rad_alloc(op,npts)
      do j=1,npts
         read(lun,*) dummy, op%f(j)
      enddo
      call rad_setup_d2(op,yp1,ypn)
      end subroutine radial_read_ascii
!
!--------------------------------------------------------------------
      subroutine radial_dump_ascii(op,lun,header)
      type(rad_func)    :: op
      integer           :: lun
      logical, intent(in), optional :: header

      integer :: j
      logical :: print_header
!
!     The standard dump is to unit "lun"
!     and includes a header with npts, delta, and cutoff
!
      print_header = .true.
      if (present(header)) then
         print_header = header
      endif
!
      if (print_header) then
         write(lun,'(i4,2g22.12,a)') op%n,
     $        op%delta, op%cutoff, " # npts, delta, cutoff"
      endif
      do j=1,op%n
         write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
      enddo

      end subroutine radial_dump_ascii
!--------------------------------------------------------------------
!
      subroutine radial_dump_xml(op,lun)

      type(rad_func)    :: op
      integer lun
      integer j

      write(lun,'(a)') '<radfunc>'
      call xml_dump_element(lun,'npts',str(op%n))
      call xml_dump_element(lun,'delta',str(op%delta))
      call xml_dump_element(lun,'cutoff',str(op%cutoff))
      write(lun,'(a)') '<data>'
      do j=1,op%n
         write(lun,'(2g22.12)') (j-1)*op%delta, op%f(j)
      enddo
      write(lun,'(a)') '</data>'
      write(lun,'(a)') '</radfunc>'
      end subroutine radial_dump_xml
!
      end module radial










