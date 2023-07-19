module flib_spline
!
! Spline interpolation, based on code in "Numerical Recipes"
!
use precision, only: dp, sp
private

public :: generate_spline
interface generate_spline
  module procedure generate_spline_sp
  module procedure generate_spline_dp
end interface generate_spline

public :: evaluate_spline
interface evaluate_spline
  module procedure evaluate_spline_sp
  module procedure evaluate_spline_dp
end interface evaluate_spline


private :: generate_spline_sp
private :: generate_spline_dp
private :: evaluate_spline_sp
private :: evaluate_spline_dp

CONTAINS  !=================

SUBROUTINE generate_spline_sp(X,Y,N,Y2,YP1,YPN,stat)
!
! Given arrays x and y of size n, this routine computes
! an array y2 holding information about the second derivative
! of the interpolating spline. YP1 and YPN, if present, are the
! requested values of the first derivative of the spline at the
! end points. If not given, a "natural" spline with zero second
! derivative at the extremes is constructed.
! The array y2 is needed to evaluate the spline, but it only needs
! to be computed once.
!
integer, intent(in)                          :: n
real(kind=sp), dimension(:), intent(in)  :: x, y
real(kind=sp), dimension(:), intent(out) :: y2
real(kind=sp), intent(in), optional      :: yp1, ypn
integer, intent(out), optional               :: stat

real(kind=sp)  ::  P,QN,SIG,UN
INTEGER            ::  I,K
integer            ::  flag

real(kind=sp), parameter  :: zero = 0.0_sp, &
                                 half = 0.5_sp, &
                                 one  = 1.0_sp, &
                                 two  = 2.0_sp, &
                                 three= 3.0_sp, &
                                 six  = 6.0_sp

real(kind=sp), dimension(n)  ::   U    ! Automatic array

!--------------------------------------------------------
! Check that the nodes are different
!
flag = 0
do i =2, n
  if (x(i) == x(i-1)) then
   flag	= -1
   if (present(stat)) stat = flag
   RETURN
  endif
enddo
	
IF (present(YP1)) THEN
    Y2(1) = -half
    U(1) = (three/ (X(2)-X(1)))* ((Y(2)-Y(1))/ (X(2)-X(1))-YP1)
ELSE
    Y2(1) = zero
    U(1) = zero
END IF

DO I = 2,N - 1
    SIG = (X(I)-X(I-1))/ (X(I+1)-X(I-1))
    P = SIG*Y2(I-1) + two
    Y2(I) = (SIG-one)/P
    U(I) = ( six * ((Y(I+1)-Y(I))/ (X(I+1)-X(I))- (Y(I)-Y(I-1))/ &
            (X(I)-X(I-1)))/ (X(I+1)-X(I-1))-SIG*U(I-1))/P
enddo
IF (present(YPN)) THEN
    QN = half
    UN = (three/(X(N)-X(N-1)))* (YPN - (Y(N)-Y(N-1))/ (X(N)-X(N-1)))
ELSE
    QN = zero
    UN = zero
END IF

Y2(N) = (UN-QN*U(N-1))/ (QN*Y2(N-1) + one)
DO K = N-1 , 1 , -1
    Y2(K) = Y2(K)*Y2(K+1) + U(K)
enddo

if (present(stat)) stat = flag

end subroutine generate_spline_sp

!..............................................
SUBROUTINE generate_spline_dp(X,Y,N,Y2,YP1,YPN,stat)
!
! Given arrays x and y of size n, this routine computes
! an array y2 holding information about the second derivative
! of the interpolating spline. YP1 and YPN, if present, are the
! requested values of the first derivative of the spline at the
! end points. If not given, a "natural" spline with zero second
! derivative at the extremes is constructed.
! The array y2 is needed to evaluate the spline, but it only needs
! to be computed once.
!
integer, intent(in)                          :: n
real(kind=dp), dimension(:), intent(in)  :: x, y
real(kind=dp), dimension(:), intent(out) :: y2
real(kind=dp), intent(in), optional      :: yp1, ypn
integer, intent(out), optional               :: stat

real(kind=dp)  ::  P,QN,SIG,UN
INTEGER            ::  I,K
integer            ::  flag

real(kind=dp), parameter  :: zero = 0.0_dp, &
                                 half = 0.5_dp, &
                                 one  = 1.0_dp, &
                                 two  = 2.0_dp, &
                                 three= 3.0_dp, &
                                 six  = 6.0_dp

real(kind=dp), dimension(n)  ::   U    ! Automatic array

!--------------------------------------------------------
! Check that the nodes are different
!
flag = 0
do i =2, n
  if (x(i) == x(i-1)) then
   flag	= -1
   if (present(stat)) stat = flag
   RETURN
  endif
enddo
	
IF (present(YP1)) THEN
    Y2(1) = -half
    U(1) = (three/ (X(2)-X(1)))* ((Y(2)-Y(1))/ (X(2)-X(1))-YP1)
ELSE
    Y2(1) = zero
    U(1) = zero
END IF

DO I = 2,N - 1
    SIG = (X(I)-X(I-1))/ (X(I+1)-X(I-1))
    P = SIG*Y2(I-1) + two
    Y2(I) = (SIG-one)/P
    U(I) = ( six * ((Y(I+1)-Y(I))/ (X(I+1)-X(I))- (Y(I)-Y(I-1))/ &
            (X(I)-X(I-1)))/ (X(I+1)-X(I-1))-SIG*U(I-1))/P
enddo
IF (present(YPN)) THEN
    QN = half
    UN = (three/(X(N)-X(N-1)))* (YPN - (Y(N)-Y(N-1))/ (X(N)-X(N-1)))
ELSE
    QN = zero
    UN = zero
END IF

Y2(N) = (UN-QN*U(N-1))/ (QN*Y2(N-1) + one)
DO K = N-1 , 1 , -1
    Y2(K) = Y2(K)*Y2(K+1) + U(K)
enddo

if (present(stat)) stat = flag

end subroutine generate_spline_dp

!..............................................
!------------------------------------------------
SUBROUTINE evaluate_spline_sp(XA,YA,Y2A,N,X,Y)
!
! Given arrays xa and ya of size n, and an array y2a computed
! previously by routine spline, this routine computes the value
! at x of the interpolating spline. The value is returned in y.
!
integer, intent(in)                          :: n
real(kind=sp), dimension(:), intent(in)  :: xa, ya, y2a
real(kind=sp), intent(in)                :: x
real(kind=sp), intent(out)               :: y


real(kind=sp)    ::  A,B,H
integer              ::  K,KHI,KLO

KLO = 1
KHI = N
do
   IF (KHI-KLO <= 1) exit
   K = (KHI+KLO)/2
   IF (XA(K) > X) THEN
      KHI = K
   ELSE
      KLO = K
   END IF
enddo

H = XA(KHI) - XA(KLO)

A = (XA(KHI)-X)/H
B = (X-XA(KLO))/H
Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO)+  &
   (B**3-B)*Y2A(KHI))* (H**2)/6.0_sp

end subroutine evaluate_spline_sp
!..............................................
SUBROUTINE evaluate_spline_dp(XA,YA,Y2A,N,X,Y)
!
! Given arrays xa and ya of size n, and an array y2a computed
! previously by routine spline, this routine computes the value
! at x of the interpolating spline. The value is returned in y.
!
integer, intent(in)                          :: n
real(kind=dp), dimension(:), intent(in)  :: xa, ya, y2a
real(kind=dp), intent(in)                :: x
real(kind=dp), intent(out)               :: y


real(kind=dp)    ::  A,B,H
integer              ::  K,KHI,KLO

KLO = 1
KHI = N
do
   IF (KHI-KLO <= 1) exit
   K = (KHI+KLO)/2
   IF (XA(K) > X) THEN
      KHI = K
   ELSE
      KLO = K
   END IF
enddo

H = XA(KHI) - XA(KLO)

A = (XA(KHI)-X)/H
B = (X-XA(KLO))/H
Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO)+  &
   (B**3-B)*Y2A(KHI))* (H**2)/6.0_dp

end subroutine evaluate_spline_dp
!..............................................
!------------------------------------------------
end module flib_spline
