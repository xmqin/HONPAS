!
! This module is not part of FFTW3, but uses the
! definitions in fftw3.f, which is part of the FFTW3 distribution
!
module fftw3_mymod

public
integer, parameter, private :: i8b = selected_int_kind(18)  !! 2^63 = 9e18

type, public :: fftw3_plan_t
  integer(kind=i8b) :: plan
end type

include "fftw3.f"

end module fftw3_mymod
