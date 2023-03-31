! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
! This module is not part of FFTW3, but uses the
! definitions in fftw3.f03, which is part of the FFTW3 distribution
!
module fftw3_mymod
use, intrinsic :: iso_c_binding
public
!integer, parameter, private :: i8b = selected_int_kind(18)  !! 2^63 = 9e18

type(C_PTR) :: plan
include "fftw3.f03"

end module fftw3_mymod
