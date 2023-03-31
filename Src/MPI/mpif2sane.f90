! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program mpif2sane

include "mpif.h"

character(len=32) :: str = "  integer, parameter, public :: "

write(unit=*,fmt="(a,/)") "module mpi__include"

write(unit=*,fmt="(2a,i5)") str, "mpi_comm_world =", mpi_comm_world
write(unit=*,fmt="(2a,i5)") str, "mpi_comm_self =", mpi_comm_self
write(unit=*,fmt="(2a,i5)") str, "mpi_sum =", mpi_sum
write(unit=*,fmt="(2a,i5)") str, "mpi_maxloc =", mpi_maxloc
write(unit=*,fmt="(2a,i5)") str, "mpi_max =", mpi_max
write(unit=*,fmt="(2a,i5)") str, "mpi_lor =", mpi_lor
write(unit=*,fmt="(2a,i5)") str, "mpi_status_size =", mpi_status_size

write(unit=*,fmt="(2a,i5)") str, "mpi_integer =", mpi_integer
write(unit=*,fmt="(2a,i5)") str, "mpi_character =", mpi_character
write(unit=*,fmt="(2a,i5)") str, "mpi_logical =", mpi_logical
write(unit=*,fmt="(2a,i5)") str, "mpi_real =", mpi_real
write(unit=*,fmt="(2a,i5)") str, "mpi_double_precision =", &
                   mpi_double_precision
write(unit=*,fmt="(2a,i5)") str, "mpi_complex =", mpi_complex
write(unit=*,fmt="(2a,i5)") str, "mpi_double_complex =", mpi_double_complex
write(unit=*,fmt="(2a,i5)") str, "mpi_2double_precision =", &
                   mpi_2double_precision

!write(unit=*,fmt="(2a,i5)") str, "mpi_2real =", mpi_2real
!write(unit=*,fmt="(2a,i5)") str, "mpi_2complex =", mpi_2complex
!write(unit=*,fmt="(2a,i5)") str, "mpi_2double_complex =", mpi_2double_complex


write(unit=*,fmt="(/,a)") "end module mpi__include"

end program mpif2sane
