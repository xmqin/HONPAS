! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module minimizer
use vars_module 

implicit none

integer, save :: job_no = 0

public :: objective_function

CONTAINS

  function objective_function(x) result(value)
    use vars_module
  real(dp), dimension(:), intent(in) :: x
  real(dp)                           :: value

  character(len=40) :: name

!$omp critical 
! This section will be executed in turns by the threads
  job_no = job_no + 1
  write(name,"(a,i4.4)") "job_", job_no
  call generate_subs_file(x,name)
!$omp end critical

  call dispatch_job(name, value)
  
end function objective_function

subroutine dispatch_job(name,value)
  character(len=*), intent(in) :: name
  real(dp), intent(out)  :: value

  real(dp) :: e
  integer  :: iu, iostat

  value = huge(0.0_dp)  ! Default value if something goes wrong

  ! Create directory
  ! apply sed file to template and copy it to directory
  call system("sh run_script.sh " // trim(name))
  call io_assign(iu)
  open(unit=iu,file=trim(name)// "/OPTIM_OUTPUT", &
        form="formatted",status="old", iostat=iostat)
  if (iostat /=0 ) then
     STOP "OPTIM_OUTPUT file not available"
  endif
  read(iu,*) e
  value = e
  close(iu)
end subroutine dispatch_job

end module minimizer

