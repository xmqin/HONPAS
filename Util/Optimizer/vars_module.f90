! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module vars_module
use precision
use parse
use sys, only: die

implicit none

type, public :: var_t
   real(dp) :: x = 0.5_dp
   real(dp) :: min = 0.0_dp
   real(dp) :: max = 1.0_dp
   real(dp) :: range = 1.0_dp
   logical  :: set_by_user = .false.
   character(len=40) :: name = "dummy"
end type var_t

integer, public  :: nvars = 0
type(var_t), dimension(:), allocatable, public :: var

public :: generate_subs_file, read_vars, constrained

CONTAINS

   subroutine generate_subs_file(x,name)
     real(dp), dimension(:), intent(in)  :: x
     character(len=*), intent(in) :: name
     
     integer :: i, iu

     call io_assign(iu)
     open(unit=iu,file=trim(name)// ".sed",form="formatted",status="replace")
     do i = 1, size(x)
        ! The f0 construct is fortran 95 only. 
        ! If your computer does not support it, you can
        ! change it to f20.7, for example, BUT then you cannot
        ! use constructs such as -$rc_s  (negative value of a variable)
        ! in the TEMPLATE file.
        write(iu,"(a,a,a,f0.7,a)") "s/$", trim(var(i)%name), "/", x(i), "/g"
     enddo
     close(iu)
   end subroutine generate_subs_file
   
   subroutine read_vars(file)
     character(len=*), intent(in) :: file

     integer iu, iostat
     character(len=132) :: line
     type(parsed_line), pointer  :: p

     logical :: filling
     
     call io_assign(iu)
     open(unit=iu,file="VARS",form="formatted",status="old", &
          iostat=iostat)
     if (iostat /=0 ) then
        STOP "VARS file not found"
     endif

     filling = .false.
     
     nvars = 0
     do
        read(iu,iostat=iostat,fmt="(a132)") line
        if (iostat /=0) exit
        p => digest(line)
        if ( match(p,"nvv")) then
           nvars = nvars + 1
           if (filling) then
              call fill_var(var(nvars),p)
           endif
        else
           call die("Wrong format in VARS file")
        endif
     enddo

     allocate(var(nvars))
     rewind(iu)
     filling = .true.

     nvars = 0
     do
        read(iu,iostat=iostat,fmt="(a132)") line
        if (iostat /=0) exit
        p => digest(line)
        if ( match(p,"nvv")) then
           nvars = nvars + 1
           if (filling) then
              call fill_var(var(nvars),p)
           endif
        else
           call die("Wrong format in VARS file")
        endif
     enddo
     close(iu)
     call print_vars()

   end subroutine read_vars

   subroutine fill_var(v,p)
     type(var_t), intent(inout) :: v
     type(parsed_line), pointer :: p

     real(dp) :: r

     v%name=trim(names(p,1))
     v%min=values(p,1)
     v%max=values(p,2)
     v%range=values(p,2) - values(p,1)
     if (v%range <= 0.0_dp) then
        call die("Negative or null interval")
     endif
     !     Set to random place in interval if not specified
     if (match(p,"nvvv")) then
        v%x=values(p,3)
        v%set_by_user = .true.
     else
        call random_number(r)
        v%x=v%min + &
             r*(v%max - v%min)
        v%set_by_user = .false.
     endif
   end subroutine fill_var

   subroutine print_vars()
     integer i
     do i = 1, nvars
        print "(a,3f10.6)", trim(var(i)%name),  &
               var(i)%min, var(i)%max, var(i)%x
     enddo
   end subroutine print_vars

   function constrained(x,i) result(y)
     real(dp), intent(in) :: x
     integer, intent(in)  :: i
     real(dp)             :: y
     
     ! Simple control of range: Stick to the nearest wall

     y = x
     if (x < var(i)%min) y = var(i)%min
     if (x > var(i)%max) y = var(i)%max

   end function constrained
     
    
end module vars_module

