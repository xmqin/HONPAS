module compare_tol_m
use m_strings

implicit none


integer, parameter ::  sp = selected_real_kind(6,30)
integer, parameter ::  dp = selected_real_kind(14,100)


public :: tol,sp,dp
public :: findTol

!Modify this value to change the default tolerance.
real(dp) :: TOL = 1.0E-4
character(len=256), public :: TOLERANCES_FILE = "tolerances.dat"

private


!List of labels/tolerances
integer,save :: n_labels
real(dp), allocatable              :: tolerances(:)
character(len=50), allocatable :: labels(:)

contains 


function findTol(label)
  
  character(len=*), intent(in) :: label
  real(dp) :: findTol

  !Internal vars.
  integer           :: i
  logical,save      :: firstTime = .true.
  character(len=50) :: lab

  if (firstTime) then
     call find_number_of_tolerance_labels()
     call find_Tolerances()
     firstTime = .false.
  endif
  !Default value.
  findTol = tol

  lab = label
  call to_lowercase(lab)
  !print *,"n_labels:",n_labels
  do i=1,n_labels
     !print *,"lab,tol: ",labels(i),tolerances(i)
     if (labels(i) == lab ) then
        findTol = tolerances(i)     
        exit
     endif
  enddo
  
end function findTol

!----------------------------------------------------

subroutine find_number_of_tolerance_labels
  logical :: fileExists
  integer :: i,io
  character(len=50) :: tmp

  inquire(file=TOLERANCES_FILE, exist=fileExists)
  if( .not. fileExists) then
     n_labels = -1
  else
     open(unit=30,file=TOLERANCES_FILE)
     do i=1,10000000
        read(30,*,iostat=io) tmp 
        if (io < 0) then 
           n_labels = i-1
           exit
        endif
     enddo
  endif
  close(30)
  

end subroutine find_number_of_tolerance_labels

!----------------------------------------------------
subroutine find_tolerances()

  logical :: fileExists
  integer :: i
  
  inquire(file=TOLERANCES_FILE, exist=fileExists)
  if( fileExists) then
     open(unit=30,file=TOLERANCES_FILE)
     
     allocate(tolerances(1:n_labels),labels(1:n_labels))
     do i=1,n_labels
        read(30,*) labels(i),tolerances(i)
        !print *, labels(i),tolerances(i)
     enddo
     do i=1,n_labels
        call to_lowercase(labels(i))
        !print *, labels(i),tolerances(i)
     enddo

     close(30)
  endif
end subroutine find_tolerances

end module compare_tol_m
