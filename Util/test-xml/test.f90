!Test the xml output of siesta 
!A xml file (called modified!) is compared against a  xml reference file.

program test_xml
use flib_dom
use compare_m, only: compare_debug, compare
use compare_m, only: STOP_ON_ERROR, MAX_NUMBER_OF_ERRORS

use compare_tol_m, only: TOLERANCES_FILE, TOL, dp

use m_getopts

implicit none

character(len=200) :: opt_arg
character(len=10)  :: opt_name 
integer            :: nargs, iostat, n_opts, nfiles, logtol
logical            :: fileExists

type(fnode), pointer     :: reference,modified 
character(len=1000) reference_file,modified_file

!     Process options
!

      n_opts = 0
      do
         call getopts('dse:g:t:x:',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
         case ('d', '+d')
            ! Debug
            compare_debug = .true.
         case ('s', '+s')
            ! Stop on any error
            STOP_ON_ERROR = .true.
         case ('e', '+e')
            ! Specify maximum number of errors before stopping
            ! 0 means no stop
            read(opt_arg,"(i10)") MAX_NUMBER_OF_ERRORS
         case ('g', '+g')
            ! Specify global tolerance (negative log10)
            ! eg: -g 4 means 1.0e-4 (the default)
            read(opt_arg,"(i2)") logtol
            TOL = 10.0_dp**(-logtol)
         case ('t', '+t')
            ! Secify tolerance file
            TOLERANCES_FILE = opt_arg
         case ('x', '+x')
            ! Future
            ! Disregard certain labels
         case ('?',':')
            write(0,*) "Invalid option: ", opt_arg(1:1)
            STOP
         end select
      enddo

      nargs = command_argument_count()
      nfiles = nargs - n_opts + 1
      if (nfiles /= 2) &
          STOP "Usage: test-xml [-d ] [ -s ] [-e NUM-ERRORS] [-g LOGTOL] [ -t TOL-FILE] file1 file2"

      call get_command_argument(n_opts,value=reference_file,status=iostat)
      if (iostat /= 0) then
          STOP "Cannot get first file"
      endif
      call get_command_argument(n_opts+1,value=modified_file,status=iostat)
      if (iostat /= 0) then
         STOP "Cannot get second file"
      endif

if (compare_debug) then
   write(0,"(a,a)") "Will try to use tolerances file: ", trim(TOLERANCES_FILE)
   inquire(file=TOLERANCES_FILE, exist=fileExists)
   if(fileExists) then
      write(6,*) " ** Tolerances file found!"
   else
      write(0,*) " ** Tolerances file not found. Will use global tol"
      write(0,"(a,g12.6)") "Using global tolerance: ", tol
   endif
  
endif
   
!Parse the files
reference => parsefile(trim(reference_file),verbose=.false.)
modified  => parsefile(trim(modified_file))

!!$print *, "Dump of Reference... before normalize"
!!$call dumpTree(reference)
!!$print *, "Dump of Modified.. before normalize"
!!$call dumpTree(modified)

call normalize(reference)
call normalize(modified)

!!$print *, "Dump of Reference... after normalize"
!!$call dumpTree(reference)
!!$print *, "Dump of Modified.. after normalize"
!!$call dumpTree(modified)

!
call compare(reference,modified)

end program 


