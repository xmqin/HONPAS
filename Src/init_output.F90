! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine init_output(IO)

  ! Subroutine to initialize the output file for SIESTA
  ! This should *only* be used to close(6) and open(6) for
  ! the new output file
  !
  ! Taken from reinit, and adapted by Nick R. Papior, 2017
  
  use files, only : stdout_file
  
  implicit none

  ! Whether this node is allowed to perform IO
  logical, intent(in) :: IO
  
  ! Quick return for non-IO node
  ! Perhaps this should be abber
  if ( .not. IO ) return

  ! First we determine whether all output should be
  ! written to stdout or to a file.
  stdout_file = cla_get_output()
  if ( stdout_file /= ' ' ) then
     ! Close default output to create a new handle
     close(unit=6)
     open(unit=6, file=trim(stdout_file), form="formatted", &
          position="rewind", action="write", status="unknown")
  end if
  
contains
  
  function cla_get_output() result(outfile)
    use files, only: label_length
    character(len=label_length) :: outfile

    ! Local variables
    integer :: narg, in, length
    character(len=128) :: line

    ! Default to no file
    outfile = ' '
    
#ifndef NO_F2003
    ! Read the output file from the command line (--out|-out|-o)
    narg = command_argument_count()
    if ( narg > 0 ) then
       in = 0
       do while ( in <= narg - 1 )

          in = in + 1
          call get_command_argument(in,line,length)

          ! If it is not an option, skip it
          if ( line(1:1) /= '-' ) cycle

          ! Truncate '-' to no '-'
          do while ( line(1:1) == '-' )
             line = line(2:)
          end do

          ! We allow these line
          if ( line(1:1) == 'o'.or.line(1:3) == 'out' ) then
             if ( in >= narg ) &
                  call die('Missing argument for output file')
             in = in + 1
             call get_command_argument(in,line,length)
             ! Capture whether the output file is too long
             if ( length > label_length ) then
                call die('Output file name is too long!')
             end if
             ! Copy over the output file
             outfile = line(1:length)

          end if
          
       end do
    end if
#endif
  end function cla_get_output
      
end subroutine init_output
