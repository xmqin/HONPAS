! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! 
! Transferred to be used in the TBTrans utility.
!
subroutine tbt_reinit( sname , slabel ) 

! Subroutine to initialise the reading of the data for SIESTA 
!
!     It uses the FDF (Flexible Data Format) package 
!     of J.M.Soler and A.Garcia
!
! Taken from redata. Written by Nick P. Andersen 2015
! **************************** OUTPUT *********************************
! character    slabel      : System Label (to name output files)
! character(len=*) sname       : System Name
! **********************************************************************

  use sys, only : bye
  use parallel, only : Node
  use fdf
  use m_verbosity

  implicit none

  character(len=*), intent(out) :: sname, slabel

!  Internal variables .................................................
  character(len=50) :: filein, fileout, string

  integer ::  count, in, length, lun, lun_tmp, iostat
  character(len=256) :: line

  logical :: debug_input, file_exists, is_pipe

! Print Welcome and Presentation .......................................
!     Non-master mpi-processes receive a copy of all the
!     pre-processed fdf input information (recursively
!     including files when necessary), and dump it on a
!     text file with name "fdf_input.<ProcessNumber>".
!     They then read from this file to construct a memory
!     image of the information.
!
!     The master process creates its memory image directly
!     from the standard fdf file (renamed as INPUT_TMP.$$,
!     as usually done for historical reasons (see below)).
!
  filein = "fdf_input" 

  if (Node.eq.0) then
     write(6,'(/a)') &
          '                           ************************ '
#ifdef TBT_PHONON
     write(6,'(a)') &
          '                           *  WELCOME TO PHtrans  * '
#else
     write(6,'(a)') &
          '                           *  WELCOME TO TBtrans  * '
#endif
     write(6,'(a)') &
          '                           ************************ '

!
!       Set name of file to read from. Done only
!       in the master node.
!
     is_pipe = .true.
#ifndef NO_F2003
     count = command_argument_count()
     if ( count > 0 ) then
        filein = ' '
        call get_command_argument(count, filein, length)
        inquire(file=filein,exist=debug_input)
        if ( debug_input ) then
          is_pipe = .false.
        else
          count = 0
        end if
     end if
#endif

!
!     Choose proper file for fdf processing
!     (INPUT_DEBUG if it exists or "standard input",
!      processed and dumped to a temporary file)
!
     inquire(file='INPUT_DEBUG',exist=debug_input)
     if ( debug_input ) then
        write(*,'(a)') 'WARNING: ' // &
             'TBTrans is reading its input from file INPUT_DEBUG'
        filein = 'INPUT_DEBUG'

#ifndef NO_F2003
     else if ( .not. is_pipe ) then

        ! Get file-name from input line (last input argument)
        filein = ' '
        call get_command_argument(count, filein, length)
        if ( length > len(filein) ) then
           call die('The argument is too long to be retrieved, please &
                &limit to 50 characters for the input file')
        end if
        inquire(file=filein, exist=debug_input)
        if ( .not. debug_input ) then
           call die('Input file '//trim(filein)//' does not exist? Have &
                &you specified the wrong file-name?')
        end if

        write(*,'(/,2a)') 'reinit: Reading from ',trim(filein)

#endif
     else
!
!          Read from standard input (dumped to a temp file)
!
        write(*,'(/a)') 'reinit: Reading from standard input'
        lun = 5
        call io_assign(lun_tmp)
        do  ! make sure we get a new file
           call system_clock( count )
           write(string,*) count
           filein = 'INPUT_TMP.'//adjustl(string)
           inquire( file=filein, exist=file_exists )
           if (.not.file_exists) exit
        end do

        open(lun_tmp,file=filein, &
             form='formatted',status='replace')
        rewind(lun_tmp)
        write(*,"(a,23('*'),a,28('*'))") &
             '***', ' Dump of input data file '

        do
           read(lun,iostat=iostat,fmt='(a)') line
           if (iostat /= 0 ) exit
           length = len_trim(line)
           if (length /= 0) then
              write(*,'(a)') line(1:length)
              if (.not. debug_input) then
                 write(lun_tmp,'(a)') line(1:length)
              endif
           endif
        enddo
        write(*,"(a,23('*'),a,29('*'))") &
             '***', ' End of input data file '
        call io_close(lun_tmp)
!
!          "filein" for fdf is now the temporary file. 
!          This was necessary historically to allow
!          the rewinds involved in fdf operation.
!
     endif
  endif

! Set up fdf ...
!
! Choose a 'unique' prefix for the log (and possible debug) fdf files
! The 5-digit sequence might be slightly different in different
! processors, depending on the system time.
  call system_clock( count )
  write(fileout,"(a,i5.5,a)") 'fdf-', mod(count,100000), ".log"

  call fdf_init(filein,trim(fileout))

#ifndef NO_F2003
  ! Read command line flags from the command line
  count = command_argument_count()
  
  ! When we are not using pipes we have to skip the last argument
  if ( .not. is_pipe ) count = count - 1
  
  if ( count > 0 ) then
    in = 0
    do while ( in <= count )
      in = in + 1
      call get_command_argument(in,line,length)
      
      ! If it is not an option, skip it
      if ( line(1:1) /= '-' ) cycle
      
      do while ( line(1:1) == '-' )
        line = line(2:)
      end do
      
      ! We allow these line
      if ( line(1:3) == 'fdf' ) then
        if ( in >= count ) &
            call die('Missing argument on command line, -fdf')
        in = in + 1
        call get_command_argument(in,line,length)
        
        ! We allow these variations:
        !  TBT.Voltage=0.1:eV
        !  TBT.Voltage:0.1:eV
        !  TBT.Voltage=0.1=eV
        line = cmd_tokenize(line)
        call fdf_overwrite(line)
        
      else if ( line(1:1) == 'V' ) then
        if ( in >= count ) &
            call die('Missing argument on command line, -V')
        in = in + 1
        call get_command_argument(in,line,length)
        line = cmd_tokenize(line)
        line = 'TBT.Voltage '//trim(line)
        call fdf_overwrite(line)
        
      else if ( line(1:1) == 'D' ) then
        if ( in >= count ) &
            call die('Missing argument on command line, -D')
        in = in + 1
        call get_command_argument(in,line,length)
        line = 'TBT.Directory '//trim(line)
        call fdf_overwrite(line)
        
      else if ( line(1:2) == 'HS' ) then
        if ( in >= count ) &
            call die('Missing argument on command line, -HS')
        in = in + 1
        call get_command_argument(in,line,length)
        line = 'TBT.HS '//trim(line)
        call fdf_overwrite(line)
        
      else if ( line(1:1) == 'L' ) then
        if ( in >= count ) &
            call die('Missing argument on command line, -L')
        in = in + 1
        call get_command_argument(in,line,length)
        line = cmd_tokenize(line)
        line = 'SystemLabel '//trim(line)
        call fdf_overwrite(line)
        
      else if ( line(1:4) == 'help' .or. line(1:1) == 'h' ) then
        write(*,'(a)') 'Help for calling the tight-binding transport code'
        write(*,'(a)') '  -out <file>'
        write(*,'(a)') '      Write all output to <file> instead of STDOUT'
        write(*,'(a)') '  -fdf <label>=<value>[:<unit>]'
        write(*,'(a)') '      Set the label to the corresponding value.'
        write(*,'(a)') '  -V <value>:<unit>'
        write(*,'(a)') '      Short-hand for setting TBT.Voltage'
        write(*,'(a)') '  -D <directory>'
        write(*,'(a)') '      Short-hand for setting TBT.Directory'
        write(*,'(a)') '  -HS <Hamiltonian>'
        write(*,'(a)') '      Short-hand for setting TBT.HS'
        write(*,'(a)') '  -L <name>'
        write(*,'(a)') '      Short-hand for setting SystemLabel'
        write(*,'(a)') '  <fdf-file>'
        write(*,'(a)') '      Use file as fdf-input, you need not to pipe it in.'
        call bye('Help-menu requested, stopping')
      end if
      
    end do
  end if
#endif

  ! Initialize the verbosity setting
  call init_verbosity('TBT.Verbosity', 5)

! Define Name of the system ...
  sname = fdf_get('SystemName',' ')
  if (Node.eq.0) then
     write(*,'(/a,71("-"))') 'reinit: '
     write(*,'(a,a)') 'reinit: System Name: ',trim(sname)
     write(*,'(a,71("-"))') 'reinit: '
  endif

! Define System Label (short name to label files) ...
  slabel = fdf_get('SystemLabel','siesta')
  if (Node.eq.0) then
     write(*,'(a,a)') 'reinit: System Label: ',trim(slabel)
     write(*,'(a,71("-"))') 'reinit: '
  endif

contains

  function cmd_tokenize(line) result(tline)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: tline

    integer :: i, n
    n = len(tline)
    tline = line
    do i = 1 , n
       if ( tline(i:i) == ':' .or. &
            tline(i:i) == '=' ) then
          tline(i:i) = ' '
       end if
    end do
  end function cmd_tokenize

end subroutine tbt_reinit
