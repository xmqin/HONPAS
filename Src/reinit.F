! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine reinit(sname) 

C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Format) package 
C     of J.M.Soler and A.Garcia
C
C Taken from redata. Writen by P.Ordejon, December'96
C **************************** OUTPUT *********************************
C character    slabel      : System Label (to name output files)
C character(len=*) sname       : System Name
C **********************************************************************

C
C  Modules
C
      use parallel,    only : Node
      use fdf
      use files,       only : slabel
      use files,       only : stdin_file, stdout_file
      use siesta_master, only: input_file  ! fdf data file
      use sys, only : bye

      implicit none

      character(len=*), intent(out) :: sname

C  Internal variables .................................................
      character(len=50) :: fileout

      integer :: narg, count
      integer :: in, length, lun, lun_tmp, iostat
      character(len=256) :: line, val

      logical debug_input, file_exists

C     Print Welcome and Presentation .......................................

      ! Default input-file
      stdin_file = "fdf_input"  
      
      if (Node.eq.0) then
      write(6,'(/a)') 
     .   '                           ***********************       '
      write(6,'(a)') 
     .   '                           *  WELCOME TO SIESTA  *       '
      write(6,'(a)')
     .   '                           ***********************       '

! ..................
!
!     Set name of file to read from.
#ifndef NO_F2003
      ! The file-input name will default to be 
      ! the last command line argument
      narg = command_argument_count()
      if ( narg > 0 ) then
         stdin_file = ' '
         call get_command_argument(narg,stdin_file,length)
         inquire(file=stdin_file,exist=file_exists)
         
         ! Force reading from stdin
         if ( .not. file_exists ) narg = 0
         if ( stdin_file == stdout_file ) narg = 0
         
      end if
#endif

c     
c     Choose proper file for fdf processing
c     (INPUT_DEBUG if it exists or "standard input",
c     processed and dumped to a temporary file)
c
      inquire(file='INPUT_DEBUG',exist=debug_input)
      if ( debug_input ) then
         write(6,'(a)') 'WARNING: ' //
     $        'Siesta is reading its input from file INPUT_DEBUG'
         stdin_file = 'INPUT_DEBUG'
         
#ifndef NO_F2003
      else if ( narg > 0 ) then
         
         ! Get file-name from input line
         stdin_file = ' '
         call get_command_argument(narg,stdin_file,length)
         if ( length > len(stdin_file) ) then
            call die ('The argument is too long to be retrieved, &
     &please limit your-self to 50 characters for the input file') 
         end if
         inquire(file=stdin_file, exist=file_exists)
         if ( .not. file_exists ) then
            call die('Input file '//trim(stdin_file)//' does not exist? &
     &Have you specified the wrong file-name?') 
         end if
           
         write(*,'(/,2a)') 'reinit: Reading from ',trim(stdin_file)

#endif
      else
         ! This is reading from std-in (or SIESTA-master)
         
         if ( input_file/=' ' ) then
!     Read data from given file
!     This is the case for running siesta under SIESTA-master
            write(6,'(a)') 'reinit: Reading from '//trim(input_file)
            call io_assign(lun)
            open(lun,file=input_file,form='formatted',status='old')
            rewind(lun)
         else                   ! Read from standard input
            write(6,'(/a)') 'reinit: Reading from standard input'
            lun = 5
         end if                 ! (input_file/=' ')
         
!     make sure we get a new file
         call io_assign(lun_tmp)
         file_exists = .true.
         do while ( file_exists )
            call system_clock( count )
            write(stdin_file,"(a,i5.5)") 'INPUT_TMP.', mod(count,100000)
            inquire(file=stdin_file, exist=file_exists )
         end do


         ! Open the input-file to dump the input data
         open(lun_tmp,file=stdin_file,
     &        form='formatted',status='replace')
         rewind(lun_tmp)
         write(6,'(a)') 'reinit: Dumped input in '//trim(stdin_file)
         write(6,"(a,23('*'),a,28('*'))")
     &        '***', ' Dump of input data file '
c
         do
            read(lun,iostat=iostat,fmt='(a)') line
            if ( iostat /= 0 ) exit
            length = len_trim(line)
            ! Truncate empty lines
            if ( length /= 0 ) then
               write(6,'(a)') line(1:length)
               if (.not. debug_input) then
                  write(lun_tmp,'(a)') line(1:length)
               end if
            end if
         end do

         write(6,"(a,23('*'),a,29('*'))")
     &        '***', ' End of input data file '
         call io_close(lun_tmp)

         if ( lun /= 5 ) then
            ! Be sure to close the input master file
            ! Otherwise we may end up opening too many files
            call io_close(lun)
         end if

!     "stdin_file" for fdf is now the temporary file. 
!     This was necessary historically to allow
!     the rewinds involved in fdf operation.
!     
      end if ! Decide how the input-file should be read
      end if ! Node .eq. 0

      
C Set up fdf ...
!
!     Choose a 'unique' prefix for the log (and possible debug) fdf files
!     The 5-digit sequence might be slightly different in different
!     processors, depending on the system time.

      call system_clock( count )
      write(fileout,"(a,i5.5,a)") 'fdf-', mod(count,100000), ".log"

      call fdf_init(stdin_file, trim(fileout))

#ifndef NO_F2003
!     Read special variables from the command line
!     This is done by ALL nodes!
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
            if ( line(1:3) == 'fdf' ) then
               if ( in >= narg ) 
     &              call die('Missing argument on command line, -fdf')
               in = in + 1
               call get_command_argument(in,line,length)
               
               ! We allow these variations:
               !  FDFLabel=0.1:eV
               !  FDFLabel:0.1:eV
               !  FDFLabel=0.1=eV
               !  "FDFLabel 0.1 eV"
               line = cmd_tokenize(line)
               call fdf_overwrite(line)
                         
            else if ( line(1:1) == 'L' ) then
               if ( in >= narg ) 
     &              call die('Missing argument on command line, -L')
               in = in + 1
               call get_command_argument(in,line,length)
               line = cmd_tokenize(line)
               line = 'SystemLabel '//trim(line)
               call fdf_overwrite(line)

            else if ( line(1:1) == 'V' ) then
               if ( in >= narg )
     &              call die('Missing argument on command line, -V')
               in = in + 1
               call get_command_argument(in,line,length)
               line = cmd_tokenize(line)
               line = 'TS.Voltage '//trim(line)
               call fdf_overwrite(line)
               
            else if ( line(1:9) == 'electrode' .or.
     &              line(1:4) == 'elec' ) then
               line = 'TS.HS.Save true'
               call fdf_overwrite(line)
               line = 'TS.DE.Save true'
               call fdf_overwrite(line)

            else if ( line(1:4) == 'help' .or.
     &              line(1:1) == 'h' ) then
      write(0,'(a)')'Help for calling SIESTA'
      write(0,'(a)')''
      write(0,'(a)')'  -out|-o <file>'
      write(0,'(a)')'      Instead of writing to stdout, write to file.'
      write(0,'(a)')'  -L <name>'
      write(0,'(a)')'      Short-hand for setting SystemLabel.'
      write(0,'(a)')'  -fdf <label>=<value>[:<unit>]'
      write(0,'(a)')'      Set the label to the corresponding value.'
      write(0,'(a)')'  -V <value>:<unit>'
      write(0,'(a)')'      Short-hand for setting TS.Voltage.'
      write(0,'(a)')'  -electrode|-elec'
      write(0,'(a)')'      Force TS.HS.Save and TS.DE.Save to true.'
      write(0,'(a)')'  <fdf-file>'
      write(0,'(a)')
     &     '      Use file as fdf-input, you need not to pipe it in.' 
      call bye('Help-menu requested, stopping')
            end if
            
         end do
      end if
#endif


C Define Name of the system ...
      sname = fdf_string('SystemName', ' ')
      if (Node.eq.0) then
        write(6,'(/a,71("-"))') 'reinit: '
        write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
        write(6,'(a,71("-"))') 'reinit: '
      endif
C ...

C Define System Label (short name to label files) ...
      slabel = fdf_string('SystemLabel', 'siesta')
      slabel = trim(slabel)
      count = len_trim(slabel)
      if ( count == 0 ) then
         call die('SystemLabel must be at least 1 character!')
      end if
      ! Check that there are no spaces in the SystemLabel
      length = index(slabel, ' ')
      if ( length > 0 .and. length < count ) then
         call die('SystemLabel must *NOT* contain any spaces!')
      end if
      if (Node.eq.0) then
        write(6,'(a,a)') 'reinit: System Label: ',trim(slabel)
        write(6,'(a,71("-"))') 'reinit: '
      endif
C ...

      contains
      
      function cmd_tokenize(line) result(tline)
      character(len=*), intent(in) :: line
      character(len=len(line)) :: tline
      
      integer :: i, n
      n = len(tline)
      tline = line
      do i = 1 , n
         if ( tline(i:i) == ':' .or. 
     &        tline(i:i) == '=' ) then
            tline(i:i) = ' '
         end if
      end do
      end function cmd_tokenize

      end
