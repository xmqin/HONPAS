! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine reinit(sname) 

C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Fromat) package 
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
      use precision
      use parallel,    only : Node
      use fdf
      use files,       only : slabel, label_length
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      character(len=*) :: sname

C  Internal variables .................................................
      character
     .  filein*20, fileout*20, line*150 

      character
     .  slabel_default*59, sname_default*20

      integer 
     .  count, length, lun, lun_tmp

      logical debug_input, file_exists
#ifdef MPI
      integer
     .  MPIerror
#endif

C Print Welcome and Presentation .......................................
      if (Node.eq.0) then
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO SIESTA  *       '
        write(6,'(a)')
     .   '                           ***********************       '
C ..................

C Dump data file to output file 
C and generate scratch file for FDF to read from
C (except if INPUT_DEBUG exists)
C

        inquire(file='INPUT_DEBUG',exist=debug_input)
        if (debug_input) then
           write(6,'(a)') 'WARNING: ' //
     $          'Siesta is reading its input from file INPUT_DEBUG'
           
           call io_assign(lun)
           filein = 'INPUT_DEBUG'
           open(lun,file='INPUT_DEBUG',form='formatted',status='old')
           rewind(lun)
        else
c
c          Read from standard input
c
           write(6,'(/a)') 'reinit: Reading from standard input'
           lun = 5
           call io_assign(lun_tmp)
           do
             call system_clock( count )
             write(filein,"(a,i5.5)") 'INPUT_TMP.', mod(count,100000)
             inquire( file=filein, exist=file_exists )
             if (.not.file_exists) exit
           end do
           open(lun_tmp,file=filein,
     $                  form='formatted',status='replace')
           rewind(lun_tmp)
        endif

        write(6,'(a,23(1h*),a,28(1h*))')
     .    '***', ' Dump of input data file '
c
 10     continue
        read(lun,err=20,end=20,fmt='(a)') line
        call chrlen(line,0,length)
        if (length .ne. 0) then
           write(6,'(a)') line(1:length)
           if (.not. debug_input) write(lun_tmp,'(a)') line(1:length)
        endif
        goto 10
 20     continue
c
        write(6,'(a,23(1h*),a,29(1h*))')
     .    '***', ' End of input data file '
c
c       Choose proper file for fdf processing
c
        if (debug_input) then
           call io_close(lun)
        else
           call io_close(lun_tmp)
        endif

C Set up fdf ...
        fileout = 'fdf.log'
        call fdf_init(filein,fileout)

C Defile Name of the system ...
        sname_default = ' '
        sname = fdf_string('SystemName',sname_default)
        write(6,'(/a,71(1h-))') 'reinit: '
        write(6,'(a,a)') 
     .   'reinit: System Name: ',trim(sname)
        write(6,'(a,71(1h-))') 'reinit: '
C ...

C Defile System Label (short name to label files) ...
        slabel_default  = 'siesta'
        slabel = fdf_string('SystemLabel',slabel_default)
        write(6,'(a,a)') 
     .   'reinit: System Label: ',slabel
        write(6,'(a,71(1h-))') 'reinit: '
C ...
      endif

C Global broadcast of strings
#ifdef MPI
      call MPI_Bcast(sname,len(sname),MPI_character,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(slabel,label_length,MPI_character,0,MPI_Comm_World,
     .  MPIerror)
#endif

      return
      end

