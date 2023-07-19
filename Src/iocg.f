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
      subroutine iocg( task, naux, cgaux, cgcntr, relaxd, found )

c**************************************************************************
c Reads/writes auxuliary arrays for conj. grag. continuation
c Written by E. Artacho. January 1999.
c**************** INPUT ***************************************************
c character task*3 : 'read' or 'write'
c integer   naux   : dimension of cgaux
c**************** INPUT OR OUTPUT (depending on task) *********************
c real*8    cgaux(naux)   : auxiliary array for CG
c real*8    cgcntr(0:20)  : same
c logical   relaxd        : whether system is relaxed or not.
c***************** OUTPUT *************************************************
c logical found : Has DM been found in disk? (Only when task='read')
c**************************************************************************

      use fdf
      use files,     only : slabel, label_length
      use precision, only : dp
      use sys,  only      : die

      implicit          none

      character(len=label_length+3) :: paste
      character         task*(*)
      logical           found, relaxd
      integer           naux
      real(dp)          cgaux(naux), cgcntr(0:20)

      external          chkdim, io_assign, io_close, paste


c Internal variables and arrays ------------------------------------------

      character(len=label_length+3) :: fname
      logical   exist1, frstme
      integer   nauxr, i, unit1

      save      frstme, fname
      data      frstme /.true./

c ------------------------------------------------------------------------

c find file name ---------------------------------------------------------

      if (frstme) then
        fname = paste(slabel,'.CG')
        frstme = .false.
      endif

c read it if it is there -------------------------------------------------

      if (task.eq.'read' .or. task.eq.'READ') then
        inquire (file=fname, exist=exist1)

        if (exist1) then
          write(6,'(/,a)') 'iocg: Reading CG continuation file'
          call io_assign(unit1)
          open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
          rewind(unit1)
          read(unit1) nauxr, relaxd
          call chkdim( 'iocg', 'cgaux', naux, nauxr, 1 )
          read(unit1) (cgcntr(i), i = 0, 20)
          read(unit1) (cgaux(i), i = 1, naux)
          call io_close(unit1)
          found = .true.
        else
          relaxd = .false.
          cgcntr(0) = 0
          cgcntr(1) = 1
          found = .false.
        endif

c write it ---------------------------------------------------------------

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        call io_assign(unit1)
        open( unit1, file=fname,
     .        form='unformatted', status='unknown' )
        rewind(unit1)
        write(unit1) naux, relaxd
        write(unit1) (cgcntr(i), i = 0, 20)
        write(unit1) (cgaux(i), i = 1, naux)
        call io_close(unit1)

      else
        call die('iocg: Incorrect task')
      endif

      return
      end

