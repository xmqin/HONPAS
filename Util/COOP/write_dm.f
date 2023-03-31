! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iodm( task, maxnd, nbasis, nspin, numd, 
     .                 listdptr, listd, dm, found )
C *******************************************************************
C Reads/writes density matrix from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*(*) : 'read' or 'write'
C integer   maxnd    : First dimension of listd and dm
C integer   nbasis   : Number of atomic orbitals
C integer   nspin    : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numd(nbasis)     : Control vector of DM matrix
C                            (number of nonzero elements of each row)
C integer listdptr(nbasis) : Control vector of DM matrix
C                            (pointer to the start of each row)
C integer listd(maxnd)     : Control vector of DM matrix
C                            (list of nonzero elements of each row)
C real*8  dm(maxnd,nspin)  : Density matrix
C ********* OUTPUT *************************************************
C logical found : Has DM been found in disk and does it have the right
C                 dimensions? (Only when task='read')
C ******************************************************************
C Behavior: If the dimensions of the DM on file are not compatible with
C what is expected, the routine returns 'found=.false.'
      
C
C  Modules
C
      use precision,    only : dp
      use sys,          only : die

      implicit  none

! -- avoid dependencies
      logical   :: ionode  = .true.
      integer   :: node    =  0
      integer   :: nodes   =  1
      integer, parameter   :: label_length    =  50
      character(len=50) :: slabel = "dummy"


      character(len=*), intent(in) :: task
      integer, intent(in) :: maxnd
      integer, intent(in) :: nbasis
      integer, intent(in) :: nspin

      integer, intent(inout) :: numd(nbasis)
      integer, intent(inout) :: listdptr(nbasis)
      integer, intent(inout) :: listd(maxnd)
      real(dp), intent(inout) :: dm(maxnd, nspin)

      logical, intent(out) :: found

C Internal variables
      logical   file_exists, okdim, dim_err
      integer   im, is, unit1, m, nb, ndmax, ns
      integer   nbasistot, ml, ndmaxg
      integer, dimension(:), allocatable, save :: numdg
      

C Saved internal variables:
      logical,           save :: frstme = .true., scndme = .false.
      character(len=label_length+3), save :: fnameu
      character(len=label_length+4), save :: fnamef
      character(len=label_length+4), save :: fnamei, fnameo
      logical,                       save :: fmti, fmto
      character(len=11),             save :: formin, formout

! Character formats for formatted I/O
! We assume that we have no integers greater than (1e10-1)
! (which is bigger than the largest 32-bit integer)
!
! and that floating point accuracy is sufficiently represented
! by 16 decimal digits of precision (which suffices for 64-bit
! IEEE floats)

      character(len=*), parameter :: intfmt = '(I11)'
      character(len=*), parameter :: floatfmt = '(ES22.14)'


! We might want to use formatted DM files in order to
! transfer them between computers. 

! However, whenever the settings are different for input/
! output, this is only the case for the first step.

! Thereafter, they must be the same; otherwise we will
! end up reading from the wrong file.

! Therefore, on all steps after the first, the 

C Find file name
      if (ionode) then
        if (frstme) then
          fmto = .false.
          fmti = fmto
          frstme = .false.
          scndme = .true.
        elseif (scndme) then
          fmti = fmto
          scndme = .false.
        endif
        fnameu = trim(slabel) // '.DM'
        fnamef = trim(slabel) // '.DMF'
        if (fmto) then
          formout = 'formatted'
          fnameo = fnamef
        else
          formout = 'unformatted'
          fnameo = fnameu
        endif
        if (fmti) then
          formin = 'formatted'
          fnamei = fnamef
        else
          formin = 'unformatted'
          fnamei = fnameu
        endif
      endif

C Find total number of basis functions over all Nodes
      nbasistot = nbasis

C Allocate local buffer array for globalised numd
      allocate(numdg(nbasistot))

      if (task.eq.'read' .or. task.eq.'READ') then
        okdim = .true.
        if (Node.eq.0) then
          inquire (file=fnamei,  exist=file_exists)
        endif


        if (file_exists) then

          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from file'
            call io_assign(unit1)
            open( unit1, file=fnamei, form=formin, status='old' )
            rewind(unit1)
            if (fmti) then
              read(unit1, intfmt) nb, ns
            else
              read(unit1) nb, ns
            endif
          endif

C Communicate the values to all Nodes and adjust to allow for
C distributed memory before checking the dimensions

C Check dimensions. It will no longer stop if wrong, but continue with fresh DM

          if ( nbasistot .ne. nb) then
             okdim = .false.
             if (Node .eq. 0)
     .       write(6,"('iodm: DM continuation file not used because ',
     .                 'of wrong dimensions.'/
     .                 'iodm:    nbasis =',i6,' nbasis on file =',
     .                 i6)") nbasistot, nb
             goto 1000
          endif

          if ( nspin .ne. ns) then
             okdim = .false.
             if (Node .eq. 0)
     .       write(6,"('iodm: DM continuation file not used because ',
     .                 'of wrong dimensions.'/
     .                 'iodm:    nspin =',i6,' nspin on file =',
     .                 i6)") nspin, ns
             goto 1000
          endif

          if (Node.eq.0) then
            if (fmti) then
              read(unit1, intfmt) (numdg(m),m=1,nbasistot)
            else
              read(unit1) (numdg(m),m=1,nbasistot)
            endif
          endif

C Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
            ml = m
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

C Check size of first dimension of dm. If we are in parallel both maxnd 
C and ndmax can vary from node to node: checking that the global maximum
C maxnd is smaller than the global minimum ndmax won't work as this may 
C miss some cases where the dimension is OK. Instead do a global .and.
C of the result of the check. However, we then need to take case on the 
C not OK path to report correctly.

          dim_err = ( maxnd .lt. ndmax )
          if ( dim_err ) then
             okdim = .false.
             write(6,"(a,a,a,i6,a,i6)")
     $            'iodm: DM continuation file not used because ',
     .            'of wrong dimensions.',
     .            'iodm:    maxnd =', maxnd,' must be at least =', ndmax
             goto 1000
          endif


          do m = 1,nbasistot
              ml = m
              if (fmti) then
                read(unit1, intfmt) 
     .               (listd(listdptr(ml)+im),im=1,numd(ml))
              else
                read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
              endif
          enddo


          do is = 1,nspin
            do m = 1,nbasistot
                ml = m
                if (fmti) then
                  read(unit1, floatfmt)
     .                 (dm(listdptr(ml)+im,is),im=1,numd(ml))
                else
                  read(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
                endif
            enddo
          enddo

 1000     if (Node.eq.0) then
            call io_close(unit1)
          endif

          found = .true.

        else

          found = .false.

        endif

C Transmit found only if the DM found had the right dimensions.
C Otherwise act as not found and start afresh
        found = found .and. okdim

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file="DMOUT", form=formout, status='unknown' )
          rewind(unit1)
          if (fmto) then
            write(unit1, intfmt) nbasistot, nspin
          else
            write(unit1) nbasistot, nspin
          endif
        endif

C Create globalised numd
        do m = 1,nbasistot
            ml = m
            numdg(m) = numd(ml)
        enddo

C Write out numd array
        if (Node.eq.0) then
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo
          if (fmto) then
            write(unit1, intfmt) (numdg(m),m=1,nbasistot)
          else
            write(unit1) (numdg(m),m=1,nbasistot)
          endif
        endif

C Write out listd array
        do m = 1,nbasistot
            ml = m
            if (fmto) then
              write(unit1, intfmt)
     .             (listd(listdptr(ml)+im),im=1,numd(ml))
            else
              write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
        enddo


C Write density matrix
        do is=1,nspin
          do m=1,nbasistot
              ml = m
              if (fmto) then
                write(unit1, floatfmt) 
     .               (dm(listdptr(ml)+im,is),im=1,numd(ml))
              else
                write(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
              endif
          enddo
        enddo

        if (Node.eq.0) then
          call io_close(unit1)
        endif

      else
        if (Node.eq.0) then
          call die('iodm: incorrect task')
        endif
      endif

C Deallocate local buffer array for globalised numd
      deallocate(numdg)

      end subroutine iodm
