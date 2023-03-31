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
program dmUnblock
!
! Converts a (possibly) blocked DM file to classic DM format
! Blocked files are written faster to disk, but they are not
! compatible with the reading routines in Siesta.
! Hence, before re-using a blocked DM or H file in Siesta, it has
! to be converted to "classic" format by this program
!
! This state of affairs will change once proper parallel-IO is implemented.
!
! The "dm" in the name refers to the DM format, but H and S, and other
! Siesta sparse matrices, can also be handled by it.
!
! Usage: 
!
!          dmUnblock [-h] -i BlockedFile -o ClassicFile'
!
!
! If the input file is actually unblocked, a simple copy is carried out.

  use m_getopts
  use f2kcli

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)

integer   ::    no_u   ! Number of atomic orbitals
integer   ::    nspin  ! Number of spins 
integer   ::    blocksize     ! blocksize

integer   ::    iog, ios, ispin, i, j
integer   ::    maxnnzbs, nblocks, norbs, n_g, norbs_g
integer   ::    base, nsize, nnzbs, nnzs_bg, ptr, nnzsi

integer, dimension(:), allocatable  :: numdg
integer, dimension(:), allocatable  :: ibuffer
real(dp), dimension(:), allocatable :: buffer

  integer  :: n_opts, iostat, nargs, nlabels
  character(len=256) :: opt_name, opt_arg
  character(len=256) :: filein, fileout
  logical :: filein_given, fileout_given, files_ok

!-----------------------------------------------------
  !
  !     Process options
  !
  filein_given = .false.
  fileout_given = .false.
  n_opts = 0
  do
     call getopts('i:o:h',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('i')
        read(opt_arg,*) filein
        filein_given = .true.
     case ('o')
        read(opt_arg,*) fileout
        fileout_given = .true.
     case ('h')
        call manual()
        STOP
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        call manual()
        STOP
     end select
  enddo

  files_ok = filein_given .and. fileout_given

!  nargs = command_argument_count()
!  nlabels = nargs - n_opts + 1
  if ((.not. files_OK)) then
     call manual()
     STOP
  endif

open(unit=1,file=filein,form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File " // trim(filein) // " cannot be opened"
  STOP
endif

open(unit=2,file=fileout,form="unformatted",status="unknown",action="write", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File " // trim(fileout) // " cannot be opened"
  STOP
endif


read(1,iostat=iostat) no_u, nspin, blocksize
if (iostat /= 0) then
   ! the reading failure is due to the absence of blocksize in the record
   rewind(1)
   read(1,iostat=iostat) no_u, nspin
   write(0,*) trim(filein) // " is a classic file with Norbs, nspin: ", &
              no_u, nspin
   write(0,*) "Siesta can read it directly"
   write(0,*) "I will just copy " // trim(filein) // " to " // trim(fileout)
   !
   ! The "copy" is actually done by setting blocksize to 1 in the general code
   !
   blocksize = 1
else
   write(0,*) trim(filein) // " is a blocked file. Norbs, nspin, blocksize: ", &
              no_u, nspin, blocksize
endif

write(2) no_u, nspin

allocate(numdg(1:no_u))

read(1) (numdg(iog),iog=1,no_u)
write(2) (numdg(iog),iog=1,no_u)

!     Find out how big the buffer has to be

         maxnnzbs = 0
         nblocks = 0
         norbs = 0
         do 
            base = nblocks*blocksize
            nsize = min(blocksize,no_u-norbs)
            nnzbs = sum(numdg(base+1:base+nsize))
            if (nnzbs > maxnnzbs) maxnnzbs = nnzbs
            norbs = norbs + nsize
            if (norbs == no_u) EXIT
            nblocks = nblocks + 1
         enddo
         !print *, "Maznnzbs = ", maxnnzbs

         allocate(buffer(maxnnzbs), ibuffer(maxnnzbs))


    n_g = 0
    do
       norbs_g = min(blocksize,no_u-n_g)
       nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
       read(1) (ibuffer(j),j=1,nnzs_bg)
       ptr = 0
       do i = 1, norbs_g
          nnzsi = numdg(n_g+i)
          write(2) (ibuffer(j),j=ptr+1,ptr+nnzsi)
          ptr = ptr + nnzsi
       enddo
       n_g = n_g + norbs_g
       if (n_g == no_u) EXIT
    enddo

    do ispin = 1, nspin        

      n_g = 0
      do
         norbs_g = min(blocksize,no_u-n_g)
         nnzs_bg = sum(numdg(n_g+1:n_g+norbs_g))
         read(1) (buffer(j),j=1,nnzs_bg)
         ptr = 0
         do i = 1, norbs_g
            nnzsi = numdg(n_g+i)
            write(2) (buffer(j),j=ptr+1,ptr+nnzsi)
            ptr = ptr + nnzsi
         enddo
         n_g = n_g + norbs_g
         if (n_g == no_u) EXIT
      enddo
    enddo

   close(1)
   close(2)


   deallocate(numdg,ibuffer,buffer)

CONTAINS
  subroutine manual()
    write(0,'(a)') ' Usage: dmUnblock [-h] -i BlockedFile -o ClassicFile'
    write(0,*) ! new line
    write(0,'(a)') ' Other options:'
    write(0,*) ! new line
    write(0,'(a)') '     -h        : print help'
  end subroutine manual


 end program dmUnblock
