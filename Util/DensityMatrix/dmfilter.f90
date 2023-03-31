! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program dmfilter

!
! Filters  a DM 
!
! Usage: dmfilter < datafile
! 
!        The input file must be named "DM" (a symbolic link would do)
!        The output file is called "DMFILTERED"
!
!        datafile contains the list of orbitals whose contribution is to be kept, one per line
!        Use a real file -- do not use with standard input
!
!        The orbital numbers can be obtained from the .stt file produced by mprop, for example.
!
implicit none

integer, parameter  :: dp = selected_real_kind(14,100)

integer   ::    norbs  ! Number of atomic orbitals
integer   ::    nnzs   ! Total number of orbital interactions
integer   ::    nspin  ! Number of spins 

integer   ::    iog, ispin, j, iostat, istart, iend, io2, ind, i

integer, dimension(:), allocatable  :: numd 
integer, dimension(:), allocatable  :: row_pointer 
integer, dimension(:), allocatable  :: column
real(dp), dimension(:,:), allocatable :: dm
logical, dimension(:), allocatable    :: mask

!-----------------------------------------------------

open(unit=1,file="DM",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File DM cannot be opened"
  STOP
endif

read(1) norbs, nspin
print *, "Norbs, nspin: ", norbs, nspin

allocate(numd(1:norbs),row_pointer(1:norbs),mask(1:norbs))

read(1) (numd(iog),iog=1,norbs)
nnzs = sum(numd(1:norbs))
!
allocate(column(1:nnzs),dm(1:nnzs,1:nspin))

 ! Compute global row pointer
   row_pointer(1) = 0
   do iog=2,norbs
      row_pointer(iog) = row_pointer(iog-1) + numd(iog-1)
   enddo

!
!   Column information
!
do iog = 1, norbs
  read(1) (column(row_pointer(iog)+j),j=1,numd(iog))
enddo

!
!  DM
!
   do ispin = 1, nspin        
      do iog = 1, norbs
            read(1) (dm(row_pointer(iog)+j,ispin),j=1,numd(iog))
      enddo
   enddo

   close(1)

!
!  Set the mask to the list of orbitals.
!  If what you want to do is exclude the orbitals in the list, you
!  can initialize mask to ".true." and set mask(i) to ".false." below.
!
mask=.false.
do
   read(5,fmt=*,iostat=iostat) i
   if (iostat /= 0) exit
   mask(i) = .true.
enddo

   do ispin = 1, nspin        
      do iog = 1, norbs
         if (.not. mask(iog))  then
            istart = row_pointer(iog) + 1
            iend = row_pointer(iog) + numd(iog)
            dm(istart:iend,ispin) = 0.0_dp
         else
            do j = 1, numd(iog)
               ind = row_pointer(iog)+j
               io2 = mod(column(ind),norbs)
               if (io2 == 0) io2 = norbs
               if (.not. mask(io2))  dm(ind,ispin) = 0.0_dp
            enddo
         endif
      enddo
   enddo

open(unit=1,file="DMFILTERED",form="unformatted",status="new",action="write", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File DMFILTERED exists"
  STOP
endif

write(1) norbs, nspin
write(1) (numd(iog),iog=1,norbs)
do iog = 1, norbs
  write(1) (column(row_pointer(iog)+j),j=1,numd(iog))
enddo
do ispin = 1, nspin        
   do iog = 1, norbs
      write(1) (dm(row_pointer(iog)+j,ispin),j=1,numd(iog))
   enddo
enddo

close(1)


end program dmfilter
