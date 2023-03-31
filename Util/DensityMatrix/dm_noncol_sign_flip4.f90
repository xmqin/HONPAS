! 
! Copyright (C) 1996-2018	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program dm_noncol_sign_flip4

!
! Changes the sign of the "4" component of a DM (useful to re-use DM files created
! with versions prior to 4.0.2 in the non-collinear case).
!
! Usage: dm_noncol_sign_flip4
! 
!        The input file must be named "DM" (a symbolic link would do)
!        The output file is called "DMSIGNFLIP4"
!
implicit none

integer, parameter  :: dp = selected_real_kind(14,100)

integer   ::    norbs  ! Number of atomic orbitals
integer   ::    nnzs   ! Total number of orbital interactions
integer   ::    nspin  ! Number of spins 

integer   ::    iog, ispin, iostat, j

integer, dimension(:), allocatable  :: numd 
integer, dimension(:), allocatable  :: row_pointer 
integer, dimension(:), allocatable  :: column
real(dp), dimension(:,:), allocatable :: dm

!-----------------------------------------------------

open(unit=1,file="DM",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File DM cannot be opened"
  STOP
endif

read(1) norbs, nspin
print *, "Norbs, nspin: ", norbs, nspin

if (nspin /= 4) then
   print *, "The DM is not from a non-collinear spin calculation"
   STOP
endif

allocate(numd(1:norbs),row_pointer(1:norbs))

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
   ! Change sign
   !
   dm(:,4) = -dm(:,4)
!

open(unit=1,file="DMSIGNFLIP4",form="unformatted",status="new",action="write", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "A file DMSIGNFLIP4 exists. Please remove or change its name"
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


end program dm_noncol_sign_flip4
