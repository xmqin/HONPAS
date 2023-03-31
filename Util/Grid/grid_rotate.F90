! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program grid_rotate
!
! Rotates a Grid (binary) file 
!
! x-->z
! y-->x
! z-->y
!
! Usage: grid_rotate
!
!        The input file must be named "GRIDFUNC" (a symbolic link would do)
!        The output file is called "GRIDFUNC_ROTATED"
!
implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

integer   ::    nspin  ! Number of spins 
integer   ::    mesh(3), nmesh(3)

integer   ::   ispin, iostat, ix, iy, iz

real(dp)  ::    cell(3,3), ncell(3,3)
real(sp), dimension(:,:,:), allocatable :: gridfunc, ngridfunc

!-----------------------------------------------------

open(unit=1,file="GRIDFUNC",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File GRIDFUNC cannot be opened"
  STOP
endif
open(unit=2,file="GRIDFUNC_ROTATED",form="unformatted",status="unknown",action="write", &
            position="rewind",iostat=iostat)

read(1) cell(1:3,1:3)
read(1) mesh(1:3), nspin

!
nmesh(1) = mesh(2)
nmesh(2) = mesh(3)
nmesh(3) = mesh(1)

ncell(1:3,1) = cell(1:3,2)
ncell(1:3,2) = cell(1:3,3)
ncell(1:3,3) = cell(1:3,1)

write(2) ncell(1:3,1:3)
write(2) nmesh(1:3), nspin

allocate(gridfunc(1:mesh(1),1:mesh(2),1:mesh(3)))
allocate(ngridfunc(1:nmesh(1),1:nmesh(2),1:nmesh(3)))

       do ispin=1,nspin

          do iz=1,mesh(3)
             do iy=1,mesh(2)
                read(1) (gridfunc(ix,iy,iz),ix=1,mesh(1))
             enddo
          enddo

          do iz=1,mesh(3)
             do iy=1,mesh(2)
                do ix=1,mesh(1)
                   ngridfunc(iy,iz,ix) = gridfunc(ix,iy,iz)
                enddo
             enddo
          enddo
          do iz=1,nmesh(3)
             do iy=1,nmesh(2)
                write(2) (ngridfunc(ix,iy,iz),ix=1,nmesh(1))
             enddo
          enddo

       enddo

       close(1)
       close(2)

     end program grid_rotate
