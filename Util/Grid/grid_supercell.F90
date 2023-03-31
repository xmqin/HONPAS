! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program grid_supercell
!
! Expands a Grid (binary) file into a supercell
!
! x-->n1*x
! y-->n2*y
! z-->n3*z
!
! Usage: grid_supercell
!
!        The input file must be named "GRIDFUNC" (a symbolic link would do)
!        The output file is called "GRIDFUNC_SUPERCELL"
!
implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

integer   ::    nspin  ! Number of spins 
integer   ::    mesh(3), nmesh(3)

integer   ::   n1, n2, n3
integer   ::   ispin, iostat, ix, iy, iz, i1, i2, i3, jx, jy, jz

real(dp)  ::    cell(3,3), ncell(3,3)
real(sp), dimension(:,:,:), allocatable :: gridfunc, ngridfunc

!-----------------------------------------------------

print *, "n1, n2, n3:"
read *, n1, n2, n3
write(0,*) "using: ", n1, n2, n3

open(unit=1,file="GRIDFUNC",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File GRIDFUNC cannot be opened"
  STOP
endif
open(unit=2,file="GRIDFUNC_SUPERCELL",form="unformatted",status="unknown",action="write", &
            position="rewind",iostat=iostat)

read(1) cell(1:3,1:3)
read(1) mesh(1:3), nspin

!
nmesh(1) = mesh(1)*n1
nmesh(2) = mesh(2)*n2
nmesh(3) = mesh(3)*n3

ncell(1:3,1) = cell(1:3,1)*n1
ncell(1:3,2) = cell(1:3,2)*n2
ncell(1:3,3) = cell(1:3,3)*n3

write(2) ncell(1:3,1:3)
write(2) nmesh(1:3), nspin

allocate(gridfunc(1:mesh(1),1:mesh(2),1:mesh(3)))
allocate(ngridfunc(1:nmesh(1),1:nmesh(2),1:nmesh(3)))

!      Can be done with less memory...

       do ispin=1,nspin

          do iz=1,mesh(3)
             do iy=1,mesh(2)
                read(1) (gridfunc(ix,iy,iz),ix=1,mesh(1))
             enddo
          enddo

	  do i1=0,n1-1
	  do i2=0,n2-1
          do i3=0,n3-1
             do iz=1,mesh(3)
	        jz = iz + i3*n3
                do iy=1,mesh(2)
                   jy = iy + i2*n2
                   do ix=1,mesh(1)
                      jx = ix + i1*n1
                      ngridfunc(jx,jy,jz) = gridfunc(ix,iy,iz)
                   enddo
                enddo
             enddo
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

     end program grid_supercell
