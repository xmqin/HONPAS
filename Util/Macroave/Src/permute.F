! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      program permute
      use f2kcli

c****************************************************************************
c
c This program reads files with info on the grid from SIESTA
c and changes the axes of the grid so that the former "x" direction
c is turned into the "z" direction.
c This is useful to run Macroave when the user has chosen "x" as the 
c direction perpendicular to a slab for performance reasons (i.e., better
c load-balancing).
c
c             OLD           NEW
c
c             x              z
c             y              x
c             z              y
c
c
c      The compiler flags should be compatible with those you used to
c      compile Siesta, as permute reads a binary file whose data
c      layout can change. This includes 32bit vs 64 bit flags, if any.
c
c Usage:
c     
c        permute gridfile 
c
c      The permutted file will be named gridfile.permutted
c
c   CAVEATS: 
c
c     -  For Files holding information for several spin components,
c        only the first is processed.
c
c****************************************************************************

      implicit none

      integer           ipt, isp, ix, iy, iz, i, j,
     .                  mesh(3), nspin, nargs, iostat,
     $                  mesh_new(3)

      character         fnamein*75
      integer, parameter  :: dp = selected_real_kind(14,100)

      real, dimension(:,:,:), allocatable  :: rho
      real(dp)  :: xfrac(3), val, fmin, fmax, cell(3,3), cell_new(3,3)

c ---------------------------------------------------------------------------


      nargs = command_argument_count()
      if (nargs /= 1)  Stop "Usage: permute filename"

      call get_command_argument(1,value=fnamein,status=iostat)

c read function from the 3D grid --------------------------------------------

      open( unit=1, file=fnamein, form="unformatted", status='old' )

      read(1) cell
  
      write(0,*) 
      write(0,*) 'Cell vectors'
      write(0,*) 
      write(0,*) cell(1,1),cell(2,1),cell(3,1)
      write(0,*) cell(1,2),cell(2,2),cell(3,2)
      write(0,*) cell(1,3),cell(2,3),cell(3,3)

      read(1) mesh, nspin

      write(0,*) 
      write(0,*) 'Grid mesh: ',mesh(1),'x',mesh(2),'x',mesh(3)
      write(0,*) 
      write(0,*) 'nspin = ',nspin
      if (nspin > 1) write(0,*)
     $     "** Only 1st spin info can be read at this point"

      allocate(rho(0:mesh(1)-1, 0:mesh(2)-1, 0:mesh(3)-1))

          do iz=0,mesh(3)-1
             do iy=0,mesh(2)-1
                read(1) (rho(ix,iy,iz),ix=0,mesh(1)-1)
             enddo
          enddo

      close(1)

      open( unit=1, file=trim(fnamein) // ".permutted",
     $                 form="unformatted", status='new' )

      cell_new(1,1) = cell(2,2)
      cell_new(2,1) = cell(3,2)
      cell_new(3,1) = cell(1,2)
      cell_new(1,2) = cell(2,3)
      cell_new(2,2) = cell(3,3)
      cell_new(3,2) = cell(1,3)
      cell_new(1,3) = cell(2,1)
      cell_new(2,3) = cell(3,1)
      cell_new(3,3) = cell(1,1)

      write(0,*) 
      write(0,*) 'New cell vectors'
      write(0,*) 
      write(0,*) cell_new(1,1),cell_new(2,1),cell_new(3,1)
      write(0,*) cell_new(1,2),cell_new(2,2),cell_new(3,2)
      write(0,*) cell_new(1,3),cell_new(2,3),cell_new(3,3)

      write(1) cell_new

      mesh_new(1) = mesh(2)
      mesh_new(2) = mesh(3)
      mesh_new(3) = mesh(1)
      nspin = 1

      write(1) mesh_new, nspin

      do iz=0,mesh_new(3)-1
         do iy=0,mesh_new(2)-1
cccc            rho_new(ix,iy,iz) = rho(iz,ix,iy)
            write(1) (rho(iz,ix,iy),ix=0,mesh_new(1)-1)
         enddo
      enddo

      end program permute
