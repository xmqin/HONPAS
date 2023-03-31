! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!     
!     This file is part of the SIESTA package.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     

      program g2c_ng

      use m_getopts
      use m_struct,   only: struct_t, read_struct
      use m_struct,   only: replicate_struct
      use m_gridfunc, only: read_gridfunc, gridfunc_t, clean_gridfunc

c****************************************************************************
c     GRID2CUBE_NG
c     Work in progress
      
c     This program transforms files with info on the grid from SIESTA
c     to the Gaussian CUBE format, that can be read by packages like MOLEKEL
c     MOLDEN, and VESTA, to plot charges, potentials, etc. on real space.

c     This version removes the restriction to orthorhombic cells, and can generate
c     data in a supercell of the original cell, but it lacks the ability to do
c     origin changes (for now).
c     It uses a STRUCT file instead of XV (cleaner, and with fractional coords)
c     
c     Written by P. Ordejon, August 2002
C     Last modification: June 2003
C     Rewritten by A. Garcia (2015)
c****************************************************************************
c     USAGE:
c     
c     This program reads files generated from SIESTA, with information of
c     physical quantities on the real space grid, such as charge densities
c     (filename.RHO, filename.DRHO, filename.LDOS) or potentials (filename.VH,
c     filename.VT) and translates them to the Gaussian CUBE format.
c     
c     The program needs two input files, specified through the command line:
c     
c     A STRUCT file, for geometry information
c     A GRIDFUNC file (LDOS, RHO, etc) containing the function data.
c     
c     Optionally, the program can handle a stride (nskip) that
c     specifies the density of grid points in the output. For each
c     direction, only one of every 'nskip' points from the input will be
c     written to the output.  nskip = 1 will yield an output with the
c     same number of grid points than the input; nskip = N will provide
c     an output N**3 less dense than the input. This option is provided
c     because some visualization programs can not handle the large memory
c     involved in large systems with very fine grids, so the grid must be
c     made coarser.
c     
c****************************************************************************

      implicit none

      integer, parameter :: dp = selected_real_kind(14,100)
      integer, parameter :: sp = selected_real_kind(5,10) ! Note

      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels

      integer           nsc(3), i1, i2
      integer           isp, ix, iy, iz, i, ii, natoms, mesh(3), nspin

      integer           idummy
      real(dp)              :: cell(3,3), oldx, celli(3,3), newx
      real(dp)              :: xfrac(3), origin(3)
      real(sp), pointer     :: rho(:,:,:,:) => null()

      character(len=256) :: fnamexv, fnamein, fnameout(2)
      
      type(struct_t) :: str, str0
      type(gridfunc_t), target :: gf
!-----------------------------------
      logical :: supercell = .false.
      logical :: debug = .false.
      logical :: translate_to_01 = .false.
      logical :: structfile_given = .false.
      logical :: gridfile_given = .false.
      integer :: sx = 1
      integer :: sy = 1
      integer :: sz = 1
      integer :: nskip = 1
!     
!     Process options
!     
      n_opts = 0
      do
      call getopts('hdts:g:x:y:z:n:',opt_name,opt_arg,n_opts,iostat)
      if (iostat /= 0) exit
      select case(opt_name)
      case ('d')
         debug = .true.
      case ('t')
         translate_to_01 = .true.
      case ('n')
         read(opt_arg,*) nskip
      case ('s')
         read(opt_arg,*) fnamexv
         structfile_given = .true.
      case ('g')
         read(opt_arg,*) fnamein
         gridfile_given = .true.
      case ('x')
         read(opt_arg,*) sx
      case ('y')
         read(opt_arg,*) sy
      case ('z')
         read(opt_arg,*) sz
      case ('h')
         call write_manual()
         STOP
      case ('?',':')
         write(0,*) "Invalid option: ", opt_arg(1:1)
         write(0,*)
         call write_manual()
         STOP
      end select
      enddo

      if (.not. structfile_given) then
         write(0,*) "You must give the structure file ( -s option)"
         write(0,*)
         call write_manual()
         STOP
      endif
      if (.not. gridfile_given) then
         write(0,*) "You must give the grid file ( -g option)"
         write(0,*)
         call write_manual()
         STOP
      endif

      write(6,*) 
      write(6,*) 'Reading grid data from file ',trim(fnamein)

! Future: enable extension recognition for dispatch to
! fortran or netcdf versions

      call read_gridfunc(fnamein,gf)

      cell = gf%cell

      write(6,*) 
      write(6,*) 'Cell vectors'
      write(6,*) 
      write(6,*) cell(1,1),cell(2,1),cell(3,1)
      write(6,*) cell(1,2),cell(2,2),cell(3,2)
      write(6,*) cell(1,3),cell(2,3),cell(3,3)

      nsc(:) = (/ sx, sy, sz /)

      mesh(:) = nsc(:) * gf%n(:)
      origin(:) = gf%origin(:)
      nspin = gf%nspin
      rho => gf%val
      if (sx*sy*sz /= 1) then
         write(0,*) " We have a supercell:", nsc(:)
         supercell = .true.
         do i = 1, 3
            if (abs(origin(i)) > 1.e-10_dp) then
               if (nsc(i) /= 1) then
                  write(0,*) " Origin reset to 0 in direction:", i
                  origin(i) = 0.0_dp
               endif
            endif
         enddo
      endif

      write(6,*) 
      write(6,*) 'Grid mesh: ',mesh(1),'x',mesh(2),'x',mesh(3)
      write(6,*) 
      write(6,*) 'nspin = ',nspin
      write(6,*) 

      if (supercell) then
         call read_struct(trim(fnamexv), str0)
         if (any(abs(cell(:,:)-str0%cell(:,:)) > 1.0e-6_dp)) then
            write(0,*) "Cell mismatch"
         endif
         call replicate_struct(str0,nsc,str)
      else
         call read_struct(trim(fnamexv), str)
         if (any(abs(cell(:,:)-str%cell(:,:)) > 1.0e-6_dp)) then
            write(0,*) "Cell mismatch"
         endif
      endif

      cell(:,:) = str%cell(:,:)
      
      natoms = str%na
      if (translate_to_01) then
         call reclat(str%cell, celli, 0)
         do i=1,natoms
            xfrac(:) = matmul(transpose(celli),str%xa(:,i))
            do ix = 1, 3
               oldx = xfrac(ix)
               newx = mod(oldx+1000_dp,1.0_dp)
               xfrac(ix) = newx
               if (debug .and. (newx /= oldx)) then
                  write(0,*) "... translating ",
     $                 oldx, " to ", newx
               endif
            enddo
            str%xa(:,i) = matmul(str%cell,xfrac(1:3))
         enddo
      endif

      if (nspin .eq. 1) then
         fnameout(1) = "Grid.cube"
      else if (nspin .eq. 2) then
         fnameout(1) = "Up.cube"
         fnameout(2) = "Down.cube"
      else 
         call die('nspin must be either 1 or 2')
      endif

      do isp=1,nspin

         write(6,*) 'Writing CUBE file ',trim(fnameout(isp))
         open( unit=2, file=fnameout(isp), form='formatted')
! Arbitrary two first lines
         write(2,*) trim(fnameout(isp))
         write(2,*) "z fastest, then y, x"
! Natoms and origin of volumetric data
         write(2,'(i5,4f12.6)') natoms, origin(:)


! Write vectors defining voxels
         do ix=1,3
            ii = mesh(ix)/nskip
            if (ii*nskip .ne. mesh(ix)) ii = ii+1
            write(2,'(i5,4f12.6)') ii,(cell(iy,ix)/ii,iy=1,3)
         enddo
         
! Atomic numbers, ?, and coordinates
         do i=1,natoms
            write(2,'(i5,4f12.6)') str%iza(i),0.0,(str%xa(ix,i),ix=1,3)
         enddo

! Volumetric data
         do ix=1,mesh(1),nskip
            i1 = rep(ix,gf%n(1))
            do iy=1,mesh(2),nskip
               i2 = rep(iy,gf%n(2))
               write(2,'(6e13.5)') (rho(i1,i2,rep(iz,gf%n(3)),isp),
     $              iz=1,mesh(3),nskip)
            enddo
         enddo

         close(2)

      enddo

      if (nspin == 2)  then

! Write diff version
         write(6,*) 'Writing diff file Up-Down.cube'

         open( unit=2, file='Up-Down.cube', form='formatted')
         write(2,*) "Up-Down"
         write(2,*) "Up-Down"
         write(2,'(i5,4f12.6)') natoms, origin(:)

         nskip = 1
         do ix=1,3
            ii = mesh(ix)/nskip
            if (ii*nskip .ne. mesh(ix)) ii = ii+1
            write(2,'(i5,4f12.6)') ii,(cell(iy,ix)/ii,iy=1,3)
         enddo

         do i=1,natoms
            write(2,'(i5,4f12.6)') str%iza(i),0.0,(str%xa(ix,i),ix=1,3)
         enddo

         do ix=1,mesh(1),nskip
            i1 = rep(ix,gf%n(1))
            do iy=1,mesh(2),nskip
               i2 = rep(iy,gf%n(2))
               write(2,'(6e13.5)') ((rho(i1,i2,rep(iz,gf%n(3)),1) -
     $                              rho(i1,i2,rep(iz,gf%n(3)),2)),
     $              iz=1,mesh(3),nskip)
            enddo
         enddo

         close(2)

! Write sum version
         write(6,*) 'Writing sum file Up+Down.cube'

         open( unit=2, file='Up+Down.cube', form='formatted')
         write(2,*) "Up+Down"
         write(2,*) "Up+Down"
         write(2,'(i5,4f12.6)') natoms, origin(:)

         nskip = 1
         do ix=1,3
            ii = mesh(ix)/nskip
            if (ii*nskip .ne. mesh(ix)) ii = ii+1
            write(2,'(i5,4f12.6)') ii,(cell(iy,ix)/ii,iy=1,3)
         enddo

         do i=1,natoms
            write(2,'(i5,4f12.6)') str%iza(i),0.0,(str%xa(ix,i),ix=1,3)
         enddo

         do ix=1,mesh(1),nskip
            i1 = rep(ix,gf%n(1))
            do iy=1,mesh(2),nskip
               i2 = rep(iy,gf%n(2))
               write(2,'(6e13.5)') ((rho(i1,i2,rep(iz,gf%n(3)),1) +
     $                              rho(i1,i2,rep(iz,gf%n(3)),2)),
     $              iz=1,mesh(3),nskip)
            enddo
         enddo

      endif
      close(2)

      write(6,*) 

      call clean_gridfunc(gf)

      CONTAINS

      ! Representative
      function rep(i,n) result(i1)
      integer, intent(in) :: i
      integer, intent(in) :: n
      integer             :: i1

      i1 = mod(i-1,n) + 1
      end function rep

      subroutine write_manual()
      
               write(0,*) "Usage: g2c_ng [-dht] [-n NSKIP]" //
     $              "  -s STRUCT_FILE -g GRIDFILE"
         write(0,*) " -n NSKIP : stride in file"
         write(0,*) " -t : translate coordinates to [0,1)"
         write(0,*) " -x : supercell factor along 1st lattice vector"
         write(0,*) " -y : supercell factor along 2nd lattice vector"
         write(0,*) " -z : supercell factor along 3rd lattice vector"
         write(0,*) " -t : translate coordinates to [0,1)"
         write(0,*) " -d : debug"

      
       end subroutine write_manual
      
      end program g2c_ng


