! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      program gen_basis

C Stand-alone program to generate the PAOs, KB projectors, Vlocal, etc
C for a set of atoms. 
C
C The input is exactly the same as that for SIESTA.
C This program generates .ion/.ion.xml (optionally .ion.nc) files readable
C by SIESTA
C
      use fdf

      use precision
      use parallel,     only : Node, Nodes, IOnode
      use basis_types, only: nsp, basis_specs_transfer, initialize
      use basis_types, only: iz, lmxkb, nkbl, 
     .           erefkb, lmxo, nzeta, rco, 
     .           lambda, filtercut,
     .           atm_label, polorb, semic, nsemic,
     .           cnfigmx, charge, smass, basistype,
     &           rinn, vcte, qcoe, qyuk, qwid, split_norm
      use basis_types, only: write_basis_specs
      use basis_types, only: basis_def_t, basis_parameters
      use basis_specs, only: read_basis_specs
      use basis_io
      use old_atmfuncs, only : nsmax
      use old_atmfuncs, only : clear_tables, allocate_old_arrays

      use chemical, only: read_chemical_types, number_of_species
      use atom,         only : atom_main, prinput
      use sys,          only : die
      use atom_options, only : get_atom_options, write_ion_plot_files
      use dftu_specs, only: read_dftu_specs
      use dftu_specs, only: dftu_proj_gen

      implicit none

      integer is
      character(len=64) :: filein
      type(basis_def_t), pointer   :: basp
      integer :: narg, length
      logical :: file_exists
      
      external :: atm_transfer

! Reading input for the pseudopotentials and atomic orbitals 
      Node =  0
      Nodes = 1
      IOnode = (Node == 0)

!     Get file-name from input line
      narg = command_argument_count()
      if ( narg > 0 ) then
        filein = ' '
        call get_command_argument(narg,filein,length)
        if ( length > len(filein) ) then
          call die ('The argument is too long to be retrieved, &
     &please limit your-self to 50 characters for the input file') 
        end if
        inquire(file=filein, exist=file_exists)
        if ( .not. file_exists ) then
          call die('Input file '//trim(filein)//' does not exist? &
     &Have you specified the wrong file-name?') 
        end if
        write(*,'(2a/)') 'gen-basis: reading from file: ', trim(filein)
      else
        filein = 'stdin'
        write(*,'(a/)') 'gen-basis: reading from std-in'
      end if
      
      call fdf_init(filein,'gen-basis_fdf.log')
      
      call get_atom_options()
      write_ion_plot_files = .true. ! get them always
      
      call read_xc_info()
      
!     Read chemical types
      call read_chemical_types(.true.)
      nsp = number_of_species()

      allocate(basis_parameters(nsp))
      do is = 1 , nsp
        call initialize(basis_parameters(is))
      end do

      call read_basis_specs()
      call basis_specs_transfer()

!     Get the parameters for the generation of the DFT+U projectors
      call read_dftu_specs()

      nsmax = nsp               ! For old_atmfuncs
      call allocate_old_arrays
      call clear_tables

      do is = 1,nsp
        call write_basis_specs(6,is)
        basp=>basis_parameters(is)
        call atom_main( iz(is), lmxkb(is),
     $      nkbl(0:,is), erefkb(1:,0:,is),lmxo(is),
     $      nzeta(0:,1:,is), rco(1:,0:,1:,is), 
     $      lambda(1:,0:,1:,is), atm_label(is),
     $      polorb(0:,1:,is), semic(is), nsemic(0:,is),
     $      cnfigmx(0:,is),charge(is),
     $      smass(is), basistype(is), is,
     $      rinn(0:,1:,is), vcte(0:,1:,is),
     &      qcoe(0:,1:,is),
     &      qyuk(0:,1:,is),
     &      qwid(0:,1:,is),
     $      split_norm(0:,1:,is),
     &      filtercut(0:,1:,is), basp)
!       Generate the projectors for the LDA+U simulations (if requested)
        call dftu_proj_gen(is)
      enddo 

      call prinput(nsp)

!     Create the new data structures
      call atm_transfer
      
      call dump_basis_ascii
      call dump_basis_netcdf
      call dump_basis_xml()

      end program gen_basis
