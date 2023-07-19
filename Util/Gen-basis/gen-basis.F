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
      program gen_basis

C Stand-alone program to generate the PAOs, KB projectors, Vlocal, etc
C for a set of atoms. 
C
C The input is exactly the same as that for SIESTA.
C This program generates .ion (optionally .ion.nc) files readable
C by SIESTA
C
C RESTRICTION: Just a single species allowed
C
!!      use f2kcli
      use fdf

      use precision
      use parallel,     only : Node, Nodes, IOnode
      use basis_types, only: nsp, basis_specs_transfer
      use basis_types, only: iz, lmxkb, nkbl, 
     .           erefkb, lmxo, nzeta, rco, 
     .           lambda,
     .           atm_label, polorb, semic, nsemic,
     .           cnfigmx, charge, smass, basistype,
     $           rinn, vcte, split_norm
      use basis_types, only: write_basis_specs
      use basis_types, only: basis_def_t, basis_parameters
      use basis_specs, only: read_basis_specs
      use basis_io
      use old_atmfuncs, only : nsmax
      use old_atmfuncs, only : clear_tables, allocate_old_arrays

      use atom,         only : atom_main, prinput
      use xcmod,        only : setXC
      use sys,          only : die
      use atom_options, only : get_atom_options, write_ion_plot_files
C
C     It has to be MPI-aware since it calls routines down the line...
C
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer is
!     integer narg
      character(len=50) filein

      type(basis_def_t), pointer   :: basp

#ifdef MPI
      integer  MPIerror
#endif
      external atm_transfer

c Reading input for the pseudopotentials and atomic orbitals 
C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

      IOnode = (Node .eq. 0)

      if (Ionode) then
!
!        Problems with mpich command-line argument handling...
!$$$         narg = command_argument_count()
!$$$         if (narg == 0) then
!$$$            filein = "stdin"
!$$$         else if (narg == 1) then
!$$$            call get_command_argument(1,filein)
!$$$         else
!$$$            stop "too many command-line arguments"
!$$$         endif
!
         filein = "stdin"
         call fdf_init(filein,'gen-basis_fdf.log')

         call get_atom_options()
         write_ion_plot_files = .true.  ! get them always

         call setXC

         call read_basis_specs()
         call basis_specs_transfer()

         nsmax = nsp            ! For old_atmfuncs
         call allocate_old_arrays
         call clear_tables

         if (nsp .gt. 1)
     $        call die("gen-basis: JUST A SINGLE SPECIES ALLOWED...")

         do is = 1,nsp
            call write_basis_specs(6,is)
            basp=>basis_parameters(is)
            call atom_main( iz(is), lmxkb(is),
     $           nkbl(0:,is), erefkb(1:,0:,is),lmxo(is),
     $           nzeta(0:,1:,is), rco(1:,0:,1:,is), 
     $           lambda(1:,0:,1:,is), atm_label(is),
     $           polorb(0:,1:,is), semic(is), nsemic(0:,is),
     $           cnfigmx(0:,is),charge(is),
     $           smass(is), basistype(is), is,
     $           rinn(0:,1:,is), vcte(0:,1:,is),
     $           split_norm(0:,1:,is), basp)
         enddo 

         call prinput(nsp)

!        Create the new data structures

         call atm_transfer

         call dump_basis_ascii
         call dump_basis_netcdf

      endif

      end program gen_basis
!






