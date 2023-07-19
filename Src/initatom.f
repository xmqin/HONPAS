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
      subroutine initatom(ns)

! Routine to initialize the Pseudopotentials and Atomic Orbitals.
! Substantially modified by Alberto Garcia (2000)
!
! The PAO and KB information can optionally be read from ASCII files
! (those produced by a standard run of Siesta or Base, but with the extension 
! renamed to '.input' instead of '.dump'), or from NetCDF files (if NetCDF
! is available). Note that there must be files for *all* species.
!
! This behavior is controlled by the FDF logicals 'user-basis' and
! 'user-basis-netcdf'.
!
! The old 'USER' basis type has been removed.
!
! This routine also outputs information about the basis specification
! determined by the routines in the 'basis_specs' modules. 
!


      use fdf
      use precision
      use basis_types, only: basis_specs_transfer, nsp
      use basis_types, only: deallocate_spec_arrays
      use basis_types, only: iz, lmxkb, nkbl, 
     .           erefkb, lmxo, nzeta, rco, 
     .           lambda,
     .           atm_label, polorb, semic, nsemic,
     .           cnfigmx, charge, smass, basistype,
     $           rinn, vcte, split_norm
      use basis_types, only: write_basis_specs
      use basis_types, only: basis_def_t, basis_parameters
      use basis_specs, only: read_basis_specs

      use basis_io, only: read_basis_ascii, read_basis_netcdf
      use basis_io, only: dump_basis_ascii, dump_basis_netcdf
      use basis_io, only: dump_basis_xml

      use old_atmfuncs, only: nsmax, allocate_old_arrays
      use old_atmfuncs, only: clear_tables, deallocate_old_arrays
      use atom, only: atom_main, prinput
      use electrostatic, only: elec_corr_setup
      use atom_options, only: get_atom_options

      implicit none

      integer, intent(out) :: ns   ! Number of species

C Internal variables ...................................................
      integer is

      type(basis_def_t), pointer   :: basp

      logical user_basis, user_basis_netcdf
      
      external atm_transfer

      external nao2gto_transfer
      
      call get_atom_options()
c
c Reading input for the pseudopotentials and atomic orbitals 

      write(6,'(/2a)') 
     .    'initatom: Reading input for the pseudopotentials ',
     .    'and atomic orbitals ----------'

      user_basis = fdf_boolean('user-basis',.false.)
      user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

      if (user_basis_netcdf) then

         write(6,'(/a)') 'Reading PAOs and KBs from NetCDF files...'
         call read_basis_netcdf(ns)
         call elec_corr_setup()

      else if (user_basis) then

         write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
         call read_basis_ascii(ns)
     
         call nao2gto_transfer()

         call elec_corr_setup()

      else
!
!     New routines in basis_specs and basis_types.
!
         call read_basis_specs()
         call basis_specs_transfer()

         nsmax = nsp             !! For old_atmfuncs
         call allocate_old_arrays()
         call clear_tables()

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

!        Create the new data structures for atmfuncs.
         call atm_transfer()
         
         call nao2gto_transfer()
        
         call deallocate_old_arrays()
         call elec_corr_setup()
         ns = nsp               ! Set number of species for main program

      endif

      call dump_basis_ascii()
      call dump_basis_netcdf()
      call dump_basis_xml()

      if (.not. user_basis .and. .not. user_basis_netcdf) then
        call deallocate_spec_arrays()
      endif

      end subroutine initatom




