! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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
      use basis_types, only: deallocate_spec_arrays, initialize
      use basis_types, only: iz, lmxkb, nkbl, 
     &           erefkb, lmxo, nzeta, rco, 
     &           lambda, filtercut,
     &           atm_label, polorb, semic, nsemic,
     &           cnfigmx, charge, smass, basistype,
     &           rinn, vcte, qcoe, qyuk, qwid, split_norm
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
      use atmparams, only: lmaxd, nkbmx, nsemx, nzetmx
      use atom_options, only: get_atom_options
      use dftu_specs, only: read_dftu_specs
      use dftu_specs, only: dftu_proj_gen

      use pseudopotential, only: pseudo_read
    
      use chemical

      use m_spin, only: SpOrb

#ifdef HAVE_LIBINT
!      use atm_types, only: nspecies
      use atomlist,  only: iaorb, indxua, iphorb, no_s
!      use siestaXC,  only: getXC
      use nao2gto_data, only: hfx_gtos, hfx_libderiv, hfx_libint,
     &                        hfx_options
      use nao2gto_io
      use siesta_options , only: hfx_wanted
#endif

      implicit none
      integer,         intent(out) :: ns   ! Number of species
!     Internal variables ...................................................
      integer                      :: is
      logical                      :: user_basis, user_basis_netcdf
      logical :: req_init_setup
      type(basis_def_t),   pointer :: basp

      external atm_transfer

      call get_atom_options()


!     Reading input for the pseudopotentials and atomic orbitals
      write(6,'(/2a)') 
     &    'initatom: Reading input for the pseudopotentials ',
     &    'and atomic orbitals ----------'

      user_basis = fdf_boolean('user-basis',.false.)
      user_basis_netcdf = fdf_boolean('user-basis-netcdf',.false.)

      ! Create list of options NOT compatible with psf/vps file
      ! reads.
      req_init_setup = fdf_defined('DFTU.proj') .or.
     &    fdf_defined('LDAU.proj')
      ! Add any other dependencies here...

      ! Initialize all basis-parameters
      call read_chemical_types()
      nsp = number_of_species()
          
      allocate(basis_parameters(nsp))
      do is = 1 , nsp
        call initialize(basis_parameters(is))
      end do
      
      ! Check that the user can perform a legal action
      req_init_setup = req_init_setup .and.
     & ( user_basis_netcdf .or. user_basis )
      if ( req_init_setup ) then
         call die('Reading PAOs and KBs from NetCDF/ascii files '//
     &'is not possible with DFTU.Proj')
      end if
      
      if (user_basis_netcdf) then
        write(6,'(/a)') 'Reading PAOs and KBs from NetCDF files...'
        call read_basis_netcdf(ns)
        call elec_corr_setup()
      else if (user_basis) then

       if ( SpOrb ) then  
          write(6,'(a)') ' initatom: Spin configuration = spin-orbit'
          do is = 1 , nsp
             basp => basis_parameters(is)
             
             basp%label = species_label(is)
             call pseudo_read(basp%label,basp%pseudopotential)
          end do
       end if
       write(6,'(/a)') 'Reading PAOs and KBs from ascii files...'
       call read_basis_ascii(ns)
       call elec_corr_setup()
      else
!       New routines in basis_specs and basis_types.
        call read_basis_specs()
        call basis_specs_transfer()

!       Get the parameters for the generation of the LDA+U projectors
        call read_dftu_specs()

        nsmax = nsp             !! For old_atmfuncs
        call allocate_old_arrays()
        call clear_tables()

        do is = 1,nsp
          call write_basis_specs(6,is)
          basp=>basis_parameters(is)
          call ATOM_MAIN( iz(is), lmxkb(is), nkbl(0:lmaxd,is),
     &                    erefkb(1:nkbmx,0:lmaxd,is), lmxo(is),
     &                    nzeta(0:lmaxd,1:nsemx,is),
     &                    rco(1:nzetmx,0:lmaxd,1:nsemx,is),
     &                    lambda(1:nzetmx,0:lmaxd,1:nsemx,is),
     &                    atm_label(is), polorb(0:lmaxd,1:nsemx,is),
     &                    semic(is), nsemic(0:lmaxd,is),
     &                    cnfigmx(0:lmaxd,is), charge(is), smass(is),
     &                    basistype(is), is, rinn(0:lmaxd,1:nsemx,is),
     &                    vcte(0:lmaxd,1:nsemx,is),
     &                    qcoe(0:lmaxd,1:nsemx,is),
     &                    qyuk(0:lmaxd,1:nsemx,is),
     &                    qwid(0:lmaxd,1:nsemx,is),
     &                    split_norm(0:lmaxd,1:nsemx,is), 
     &                    filtercut(0:lmaxd,1:nsemx,is), basp)
!         Generate the projectors for the LDA+U simulations (if requested)
          call dftu_proj_gen(is)
        enddo 

        call prinput(nsp)

!       Create the new data structures for atmfuncs.
        call atm_transfer()
        call deallocate_old_arrays()
        call elec_corr_setup()
        ns = nsp               ! Set number of species for main program

      ! Prepare hybrid XC calculations
      ! Read options and Gaussian fitting parameters from FDF
!      if ( hfx_wanted ) then
        allocate(hfx_gtos(ns))
        call nao2gto_transfer(hfx_gtos, hfx_options)
!      endif

      endif

      call dump_basis_ascii()
      call dump_basis_netcdf()
      call dump_basis_xml()

      if (.not. user_basis .and. .not. user_basis_netcdf) then
        call deallocate_spec_arrays()
      endif

      end subroutine initatom

