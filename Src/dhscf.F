! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!
!  The old dhscf has been split in two parts: an initialization routine
!  (dhscf_init), which is called after every geometry change but before
!  the main scf loop, and a dhscf proper, which is called at every step
!  of the scf loop
!
!  NOTE that the mesh initialization part is now done *unconditionally*
!  in dhscf_init, i.e, after *every* geometry change, even if the
!  change does not involve a cell change. The reason is to avoid 
!  complexity, since now the mesh parallel distributions will depend on 
!  the detailed atomic positions even if the cell does not change.
!
!  Besides, the relative cost of a "mesh only" initialization is negligible.
!  The only real observable effect would be a printout of "initmesh" data
!  at every geometry iteration.
!
      module m_dhscf

!
!     To facilitate the communication among dhscf_init and dhscf,
!     some arrays that hold data which do not change during the SCF loop
!     have been made into module variables

!     Some others are scratch, such as nmpl, ntpl, etc

      use precision,      only : dp, grid_p
      use m_dfscf,        only : dfscf
      implicit none

      real(grid_p),   pointer :: rhopcc(:), rhoatm(:), Vna(:)
      real(dp)                :: Uharrs     ! Harris energy
      logical                 :: IsDiag
      integer                 :: nml(3), ntml(3), npcc, nmpl, ntpl
      integer                 :: nbcell
      real(dp)                :: bcell(3,3), cell(3,3),
     &                           dvol, field(3), rmax, scell(3,3)
      real(dp)                :: G2mesh =  0.0_dp
      logical :: debug_dhscf = .false.
      logical :: use_bsc_cellxc
      
      character(len=*), parameter :: debug_fmt =
     &     '('' -- dhscf : Node '',i3,tr1,a,4(tr1,e20.7))'

      public :: dhscf_init, dhscf

      CONTAINS

      subroutine dhscf_init(nspin, norb, iaorb, iphorb,
     &                      nuo, nuotot, nua, na,
     &                      isa, xa, indxua, ucell,
     &                      mscell, G2max, ntm,
     &                      maxnd, numd, listdptr, listd, datm,
     &                      Fal, stressl)

      use precision,      only : dp, grid_p
      use parallel,       only : Node, Nodes
      use atmfuncs,       only : rcut, rcore
      use fdf
      use sys,            only : die
      use mesh,           only : xdsp, nsm, nsp, meshLim
      use parsing

      use alloc,          only : re_alloc, de_alloc
      use siesta_options, only : harrisfun
      use meshsubs,       only : PartialCoreOnMesh
      use meshsubs,       only : NeutralAtomOnMesh
      use meshsubs,       only : PhiOnMesh
      use meshsubs,       only : InitMesh
      use meshsubs,       only : InitAtomMesh
      use meshsubs,       only : setupExtMesh
      use meshsubs,       only : distriPhiOnMesh

      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use moreMeshSubs,   only : UNIFORM, QUADRATIC, LINEAR
      use moreMeshSubs,   only : TO_SEQUENTIAL, TO_CLUSTER, KEEP

      use meshdscf,       only : createLocalDscfPointers

      use iogrid_netcdf, only: set_box_limits
#ifdef NCDF_4
      use files,         only  : slabel
      use m_ncdf_siesta, only : cdf_init_grid
      use m_ncdf_io, only : cdf_init_mesh
#endif
      use m_efield,       only : initialize_efield, acting_efield
      use m_efield,       only : user_specified_field
      use dipole_m,       only : init_dipole_correction

      use m_doping_uniform,       only: initialize_doping_uniform
      use m_doping_uniform,       only: compute_doping_structs_uniform,
     $                                  doping_active
      use m_rhog,                 only: rhog, rhog_in
      use m_rhog,                 only: order_rhog
      use siesta_options,         only: mix_charge
#ifdef NCDF_4
      use siesta_options,         only: write_cdf
#endif
#ifdef MPI
      use mpi_siesta
#endif

      use m_mesh_node,   only: init_mesh_node
      use m_charge_add,  only: init_charge_add
      use m_hartree_add, only: init_hartree_add

      use m_ts_global_vars,only: TSmode
      use m_ts_options,   only : IsVolt, N_Elec, Elecs
      use m_ts_electype,  only : Elec_ortho2semi_inf
      use m_ts_voltage,   only : ts_init_voltage

      implicit none
      integer, intent(in)     :: nspin, norb, iaorb(norb), iphorb(norb),
     &                           nuo, nuotot, nua, na, isa(na),
     &                           indxua(na), mscell(3,3), maxnd,
     &                           numd(nuo), listdptr(nuo), listd(maxnd)
      real(dp), intent(in)    :: xa(3,na), ucell(3,3), datm(norb)
      real(dp), intent(inout) :: g2max
      integer, intent(inout)  :: ntm(3)
      real(dp), intent(inout) :: Fal(3,nua), stressl(3,3)

      real(dp), parameter     :: tiny  = 1.e-12_dp 
      integer                 :: io, ia, iphi, is, n, i, j
      integer                 :: nsc(3), nsd
      real(dp)                :: DStres(3,3), volume
      real(dp), external      :: volcel, ddot
      real(grid_p)            :: dummy_Drho(1,1), dummy_Vaux(1),
     &                           dummy_Vscf(1)
      logical, save           :: frstme = .true.   ! Keeps state
      real(grid_p),   pointer :: Vscf(:,:), rhoatm_par(:)
      integer,        pointer :: numphi(:), numphi_par(:)
      character(len=10)       :: shape

      integer                 :: nm(3)   ! For call to initMesh
      logical                 :: need_gradients_in_xc, is_vdw
      
      ! Transport direction (unit-cell aligned)
      integer                 :: iE
      logical :: is_ortho
!--------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '    PRE dhscf_init' )
#endif
C ----------------------------------------------------------------------
C     General initialisation
C ----------------------------------------------------------------------
C     Start time counter
      call timer( 'DHSCF_Init', 1 )

      nsd = min(nspin,2)
      nullify(Vscf,rhoatm_par)

      if (frstme) then
        debug_dhscf = fdf_get('Debug.DHSCF', .false.)
         
        nullify( xdsp, rhopcc, Vna, rhoatm )
C       nsm lives in module m_dhscf now    !! AG**
        nsm = fdf_integer( 'MeshSubDivisions', 2 )
        nsm = max( nsm, 1 )

C       Set mesh sub-division variables & perform one off allocation
        nsp = nsm*nsm*nsm
        call re_alloc( xdsp, 1, 3, 1, nsp, 'xdsp', 'dhscf_init' )

      endif   ! First time


C ----------------------------------------------------------------------
C     Orbital initialisation : part 1
C ----------------------------------------------------------------------

C     Find the maximum orbital radius
      rmax = 0.0_dp
      do io = 1, norb
        ia   = iaorb(io)    ! Atomic index of each orbital
        iphi = iphorb(io)   ! Orbital index of each  orbital in its atom
        is   = isa(ia)      ! Species index of each atom
        rmax = max( rmax, rcut(is,iphi) )
      enddo

C     Start time counter for mesh initialization
      call timer( 'DHSCF1', 1 )

C ----------------------------------------------------------------------
C     Unit cell handling
C ----------------------------------------------------------------------
C     Find diagonal unit cell and supercell
      call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
      if (.not.IsDiag) then
        if (Node.eq.0) then
          write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     &      'DHSCF: WARNING: New shape of unit cell and supercell:',
     &      ('DHSCF:',(cell(i,j),i=1,3),'   x',nsc(j),j=1,3)
        endif
      endif

C     Find the system shape
      call shaper( cell, nua, isa, xa, shape, nbcell, bcell )

C     Find system volume
      volume = volcel( cell )

C ----------------------------------------------------------------------
C     Mesh initialization 
C ----------------------------------------------------------------------
      call InitMesh( na, cell, norb, iaorb, iphorb, isa, rmax, 
     &               G2max, G2mesh, nsc, nmpl, nm,
     &               nml, ntm, ntml, ntpl, dvol )

!     Setup box descriptors for each processor,
!     held in module iogrid_netcdf
      call set_box_limits( ntm, nsm )

      ! Initialize information on local mesh for each node
      call init_mesh_node( cell, ntm , meshLim , nsm )

      ! Setup charge additions in the mesh
      call init_charge_add(cell, ntm)

      ! Setup Hartree additions in the mesh
      call init_hartree_add(cell, ntm)

#ifdef NCDF_4
      ! Initialize the box for each node...
      call cdf_init_mesh( ntm, nsm )
      if ( write_cdf ) then
        call cdf_init_grid( trim(slabel)//'.nc', ntm)
      end if
#endif

C     Stop time counter for mesh initialization
      call timer( 'DHSCF1', 2 )
C ----------------------------------------------------------------------
C     End of mesh initialization 
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C     Initialize atomic orbitals, density and potential
C ----------------------------------------------------------------------
C     Start time counter for atomic initializations
      call timer( 'DHSCF2', 1 )

      if (nodes.gt.1) then
        call setMeshDistr( UNIFORM, nsm, nsp,
     &                     nml, nmpl, ntml, ntpl )
      endif

C     Initialise quantities relating to the atom-mesh positioning
      call InitAtomMesh( UNIFORM, na, xa )

!     Check if we need stencils in cellxc, and detect incompatibility
!     with Harris forces 
      call xc_setup(need_gradients_in_xc, is_vdw)

      if (need_gradients_in_xc .and. harrisfun) then
         call die("** Harris forces not implemented for non-LDA XC")
      endif

      ! Give the opportunity to use BSC's version,
      ! unless we have a vdw functional
      use_bsc_cellxc = fdf_get("XC.Use.BSC.CellXC", .false.)
      if (is_vdw .and. use_bsc_cellxc) then
         call message("INFO",
     $        "BSC's cellxc cannot be used with vdW functionals")
         use_bsc_cellxc = .false.
      endif
      
C     Compute the number of orbitals on the mesh and recompute the
C     partions for every processor in order to have a similar load
C     in each of them.
      nullify( numphi )
      call re_alloc( numphi, 1, nmpl, 'numphi', 'dhscf_init' )
!$OMP parallel do default(shared), private(i)
      do i= 1, nmpl
        numphi(i) = 0
      enddo
!$OMP end parallel do

      call distriPhiOnMesh( nm, nmpl, norb, iaorb, iphorb,
     &                      isa, numphi, need_gradients_in_xc )

C     Find if there are partial-core-corrections for any atom
      npcc = 0
      do ia = 1,na
        if (rcore(isa(ia)) .gt. tiny) npcc = 1
      enddo

C     Find partial-core-correction energy density
C     Vscf and Vaux are not used here
      call re_alloc( rhopcc, 1, ntpl*npcc+1, 'rhopcc', 'dhscf_init' )
      if (npcc .eq. 1) then
        call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     &    nsd, dvol, volume, dummy_Vscf, dummy_Vaux, Fal, stressl,
     &    .false., .false. )
        call reord( rhopcc, rhopcc, nml, nsm, TO_SEQUENTIAL )
        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'rhopcc',sum(rhopcc)*dvol
        end if
      endif

C     Find neutral-atom potential
C     Drho is not used here
      call re_alloc( Vna, 1, ntpl, 'Vna', 'dhscf_init' )
      call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                        volume, dummy_DRho, Fal, stressl, 
     &                        .false., .false. )
      call reord( Vna, Vna, nml, nsm, TO_SEQUENTIAL )
      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'Vna',sqrt(sum(Vna**2))
      end if

      if (nodes.gt.1) then
        if (node .eq. 0) then
           write(6,"(a)") "Setting up quadratic distribution..."
        endif
        call setMeshDistr( QUADRATIC, nsm, nsp,
     &                     nml, nmpl, ntml, ntpl )

C       Create extended mesh arrays for the second data distribution
        call setupExtMesh( QUADRATIC, rmax )

C       Compute atom positions for the second data distribution 
        call InitAtomMesh( QUADRATIC, na, xa )
      endif

C     Calculate orbital values on mesh
!     numphi has already been computed in distriPhiOnMesh
!     in the UNIFORM distribution 
      if (nodes.eq.1) then
        numphi_par => numphi
      else
        nullify(numphi_par)
        call re_alloc( numphi_par, 1, nmpl, 'numphi_par',
     &                 'dhscf_init' )
        call distMeshData( UNIFORM, numphi, QUADRATIC,
     &                     numphi_par, KEEP )
      endif

      call PhiOnMesh( nmpl, norb, iaorb, iphorb, isa, numphi_par )

      if (nodes.gt.1) then
        call de_alloc( numphi_par, 'numphi_par', 'dhscf_init' )
      endif
      call de_alloc( numphi, 'numphi', 'dhscf_init' )

C ----------------------------------------------------------------------
C       Create sparse indexing for Dscf as needed for local mesh
C       Note that this is done in the QUADRATIC distribution
C       since 'endpht' (computed finally in PhiOnMesh and stored in
C       meshphi module) is in that distribution.
C ----------------------------------------------------------------------
      if (Nodes.gt.1) then
        call CreateLocalDscfPointers( nmpl, nuotot, numd, listdptr, 
     &                                listd )
      endif

C ----------------------------------------------------------------------
C     Calculate terms relating to the neutral atoms on the mesh
C ----------------------------------------------------------------------
C     Find Harris (sum of atomic) electron density
      call re_alloc( rhoatm_par, 1, ntpl, 'rhoatm_par', 'dhscf_init' )
      call rhooda( norb, nmpl, datm, rhoatm_par, iaorb, iphorb, isa )
!     rhoatm_par comes out of here in clustered form in QUADRATIC dist

C     Routine Poison should use the uniform data distribution
      if (nodes.gt.1) then
         call setMeshDistr( UNIFORM, nsm, nsp,
     &                      nml, nmpl, ntml, ntpl )
      endif

C     Create Rhoatm using UNIFORM distr, in sequential form
      call re_alloc( rhoatm, 1, ntpl, 'rhoatm', 'dhscf_init' )
      call distMeshData( QUADRATIC, rhoatm_par,
     &     UNIFORM, rhoatm, TO_SEQUENTIAL )

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'rhoatm', sum(rhoatm)*dvol
      end if

!
!  AG: The initialization of doping structs could be done here now,
!      in the uniform distribution, and with a simple loop over
!      rhoatm.

        if (frstme) call initialize_doping_uniform()
        if (doping_active) then  
           call compute_doping_structs_uniform(ntpl,rhoatm,nsd)
           ! Will get the global number of hit points
           ! Then, the doping density to be added can be simply computed
        endif

C     Allocate Temporal array
      call re_alloc( Vscf, 1, ntpl, 1, 1, 'Vscf', 'dhscf_init' )

C     Vscf is filled here but not used later
C     Uharrs is computed (and saved)
C     DStres is computed but not used later
      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, rhoatm,
     &             Uharrs, Vscf, DStres, nsm )
      call de_alloc( Vscf, 'Vscf', 'dhscf_init' )

!     Always deallocate rhoatm_par, as it was used even if nodes=1
      call de_alloc( rhoatm_par, 'rhoatm_par', 'dhscf_init' )

      if (mix_charge) then
         call re_alloc( rhog, 1, 2, 1, ntpl, 1, nspin,
     $                 'Rhog', 'dhscf_init' )
         call re_alloc( rhog_in, 1, 2, 1, ntpl, 1, nspin,
     $                 'Rhog_in', 'dhscf_init' )
         call order_rhog(cell, ntml(1), ntml(2), ntml(3), ntm, nsm)
      endif

C     Stop time counter for atomic initializations
      call timer( 'DHSCF2', 2 )

C ----------------------------------------------------------------------
C     At the end of initializations:
!     We are in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna are in UNIFORM dist, sequential form
!     The index array endpht is in the QUADRATIC distribution
C ----------------------------------------------------------------------
      if (frstme) then
        call initialize_efield()
        call init_dipole_correction(nbcell, bcell, acting_efield)
      end if

      ! Check if we need to add the potential 
      ! corresponding to the voltage-drop.
      if ( TSmode ) then
        ! These routines are important if there are cell-changes
        if ( IsVolt ) then
           call ts_init_voltage( cell, nua, xa, ntm)
        end if

        if ( acting_efield ) then
         ! We do not allow the electric field for 
         ! transiesta runs with V = 0, either.
         ! It does not make sense, only for fields perpendicular
         ! to the applied bias.
           
         ! We need to check that the e-field is perpendicular
         ! to the transport direction, and that the system is
         ! either a chain, or a slab.
         ! However, due to the allowance of a dipole correction
         ! along the transport direction for buffer calculations
         ! we have to allow all shapes. (atom is not transiesta
         ! compatible anyway)
            
         ! check that we do not violate the periodicity
         if ( Node .eq. 0 ) then
            write(*,'(/,2(2a,/))') 'ts-WARNING: ',
     &           'E-field/dipole-correction! ',
     &           'ts-WARNING: ',
     &           'I hope you know what you are doing!'
         end if

         ! This is either dipole or user, or both
         is_ortho = .true.
         if ( N_Elec == 1 ) then
           ! When using a single electrode
           ! one may always have an electric field.
           ! For fields along the semi-infinite direction
           ! it corresponds to a capacitor setup
           ! For other fields it is simply a gate-like field.
           ! However, when the electrode uses real space self-energies
           ! then we do not allow *any* fields along the semi-inf
           ! directions since that would mean an infinite field.
           select case ( Elecs(1)%t_dir )
           case ( 4 , 5 , 6, 7 )
             is_ortho = Elec_ortho2semi_inf(Elecs(1),
     &           user_specified_field)
           end select
         else
           ! For more than 1 electrode we do not allow any electric
           ! fields along the semi-infinite directions.
           do iE = 1 , N_Elec
             is_ortho = is_ortho .and.
     &           Elec_ortho2semi_inf(Elecs(iE), user_specified_field)
           end do
         end if
         if ( .not. is_ortho ) then
           call die('User defined E-field must be
     &perpendicular to semi-infinite directions.')
         end if
        end if ! acting_efield
         
        ! We know that we currently allow people to do more than
        ! they probably should be allowed. However, there are many
        ! corner cases that may require dipole corrections, or 
        ! electric fields to "correct" an intrinsic dipole.
        ! For instance, what should we do with a dipole in a transiesta
        ! calculation?
        ! Should we apply a field to counter act it in a device
        ! calculation?
         
      end if

      frstme = .false.

      call timer( 'DHSCF_Init', 2 )
#ifdef DEBUG
      call write_debug( '    POS dhscf_init' )
#endif
!------------------------------------------------------------------------- END
      end subroutine dhscf_init

      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo, 
     &                  nuotot, nua, na, isa, xa, indxua,
     &                  ntm, ifa, istr, iHmat,
     &                  filesOut, maxnd, numd,
     &                  listdptr, listd, Dscf, datm, maxnh, Hmat,
     &                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, 
     &                  Exc, Dxc, dipol, stress, Fal, stressl,
     &                  use_rhog_in, charge_density_only )

C
C Calculates the self-consistent field contributions to Hamiltonian
C matrix elements, total energy and atomic forces.
C Coded by J.M. Soler, August 1996. July 1997.
C Modified by J.D. Gale, February 2000.
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nspin         : Number of different spin polarisations
C                         nspin=1 => Unpolarized, nspin=2 => polarized
C                         nspin=4 => Noncollinear spin or spin-orbit
C integer norb          : Total number of basis orbitals in supercell
C integer iaorb(norb)   : Atom to which each orbital belongs
C integer iphorb(norb)  : Orbital index (within atom) of each orbital
C integer nuo           : Number of orbitals in a unit cell in this node
C integer nuotot        : Number of orbitals in a unit cell
C integer nua           : Number of atoms in unit cell
C integer na            : Number of atoms in supercell
C integer isa(na)       : Species index of all atoms in supercell
C real*8  xa(3,na)      : Atomic positions of all atoms in supercell
C integer indxua        : Index of equivalent atom in unit cell
C integer ifa           : Switch which fixes whether the SCF contrib.
C                         to atomic forces is calculated and added to fa
C integer istr          : Switch which fixes whether the SCF contrib.
C                         to stress is calculated and added to stress
C integer iHmat         : Switch which fixes whether the Hmat matrix
C                         elements are calculated or not.
C type(filesOut_t) filesOut : output file names (If blank => not saved)
C integer maxnd             : First dimension of listd and Dscf
C integer numd(nuo)         : Number of nonzero density-matrix
C                             elements for each matrix row
C integer listdptr(nuo)     : Pointer to start of rows of density-matrix
C integer listd(maxnd)      : Nonzero-density-matrix-element column
C                             indexes for each matrix row
C real*8  Dscf(maxnd,h_spin_dim): SCF density-matrix elements
C real*8  datm(norb)        : Harris density-matrix diagonal elements
C                             (atomic occupation charges of orbitals)
C integer maxnh             : First dimension of listh and Hmat
C real*8  Hmat(maxnh,h_spin_dim) : Hamiltonian matrix in sparse form,
C                             to which are added the matrix elements
C                                 <ORB_I | DeltaV | ORB_J>, where
C                             DeltaV = Vna + Vxc(SCF) + 
C                                      Vhartree(RhoSCF-RhoHarris)
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C integer ntm(3) : Number of mesh divisions of each cell
C                  vector, including subgrid.
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C real*8  Enaatm : Integral of Vna * rhoatm
C real*8  Enascf : Integral of Vna * rhoscf
C real*8  Uatm   : Harris hartree electron-interaction energy
C real*8  Uscf   : SCF hartree electron-interaction energy
C real*8  DUscf  : Electrostatic (Hartree) energy of
C                    (rhoscf - rhoatm) density
C real*8  DUext  : Interaction energy with external electric field
C real*8  Exc    : SCF exchange-correlation energy
C real*8  Dxc    : SCF double-counting correction to Exc
C                    Dxc = integral of ( (epsxc - Vxc) * Rho )
C                    All energies in Rydbergs
C real*8  dipol(3): Electric dipole (in a.u.)
C                   only when the system is a molecule/slab
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C real*8  Fal(3,nua) : Atomic forces, to which the SCF contribution
C                        is added by this routine when ifa=1.
C                        the SCF contribution is minus the derivative
C                        of ( Enascf - Enaatm + DUscf + Exc ) with
C                        respect to atomic positions, in  Ry/Bohr.
C                        contributions local to this node
C real*8 stressl(3,3): Stress tensor, to which the SCF contribution
C                      is added by this routine when ifa=1.
C                      the SCF contribution is minus the derivative of
C                         ( Enascf - Enaatm + DUscf + Exc ) / volume
C                      with respect to the strain tensor, in Ry.
C                        contributions local to this node
C ----------------------------------------------------------------------
C Units :
C ----------------------------------------------------------------------
C Energies in Rydbergs
C Distances in Bohr
C ----------------------------------------------------------------------
C     Modules
      use precision,     only  : dp, grid_p

      use parallel,      only  : Node, Nodes
      use atmfuncs,      only  : rcut, rcore
      use units,         only  : Debye, eV, Ang
      use fdf
      use sys,           only  : die, bye
      use mesh,          only  : nsm, nsp, meshLim
      use parsing
      use m_iorho,       only  : write_rho
      use m_forhar,      only  : forhar
      use alloc,         only  : re_alloc, de_alloc
      use files,         only  : slabel
      use files,         only  : filesOut_t ! derived type for output file names
      use siesta_options, only : harrisfun, save_initial_charge_density
      use siesta_options, only : analyze_charge_density_only
      use meshsubs,       only : LocalChargeOnMesh
      use meshsubs,       only : PartialCoreOnMesh
      use meshsubs,       only : NeutralAtomOnMesh

      use moreMeshSubs,   only : setMeshDistr, distMeshData
      use moreMeshSubs,   only : UNIFORM, QUADRATIC, LINEAR
      use moreMeshSubs,   only : TO_SEQUENTIAL, TO_CLUSTER, KEEP
      use m_partial_charges, only: compute_partial_charges
      use m_partial_charges, only: want_partial_charges

      use siestaXC, only : cellXC ! Finds xc energy and potential
      use siestaXC, only : myMeshBox    ! Returns my processor mesh box
      use siestaXC, only : jms_setMeshDistr => setMeshDistr
                                        ! Sets a distribution of mesh
                                        ! points over parallel processors
      use m_vmat, only: vmat
      use m_vmatsp, only: vmatsp
      use m_rhoofd, only: rhoofd
      use m_rhoofdsp, only: rhoofdsp
#ifdef MPI
      use mpi_siesta
#endif
      use iogrid_netcdf, only: write_grid_netcdf
      use iogrid_netcdf, only: read_grid_netcdf
      use siesta_options, only: read_charge_cdf
      use siesta_options, only: savebader
      use siesta_options, only: read_deformation_charge_cdf
      use siesta_options, only: mix_charge

      use m_efield,       only: add_potential_from_field
      use m_efield,       only: user_specified_field, acting_efield
      use dipole_m, only: get_dipole_origin
      use dipole_m, only: get_field_from_dipole
      use dipole_m, only: dipole_charge, dipole_potential
      use dipole_m, only: dip_vacuum
      use dipole_m, only: SLABDIPOLECORR
      use dipole_m, only: SLABDIPOLECORR_NONE
      use dipole_m, only: SLABDIPOLECORR_CHARGE, SLABDIPOLECORR_VACUUM

      use m_doping_uniform, only: doping_active, doping_uniform

      use m_charge_add,   only: charge_add
      use m_hartree_add,  only: hartree_add

#ifdef NCDF_4
      use siesta_options, only: write_cdf
      use m_ncdf_siesta, only: cdf_save_grid
#endif
      use m_rhofft,       only: rhofft, FORWARD, BACKWARD
      use m_rhog,         only: rhog_in, rhog
      use m_spin,         only: spin
      use m_spin,         only: Spiral, qSpiral

      use m_ts_global_vars,only: TSmode, TSrun
      use m_ts_options, only: IsVolt, Elecs, N_elec
      use m_ts_voltage, only: ts_voltage
      use m_ts_hartree, only: ts_hartree_fix

      implicit none

      integer
     &  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot, 
     &  iaorb(norb), ifa, iHmat,
     &  indxua(na), iphorb(norb), isa(na), 
     &  istr, listd(*), listdptr(nuo), ntm(3), numd(nuo)

      real(dp)
     &  datm(norb), dipol(3), Dscf(:,:),
     &  DUext, DUscf, Dxc, Enaatm, Enascf, Exc,
     &  Hmat(:,:), Uatm, Uscf, xa(3,na),
     &  stress(3,3), Fal(3,nua), stressl(3,3)

      type(filesOut_t) filesOut

      logical, intent(in), optional :: use_rhog_in
      logical, intent(in), optional :: charge_density_only

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        cellxc(...)    : Finds total exch-corr energy and potential
C        cross(a,b,c)   : Finds the cross product of two vectors
C        dfscf(...)     : Finds SCF contribution to atomic forces
C        dipole_charge  : Finds electric dipole moment using the charges
C        doping(...)    : Adds a background charge for doped systems
C        write_rho(...)     : Saves electron density on a file
C        poison(...)    : Solves Poisson equation
C        reord(...)     : Reorders electron density and potential arrays
C        rhooda(...)    : Finds Harris electron density in the mesh
C        rhoofd(...)    : Finds SCF electron density in the mesh
C        rhoofdsp(...)  : Finds SCF electron density in the mesh for
C                         spiral arrangement of spins
C        timer(...)     : Finds CPU times
C        vmat(...)      : Finds matrix elements of SCF potential
C        vmatsp(...)    : Finds matrix elements of SCF potential for
C                         spiral arrangement of spins
C        delk(...)      : Finds matrix elements of exp(i \vec{k} \cdot \vec{r})
C real*8 volcel( cell ) : Returns volume of unit cell
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C integer nbcell        : Number of bulk lattice vectors
C real*8  bcell(3,3)    : Bulk lattice vectors
C real*8  cell(3,3)     : Auxiliary lattice vectors (same as ucell)
C real*8  const         : Auxiliary variable (constant within a loop)
C real*8  DEc           : Auxiliary variable to call cellxc
C real*8  DEx           : Auxiliary variable to call cellxc
C real*8  dvol          : Mesh-cell volume
C real*8  Ec            : Correlation energy
C real*8  Ex            : Exchange energy
C real*8  field(3)      : External electric field
C integer i             : General-purpose index
C integer ia            : Atom index
C integer io            : Orbital index
C integer ip            : Point index
C integer is            : Species index
C logical IsDiag        : Is supercell diagonal?
C integer ispin         : Spin index
C integer j             : General-purpose index

C integer myBox(2,3)    : My processor's mesh box 

C integer nbcell        : Number of independent bulk lattice vectors
C integer npcc          : Partial core corrections? (0=no, 1=yes)
C integer nsd           : Number of diagonal spin values (1 or 2)
C integer ntpl          : Number of mesh Total Points in unit cell
C                         (including subpoints) locally
C real*4  rhoatm(ntpl)  : Harris electron density
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C real*4  DRho(ntpl)    : Selfconsistent electron density difference
C real*8  rhotot        : Total density at one point
C real*8  rmax          : Maximum orbital radius
C real*8  scell(3,3)    : Supercell vectors
C real*4  Vaux(ntpl)    : Auxiliary potential array
C real*4  Vna(ntpl)     : Sum of neutral-atom potentials
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of selfconsistent density
C real*8  x0(3)         : Center of molecule
C logical harrisfun     : Harris functional or Kohn-Sham?
C ----------------------------------------------------------------------

C     Local variables
      integer  :: i, ia, ip, ispin, nsd, np_vac
      real(dp) :: b1Xb2(3), const, DEc, DEx, DStres(3,3),
     &            Ec, Ex, rhotot, x0(3), volume, Vmax_vac, Vmean_vac

      logical :: use_rhog
      
      real(dp), external :: volcel, ddot

      external
     &  cross,
     &  poison, 
     &  reord, rhooda,
     &  timer

!     Dummy arrays for bsc_cellxc call
      real(grid_p) :: aux3(3,1)
      real(grid_p) :: dummy_DVxcdn(1,1,1)
      real(grid_p), pointer :: Vscf_gga(:,:), DRho_gga(:,:)
      external bsc_cellxc

!     Interface to JMS's SiestaXC
      integer       :: myBox(2,3)
      real(dp) :: stressXC(3,3)

C     Work arrays
      real(grid_p), pointer :: Vscf(:,:), Vscf_par(:,:),
     &                         DRho(:,:), DRho_par(:,:),
     &                         Vaux(:), Vaux_par(:), Chlocal(:),
     &                         Totchar(:), fsrc(:), fdst(:),
     &                         rhoatm_quad(:) => null(),
     &                         DRho_quad(:,:) => null()
      ! Temporary reciprocal spin quantity
      real(grid_p) :: rnsd

#ifdef MPI
      integer  :: MPIerror
      real(dp) :: sbuffer(7), rbuffer(7)
#endif

#ifdef DEBUG
      call write_debug( '    PRE DHSCF' )
#endif

      if ( spin%H /= size(Dscf,dim=2) ) then
         call die('Spin components is not equal to options.')
      end if

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'DM',
     &        (sqrt(sum(Dscf(:,ispin)**2)),ispin=1,spin%H)
         write(*,debug_fmt) Node,'H',
     &        (sqrt(sum(Hmat(:,ispin)**2)),ispin=1,spin%H)
      end if
      
!-------------------------------------------------------------------- BEGIN
C ----------------------------------------------------------------------
C Start of SCF iteration part
C ----------------------------------------------------------------------

C ----------------------------------------------------------------------
C     At the end of DHSCF_INIT, and also at the end of any previous
!     call to dhscf, we were in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna were in UNIFORM dist, sequential form
!     The index array endpht was in the QUADRATIC distribution
C ----------------------------------------------------------------------

#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_restart( )
#endif
      call timer( 'DHSCF', 1 )
      call timer( 'DHSCF3', 1 )

      nullify( Vscf, Vscf_par, DRho, DRho_par,
     &         Vaux, Vaux_par, Chlocal, Totchar )
      nullify( Vscf_gga, DRho_gga)

      volume = volcel(cell)

C-------------------------------------------------------------------------
      
      if (analyze_charge_density_only) then
         !! Use the functionality in the first block
         !! of the routine to get charge files and partial charges
         call setup_analysis_options()
      endif
      
      if (filesOut%vna .ne. ' ') then
        ! Uniform dist, sequential form
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vna',1,ntml,Vna)
        else
           call write_rho( filesOut%vna, 
     &          cell, ntm, nsm, ntpl, 1, Vna)
           call write_grid_netcdf( cell, ntm, 1, ntpl, Vna, "Vna")
        end if
#else
        call write_rho( filesOut%vna,
     &       cell, ntm, nsm, ntpl, 1, Vna )
        call write_grid_netcdf( cell, ntm, 1, ntpl, Vna, "Vna")
#endif
      endif

C     Allocate memory for DRho using the UNIFORM data distribution
      call re_alloc( DRho, 1, ntpl, 1, nspin, 'DRho','dhscf' )

C Find number of diagonal spin values
      nsd  = min( nspin, 2 )
      if ( nsd == 1 ) then
         rnsd = 1._grid_p
      else
         rnsd = 1._grid_p / nsd
      end if

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------
!
!     The reading routine works in the uniform distribution, in
!     sequential form
!
      if (present(use_rhog_in)) then
            use_rhog = use_rhog_in
         else
            use_rhog = .false.
      endif
      if (use_rhog) then
         ! fourier transform back into drho
         call rhofft(cell, ntml(1), ntml(2), ntml(3), ntm, nspin,
     $                  DRho,rhog_in,BACKWARD)

      else if (read_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"Rho")
         read_charge_cdf = .false.
      else if (read_deformation_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"DeltaRho")
         ! Add to diagonal components only
         do ispin = 1,nsd
            do ip= 1, ntpl
C             rhoatm and Drho are in sequential mode
              DRho(ip,ispin) = DRho(ip,ispin) + rhoatm(ip) * rnsd
            enddo
         enddo
         read_deformation_charge_cdf = .false.
      else
        ! Set the QUADRATIC distribution and allocate memory for DRho_par
        ! since the construction of the density from the DM and orbital
        ! data needs that distribution
        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
        call re_alloc( DRho_par, 1, ntpl, 1, nspin,
     &                 'DRho_par','dhscf' )       
        if (Spiral) then
          call rhoofdsp( norb, nml, nmpl, maxnd, numd, listdptr, listd,
     &                   nspin, Dscf, DRho_par, nuo, nuotot, iaorb,
     &                   iphorb, isa, qspiral )
        else
          call rhoofd( norb, nmpl, maxnd, numd, listdptr, listd,
     &                 nspin, Dscf, DRho_par, 
     &                 nuo, nuotot, iaorb, iphorb, isa )
        endif
        ! DRHO_par is here in QUADRATIC, clustered form

C       Set the UNIFORM distribution again and copy DRho to it
        if (nodes.gt.1) then
           call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
        endif

        do ispin = 1, nspin
          fsrc => DRho_par(:,ispin)
          fdst => DRho(:,ispin)
         ! Sequential to be able to write it out
         ! if nodes==1, this call will just reorder
          call distMeshData( QUADRATIC, fsrc,
     &                       UNIFORM, fdst, TO_SEQUENTIAL )
        enddo
        call de_alloc( DRho_par, 'DRho_par','dhscf' )

        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'Rho',
     &          (sum(DRho(:,ispin))*dvol,ispin=1,nspin)
        end if

        if (save_initial_charge_density) then
           ! This section is to be deprecated in favor
           ! of "analyze_charge_density_only"
           ! (except for the special name for the .nc file)
#ifdef NCDF_4
          if ( write_cdf ) then
             call cdf_save_grid(trim(slabel)//'.nc','RhoInit',nspin,
     &            ntml,DRho)
          else
             call write_rho( "RHO_INIT", cell, ntm, nsm, ntpl,
     $            nspin, DRho)
             call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho,
     $            "RhoInit")
          end if
#else
          call write_rho( "RHO_INIT", cell, ntm, nsm, ntpl,
     $         nspin, DRho)
          call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho,
     $         "RhoInit")
#endif
          call timer('DHSCF3',2)
          call timer('DHSCF',2)
          call bye("STOP after producing RHO_INIT from input DM")
        endif

      endif

      if (mix_charge) then
         ! Save fourier transform of charge density
         call rhofft(cell, ntml(1), ntml(2), ntml(3), ntm, nspin,
     $        DRho,rhog,FORWARD)
      endif
!
!     Proper place to integrate Hirshfeld and Voronoi code,
!     since we have just computed rhoatm and Rho.

      if (want_partial_charges) then
        ! The endpht array is in the quadratic distribution, so
        ! we need to use it for this...
        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
        call re_alloc( DRho_quad, 1, ntpl, 1, nspin,
     &                 'DRho_quad','dhscf' )       
        call re_alloc( rhoatm_quad, 1, ntpl,
     &                 'rhoatm_quad','dhscf' )
        ! Redistribute grid-density
        do ispin = 1, nspin
          fsrc => DRho(:,ispin)
          fdst => DRho_quad(:,ispin)
         ! if nodes==1, this call will just reorder
          call distMeshData( UNIFORM, fsrc,
     &                       QUADRATIC, fdst, TO_CLUSTER )
        enddo
        call distMeshData( UNIFORM, rhoatm,
     &                     QUADRATIC, rhoatm_quad, TO_CLUSTER )

        call compute_partial_charges(DRho_quad,rhoatm_quad,
     .                  nspin, iaorb, iphorb, 
     .                  isa, nmpl,dvol)

        call de_alloc(rhoatm_quad,'rhoatm_quad','dhscf')
        call de_alloc(Drho_quad,'DRho_quad','dhscf')
        if (nodes.gt.1) then
           call setMeshDistr( UNIFORM, nsm, nsp,
     &                        nml, nmpl, ntml, ntpl )
        endif
      endif

C ----------------------------------------------------------------------
C Save electron density
C ----------------------------------------------------------------------
      if (filesOut%rho .ne. ' ') then
        !  DRho is already using a uniform, sequential form
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Rho',nspin,ntml,
     &          DRho)
        else
           call write_rho( filesOut%rho, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "Rho")
        end if
#else
        call write_rho( filesOut%rho, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "Rho")
#endif
      endif
C ----------------------------------------------------------------------
C Save the diffuse ionic charge and/or the total (ionic+electronic) charge
C ----------------------------------------------------------------------
      if (filesOut%psch .ne. ' ' .or. filesOut%toch .ne. ' ') then
C       Find diffuse ionic charge on mesh
        ! Note that the *OnMesh routines, except PhiOnMesh,
        ! work with any distribution, thanks to the fact that
        ! the ipa, idop, and indexp arrays are distro-specific
        call re_alloc( Chlocal, 1, ntpl, 'Chlocal', 'dhscf' )
        call LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
        ! Chlocal comes out in clustered form, so we convert it
        call reord( Chlocal, Chlocal, nml, nsm, TO_SEQUENTIAL )

        if ( debug_dhscf ) then
           write(*,debug_fmt) Node,'Chlocal',sum(Chlocal)*dvol
        end if

C       Save diffuse ionic charge 
        if (filesOut%psch .ne. ' ') then
#ifdef NCDF_4
          if ( write_cdf ) then
             call cdf_save_grid(trim(slabel)//'.nc','Chlocal',1,ntml,
     &            Chlocal)
          else
             call write_rho( filesOut%psch, cell, ntm, nsm, ntpl, 1,
     &            Chlocal)
             call write_grid_netcdf( cell, ntm, 1, ntpl, Chlocal,
     &            'Chlocal' )
          end if
#else
          call write_rho( filesOut%psch, cell, ntm, nsm, ntpl, 1,
     &         Chlocal)
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     &         Chlocal, 'Chlocal' )
#endif
        endif

C       Save total (ionic+electronic) charge 
        if ( filesOut%toch .ne. ' ') then
           ! *****************
           ! **  IMPORTANT  **
           ! The Chlocal array is re-used to minimize memory
           ! usage. In the this small snippet the Chlocal
           ! array will contain the total charge, and
           ! if the logic should change, (i.e. should Chlocal
           ! be retained) is the Totchar needed to be re-instantiated.
           ! *****************

!$OMP parallel default(shared), private(ispin,ip)
           do ispin = 1, nsd
!$OMP do 
           do ip = 1, ntpl
              Chlocal(ip) = Chlocal(ip) + DRho(ip,ispin)
           end do
!$OMP end do 
           end do
!$OMP end parallel

          ! See note above
#ifdef NCDF_4
           if ( write_cdf ) then
              call cdf_save_grid(trim(slabel)//'.nc','RhoTot',1,ntml,
     &             Chlocal)
           else
              call write_rho(filesOut%toch,cell,ntm,nsm,ntpl,1,Chlocal)
              call write_grid_netcdf( cell, ntm, 1, ntpl, Chlocal, 
     &             "TotalCharge")
           end if
#else
           call write_rho( filesOut%toch, cell, ntm, nsm, ntpl, 1,
     &          Chlocal )
           call write_grid_netcdf( cell, ntm, 1, ntpl,
     &          Chlocal, "TotalCharge")
#endif
        end if 
        call de_alloc( Chlocal, 'Chlocal', 'dhscf' )
      endif

C ----------------------------------------------------------------------
C Save the total charge (model core + valence) for Bader analysis
C ----------------------------------------------------------------------
      
      ! The test for toch guarantees that we are in "analysis mode"
      if (filesOut%toch .ne. ' ' .and. savebader) then
        call save_bader_charge()
      endif

C Find difference between selfconsistent and atomic densities

      !Both DRho and rhoatm are using a UNIFORM, sequential form
!$OMP parallel default(shared), private(ispin,ip)
      do ispin = 1,nsd
!$OMP do
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) * rnsd
        enddo
!$OMP end do nowait
      enddo
!$OMP end parallel

C ----------------------------------------------------------------------
C Save electron density difference
C ----------------------------------------------------------------------
      if (filesOut%drho .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','RhoDelta',nspin,ntml,
     &          DRho)
        else
           call write_rho( filesOut%drho, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &          "DeltaRho")
        end if
#else
        call write_rho( filesOut%drho, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     &       DRho, "DeltaRho")
#endif
      endif

      if (present(charge_density_only)) then 
         if (charge_density_only) then
            call timer('DHSCF3',2)
            call timer('DHSCF',2)
            call de_alloc( DRho, 'DRho', 'dhscf' )
            RETURN
         endif
      endif

! End of analysis section
! Can exit now, if requested
      
        if (analyze_charge_density_only) then 
            call timer('DHSCF3',2)
            call timer('DHSCF',2)
           call bye("STOP after analyzing charge from input DM")
        endif

!-------------------------------------------------------------
C     Transform spin density into sum and difference

      ! TODO Check for NC/SO:
      ! Should we diagonalize locally at every point first?
      if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(rhotot,ip)
        do ip = 1,ntpl
          rhotot     = DRho(ip,1) + DRho(ip,2)
          DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
          DRho(ip,1) = rhotot
        enddo
!$OMP end parallel do
      endif

      if ( debug_dhscf ) then
        write(*,debug_fmt) Node,'DRho-clean',
     &      (sum(DRho(:,ispin))*dvol,ispin=1,nspin)
      end if

C Add a background charge to neutralize the net charge, to
C model doped systems. It only adds the charge at points
C where there are atoms (i.e., not in vacuum).          
C First, call with 'task=0' to add background charge    
      if (doping_active) call doping_uniform(cell,ntpl,0,
     $     DRho(:,1),rhoatm)
      
      ! Add doping in cell (from ChargeGeometries/Geometry.Charge)
      ! Note that this routine will return immediately if no dopant is present
      call charge_add('+',cell, ntpl, DRho(:,1) )

      if ( debug_dhscf ) then
        write(*,debug_fmt) Node,'DRho-doped',
     &      (sum(DRho(:,ispin))*dvol,ispin=1,nspin)
      end if

C ----------------------------------------------------------------------
C Calculate the dipole moment
C ----------------------------------------------------------------------
      dipol(1:3) = 0.0_dp
      if ( nbcell /= 3 ) then
        ! Find dipole in system

        select case ( SLABDIPOLECORR )
        case ( SLABDIPOLECORR_VACUUM )
          ! we should not calculate the dipole here
          ! User requests the vacuum dipole calculation
        case default
          ! Regardless of whether the user requests a
          ! dipole correction or not, then we will calculate the dipole
          ! using the charge method.

          ! Get dipole origin
          call get_dipole_origin(cell, nua, xa, x0)

          ! This routine is distribution-blind
          ! and will reduce over all processors.
          call dipole_charge(cell, dvol, ntm, ntml, nsm, DRho,
     &        x0, dipol)

          ! Orthogonalize dipole to bulk directions
          if ( nbcell == 1 ) then ! chain
            const = ddot(3,dipol,1,bcell,1) / ddot(3,bcell,1,bcell,1)
            dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
          else if ( nbcell == 2 ) then ! slab
            call cross( bcell(1,1), bcell(1,2), b1Xb2 )
            const = ddot(3,dipol,1,b1Xb2,1) / ddot(3,b1Xb2,1,b1Xb2,1)
            dipol(1:3) = const * b1Xb2(1:3)
          end if

        end select

      endif

C ----------------------------------------------------------------------
C     Find Hartree potential of DRho = rhoscf-rhoatm. Store it in Vaux
C ----------------------------------------------------------------------
!     Solve Poisson's equation
      call re_alloc( Vaux, 1, ntpl, 'Vaux', 'dhscf' )

      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, DRho,
     &             DUscf, Vaux, DStres, nsm )

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'Poisson',sqrt(sum(Vaux(:)**2))
      end if

      ! Vscf is in the UNIFORM, sequential form, and only using
      ! the first spin index

      ! We require that even the SIESTA potential is "fixed"
      ! NOTE, this will only do something if
      !   TS.Hartree.Fix is set
      call ts_hartree_fix( ntm, ntml, Vaux)
      
C Add contribution to stress from electrostatic energy of rhoscf-rhoatm
      if (istr .eq. 1) then
        stressl(1:3,1:3) = stressl(1:3,1:3) + DStres(1:3,1:3)
      endif

C ----------------------------------------------------------------------
C     Find electrostatic (Hartree) energy of full SCF electron density
C     using the original data distribution 
C ----------------------------------------------------------------------
      Uatm = Uharrs
      Uscf = 0._dp
!$OMP parallel do default(shared), private(ip), 
!$OMP&reduction(+:Uscf)
      do ip = 1, ntpl
        Uscf = Uscf + Vaux(ip) * rhoatm(ip)
      enddo
!$OMP end parallel do
      Uscf = Uscf * dVol + Uatm + DUscf

C Call doping with 'task=1' to remove background charge added previously
C The extra charge thus only affects the Hartree energy and potential,
C but not the contribution to Enascf ( = \Int_{Vna*\rho})
      if (doping_active) call doping_uniform(cell,ntpl,1,
     $     DRho(:,1),rhoatm)

      ! Remove doping in cell (from ChargeGeometries/Geometry.Charge)
      ! Note that this routine will return immediately if no dopant is present
      call charge_add('-',cell, ntpl, DRho(:,1) )

C ----------------------------------------------------------------------
C Add neutral-atom potential to Vaux
C ----------------------------------------------------------------------
      Enaatm = 0.0_dp
      Enascf = 0.0_dp
!$OMP parallel do default(shared), private(ip),
!$OMP&reduction(+:Enaatm,Enascf)
      do ip = 1, ntpl
        Enaatm   = Enaatm + Vna(ip) * rhoatm(ip)
        Enascf   = Enascf + Vna(ip) * DRho(ip,1)
        Vaux(ip) = Vaux(ip) + Vna(ip)
      enddo
!$OMP end parallel do
      Enaatm = Enaatm * dVol
      Enascf = Enaatm + Enascf * dVol


      if ( SLABDIPOLECORR == SLABDIPOLECORR_VACUUM ) then
        ! In case the user requested a vacuum dipole, calculate it here
        ! REMARK
        ! This method will only ever work for 3D-FFT poisson solved
        ! potential
        ! Note we do this on the potential which will be written
        ! to ElectrostaticPotential
        call dipole_potential(cell, nua, isa, xa, ntm, ntml, nsm,
     &      dip_vacuum, Vaux, dipol)

        ! Orthogonalize dipole to bulk directions
        if ( nbcell == 1 ) then ! chain
          const = ddot(3,dipol,1,bcell,1) / ddot(3,bcell,1,bcell,1)
          dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
        else if ( nbcell == 2 ) then ! slab
          call cross( bcell(1,1), bcell(1,2), b1Xb2 )
          const = ddot(3,dipol,1,b1Xb2,1) / ddot(3,b1Xb2,1,b1Xb2,1)
          dipol(1:3) = const * b1Xb2(1:3)
        end if

      end if

C ----------------------------------------------------------------------
C Add potential from external electric field (if present)
C ----------------------------------------------------------------------
      if ( acting_efield ) then
        select case ( SLABDIPOLECORR )
        case ( SLABDIPOLECORR_NONE )
          ! No dipole correction
          field(:) = 0._dp
          DUext = 0._dp

        case default
          ! some kind of dipole correction
          field(:) = get_field_from_dipole(dipol, cell)

          if ( Node == 0 ) then
            write(6,'(a,3(tr1,f12.4),a)')
     $          'Dipole moment in unit cell =', dipol/Debye, ' D'
            write(6,'(a,3(tr1,f12.6),a)')
     $          'Electric field for dipole correction =',
     $          field/eV*Ang, ' eV/Ang/e'
          end if

          ! The dipole correction energy has an extra factor
          ! of one half because the field involved is internal.
          ! See the paper by Bengtsson DOI:10.1103/PhysRevB.59.12301
          ! Hence we compute this part separately
          DUext = -0.5_dp * ddot(3,field,1,dipol,1)
        end select

        ! Add the external electric field
        field(:) = field(:) + user_specified_field(:)

        ! This routine expects a sequential array,
        ! but it is distribution-blind
        ! Note that we do not need to specify a center
        ! Since we *just* need the saw-tooth to be located
        ! in the vacuum region (beyond orbitals)
        call add_potential_from_field( field, cell, nua, isa, xa,
     &      ntm, ntml, nsm, Vaux )

        ! Add energy of external electric field
        DUext = DUext - ddot(3,user_specified_field,1,dipol,1)

      end if

! ---------------------------------------------------------------------
!     Transiesta:
!     add the potential corresponding to the (possible) voltage-drop.
!     note that ts_voltage is not sharing the reord wih efield since
!     we should not encounter both at the same time.
! ---------------------------------------------------------------------
      if (TSmode.and.IsVolt.and.TSrun) then
         ! This routine expects a sequential array,
         ! in whatever distribution
         call ts_voltage(cell, ntm, ntml, Vaux)
      endif

! ----------------------------------------------------------------------
! Add potential from user defined geometries (if present)
! ----------------------------------------------------------------------
      call hartree_add( cell, ntpl, Vaux )

C ----------------------------------------------------------------------
C     Save electrostatic potential
C ----------------------------------------------------------------------
      if (filesOut%vh .ne. ' ') then
        ! Note that only the first spin component is used
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vh',1,ntml,
     &          Vaux)
        else
           call write_rho( filesOut%vh, cell, ntm, nsm, ntpl, 1, Vaux )
           call write_grid_netcdf( cell, ntm, 1, ntpl, Vaux, 
     &          "ElectrostaticPotential")
        end if
#else
        call write_rho( filesOut%vh, cell, ntm, nsm, ntpl, 1, Vaux )
        call write_grid_netcdf( cell, ntm, 1, ntpl,
     &       Vaux, "ElectrostaticPotential")
#endif
      endif

C     Get back spin density from sum and difference
      ! TODO Check for NC/SO:
      ! Should we diagonalize locally at every point first?
      if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(ip,rhotot)
        do ip = 1, ntpl
          rhotot     = DRho(ip,1)
          DRho(ip,1) = 0.5_dp * (rhotot - DRho(ip,2))
          DRho(ip,2) = 0.5_dp * (rhotot + DRho(ip,2))
        enddo
!$OMP end parallel do
      endif

C Exchange-correlation energy
C ----------------------------------------------------------------------
      call re_alloc( Vscf, 1, ntpl, 1, nspin, 'Vscf', 'dhscf' )

      if (npcc .eq. 1) then
         
!$OMP parallel default(shared), private(ip,ispin)
       do ispin = 1,nsd
!$OMP do
          do ip= 1, ntpl
             DRho(ip,ispin) = DRho(ip,ispin) +
     &            (rhopcc(ip)+rhoatm(ip)) * rnsd
          enddo
!$OMP end do nowait
       enddo
!$OMP end parallel

      else

!$OMP parallel default(shared), private(ip,ispin)
       do ispin = 1,nsd
!$OMP do
          do ip= 1, ntpl
             DRho(ip,ispin) = DRho(ip,ispin) +
     &            rhoatm(ip) * rnsd
          enddo
!$OMP end do nowait
       enddo
!$OMP end parallel

      end if

      ! Write the electron density used by cellxc
      if (filesOut%rhoxc .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','RhoXC',nspin,ntml,
     &          DRho )
        else
           call write_rho( filesOut%rhoxc, cell, ntm, nsm, ntpl, nspin,
     &          DRho )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &          "RhoXC")
        end if
#else
        call write_rho( filesOut%rhoxc, cell, ntm, nsm, ntpl, nspin,
     &       DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, 
     &       "RhoXC")
#endif
      endif

!     Everything now is in UNIFORM, sequential form

         call timer("XC",1)

         ! Switch to "zero/not-zero rho" distribution (miscalled 'LINEAR')
         if (nodes.gt.1) then
            call setMeshDistr( LINEAR, nsm, nsp, nml, nmpl, ntml, ntpl )
         endif

         call re_alloc( Vscf_gga, 1, ntpl, 1, nspin, 'Vscf_gga','dhscf')
         call re_alloc( DRho_gga, 1, ntpl, 1, nspin, 'DRho_gga','dhscf')

         ! Redistribute all spin densities
         do ispin = 1, nspin
            fsrc => DRho(:,ispin)
            fdst => DRho_gga(:,ispin)
            call distMeshData( UNIFORM, fsrc, LINEAR, fdst, KEEP )
         enddo

      if (use_bsc_cellxc) then

         call timer("BSC-CellXC",1)
         call bsc_cellxc( 0, 0, cell, ntml, ntml, ntpl, 0, aux3, nspin,
     &                    DRho_gga, Ex, Ec, DEx, DEc, Vscf_gga,
     &                    dummy_DVxcdn, stressl )

         call timer("BSC-CellXC",2)

      else
         
         call timer("GXC-CellXC",1)

         ! Note that RG's meshLim is in blocks of nsm points, and 1-based
         ! Example.
         ! 1:16 17:32  (base 1)
         ! 0:15 16:31  (base 0)
         ! Times nsm=2:
         ! 0:30 32:62
         ! Add 1 to lb, and nsm to ub:
         ! 1:32 33:64  (base 1, in units of small-points)

         myBox(1,:) = (meshLim(1,:)-1)*nsm + 1
         myBox(2,:) = (meshLim(2,:)-1)*nsm + nsm
         
         call cellXC( 0, cell, ntm, myBox(1,1), myBox(2,1),
     &                myBox(1,2), myBox(2,2),
     &                myBox(1,3), myBox(2,3), nspin,
     &                DRho_gga, Ex, Ec, DEx, DEc, stressXC, Vscf_gga,
     &                keep_input_distribution = .true. )
   !     Vscf is still sequential after the call to JMS's cellxc
         stress = stress + stressXC
         call timer("GXC-CellXC",2)
         
      endif  ! use_bsc_cellxc
      
         ! Go back to uniform distribution
         if (nodes.gt.1) then
            call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl)
         endif

         ! Redistribute to the Vxc array
         do ispin = 1,nspin
            fsrc => Vscf_gga(:,ispin)
            fdst => Vscf(:,ispin)
            call distMeshData( LINEAR, fsrc, UNIFORM, fdst, KEEP )
         enddo
         call de_alloc( DRho_gga, 'DRho_gga', 'dhscf' )
         call de_alloc( Vscf_gga, 'Vscf_gga', 'dhscf' )

         call timer("XC",2)

      if ( debug_dhscf ) then
         write(*,debug_fmt) Node,'XC',
     &        (sqrt(sum(Vscf(:,ispin)**2)),ispin=1,nspin)
      end if
      
      Exc =  Ex +  Ec
      Dxc = DEx + DEc

!     Vscf contains only Vxc, and is UNIFORM and sequential
!     Now we add up the other contributions to it, at 
!     the same time that we get DRho back to true DeltaRho form
!$OMP parallel default(shared), private(ip,ispin)

      ! Hartree potential only has diagonal components
      do ispin = 1,nsd
        if (npcc .eq. 1) then
!$OMP do
           do ip = 1,ntpl
              DRho(ip,ispin) = DRho(ip,ispin) - 
     &             (rhoatm(ip)+rhopcc(ip)) * rnsd
              Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
           enddo
!$OMP end do
        else
!$OMP do
           do ip = 1,ntpl
              DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) * rnsd
              Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
           enddo
!$OMP end do
        endif
      enddo
!$OMP end parallel

C ----------------------------------------------------------------------
C     Save total potential
C ----------------------------------------------------------------------
      if (filesOut%vt .ne. ' ') then
#ifdef NCDF_4
        if ( write_cdf ) then
           call cdf_save_grid(trim(slabel)//'.nc','Vt',nspin,ntml,
     &          Vscf)
        else 
           call write_rho( filesOut%vt, cell, ntm, nsm, ntpl, nspin,
     &          Vscf )
           call write_grid_netcdf( cell, ntm, nspin, ntpl, Vscf, 
     &          "TotalPotential")
        end if
#else
        call write_rho( filesOut%vt, cell, ntm, nsm, ntpl, nspin, Vscf)
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     &       Vscf, "TotalPotential")
#endif
      endif

C ----------------------------------------------------------------------
C Print vacuum level
C ----------------------------------------------------------------------

      if (filesOut%vt/=' ' .or. filesOut%vh/=' ') then
        forall(ispin=1:nsd) 
     .    DRho(:,ispin) = DRho(:,ispin) + rhoatm(:) * rnsd
        call vacuum_level( ntpl, nspin, DRho, Vscf, 
     .                     np_vac, Vmax_vac, Vmean_vac )
        forall(ispin=1:nsd) 
     .    DRho(:,ispin) = DRho(:,ispin) - rhoatm(:) * rnsd
        if (np_vac>0 .and. Node==0) print'(/,a,2f12.6,a)', 
     .    'dhscf: Vacuum level (max, mean) =', 
     .    Vmax_vac/eV, Vmean_vac/eV, ' eV'
      endif

C ----------------------------------------------------------------------
C     Find SCF contribution to hamiltonian matrix elements
C ----------------------------------------------------------------------
      if (iHmat .eq. 1) then
         if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
         endif

!        This is a work array, to which we copy Vscf
         call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                 'Vscf_par', 'dhscf' )

         do ispin = 1, nspin
           fsrc => Vscf(:,ispin)
           fdst => Vscf_par(:,ispin)
           call distMeshData( UNIFORM, fsrc,
     &                        QUADRATIC, fdst, TO_CLUSTER )
         enddo

         if (Spiral) then
            call vmatsp( norb, nml, nmpl, dvol, nspin, Vscf_par, maxnd,
     &           numd, listdptr, listd, Hmat, nuo,
     &           nuotot, iaorb, iphorb, isa, qspiral )
         else
            call vmat( norb, nmpl, dvol, spin, Vscf_par, maxnd,
     &           numd, listdptr, listd, Hmat, nuo,
     &           nuotot, iaorb, iphorb, isa )
         endif

         call de_alloc( Vscf_par,  'Vscf_par', 'dhscf' )
         if (nodes.gt.1) then
!          Everything back to UNIFORM, sequential
           call setMeshDistr( UNIFORM, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
         endif
      endif


#ifdef MPI
C     Global reduction of Uscf/DUscf/Uatm/Enaatm/Enascf
      sbuffer(1) = Uscf
      sbuffer(2) = DUscf
      sbuffer(3) = Uatm
      sbuffer(4) = Enaatm
      sbuffer(5) = Enascf
      if (use_bsc_cellxc) then
         sbuffer(6) = Exc
         sbuffer(7) = Dxc
      endif

      call MPI_AllReduce( sbuffer, rbuffer, 7, MPI_double_precision,
     &                     MPI_Sum, MPI_Comm_World, MPIerror )
      Uscf   = rbuffer(1)
      DUscf  = rbuffer(2)
      Uatm   = rbuffer(3)
      Enaatm = rbuffer(4)
      Enascf = rbuffer(5)
      if (use_bsc_cellxc) then
         Exc    = rbuffer(6)
         Dxc    = rbuffer(7)
      endif
#endif /* MPI */

C     Add contribution to stress from the derivative of the Jacobian of ---
C     r->r' (strained r) in the integral of Vna*(rhoscf-rhoatm)
      if (istr .eq. 1) then
        do i = 1,3
          stress(i,i) = stress(i,i) + ( Enascf - Enaatm ) / volume
        enddo
      endif

C     Stop time counter for SCF iteration part
      call timer( 'DHSCF3', 2 )

C ----------------------------------------------------------------------
C     End of SCF iteration part
C ----------------------------------------------------------------------

      if (ifa.eq.1 .or. istr.eq.1) then
C ----------------------------------------------------------------------
C Forces and stress : SCF contribution
C ----------------------------------------------------------------------
C       Start time counter for force calculation part
        call timer( 'DHSCF4', 1 )

C       Find contribution of partial-core-correction
        if (npcc .eq. 1) then
          call reord( rhopcc, rhopcc, nml, nsm, TO_CLUSTER )
          call reord( Vaux, Vaux, nml, nsm, TO_CLUSTER )
          ! The partial core calculation only acts on
          ! the diagonal spin-components (no need to
          ! redistribute un-used elements)
          do ispin = 1, nsd
            call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                  nml, nsm, TO_CLUSTER )
          enddo

          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua, nsd,
     &                            dvol, volume, Vscf, Vaux, Fal,
     &                            stressl, ifa.ne.0, istr.ne.0 )
   
          call reord( rhopcc, rhopcc, nml, nsm, TO_SEQUENTIAL )
          call reord( Vaux, Vaux, nml, nsm, TO_SEQUENTIAL )
          ! ** see above
          do ispin = 1, nsd
            call reord( Vscf(:,ispin), Vscf(:,ispin),
     &                  nml, nsm, TO_SEQUENTIAL )
          enddo

          if ( debug_dhscf ) then
             write(*,debug_fmt) Node,'PartialCore',
     &            (sqrt(sum(Vscf(:,ispin)**2)),ispin=1,nsd)
          end if
        endif

        if ( harrisfun) then
!         Forhar deals internally with its own needs
!         for distribution changes
!         Upon entry, everything is UNIFORM, sequential form
          call forhar( ntpl, nspin, nml, ntml, ntm, npcc, cell, 
     &                 rhoatm, rhopcc, Vna, DRho, Vscf, Vaux )
!         Upon return, everything is UNIFORM, sequential form
        endif

C     Transform spin density into sum and difference
        ! TODO NC/SO
        ! Should we perform local diagonalization?
        if (nsd .eq. 2) then
!$OMP parallel do default(shared), private(rhotot,ip)
          do ip = 1,ntpl
            rhotot     = DRho(ip,1) + DRho(ip,2)
            DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
            DRho(ip,1) = rhotot
          enddo
!$OMP end parallel do
        endif

C       Find contribution of neutral-atom potential
        call reord( Vna, Vna, nml, nsm, TO_CLUSTER )
        call reord( DRho, DRho, nml, nsm, TO_CLUSTER )
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     &                          volume, DRho, Fal, stressl, 
     &                          ifa.ne.0, istr.ne.0 )
        call reord( DRho, DRho, nml, nsm, TO_SEQUENTIAL )
        call reord( Vna, Vna, nml, nsm, TO_SEQUENTIAL )

        if (nodes.gt.1) then
           call setMeshDistr( QUADRATIC, nsm, nsp,
     &                       nml, nmpl, ntml, ntpl )
        endif

        call re_alloc( Vscf_par, 1, ntpl, 1, nspin,
     &                 'Vscf_par', 'dhscf' )
        do ispin = 1, nspin
          fsrc => Vscf(:,ispin)
          fdst => Vscf_par(:,ispin)
          call distMeshData( UNIFORM, fsrc,
     &                       QUADRATIC, fdst, TO_CLUSTER )
        enddo

!       Remember that Vaux contains everything except Vxc
        call re_alloc( Vaux_par, 1, ntpl, 'Vaux_par', 'dhscf' )
        call distMeshData( UNIFORM, Vaux,
     &                     QUADRATIC, Vaux_par, TO_CLUSTER )

        call dfscf( ifa, istr, na, norb, nuo, nuotot, nmpl, nspin,
     &              indxua, isa, iaorb, iphorb, 
     &              maxnd, numd, listdptr, listd, Dscf, datm,
     &              Vscf_par, Vaux_par, dvol, volume, Fal, stressl )

        call de_alloc( Vaux_par, 'Vaux_par', 'dhscf' )
        call de_alloc( Vscf_par, 'Vscf_par', 'dhscf' )

        if (nodes.gt.1) then
          call setMeshDistr( UNIFORM, nsm, nsp, nml, nmpl, ntml, ntpl )
        endif


C       Stop time counter for force calculation part
        call timer( 'DHSCF4', 2 )
C ----------------------------------------------------------------------
C       End of force and stress calculation
C ----------------------------------------------------------------------
      endif

!     We are in the UNIFORM distribution
!     Rhoatm, Rhopcc and Vna are in UNIFORM dist, sequential form
!     The index array endpht is in the QUADRATIC distribution

C     Stop time counter
      call timer( 'DHSCF', 2 )

C ----------------------------------------------------------------------
C     Free locally allocated memory
C ----------------------------------------------------------------------
      call de_alloc( Vaux, 'Vaux', 'dhscf' )
      call de_alloc( Vscf, 'Vscf', 'dhscf' )
      call de_alloc( DRho, 'DRho', 'dhscf' )

#ifdef DEBUG
      call write_debug( '    POS DHSCF' )
#endif
!------------------------------------------------------------------------ END
      CONTAINS

      subroutine save_bader_charge()
      use meshsubs, only: ModelCoreChargeOnMesh
#ifdef NCDF_4
      use siesta_options, only: write_cdf
      use m_ncdf_siesta, only: cdf_save_grid
#endif
      ! Auxiliary routine to output the Bader Charge
      !
      real(grid_p), pointer :: BaderCharge(:) => null()

      call re_alloc( BaderCharge, 1, ntpl, name='BaderCharge',
     &                 routine='dhscf' )

      ! Find a model core charge by re-scaling the local charge
      call ModelCoreChargeOnMesh( na, isa, ntpl, BaderCharge, indxua )
      ! It comes out in clustered form, so we convert it
      call reord( BaderCharge, BaderCharge, nml, nsm, TO_SEQUENTIAL)
      do ispin = 1,nsd
         BaderCharge(1:ntpl) = BaderCharge(1:ntpl) + DRho(1:ntpl,ispin)
      enddo
#ifdef NCDF_4
      if ( write_cdf ) then
         call cdf_save_grid(trim(slabel)//'.nc','RhoBader',1,ntml,
     &        BaderCharge)
      else
         call write_rho( trim(slabel)// ".BADER", cell,
     $        ntm, nsm, ntpl, 1, BaderCharge )
         call write_grid_netcdf( cell, ntm, 1, ntpl,
     $        BaderCharge, "BaderCharge")
      end if
#else
      call write_rho( trim(slabel)// ".BADER", cell,
     $     ntm, nsm, ntpl, 1, BaderCharge )
      call write_grid_netcdf( cell, ntm, 1, ntpl,
     $     BaderCharge, "BaderCharge")
#endif

      call de_alloc( BaderCharge, name='BaderCharge' )
      end subroutine save_bader_charge

      subroutine setup_analysis_options()
      !! For the analyze-charge-density-only case,
      !! avoiding any diagonalization
      
      use siesta_options, only: hirshpop, voropop
      use siesta_options, only: saverho, savedrho, saverhoxc
      use siesta_options, only: savevh, savevt, savevna
      use siesta_options, only: savepsch, savetoch
      
      want_partial_charges = (hirshpop .or. voropop) 

      if (saverho)   filesOut%rho   = trim(slabel)//'.RHO' 
      if (savedrho)  filesOut%drho  = trim(slabel)//'.DRHO'
      if (saverhoxc) filesOut%rhoxc = trim(slabel)//'.RHOXC'
      if (savevh)    filesOut%vh    = trim(slabel)//'.VH'  
      if (savevt)    filesOut%vt    = trim(slabel)//'.VT'  
      if (savevna)   filesOut%vna   = trim(slabel)//'.VNA' 
      if (savepsch)  filesOut%psch  = trim(slabel)//'.IOCH'
      if (savetoch)  filesOut%toch  = trim(slabel)//'.TOCH'
      
      end subroutine setup_analysis_options
      
      end subroutine dhscf

      subroutine delk_wrapper(isigneikr, norb, maxnd,
     &                        numd, listdptr, listd,
     &                        nuo,  nuotot, iaorb, iphorb, isa )

      use m_delk,  only  : delk  ! The real workhorse, similar to vmat

      use moreMeshSubs,   only : setMeshDistr
      use moreMeshSubs,   only : UNIFORM, QUADRATIC
      use parallel,       only : Nodes
      use mesh,           only : nsm, nsp

!
!     This is a wrapper to call delk, using some of the module
!     variables of m_dhscf, but from outside dhscf itself.
!
      integer                  :: isigneikr, 
     &                            norb, nuo, nuotot, maxnd,
     &                            iaorb(*), iphorb(*), isa(*), 
     &                            numd(nuo),     
     &                            listdptr(nuo), listd(maxnd)

! The dhscf module variables used are:
! nmpl
! dvol
! nml
! nmpl
! ntml
! ntpl
!
! Some of them might be put somewhere else (mesh?) to allow some
! of the kitchen-sink functionality of dhscf to be made more modular.
! For example, this wrapper might live independently if enough mesh
! information is made available to it.

C ----------------------------------------------------------------------
C Calculate matrix elements of exp(i \vec{k} \cdot \vec{r})
C ----------------------------------------------------------------------
        if (isigneikr .eq. 1 .or. isigneikr .eq. -1) then

          if (nodes.gt.1) then
            call setMeshDistr( QUADRATIC, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
          endif

          call delk( isigneikr, norb, nmpl, dvol, maxnd,
     &               numd, listdptr, listd,
     &               nuo,  nuotot, iaorb, iphorb, isa )

          if (nodes.gt.1) then
!           Everything back to UNIFORM, sequential
            call setMeshDistr( UNIFORM, nsm, nsp,
     &                         nml, nmpl, ntml, ntpl )
          endif

        endif
      end subroutine delk_wrapper

      !> Chech whether we need to initialize stencils for
      !> BSC's version of cellxc, and whether we have a vdw
      !> functional
      subroutine xc_setup(need_grads,is_vdw)
      use siestaxc, only: getxc

      logical, intent(out) :: need_grads
      logical, intent(out) :: is_vdw

      integer         :: nf, nXCfunc
      character(len=20):: XCauth(10), XCfunc(10)

      need_grads = .false.
      is_vdw = .false.
      call getXC( nXCfunc, XCfunc, XCauth )
      do nf = 1,nXCfunc
        if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
           need_grads = .true.
        endif
        if ( XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
           need_grads = .true.
           is_vdw = .true.
        endif
      enddo
      end subroutine xc_setup
      
      end module m_dhscf
