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
      subroutine dhscf( nspin, norb, iaorb, iphorb, nuo,
     .                  nuotot, nua, na, isa, xa, indxua, ucell, mscell,
     .                  G2max, ntm, ifa, istr, iHmat, filrho, fildrh, 
     .                  filevh, filevt, filepsch, filetoch, maxnd, numd,
     .                  listdptr, listd, 
     .                  Dscf, datm, maxnh, numh, listhptr, listh, Hmat,
     .                  Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, 
     .                  Exc, Dxc, dipol, Fa, stress, Fal, stressl )

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
C                         nspin=4 => Noncollinear spin
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
C real*8  ucell(3,3)    : Unit cell lattice vectors: ucell(ixyz,ivect)
C integer mscell(3,3)   : Supercell vectors in units of ucell:
C                         Supercell(ix,iv) = Sum_jv( ucell(ix,jv) *
C                                                   mscell(jv,iv) )
C integer ifa           : Switch which fixes whether the SCF contrib.
C                         to atomic forces is calculated and added to fa
C integer istr          : Switch which fixes whether the SCF contrib.
C                         to stress is calculated and added to stress
C integer iHmat         : Switch which fixes whether the Hmat matrix
C                         elements are calculated or not.
C character*(*) filrho  : Name of file to saving the electron density
C                         (If blank => not saved)
C character*(*) fildrh  : Name of file to saving Delta Rho (Rho - Rho_atoms)
C                         (If blank => not saved)
C character*(*) filevh  : Name of file to save electrostatic potential
C                         (If blank => not saved)
C character*(*) filevt  : Name of file to save total potential
C                         (If blank => not saved)
C character*(*) filepsch: Name of file to save the soft diffuse ionic charge
C                         (If blank => not saved)
C character*(*) filetoch: Name of file to save the total ionic+electronic charge
C                         (If blank => not saved)
C integer maxnd             : First dimension of listd and Dscf
C integer numd(nuo)         : Number of nonzero density-matrix
C                             elements for each matrix row
C integer listdptr(nuo)     : Pointer to start of rows of density-matrix
C integer listd(maxnd)      : Nonzero-density-matrix-element column
C                             indexes for each matrix row
C real*8  Dscf(maxnd,nspin): SCF density-matrix elements
C real*8  datm(norb)        : Harris density-matrix diagonal elements
C                             (atomic occupation charges of orbitals)
C integer maxnh             : First dimension of listh and Hmat
C integer numh(nuo)         : Number of nonzero Hamiltonian-matrix
C                             elements for each matrix row
C integer listhptr(nuo)     : Pointer to start of rows of H-matrix
C integer listh(maxnh)      : Nonzero-Hamiltonian-matrix-element column
C                             indexes for each matrix row
C real*8  Hmat(maxnh,nspin) : Hamiltonian matrix in sparse form,
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
C                   only when the system is a molecule
C ----------------------------------------------------------------------
C Input/output :
C ----------------------------------------------------------------------
C real*8  G2max      : Effective planewave cutoff in Ry. determines
C                        the mesh density and the precision of integrals
C                          on input : value required
C                          on output: value used, which may be larger
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

C
C  Modules
C
      use precision
      use parallel,      only : Node, Nodes
      use atmfuncs,      only : rcut, rcore
      use units,         only : Debye
      use fdf
      use sys,           only : die, bye
      use mesh,          only : xdsp, nsm, nsp
      use parsing
      use m_iorho,       only : write_rho
      use m_forhar,      only : forhar
      use xcmod,         only : nXCfunc, XCauth
      use alloc,         only : re_alloc, de_alloc
      use files,         only : slabel
#ifdef MPI
      use mpi_siesta
#endif
      use iogrid_netcdf, only: write_grid_netcdf
      use iogrid_netcdf, only: read_grid_netcdf
      use iogrid_netcdf, only: set_box_limits
      use siesta_options, only: read_charge_cdf
      use siesta_options, only: savevna, savebader
      use siesta_options, only: read_deformation_charge_cdf
      use siesta_options, only: save_initial_charge_density
!
      use m_efield,       only: initialize_efield, acting_efield
      use m_efield,       only: get_field_from_dipole, dipole_correction
      use m_efield,       only: add_potential_from_field
      use m_efield,       only: user_specified_field

      use m_doping,       only: initialize_doping, doping
      use m_doping,       only: compute_doping_structs, doping_active
      use m_doping,       only: ADD_BACKGROUND, REMOVE_BACKGROUND
      use m_partial_charges, only: compute_partial_charges
      use m_partial_charges, only: want_partial_charges
      
#ifdef TRANSIESTA
      use siesta_options, only: isolve
      use m_ts_options, only: UseVFix
      use m_ts_in_siesta, only: TSVHFix, TSVoltage2
      use m_ts_global_vars, only: TSrun
#endif /* TRANSIESTA */

      implicit none

      integer
     .  maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,
     .  iaorb(norb), ifa, iHmat, indxua(na), iphorb(norb), isa(na), 
     .  istr, listd(*), listh(maxnh), listdptr(nuo), listhptr(nuo),
     .  mscell(3,3), ntm(3), numd(nuo), numh(nuo)

      real(dp)
     .  datm(norb), dipol(3), Dscf(maxnd,nspin),
     .  DUext, DUscf, Dxc, Enaatm, Enascf, Exc, G2max, G2mesh,
     .  Hmat(maxnh,nspin), Uatm, ucell(3,3), Uscf, xa(3,na),
     .  Fa(3,nua), stress(3,3), Fal(3,nua), stressl(3,3)

      character
     .  filevh*(*), filevt*(*), fildrh*(*), filrho*(*), 
     .  filepsch*(*), filetoch*(*)

C ----------------------------------------------------------------------
C Routines called internally:
C ----------------------------------------------------------------------
C        cellxc(...)    : Finds total exch-corr energy and potential
C        cross(a,b,c)   : Finds the cross product of two vectors
C        dfscf(...)     : Finds SCF contribution to atomic forces
C        dipole(...)    : Finds electric dipole moment
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
C real*8 volcel( cell ) : Returns volume of unit cell
C ----------------------------------------------------------------------
C Internal variables and arrays:
C ----------------------------------------------------------------------
C real*8  bcell(3,3)    : Bulk lattice vectors
C real*8  cell(3,3)     : Auxiliary lattice vectors (same as ucell)
C real*8  const         : Auxiliary variable (constant within a loop)
C real*8  DEc           : Auxiliary variable to call cellxc
C real*8  DEx           : Auxiliary variable to call cellxc
C real*8  dvol          : Mesh-cell volume
C real*8  Ec            : Correlation energy
C real*8  Ex            : Exchange energy
C logical frstme        : First time this routine is called?
C integer i             : General-purpose index
C integer ia            : Atom index
C integer io            : Orbital index
C integer ip            : Point index
C integer is            : Species index
C logical IsDiag        : Is supercell diagonal?
C integer ispin         : Spin index
C integer j             : General-purpose index
C integer nbcell        : Number of independent bulk lattice vectors
C integer nm(3)         : Number of Mesh divisions of each cell vector
C integer npcc          : Partial core corrections? (0=no, 1=yes)
C integer nsd           : Number of diagonal spin values (1 or 2)
C integer ntpl          : Number of mesh Total Points in unit cell
C                         (including subpoints) locally
C integer ntopl         : Total number of nonzero orbital points locally
C real*4  rhoatm(ntpl)  : Harris electron density
C real*4  rhopcc(ntpl)  : Partial-core-correction density for xc
C real*4  DRho(ntpl)    : Selfconsistent electron density difference
C real*8  rhotot        : Total density at one point
C real*8  rmax          : Maximum orbital radius
C logical samesh        : Same mesh of last call?
C logical samexa        : Same atomic positions of last call?
C real*8  scell(3,3)    : Supercell vectors
C character shape*10    : Name of system shape
C real*8  tiny          : A small constant
C real*8  Uharrs        : Hartree energy of Harris electron density
C real*4  Vaux(ntpl)    : Auxiliary potential array
C real*4  Vna(ntpl)     : Sum of neutral-atom potentials
C real*8  volume        : Unit cell volume
C real*4  Vscf(ntpl)    : Hartree potential of selfconsistent density
C real*8  x0(3)         : Center of molecule
C logical harrisfun     : Harris functional or Kohn-Sham?
C ----------------------------------------------------------------------

C
C Internal variables
C
      integer
     .  i, ia, io, ip, iphi, is, ispin, j, ifor,
     .  n, nbcell, nm(3), npcc, nsd, nv, nzero(3), 
     .  ntopl, nmpl, nsc(3), ntpl, nml(3), ntml(3)

#ifdef MPI
      integer
     .  MPIerror
#endif

      real(grid_p), dimension(:), pointer, save :: rhopcc, rhoatm, Vna
      real(grid_p), dimension(:),     pointer, save   :: Vaux
      real(grid_p), dimension(:,:),   pointer, save   :: DRho, Vscf
      real(grid_p), dimension(:), pointer, save   :: Chlocal, Totchar
!
!     Dummy arrays for cellxc call
!
      real(grid_p) ::  aux3(3,1)
      real(grid_p), dimension(1,1,1)  :: DVxcdn_dummy

      real(dp)
     .  bcell(3,3), b1Xb2(3), cell(3,3), const,
     .  DEc, DEx, ddot, DStres(3,3), dvol, Ec, Ex, field(3),
     .  rhotot, rmax, scell(3,3), tiny, 
     .  Uharrs, volume, volcel, x0(3)

      real(dp), save ::  qspiral(3)

#ifdef MPI
      real(dp)  ::  Eloc
#endif

      logical :: leqi

      logical, save :: 
     .  frstme, harrisfun, IsDiag, samesh, samexa, spiral

      character(len=10), save ::
     .  shape

      external
     .  cellxc, cross, dfscf, dipole, ddot, 
     .  poison, 
     .  reord, rhooda, rhoofd, rhoofdsp, 
     .  timer, vmat, vmatsp, volcel,
     .  readsp
     .  , update_hfx, update_hfx_force

      save
     .  bcell, cell, dvol, G2mesh,
     .  nm, nml, ntml, npcc, nzero, nmpl, ntpl, ntopl,
     .  rmax, scell, tiny, Uharrs, volume

      data
     .  G2mesh / 0.0_dp /,
     .  nzero / 3*0 /,
     .  tiny  / 1.e-12_dp /,
     .  frstme /.true./
#ifdef TRANSIESTA
! TSS Begin
      character(6) :: method
      logical, save :: transmode=.false., isvolt
! TSS End

#endif /* TRANSIESTA */

C ----------------------------------------------------------------------
C General initialisation
C ----------------------------------------------------------------------

C Start time counter
      call timer( 'DHSCF', 1 )

      if (frstme) then
#ifdef TRANSIESTA
! TSS Begin
        if(isolve==2) transmode=.true.
! TSS End
#endif /* TRANSIESTA */
        nullify(Vna,rhoatm,rhopcc)
        nullify(vaux,vscf,DRho)
        nullify(chlocal,totchar)
        if (Node.eq.0) then
          harrisfun = fdf_boolean('Harris_functional',.false.)
          nsm    = fdf_integer('MeshSubDivisions',2)
          nsm    = max(nsm,1)
        endif

#ifdef MPI
        call MPI_Bcast(harrisfun,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
        call MPI_Bcast(nsm,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
#ifdef TRANSIESTA
#ifdef MPI
! TSS Begin
        call MPI_Bcast(transmode,1,MPI_logical,0,MPI_Comm_World,
     .       MPIerror)
! TSS End
#endif
#endif /* TRANSIESTA */

C Set mesh sub-division variables & perform one off allocation
        nsp = nsm*nsm*nsm
        nullify( xdsp )
        call re_alloc( xdsp, 1, 3, 1, nsp, name='xdsp',
     &                 routine='dhscf' )

C Read spin-spiral wavevector (if defined)
        call readsp( qspiral, spiral )
        if (spiral .and. nspin.lt.4)
     $    call die('dhscf: ERROR: spiral defined but nspin < 4')

      endif

      if (.not.harrisfun) then
        ifor = 0
      else if (harrisfun) then
        ifor = 1
        do n = 1,nXCfunc
          if (.not.(leqi(XCauth(n),'PZ').or.leqi(XCauth(n),'CA'))) then
            if (Node .eq. 0) then
              write(6,100)
              write(6,101)
              write(6,'(a)')
     .         'dhscf:   Harris forces for a XC functional different of'
              write(6,'(a)')
     .         'dhscf:   Perdew-Zunger has not been implemented yet'
              write(6,102)
            endif
            call die()
          endif
        enddo
      endif

C ----------------------------------------------------------------------
C Orbital initialisation : part 1
C ----------------------------------------------------------------------

C Find the maximum orbital radius
      rmax = 0.0_dp
      do io = 1,norb
        ia = iaorb(io)
        iphi = iphorb(io)
        is = isa(ia)
        rmax = max( rmax, rcut(is,iphi) )
      enddo

C ----------------------------------------------------------------------
C Check whether mesh initialisation is needed
C ----------------------------------------------------------------------
      call SameMeshAndAtoms( na, xa, ucell, rmax, G2max, G2mesh,
     .  samesh, samexa )

C ----------------------------------------------------------------------
C Initialization due to change of mesh
C ----------------------------------------------------------------------
      if ( .not.samesh ) then

C Start time counter for mesh initialization
        call timer( 'DHSCF1', 1 )


C ----------------------------------------------------------------------
C Unit cell handling
C ----------------------------------------------------------------------
C Find diagonal unit cell and supercell
        call digcel( ucell, mscell, cell, scell, nsc, IsDiag )
        if (.not.IsDiag) then
          if (Node.eq.0) then
            write(6,'(/,a,3(/,a,3f12.6,a,i6))')
     .      'DHSCF: WARNING: New shape of unit cell and supercell:',
     .     ('DHSCF:',(cell(i,j),i=1,3),'   x',nsc(j),j=1,3)
          endif
        endif

C Find the system shape
        call shaper( cell, nua, isa, xa, shape, nbcell, bcell )

C Find system volume
        volume = volcel( cell )

C ----------------------------------------------------------------------
C Mesh initialization 
C ----------------------------------------------------------------------
        call InitMesh( na, cell, norb, iaorb, iphorb, isa, rmax, 
     .                 G2max, G2mesh, nsc, nmpl, nm,
     .                 nml, ntm, ntml, ntpl, ntopl, dvol)

        ! Setup box descriptors for each processor,
        ! held in module iogrid_netcdf
        call set_box_limits(ntm,nsm)

C ----------------------------------------------------------------------
C Allocate memory that depends on the mesh
C ----------------------------------------------------------------------

C Allocate memory related to ntpl
        call re_alloc( rhoatm, 1, ntpl, name='rhoatm', routine='dhscf' )
        call re_alloc( Vna, 1, ntpl, name='Vna', routine='dhscf' )

C Find if there are partial-core-corrections for any atom
        npcc = 0
        do ia = 1,na
          if (rcore(isa(ia)) .gt. tiny) npcc = 1
        enddo

C Allocate memory for PCC density
        call re_alloc( rhopcc, 1, ntpl*npcc+1, name='rhopcc',
     &                 routine='dhscf' )

C Stop time counter for mesh initialization
        call timer( 'DHSCF1', 2 )

      endif
C ----------------------------------------------------------------------
C End of mesh initialization 
C ----------------------------------------------------------------------

C Allocate local memory
      call re_alloc( DRho, 1, ntpl, 1, nspin, name='DRho',
     &               routine='dhscf' )
      call re_alloc( Vscf, 1, ntpl, 1, nspin, name='Vscf',
     &               routine='dhscf' )
      call re_alloc( Vaux, 1, ntpl, name='Vaux', routine='dhscf' )

      if (.not.samexa) then
C ----------------------------------------------------------------------
C Initialize atomic orbitals, density and potential
C ----------------------------------------------------------------------

C Start time counter for atomic initializations
        call timer( 'DHSCF2', 1 )

C Initialise quantities relating to the atom-mesh positioning
        call InitAtomMesh( na, xa)

C Find partial-core-correction energy density
        if (npcc .eq. 1) then
          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     .      nsd, dvol, volume, Vscf, Vaux, Fal, stressl, .false.,
     .      .false. )
        endif

C Find neutral-atom potential
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     .                          volume, DRho, Fal, stressl, 
     .                          .false., .false. )
        if (savevna) then
           call reord( Vna, Vna, nml, nsm, +1)
           call write_rho( trim(slabel)// ".VNA",
     $                     cell, ntm, nsm, ntpl, 1, Vna)
           call write_grid_netcdf( cell, ntm, 1, ntpl, Vna, "Vna")
           call reord( Vna, Vna, nml, nsm, -1)
        endif

C Calculate orbital values on mesh
        call PhiOnMesh( nmpl, norb, iaorb, iphorb, isa )

C ----------------------------------------------------------------------
C Create sparse indexing for Dscf as needed for local mesh
C ----------------------------------------------------------------------
        if (Nodes.gt.1) then
          call CreateLocalDscfPointers( nmpl, nuotot, numd, listdptr, 
     .                                  listd )
        endif

C ----------------------------------------------------------------------
C Calculate terms relating to the neutral atoms on the mesh
C ----------------------------------------------------------------------

C Find Harris (sum of atomic) electron density
        call rhooda( norb, nmpl, datm, rhoatm, iaorb, iphorb, isa )

C Solve Poisson's equation to find Hartree energy of rhoatm
        call reord( rhoatm, rhoatm, nml, nsm, +1)
        call poison( cell, ntml(1), ntml(2), ntml(3), ntm, rhoatm,
     .               Uharrs, Vscf, DStres, nsm )
        call reord( rhoatm, rhoatm, nml, nsm, -1)

C Stop time counter for atomic initializations
        call timer( 'DHSCF2', 2 )

C ----------------------------------------------------------------------
C End of atomic initializations
C ----------------------------------------------------------------------
      endif

C ----------------------------------------------------------------------
C Start of SCF iteration part
C ----------------------------------------------------------------------
      call timer( 'DHSCF3', 1 )

C Find number of diagonal spin values
      nsd = min( nspin, 2 )

C ----------------------------------------------------------------------
C Find SCF electron density at mesh points. Store it in array DRho
C ----------------------------------------------------------------------

      if (read_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"Rho")
         do ispin = 1,nspin
            call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
         enddo
         read_charge_cdf = .false.
      else if (read_deformation_charge_cdf) then
         call read_grid_netcdf(ntm(1:3),nspin,ntpl,DRho,"DeltaRho")
         do ispin = 1,nspin
            call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
            DRho(:,ispin) = DRho(:,ispin) + rhoatm(:)
         enddo
         read_deformation_charge_cdf = .false.
      else
         if (spiral) then
            call rhoofdsp(norb, nmpl, maxnd, numd, listdptr, listd,
     .                nspin, Dscf, DRho, nuo, nuotot, iaorb, iphorb,
     .                isa, qspiral)
         else
            call rhoofd(norb, nmpl, maxnd, numd, listdptr, listd,
     .              nspin, Dscf, DRho, nuo, nuotot, iaorb, iphorb, isa)
!
            if (save_initial_charge_density) then
               do ispin = 1,nspin
                  call reord( DRho(1,ispin), DRho(1,ispin),
     $                 nml, nsm, +1)
               enddo
               call write_rho( "RHO_INIT", cell, ntm, nsm, ntpl,
     $              nspin, DRho)
               call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho,
     $               "RhoInit")
               do ispin = 1,nspin
                  call reord( DRho(1,ispin), DRho(1,ispin),
     $                 nml, nsm, -1)
               enddo
              call bye("STOP after producing RHO_INIT")
            endif
         endif
      endif
!
!     Proper place to integrate Hirshfeld and Voronoi code,
!     since we have just computed rhoatm and Rho.
      if (want_partial_charges) then
         call compute_partial_charges(DRho,rhoatm,
     .                  nspin, iaorb, iphorb, 
     .                  isa, nmpl,dvol)
      endif

C ----------------------------------------------------------------------
C Save electron density
C ----------------------------------------------------------------------
      if (filrho .ne. ' ') then
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, +1)
        enddo
        call write_rho( filrho, cell, ntm, nsm, ntpl, nspin, DRho)
        call write_grid_netcdf( cell, ntm, nspin, ntpl, DRho, "Rho")
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
        enddo
      endif

C ----------------------------------------------------------------------
C Save the diffuse ionic charge and/or the total (ionic+electronic) charge
C ----------------------------------------------------------------------

      if (filepsch .ne. ' ' .or. filetoch .ne. ' ') then

        call re_alloc( Chlocal, 1, ntpl, name='Chlocal',
     &                 routine='dhscf' )

C Find diffuse ionic charge on mesh
        call LocalChargeOnMesh( na, isa, ntpl, Chlocal, indxua )
      
C Save diffuse ionic charge 
        if (filepsch .ne. ' ') then
          call reord( Chlocal, Chlocal, nml, nsm, +1)
          call write_rho( filepsch, cell, ntm, nsm, ntpl, 1, Chlocal)
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                            Chlocal, "Chlocal")
          call reord( Chlocal, Chlocal, nml, nsm, -1)
        endif

C Save total (ionic+electronic) charge 
        if ( filetoch .ne. ' ') then
           call re_alloc( Totchar, 1, ntpl, name='Totchar',
     &          routine='dhscf' )
        
          Totchar(:) = 0.0_grid_p
          do ispin = 1,nsd
            Totchar(1:ntpl) = Totchar(1:ntpl) + 
     .                        Chlocal(1:ntpl)/nsd + DRho(1:ntpl,ispin)
          enddo
          call reord( Totchar, Totchar, nml, nsm, +1)
          call write_rho( filetoch, cell, ntm, nsm, ntpl, 1, Totchar )
          call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                            Totchar, "TotalCharge")
          call reord( Totchar, Totchar, nml, nsm, -1)

          call de_alloc( Totchar, name='Totchar' )
        endif 

        call de_alloc( Chlocal, name='Chlocal' )
      endif

C ----------------------------------------------------------------------
C Save the total charge (model core + valence) for Bader analysis
C ----------------------------------------------------------------------

      if (filetoch .ne. ' ' .and. savebader) then
        call save_bader_charge()
      endif


C Find difference between selfconsistent and atomic densities
      do ispin = 1,nsd
        DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) - rhoatm(1:ntpl)/nsd
      enddo

C ----------------------------------------------------------------------
C Save electron density difference
C ----------------------------------------------------------------------
      if (fildrh .ne. ' ') then
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, +1)
        enddo
        call write_rho( fildrh, cell, ntm, nsm, ntpl, nspin, DRho )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     $                          DRho, "DeltaRho")
        do ispin = 1,nspin
          call reord( DRho(1,ispin), DRho(1,ispin), nml, nsm, -1)
        enddo
      endif

C Transform spin density into sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot = DRho(ip,1) + DRho(ip,2)
          DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
          DRho(ip,1) = rhotot
        enddo
      endif

C Add a background charge to neutralize the net charge, to
C model doped systems. It only adds the charge at points
C where there are atoms (i.e., not in vacuum).          

C NOTE: This section was formerly reached only for non-bulk systems,
C       while the corresponding removal was done unconditionally.
C       In the interest of generality, even if it might not make sense
C       for bulk systems, it can now be used for all systems.

        if (frstme) call initialize_doping()
        if (doping_active .and. (.not. (samexa.and.samesh))) then
           call compute_doping_structs(nmpl,rhoatm,nsd)
        endif
        if (doping_active) call doping(cell,nmpl,ADD_BACKGROUND,DRho)

C ----------------------------------------------------------------------
C Calculate the dipole moment
C ----------------------------------------------------------------------
      dipol(1:3) = 0.0_dp
      if (shape .ne. 'bulk') then

C Find center of system
        x0(1:3) = 0.0_dp
        do ia = 1,nua
          x0(1:3) = x0(1:3) + xa(1:3,ia) / nua
        enddo

C Find dipole
        call reord( DRho, DRho, nml, nsm, +1)
        call dipole( cell, ntm, ntml(1), ntml(2), ntml(3), nsm,
     .    DRho, x0, dipol )
        call reord( DRho, DRho, nml, nsm, -1)

C Orthogonalize dipole to bulk directions
        if (shape .eq. 'chain') then
          const = ddot(3,dipol,1,bcell,1) / ddot(3,bcell,1,bcell,1)
          dipol(1:3) = dipol(1:3) - const * bcell(1:3,1)
        elseif (shape .eq. 'slab') then
          call cross( bcell(1,1), bcell(1,2), b1Xb2 )
          const = ddot(3,dipol,1,b1Xb2,1) / ddot(3,b1Xb2,1,b1Xb2,1)
          dipol(1:3) = const * b1Xb2(1:3)
        endif
      endif

C ----------------------------------------------------------------------
C Find Hartree potential of DRho = rhoscf-rhoatm. Store it in Vscf
C ----------------------------------------------------------------------

C Initialise Vscf
      Vscf(1:ntpl,1) = 0.0_grid_p

C Solve Poisson's equation
      call reord( DRho, DRho, nml, nsm, +1)
      call poison( cell, ntml(1), ntml(2), ntml(3), ntm, DRho,
     .             DUscf, Vscf, DStres, nsm )
#ifdef TRANSIESTA

! TSS BEgin
      if(transmode .and. UseVFix)
     .     call TSVHFix( ntm, nsm, vscf )
! TSS End
#endif /* TRANSIESTA */
      call reord( DRho, DRho, nml, nsm, -1)
      call reord( Vscf, Vscf, nml, nsm, -1)

C Add contribution to stress from electrostatic energy of rhoscf-rhoatm
      if (istr .eq. 1) then
        stressl(1:3,1:3) = stressl(1:3,1:3) + DStres(1:3,1:3)
      endif

C ----------------------------------------------------------------------
C Find electrostatic (Hartree) energy of full SCF electron density
C ----------------------------------------------------------------------
      Uatm = Uharrs
      Uscf = Uatm + DUscf
      do ip = 1,ntpl
        Uscf = Uscf + dvol * Vscf(ip,1) * rhoatm(ip)
      enddo

#ifdef MPI
C Global reduction of Uscf/DUscf/Uatm
      call MPI_AllReduce(Uscf,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
      Uscf = Eloc
      call MPI_AllReduce(DUscf,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
      DUscf = Eloc
      call MPI_AllReduce(Uatm,Eloc,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
      Uatm = Eloc
#endif

C Call doping again to remove background charge added previously
C The extra charge thus only affects the Hartree energy and potential,
C but not the contribution to Enascf ( = \Int_{Vna*\rho})

      if (doping_active) call doping(cell,nmpl,REMOVE_BACKGROUND,DRho)

C ----------------------------------------------------------------------
C Add neutral-atom potential to Vscf
C ----------------------------------------------------------------------

      Enaatm = 0.0_dp
      Enascf = 0.0_dp
      do ip = 1,ntpl
        Enaatm = Enaatm + dvol * Vna(ip) * rhoatm(ip)
        Enascf = Enascf + dvol * Vna(ip) * ( rhoatm(ip) + DRho(ip,1) )
        Vscf(ip,1) = Vscf(ip,1) + Vna(ip)
      enddo

#ifdef MPI
C Global reduction of Enaatm and Enascf
      Eloc = Enaatm
      call MPI_AllReduce(Eloc,Enaatm,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
      Eloc = Enascf
      call MPI_AllReduce(Eloc,Enascf,1,MPI_double_precision,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
#endif

C ----------------------------------------------------------------------
C Add potential from external electric field (if present)
C ----------------------------------------------------------------------

      if (frstme) call initialize_efield()

      if (acting_efield) then
         if (dipole_correction) then
            field = get_field_from_dipole(dipol,cell)
            if (Node == 0) then
	       write(6,'(a,3f12.4,a)')
     $              'Dipole moment in unit cell   =', dipol/Debye, ' D'
	       write(6,'(a,3f12.6,a)')
     $              'Electric field for dipole correction =',
     $              field, ' Ry/Bohr/e'
            endif
         else
            field = user_specified_field
         endif
         call reord( Vscf, Vscf, nml, nsm, +1)
         call add_potential_from_field( field, cell, nua, isa, xa,
     $                                  ntm, nsm, Vscf )
         call reord( Vscf, Vscf, nml, nsm, -1)
         DUext = - ddot(3,field,1,dipol,1)
      endif

#ifdef TRANSIESTA
c TSS Begin
c ---------------------------------------------------------------------
c transiesta:
c add the potential corresponding to the (possible) voltage-drop.
c note that tsvoltage is not sharing the reord wih efield since
c we should not encounter both at the same time.
c ---------------------------------------------------------------------
      if (transmode) then

      if (frstme) then

         call tsvoltage2(cell, nzero, nsm, Vscf, isvolt )

         if(isvolt .and. acting_efield) then
           if(Node.eq.0) then
             write(6,*) 'WARNING: both ext. voltage and field !!'
             write(6,*) 'WARNING: will neglect voltage !!       '
           end if
           isvolt = .false.
         end if
      endif

      if(isvolt ) then
        if(TSrun) then
         call reord( Vscf, Vscf, nml, nsm, +1)
         call tsvoltage2(cell, ntm, nsm, Vscf, isvolt )
         call reord( Vscf, Vscf, nml, nsm, -1)
        endif
      endif
      endif !TSS
c TSS End


#endif /* TRANSIESTA */

C ----------------------------------------------------------------------
C Save electrostatic potential
C ----------------------------------------------------------------------

      if (filevh .ne. ' ') then
        call reord( Vscf, Vscf, nml, nsm, +1)
        call write_rho( filevh, cell, ntm, nsm, ntpl, 1, Vscf)
        call write_grid_netcdf( cell, ntm, 1, ntpl,
     $                          Vscf, "ElectrostaticPotential")
        call reord( Vscf, Vscf, nml, nsm, -1)
      endif

C Add contribution to stress from the derivative of the Jacobian of ---
C r->r' (strained r) in the integral of Vna*(rhoscf-rhoatm)
        if (istr .eq. 1) then
          do i = 1,3
            stress(i,i) = stress(i,i) + ( Enascf - Enaatm ) / volume
          enddo
        endif

C Get back spin density from sum and difference
      if (nsd .eq. 2) then
        do ip = 1,ntpl
          rhotot = DRho(ip,1)
          DRho(ip,1) = 0.5_dp * ( rhotot - DRho(ip,2) )
          DRho(ip,2) = 0.5_dp * ( rhotot + DRho(ip,2) )
        enddo
      endif

C ----------------------------------------------------------------------
C Exchange-correlation energy
C ----------------------------------------------------------------------

      Vaux(1:ntpl) = Vscf(1:ntpl,1)
      do ispin = 1,nsd
        DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) + rhoatm(1:ntpl) / nsd
        if (npcc .eq. 1) 
     .    DRho(1:ntpl,ispin) = DRho(1:ntpl,ispin) + rhopcc(1:ntpl)/nsd
      enddo

      do ispin = 1,nspin
        call reord(DRho(1,ispin),DRho(1,ispin),nml,nsm,+1)
      enddo

      call cellxc( 0, 0, cell, ntml, ntml, ntpl, 0, aux3, nspin, DRho, 
     .             Ex, Ec, DEx, DEc, Vscf, DVxcdn_dummy, stressl )

      do ispin = 1,nspin
        call reord(DRho(1,ispin),DRho(1,ispin),nml,nsm,-1)
        call reord(Vscf(1,ispin),Vscf(1,ispin),nml,nsm,-1)
      enddo

      Exc = Ex + Ec
      Dxc = DEx + DEc
      do ispin = 1,nsd
        do ip = 1,ntpl
          DRho(ip,ispin) = DRho(ip,ispin) - rhoatm(ip) / nsd
          if (npcc .eq. 1) 
     .      DRho(ip,ispin) = DRho(ip,ispin) - rhopcc(ip) / nsd
          Vscf(ip,ispin) = Vscf(ip,ispin) + Vaux(ip)
        enddo
      enddo


C ----------------------------------------------------------------------
C Save total potential
C ----------------------------------------------------------------------

      if (filevt .ne. ' ') then
        do ispin = 1,nspin
          call reord( Vscf(1,ispin), Vscf(1,ispin), nml, nsm, +1)
        enddo
        call write_rho( filevt, cell, ntm, nsm, ntpl, nspin, Vscf )
        call write_grid_netcdf( cell, ntm, nspin, ntpl,
     $                          Vscf, "TotalPotential")
        do ispin = 1,nspin
          call reord( Vscf(1,ispin), Vscf(1,ispin), nml, nsm, -1)
        enddo
      endif

C ----------------------------------------------------------------------
C Find SCF contribution to hamiltonian matrix elements
C ----------------------------------------------------------------------

      if (iHmat .eq. 1) then
        if (spiral) then
          call vmatsp( norb, nmpl, dvol, nspin, Vscf, maxnd,
     .               numd, listdptr, listd, Hmat, nuo,
     .               nuotot, iaorb, iphorb, isa, qspiral )
        else
          call vmat( norb, nmpl, dvol, nspin, Vscf, maxnd,
     .               numd, listdptr, listd, Hmat, nuo,
     .               nuotot, iaorb, iphorb, isa )
        endif
      endif

C Stop time counter for SCF iteration part
      call timer( 'DHSCF3', 2 )

C ----------------------------------------------------------------------
C End of SCF iteration part
C ----------------------------------------------------------------------

C-----------------------------------------------------------------------
C add Hartree-Fock exchange potential to Hamiltonian (nuo, norb)
C xmqin, January 2014
C----------------------------------------------------------------------
      call timer( 'DHFX', 1 )
 
C      if(node.eq.0) write(6,*) 'Hybrid DFT calculations'
 
      do n = 1,nXCfunc
         if (leqi(XCauth(n),'HSE06') .or. leqi(XCauth(n),'hse06')
     .       .or. leqi(XCauth(n),'PBE0') .or . leqi(XCauth(n),'pbe0'))
     .   then
!! Both cell and supercell have been diagonalized.
      call update_hfx( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,
     .                 na, isa, xa, indxua, cell, nsc, maxnd, numd, 
     .                 listdptr, listd, Dscf, maxnh, numh, listhptr, 
     .                 listh, samexa, Hmat, Exc)
         endif
      enddo

      call timer( 'DHFX', 2 )
C------------------------------------------------------------------------

      if (ifa.eq.1 .or. istr.eq.1) then
C ----------------------------------------------------------------------
C Forces and stress : SCF contribution
C ----------------------------------------------------------------------

C Start time counter for force calculation part
        call timer( 'DHSCF4', 1 )

C Find contribution of partial-core-correction
        if (npcc .eq. 1) then
          call PartialCoreOnMesh( na, isa, ntpl, rhopcc, indxua,
     .      nsd, dvol, volume, Vscf, Vaux, Fal, stressl, ifa.ne.0,
     .      istr.ne.0 )
        endif

        if ( ifor .eq. 1)
     .    call forhar( ntpl, nspin, nml, ntml, ntm, npcc, cell, 
     .                 rhoatm, rhopcc, Vna, DRho, Vscf, Vaux )

C Transform spin density into sum and difference
        if (nsd .eq. 2) then
          do ip = 1,ntpl
            rhotot = DRho(ip,1) + DRho(ip,2)
            DRho(ip,2) = DRho(ip,2) - DRho(ip,1)
            DRho(ip,1) = rhotot
          enddo
        endif

C Find contribution of neutral-atom potential
        call NeutralAtomOnMesh( na, isa, ntpl, Vna, indxua, dvol, 
     .                          volume, DRho, Fal, stressl, 
     .                          ifa.ne.0, istr.ne.0 )

C Vaux is (minus) the potential which multiplies rhoatm
*       if (nsd .eq. 2) then
*         Vaux(1:ntpl) = 0.5_grid_p * Vaux(1:ntpl)
*       endif

        call dfscf( ifa, istr, na, norb, nuo, nuotot, nmpl, nspin,
     .              indxua, isa, iaorb, iphorb,
     .              maxnd, numd, listdptr, listd, Dscf, datm,
     .              Vscf, Vaux, dvol, volume, Fal, stressl )


C Stop time counter for force calculation part
        call timer( 'DHSCF4', 2 )


C-----------------------------------------------------------------------
C add Hartree-Fock exchange gradients to force (1:3, nua)
C xmqin, January 2014
C----------------------------------------------------------------------

      call timer( 'DFHFX', 1 )
      if(node.eq.0) write(6,*) "To calculate HFX gradient"
      do n = 1,nXCfunc
         if (leqi(XCauth(n),'HSE06') .or. leqi(XCauth(n),'hse06')
     .       .or. leqi(XCauth(n),'PBE0') .or . leqi(XCauth(n),'pbe0'))
     .   then

      call update_hfx_force( nspin, norb, iaorb, iphorb, nuo, nuotot,
     .               nua, na, isa, xa, indxua, cell, nsc, maxnd, numd,
     .               listdptr, listd, Dscf, maxnh, numh, listhptr,
     .               listh, Fal)

         endif
       enddo
      call timer( 'DFHFX', 2 )

C ----------------------------------------------------------------------
C End of force and stress calculation
C ----------------------------------------------------------------------
      endif


C Stop time counter
      call timer( 'DHSCF', 2 )

C ----------------------------------------------------------------------
C Free locally allocated memory
C ----------------------------------------------------------------------
      call de_alloc( Vaux, name='Vaux' )
      call de_alloc( Vscf, name='Vscf' )
      call de_alloc( DRho, name='DRho' )

      frstme = .false.

100   format(/,'dhscf: ',71(1h*))
101   format('dhscf:                  INPUT ERROR')
102   format('dhscf: ',71(1h*))

      CONTAINS

      subroutine save_bader_charge()
      ! Auxiliary routine to output the Bader Charge
      !
      real(grid_p), pointer :: BaderCharge(:) => null()

      call re_alloc( BaderCharge, 1, ntpl, name='BaderCharge',
     &                 routine='dhscf' )

      ! Find a model core charge by re-scaling the local charge
      call ModelCoreChargeOnMesh( na, isa, ntpl, BaderCharge, indxua )
      
      do ispin = 1,nsd
         BaderCharge(1:ntpl) = BaderCharge(1:ntpl) + DRho(1:ntpl,ispin)
      enddo
      call reord( BaderCharge, BaderCharge, nml, nsm, +1)
      call write_rho( trim(slabel)// ".BADER", cell,
     $     ntm, nsm, ntpl, 1, BaderCharge )
      call write_grid_netcdf( cell, ntm, 1, ntpl,
     $     BaderCharge, "BaderCharge")

      call de_alloc( BaderCharge, name='BaderCharge' )
      end subroutine save_bader_charge

      end subroutine dhscf

