      MODULE m_siesta_init
      private
      public :: siesta_init

      CONTAINS

      subroutine siesta_init()
      USE Kpoint_grid, only: setup_Kpoint_grid, gamma_scf, maxk
      USE Kpoint_pdos, only: gamma_pdos
      USE Band, only: gamma_bands, setup_bands
      USE m_ksvinit, only: gamma_polarization, estimate_pol_kpoints
      USE m_struct_init, only: struct_init
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: no_u, rmaxkb, amass, lasto, qtot, iza, rmaxo,
     .                    zvaltot, superc, initatomlists, no_l
      use m_fdf_global, only: fdf_global_get
      use sys, only: die, bye
      use xcmod, only: setXC
      use molecularmechanics, only : inittwobody
      use metaforce, only: initmeta
      use m_mpi_utils, only : broadcast
      use alloc, only: re_alloc, alloc_report
      use phonon, only: phonon_num_disps, phonon_setup
      use parallelsubs, only: getnodeorbs
      use m_iostruct,   only: write_struct, read_struct
      use zmatrix,      only: lUseZmatrix
      use zmatrix,      only: write_canonical_ucell_and_Zmatrix
      use m_check_supercell, only: check_sc_factors
      use files,       only : slabel
      use siesta_cmlsubs, only: siesta_cml_init
      use m_timestamp, only : timestamp
      use m_wallclock, only : wallclock
      use parallel, only: Node, Nodes, IOnode
      USE densematrix, only: Haux, Saux, Psi
      use m_energies
      use m_steps
      use m_spin
      use m_rmaxh
      use m_forces 
      use m_eo
      use m_gamma
      use m_ioxv, only: xv_file_read
      use m_projected_DOS, only: init_projected_DOS
      USE writewave, only: gamma_wavefunctions, setup_wf_kpoints
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
#endif
#ifdef TRANSIESTA
! TSS Begin
      use m_ts_options, only: read_ts_options
      use m_ts_in_siesta, only: setup_ts_kpoint_grid
! TSS End
#endif /* TRANSIESTA */

      implicit none

#ifdef MPI
      integer :: MPIerror  ! Return error code in MPI routines
#endif

      real(dp):: veclen     ! Length of a unit-cell vector


      integer :: i, is, ispin
      integer :: neigmin     ! Min. number of eigenstates (per k point)
      integer :: ns          ! Number of species

      logical :: user_wants_supercell

!------------------------------------------------------------------------- BEGIN
! Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#endif

      IOnode = (Node .eq. 0)

! Print version information ...........................................
      if (IOnode) then
        call prversion()
#ifdef MPI
        if (Nodes.gt.1) then
          write(6,'(/,a,i4,a)')
     .        '* Running on ', Nodes, ' nodes in parallel'
        else
          write(6,'(/,a,i4,a)')
     .        '* Running in serial mode with MPI'
        endif
#else
        write(6,'(/,a,i4,a)')
     .        '* Running in serial mode'
#endif
         call timestamp('Start of run')
         call wallclock('Start of run')
      endif
! ..................

! Start time counter ..................................................
      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

! Nullify arrays
      nullify(Haux,Saux,psi)
      nullify(listhptr)
      nullify(numh)
      nullify(Dscf,Dold,Eold,Escf)
      nullify(listh)
      nullify(H,H0,S)
      nullify(xijo)

! Initialize some variables
      DUext = 0.0_dp
      Eharrs = 0.0_dp
      Eharrs1 = 0.0_dp
      Eions = 0.0_dp
      Ekinion = 0.0_dp
      Elast = 0.0_dp
      Emad = 0.0_dp
      Emeta = 0.0_dp
      Emm = 0.0_dp
      Entrop = 0.0_dp
      Entropy = 0.0_dp
      FreeE = 0.0_dp

! Initialise read .....................................................
      call reinit(sname)

      call siesta_cml_init() ! Initialize CML (relies on reinit)

! Set allocation report level .........................................
      call fdf_global_get(level, 'alloc_report_level', 0 )
      call alloc_report( level=level, file=trim(slabel)//'.alloc' )
! ..................

! Initialise exchange-correlation functional information
      call setXC()

! Initialise force field component
      call inittwobody()

! Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom(ns)
      call broadcast(ns)
      call broadcast_basis()

      call fdf_global_get(atmonly,'Atom-Setup-Only',.false.)
      if (atmonly) call bye("End of atom setup")

! Read geometry
      call struct_init()      ! Sets na_u, isa, ucell

! Initialize atom lists

      call initatomlists()    ! Sets iza

! early exit if only checking the structure

      call fdf_global_get(struct_only,'Output-Structure-Only',.false.)
      if (IONode) then
         call write_struct( ucell, na_u, isa, iza, xa )
         if (lUseZmatrix) then
            call write_canonical_ucell_and_Zmatrix(
     $                        filename="OUT.UCELL.ZMATRIX.INITIAL")
         endif
      endif
      if (struct_only) then
         call bye("End of structure processing")
      endif
      
!  --   End of Initial Structure Processing --

      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     .    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     .                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     .                           'used, including default values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

! Allocate other arrays based on read sizes ................

      nullify(fa,cfa)
      call re_alloc(fa,1,3,1,na_u,name="fa",
     $                routine="siesta_init")
      call re_alloc(cfa,1,3,1,na_u,name="cfa",
     $                routine="siesta_init")

      call init_spin()

! Read simulation data ................................................
      call read_options( na_u, ns, nspin )

#ifdef TRANSIESTA
! Read TS Input options
      call read_ts_options()
! TSS End
#endif /* TRANSIESTA */


      qtot = qtot - charnet     ! qtot set in initatomlists
                                ! charnet set in redata
      if (IOnode) then
         write(6,fmt="(a,f12.6)") "Total number of electrons: ", qtot
         write(6,fmt="(a,f12.6)") "Total ionic charge: ", zvaltot
      endif
!
!     Warn the user: if not doing a direct optimization, the Zmatrix
!     coordinates are no longer updated. Only coordinates are treated.
!
      if (lUseZmatrix) then
         if (idyn .ne. 0) then
            write(6,"(a)")
     $      " *** WARNING: Zmatrix form will be used only for input !!"
            write(0,"(a)")
     $      " *** WARNING: Zmatrix form will be used only for input !!"
         endif
      endif

! Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (nspin .ne. 2)
     $   call die('siesta: ERROR: ' //
     $        'You can only fix the spin of the system' //
     $        ' for collinear spin polarized calculations.')
        do i = 1,2
          qs(i) = (qtot + (3-2*i)*ts) / 2.0_dp
        enddo
      else
        qs(1:nspin) = 0.0_dp
        if (nspin .le. 2) then
          do ispin = 1,nspin
            qs(ispin) = qtot/nspin
          enddo
        endif
      endif
! ..................

! Find maximum interaction range ......................................
      if (negl) then
        rmaxh = 2.0_dp*rmaxo
      else
        rmaxh = 2.0_dp*rmaxo + 2.0_dp*rmaxkb
      endif
! ......................

! Madelung correction for charged systems .............................
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

! Parallel initialisation
      call initparallel(no_u,na_u,lasto,xa,ucell,rmaxh,rcoor,isolve)
      ! if (IOnode) call show_distribution()

! Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

! Initialise arrays
      call re_alloc(listhptr,1,no_l,name='listhptr',routine='siesta',
     .              copy=.false.)
      call re_alloc(numh,1,no_l,name='numh',routine='siesta',
     .              copy=.false.)
      listhptr(1:no_l) = 0
      numh(:) = 0

! Get number of eigenstates that need to be calculated

      call fdf_global_get(neigwanted,'NumberOfEigenStates',no_u)


! Check number of eigenstates - cannot be larger than number of
! basis functions or smaller than number of occupied states + 1
! so that the Fermi level can be estimated
      do is = 1,nspin
        neigmin = nint(qs(is)/real(3 - min(nspin,2), kind=dp)) + 1
        neigwanted = max(neigwanted,neigmin)
      enddo
      neigwanted = min(neigwanted,no_u)

      ! Find k-grid for Brillouin zone integration 
      ! NOTE: We need to know whether gamma is .true. or
      ! not early, in order to decide whether to use an 
      ! auxiliary supercell for the calculation of matrix elements.

      call setup_Kpoint_grid( ucell )
      gamma = gamma_scf

      ! Call initialisation of PDOS here since we need to check if 
      ! the auxiliary supercell is needed for a non-gamma calculation
      call init_projected_DOS()
      if (do_pdos) then
        gamma = gamma .and. gamma_pdos
      endif

#ifdef TRANSIESTA
! TSS Begin
      call setup_ts_kpoint_grid( ucell )
! TSS End
#endif /* TRANSIESTA */

      nullify(eo,qo)
      call re_alloc(eo,1,no_u,1,nspin,1,maxk,name="eo",
     $              routine="state_init")
      call re_alloc(qo,1,no_u,1,nspin,1,maxk,name="qo",
     $              routine="state_init")

      call setup_bands()
      gamma = gamma .and. gamma_bands

      call setup_wf_kpoints()
      gamma = gamma .and. gamma_wavefunctions

      call estimate_pol_kpoints(ucell)
      gamma = gamma .and. gamma_polarization
!!      print *, "gamma, gamma_pol: ", gamma, gamma_polarization
!
!     User can request that the calculation is done with an explicit
!     auxiliary supercell and hermitian version, even if using only
!     the gamma point
!
      call fdf_global_get(user_wants_supercell,'ForceAuxCell',.false.)
      if (user_wants_supercell) gamma = .false.

! Find required supercell
! 2*rmaxh is used to guarantee that two given orbitals in the
! supercell can only overlap once

      if (gamma) then
         nsc(1:3) = 1
      else
         do i=1,3
            veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
            nsc(i) = ceiling( 2 * rmaxh / veclen )
         enddo
         if (.not. naiveauxcell)
     $        call check_sc_factors(ucell,nsc,2*rmaxh)
      endif

      mscell = 0.0_dp
      do i = 1, 3
         mscell(i,i) = nsc(i)
         nscold(i) = nsc(i)
      enddo


! Find auxiliary supercell (required only for k sampling) ............
      call superc( ucell, scell, nsc)

! Initialise metadynamic forces if required
      call initmeta()

      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .ge. 1 .and. idyn .le. 5) then
        inicoor = istart
        fincoor = ifinal
      else if (idyn .eq. 6) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else if (idyn .eq. 7) then
        call phonon_setup()
        inicoor = 1
        fincoor = phonon_num_disps
      else if (idyn .eq. 8) then
        inicoor = 0
        fincoor = huge(1)
      else
         call die('siesta: wrong idyn')
      endif

! Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. xv_file_read)) 
     .    call vmb(na_u,tempinit,amass,xa,isa,va)
! ..................

      istp = 0
      call timer( 'Setup', 2 )

! Output memory use before main loop
      call printmemory( 6, 0 )

! Initialization now complete. Flush stdout.
      if (ionode) call pxfflush(6)
!--------------------------------------------------------------------------- END

      END subroutine siesta_init

      End MODULE m_siesta_init
