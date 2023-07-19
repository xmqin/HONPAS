      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      USE Kpoint_grid, only: setup_Kpoint_grid, maxk
      USE m_new_dm, only: new_dm
      USE m_proximity_check, only: proximity_check
      USE siesta_options
      use sparse_matrices, only: nh, maxnh, numh, listh, listhptr
      use sparse_matrices, only: Dscf, Dold, Escf, Eold
      use sparse_matrices, only: xijo, H, H0, S
      use siesta_geom
      use atomlist, only: iphorb, iphkb, rco, rckb, indxua,no_s, rmaxkb,
     .                    rmaxo, no_u, lastkb, lasto, superc, rmaxv
      use atomlist, only: no_s, indxuo
      use atomlist, only: no_l
      use alloc, only: re_alloc, alloc_report, de_alloc
      use m_hsparse,    only: hsparse
      use m_overlap,    only: overlap
      use m_check_supercell, only: check_sc_factors
      use siesta_cml,   only: cml_p, cmlStartStep, mainXML
      use zmatrix,      only: lUseZmatrix, write_zmatrix
      use m_energies, only: Emad
      use write_subs
      use m_steps
      use parallel, only: IOnode, Nodes
      use m_spin, only: nspin
      use m_rmaxh
      use m_pulay, only: init_pulay_arrays
      use m_eo
      use m_gamma
      use files,          only : slabel
      use m_mpi_utils, only: globalize_or
      use parallel, only: Node
      use domain_decom, only: domainDecom, use_dd, use_dd_perm

#ifdef TRANSIESTA
      use m_ts_options, only : onlyS, mixH
      use m_ts_global_vars, only : Hold
      use sys, only : bye
      use m_ts_io, only : ts_iohs
      use files, only: slabel, label_length
#endif
#ifdef CDF
      use iodm_netcdf, only: setup_dm_netcdf_file
      use iodmhs_netcdf, only: setup_dmhs_netcdf_file
#endif

      implicit none

      integer :: istep

      real(dp):: veclen     ! Length of a unit-cell vector
      real(dp):: rmax
      logical :: cell_can_change

      integer :: i, ix, iadispl, ixdispl

      logical:: auxchanged ! Has auxiliary supercell changed?
      logical :: folding, folding1
      external ::  madelung, timer

      real(dp), external :: volcel
#ifdef TRANSIESTA
      integer, dimension(3,3) :: ts_kscell_file = 0
      real(dp), dimension(3) :: ts_kdispl_file = 0.0
      logical :: ts_gamma_scf_file = .true.
      character(len=label_length+6) :: fname
      integer :: fnlength
      real(dp) :: dummyef=0.0, dummyqtot=0.0
      real(dp), pointer :: dummyH(:,:)
!     SIESTA External
      character(len=label_length+5) PASTE
#endif

!------------------------------------------------------------------------- BEGIN
      call timer( 'IterMD', 1 )
      istp = istp + 1
      if (IOnode) then
        write(6,'(/2a)') 'siesta:                 ',
     .                    '=============================='
        select case (idyn)
        case (0)
          write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='CG', index=istp)
          ! Print Z-matrix coordinates
          if (lUseZmatrix) then
             call write_Zmatrix()
          endif
        case (1:5)
          write(6,'(28(" "),a,i6)') 'Begin MD step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='MD', index=istep)
        case (6)
          write(6,'(28(" "),a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC', index=istep)

          if (istep .eq. 0) then
            write(6,'(28(" "),a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(28(" "),a,i6)') 'displace atom   ', iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ', ixdispl
            write(6,'(28(" "),a,f8.4,a)') 'by       ', dx, ' Bohr'
          endif

        case (8)
          write(6,'(28(" "),a,i6)') 'Begin Server step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)
          ! Find interatomic distances (output in file BONDS)
          call bonds( ucell, na_u, isa, xa,
     $                rmax_bonds, trim(slabel)// ".BONDS" )

        end select
        write(6,'(2a)') '                        ',
     .                    '=============================='


        ! Print atomic coordinates
        call outcoor( ucell, xa, na_u, ' ', writec )
        call siesta_write_positions()
      endif

      ! Actualize things if variable cell 
      auxchanged = .false.
      cell_can_change = ( varcel .or.
     $                    (idyn .eq. 8)  ! Force/stress evaluation
     $                  )
      if (change_kgrid_in_md) then
         cell_can_change = cell_can_change .or.
     $        (idyn .eq. 3)      ! Parrinello-Rahman
     $        .or. (idyn .eq. 4) ! Nose-Parrinello-Rahman
     $        .or. (idyn .eq. 5) ! Anneal
      endif

      if ( cell_can_change
     $     .and. (istep.ne.inicoor) .and. (.not.gamma) ) then

        ! Will print k-points also
        call setup_Kpoint_grid( ucell )

        call re_alloc(eo,1,no_u,1,nspin,1,maxk,name="eo",
     $                routine="state_init")
        call re_alloc(qo,1,no_u,1,nspin,1,maxk,name="qo",
     $                routine="state_init")

        ! Find required supercell
        if (gamma) then
           nsc(1:3) = 1
        else if (fixauxcell) then
           nsc(1:3) = nscold(1:3)
        else
           do i=1,3
              veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
              nsc(i) = ceiling( 2 * rmaxh / veclen )
           enddo
           if (.not. naiveauxcell)
     $         call check_sc_factors(ucell,nsc,2*rmaxh)
        endif

        mscell = 0.0_dp
        do i = 1, 3
           mscell(i,i) = nsc(i)
           if (nsc(i).ne.nscold(i)) auxchanged = .true.
           nscold(i) = nsc(i)
        enddo

        ! Madelung correction for charged systems 
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        endif

      endif
! End variable cell actualization

! Auxiliary supercell
      ! Do not move from here, as the coordinates might have changed
      ! even if not the unit cell
      call superc(ucell, scell, nsc)

! Print unit cell and compute cell volume
! Possible BUG: 
! Note that this volume is later used in write_subs and the md output
! routines, even if the cell later changes.
      if (IOnode) call outcell(ucell)
      volume_of_some_cell = volcel(ucell)
! ...................

      ! Use largest possible range in program, except hsparse...
      ! 2 * rmaxv: Vna overlap
      ! 2 * rmaxo: orbital overlap
      ! rmaxo + rmaxkb: Non-local KB action
      ! 2.0_dp * (rmaxo+rmaxkb) : Orbital interaction through KB projectors

      rmax = max( 2._dp*rmaxv, 2._dp*rmaxo, rmaxo+rmaxkb)

      if (.not. negl) then
         rmax = max(rmax, 2.0_dp * (rmaxo+rmaxkb) )
      endif

      ! Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      ! List of nonzero Hamiltonian matrix elements
      ! and, if applicable,  vectors between orbital centers

      ! Listh and xijo are allocated inside hsparse
      ! nh is set by the routine also
      ! Note: We always generate xijo now, for COOP and other
      !       analyses.

      nh = maxnh
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto, lastkb, 
     .                 iphorb, iphKB, nh, numh, listhptr, listh,
     $                 set_xijo=.true., xijo=xijo, gamma=gamma,
     $                 folding=folding1)
      maxnh = nh
!
      call globalize_or(folding1,folding)
      if (folding) then
         if (IOnode) then
            print *, "Folding of H and S is implicitly performed"
         endif
      endif
      !
      ! If using domain decomposition, redistribute orbitals
      ! for this geometry, based on the hsparse info. 
      ! The first time round, the initial distribution is a
      ! simple block one (given by preSetOrbitLimits).
      !
      ! Any DM, etc, read from file will be redistributed according
      ! to the new pattern. 
      ! Inherited DMs from a previous geometry cannot be used if the
      ! orbital distribution changes. For now, we avoid changing the
      ! distribution (the variable use_dd_perm is .true. if domain
      ! decomposition is in effect). Names should be changed...

      if (use_dd .and. (.not. use_dd_perm)) then
         call domainDecom( no_u, no_l, maxnh )  ! maxnh intent(in) here
         maxnh = sum(numh(1:no_l))
         nh = maxnh
         ! We still need to re-create Julian Gale's
         ! indexing for O(N) in parallel.
         print "(a5,i3,a20,3i8)",
     $         "Node: ", Node, "no_u, no_l, maxnh: ", no_u, no_l, maxnh
         call setup_ordern_indexes(no_l, no_u, Nodes)
      endif

      ! Note that the resizing of Dscf is done inside new_dm below
      call re_alloc(Dold,1,maxnh,1,nspin,name='Dold',
     .     routine='state_init',copy=.false.,shrink=.true.)
      call re_alloc(Eold,1,maxnh,1,nspin,name='Eold',
     .     routine='state_init',copy=.false.,shrink=.true.)
      call re_alloc(Escf,1,maxnh,1,nspin,name='Escf',
     .     routine='state_init',copy=.false.,shrink=.true.)
#ifdef TRANSIESTA
      if(mixH) then
         call re_alloc(Hold,1,maxnh,1,nspin,name='Hold',
     .     routine='state_init',copy=.false.,shrink=.true.)
      end if
#endif /* TRANSIESTA */
      
      ! Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      call re_alloc(H,1,maxnh,1,nspin,name='H',
     .              routine='state_init',shrink=.true.,copy=.false.)
      call re_alloc(H0,1,maxnh,name='H0',routine='state_init',
     .              shrink=.true.,copy=.false.)
      call re_alloc(S,1,maxnh,name='S',routine='state_init',
     .              copy=.false.,shrink=.true.)

      ! Find overlap matrix 
      call overlap(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     .             maxnh, lasto, iphorb, isa, 
     .             numh, listhptr, listh, S)

      ! Initialize density matrix
      ! The resizing of Dscf is done inside new_dm
      call new_dm( auxchanged, numh, listhptr, listh, Dscf)

      ! Check for size of Pulay auxiliary matrices
      call init_pulay_arrays()


#ifdef TRANSIESTA  
      !If onlyS, Save overlap matrix and exit                      ! MPN 
      if (onlyS) then                                              ! MPN    
         fname = paste( slabel, '.onlyS' )
         fnlength = label_length+6
         nullify(dummyH)
         call re_alloc( dummyH, 1, maxnh, 1, 1, name='dummyH', 
     .        routine='state_init' )

         call ts_iohs('write   ', gamma, .true., no_u, no_s, 
     .        1, indxuo, maxnh, numh, listhptr, listh, dummyH, 
     .        S, dummyqtot, temp, xijo, fnlength, fname, na_u, 
     .        lasto, isa, dummyef, 
     .        ucell, ts_kscell_file, ts_kdispl_file,            
     .        ts_gamma_scf_file, xa, 0, 0)
         call de_alloc(dummyH, name='dummyH')
         call bye('Save overlap matrix and exit') ! Exit siesta    ! MPN 
      endif                                                        ! MPN    
#endif /* TSS     ! MPN    */

#ifdef CDF
      if (writedm_cdf) then
         call setup_dm_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh)
      endif
      if (writedm_cdf_history) then
         call setup_dm_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh,
     $                              istep)
      endif
      if (writedmhs_cdf) then
         call setup_dmhs_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh,
     $                              s)
      endif
      if(writedmhs_cdf_history) then
         call setup_dmhs_netcdf_file( maxnh, no_l, nspin,
     $                              no_s, indxuo,
     $                              numh, listhptr, listh,
     $                              s,
     $                              istep)
      endif
#endif

!--------------------------------------------------------------------------- END
      END subroutine state_init
      END module m_state_init





