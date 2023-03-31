! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_state_init

      private
      public :: state_init

      CONTAINS

      subroutine state_init( istep )
      use Kpoint_grid,       only: setup_Kpoint_grid, nkpnt, gamma_scf
      use m_os,              only: file_exist
      use m_new_dm,          only: new_dm
      use m_proximity_check, only: proximity_check
      use siesta_options
      use units, only: Ang
      use sparse_matrices, only: maxnh, numh, listh, listhptr
      use sparse_matrices, only: Dold, Dscf, DM_2D
      use sparse_matrices, only: Eold, Escf, EDM_2D
      use sparse_matrices, only: Hold, H, H_2D
      use sparse_matrices, only: xijo, xij_2D
      use sparse_matrices, only: S, S_1D

      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_dftu_2D, H_so_2D

      use sparse_matrices, only: sparse_pattern
      use sparse_matrices, only: block_dist, single_dist
      use sparse_matrices, only: DM_history


      use create_Sparsity_SC, only: crtSparsity_SC
      use m_sparsity_handling, only: SpOrb_to_SpAtom
      use m_sparsity_handling, only: Sp_to_Spglobal
      use m_pivot_methods, only: sp2graphviz

      use siesta_geom
      use atomlist,          only: iphorb, iphkb, indxua, iaorb,
     &                             rmaxo, rmaxkb, rmaxv, rmaxdftu,
     &                             lastkb, lasto, superc, indxuo,
     &                             no_u, no_s, no_l, iza
      use alloc,             only: re_alloc, de_alloc, alloc_report
      use m_hsparse,         only: hsparse
      use m_overlap,         only: overlap
      use m_supercell,       only: exact_sc_ag
      use siesta_cml,        only: cml_p, cmlStartStep, mainXML
      use siesta_cml,        only: cmlStartPropertyList
      use siesta_cml,        only: cmlEndPropertyList
      use siesta_cml,        only: cmlAddProperty
      use zmatrix,           only: lUseZmatrix, write_zmatrix
      use m_energies,        only: Emad
      use write_subs
      use m_ioxv,            only: ioxv
      use m_steps
      use parallel,          only: IOnode, node, nodes
      use m_spin,            only: spin
      use m_rmaxh
      use m_mixing,          only: mixers_history_init
      use m_mixing_scf,      only: scf_mixs, scf_mix

      use m_normalize_dm, only: normalize_dm

      use m_eo
      use files,             only: slabel
      use m_mpi_utils,       only: globalize_or
      use m_mpi_utils,       only: globalize_sum
      use domain_decom,      only: domainDecom, use_dd, use_dd_perm
      use dftu_specs,        only: switch_dftu, dftu_init
      use fdf,               only: fdf_get
      use sys,               only: message, die, bye
      use m_sparse, only : xij_offset

      use m_ts_kpoints, only: setup_ts_kpoint_grid
      use ts_dq_m, only : TS_DQ_METHOD, TS_DQ_METHOD_FERMI
      use m_ts_options, only : BTD_method
      use m_ts_options, only : TS_Analyze
      use m_ts_options, only : N_Elec, Elecs, IsVolt
      use m_ts_electype
      use m_ts_global_vars, only: TSrun, TSmode, onlyS
      use m_ts_io, only : fname_TSHS, ts_write_tshs
      use m_ts_sparse, only : ts_sparse_init
      use m_ts_tri_init, only : ts_tri_init, ts_tri_analyze
      use files, only: slabel, label_length
#ifdef CDF
      use iodm_netcdf, only: setup_dm_netcdf_file
      use iodmhs_netcdf, only: setup_dmhs_netcdf_file
#endif
#ifdef NCDF_4
      use dictionary
      use m_ncdf_siesta, only : cdf_init_file, cdf_save_settings
      use m_ncdf_siesta, only : cdf_save_state, cdf_save_basis
#ifdef MPI
      use mpi_siesta
#endif
#endif

      use class_Sparsity
      use class_dSpData1D
      use class_dSpData2D
      use class_dData2D
      use class_Pair_Geometry_dSpData2D
      use class_Fstack_Pair_Geometry_dSpData2D

#ifdef TEST_IO
      use m_test_io
#endif
#ifdef SIESTA__FLOOK
      use siesta_dicts, only : dict_repopulate_MD
      use siesta_dicts, only : dict_repopulate_sparse_matrices
#endif

      use m_handle_sparse, only: correct_supercell_SpD
      use m_restruct_SpData2D, only: restruct_dSpData2D

      implicit none

      integer            :: istep, nnz
      real(dp)           :: veclen      ! Length of a unit-cell vector
      real(dp)           :: rmax
      logical            :: cell_can_change
      integer            :: i, ix, iadispl, ixdispl
      logical            :: auxchanged   ! Auxiliary supercell changed?
      logical            :: folding, folding1
      logical            :: diag_folding, diag_folding1
      logical            :: foundxv  ! dummy for call to ioxv

      external           ::  madelung, timer
      real(dp), external :: volcel

      integer                       :: ts_kscell_file(3,3) = 0
      real(dp)                      :: ts_kdispl_file(3) = 0.0
      logical                       :: ts_Gamma_file = .true.
      character(len=label_length+6) :: fname
      real(dp)                      :: dummyef=0.0, dummyqtot=0.0
      type(Sparsity) :: g_Sp
#ifdef NCDF_4
      type(dictionary_t) :: d_sav
#ifdef MPI
      integer :: MPIerror
#endif
#endif
      
      character(len=256)            :: oname

      type(dData2D) :: tmp_2D

!     Variables required to correct the DM in the history
      type(dSpData2D), pointer :: tmp_Sp2D
      type(Pair_Geometry_dSpData2D), pointer :: pair

      real(dp) :: dummy_qspin(8)

!------------------------------------------------------------------- BEGIN
      call timer( 'IterGeom', 1 )

      call timer( 'state_init', 1 )

      istp = istp + 1

      if (IOnode) then
      write(6,'(/,t22,a)') repeat('=',36)
      select case (idyn)
      case (0)
         if (nmove == 0) then
            write(6,'(t25,a)') 'Single-point calculation'
            if (cml_p) call cmlStartStep(mainXML, type='Single-Point',
     $           index=istp)
         else
            if (broyden_optim) then
               write(6,'(t25,a,i6)') 'Begin Broyden opt. move = ',
     $              istep
            else if (fire_optim) then
               write(6,'(t25,a,i6)') 'Begin FIRE opt. move = ',
     $              istep
            else
               write(6,'(t25,a,i6)') 'Begin CG opt. move = ',
     $              istep
            end if
            if (cml_p) call cmlStartStep(mainXML, type='Geom. Optim',
     $           index=istp)
         endif
         
!        Print Z-matrix coordinates
         if (lUseZmatrix) then
            call write_Zmatrix()
         endif
      case (1, 3)
         if (iquench > 0 ) then
            write(6,'(t25,a,i6)') 'Begin MD quenched step = ',
     $           istep
            if (cml_p) call cmlStartStep(mainXML, type='MD-quenched',
     $           index=istep)
         else
            write(6,'(t25,a,i6)') 'Begin MD step = ',
     $           istep
            if (cml_p) call cmlStartStep(mainXML, type='MD',
     $             index=istep)
         endif
      case (2,4,5)
         write(6,'(t25,a,i6)') 'Begin MD step = ', istep
         if (cml_p) call cmlStartStep(mainXML, type='MD', index=istep)
      case (6)
          write(6,'(t25,a,i6)') 'Begin FC step = ',istep
          if (cml_p) call cmlStartStep(mainXML, type='FC', index=istep)
          
          if (istep .eq. 0) then
            write(6,'(t25,a)') 'Undisplaced coordinates'
         else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(t26,a,i0,/,t26,a,i1,a,f10.6,a)') 'displace atom ',
     &           iadispl,'in direction ',ixdispl,' by', dx/Ang,' Ang'
         endif
         
      case (8)
         write(6,'(t25,a,i6)') 'Begin Server step = ',istep
         if (cml_p) call cmlStartStep(mainXML, type='FS', index=istep)
         
      case (9)
         if ( istep == 0 ) then
            write(6,'(t25,a,i7)')'Explicit coord. initialization'
         else
            write(6,'(t25,a,i7)')'Explicit coord. step =',istep
         end if
         if (cml_p) call cmlStartStep(mainXML, type='ECS', index=istep)
         
      case (10)
         write(6,'(t25,a,i7)')'LUA coord. step =',istep
         if (cml_p) call cmlStartStep(mainXML, type='LUA', index=istep)
         
      end select
      
      write(6,'(t22,a)') repeat('=',36)
      
!     Print atomic coordinates
      call outcoor( ucell, xa, na_u, ' ', writec )

!     Save structural information in crystallographic format
!     (in file SystemLabel.STRUCT_OUT),
!     canonical Zmatrix (if applicable), and CML record
      
      call siesta_write_positions(moved=.false.)
       
      endif ! IONode

      ! Write the XV file for single-point calculations, so that
      ! it is there at the end for those users who rely on it
      call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &           foundxv)

!     Actualize things if variable cell
!     These checks are made to ensure that the k-points
!     and Madelung terms are corrected in case the cell is changed.
      auxchanged = .false.
      cell_can_change = ( varcel .or.
     &                    (idyn .eq. 8)  ! Force/stress evaluation
     &                  )
      if (change_kgrid_in_md) then
         cell_can_change = cell_can_change .or.
     &                     (idyn .eq. 3)   .or. ! Parrinello-Rahman
     &                     (idyn .eq. 4)   .or. ! Nose-Parrinello-Rahman
     &                     (idyn .eq. 5)        ! Anneal
      endif

      if ( cell_can_change .and. istep /= inicoor ) then

!       Madelung correction for charged systems 
        if (charnet .ne. 0.0_dp) then
          call madelung(ucell, shape, charnet, Emad)
        end if

        if ( .not. gamma_scf ) then

!         Will print k-points also
          call setup_Kpoint_grid( ucell )
          
          call setup_ts_kpoint_grid( ucell )
          
          call re_alloc( eo, 1, no_u, 1, spin%spinor, 1, nkpnt, 'eo',
     &        'state_init')
          call re_alloc( qo, 1, no_u, 1, spin%spinor, 1, nkpnt, 'qo',
     &        'state_init' )
        end if
        
      end if
!     End variable cell actualization

!     Always (in case of auxiliary cell) check the new sparse pattern
      if ( use_aux_cell ) then
!       Determine the tightest auxiliary supercell using
!       also the atomic positions 
        call exact_sc_ag(negl,ucell,na_u,isa,xa,nsc)
        
        mscell = 0.0_dp
        do i = 1, 3
          mscell(i,i) = nsc(i)
          if (nsc(i) /= nsc_old(i)) auxchanged = .true.
        end do
        
      end if
      
!     Auxiliary supercell
!     Do not move from here, as the coordinates might have changed
!     even if not the unit cell
      call superc(ucell, scell, nsc)
#ifdef SIESTA__FLOOK
      call dict_repopulate_MD()
#endif

!     Print unit cell and compute cell volume
!     Possible BUG: 
!     Note that this volume is later used in write_subs and the md output
!     routines, even if the cell later changes.
      if (IOnode) call outcell(ucell)
      volume_of_some_cell = volcel(ucell)

!     Use largest possible range in program, except hsparse...
!     2 * rmaxv: Vna overlap
!     rmaxo + rmaxkb: Non-local KB action
!     2 * (rmaxo + rmaxdftu): Interaction through DFTU projector
!     2.0_dp * (rmaxo+rmaxkb) : Orbital interaction through KB projectors
      rmax = max( 2._dp*rmaxv, 2._dp*(rmaxo+rmaxdftu), rmaxo+rmaxkb)

      if ( .not. negl ) then
        rmax = max(rmax, 2.0_dp * (rmaxo+rmaxkb) )
      endif

!     Check if any two atoms are unreasonably close
      call proximity_check(rmax)

      ! Clear history of mixing parameters
      call mixers_history_init( scf_mixs )
      scf_mix => scf_mixs(1)

      ! Ensure sparsity pattern is empty
      call delete(sparse_pattern)
      ! sadly deleting the sparse pattern does not necessarily
      ! mean that the arrays are de-associated.
      ! Remember that the reference counter could (in MD)
      ! be higher than 1, hence we need to create "fake"
      ! containers and let the new<class> delete the old
      ! sparsity pattern
      nullify(numh,listhptr,listh)
      allocate(numh(no_l),listhptr(no_l))
      ! We do not need to allocate listh
      ! that will be allocated in hsparse

!     List of nonzero Hamiltonian matrix elements
!     and, if applicable,  vectors between orbital centers

!     Listh and xijo are allocated inside hsparse
!     Note: We always generate xijo now, for COOP and other
!           analyses.
      call delete(xij_2D) ! as xijo will be reallocated
      nullify(xijo)
      call hsparse( negl, scell, nsc, na_s, isa, xa, lasto,
     &              lastkb, iphorb, iphKB, maxnh, .not. use_aux_cell,
     $              set_xijo=.true., folding=folding1,
     $              diagonal_folding=diag_folding1,
     $              debug_folding=fdf_get('debug-folding',.false.))
!
      call globalize_or(diag_folding1,diag_folding)
      call globalize_or(folding1,folding)
      if (diag_folding .and. .not. use_aux_cell ) then
         call message("WARNING","Gamma-point calculation " //
     $         "with interaction between periodic images")
         call message("WARNING",
     $           "Some features might not work optimally:")
         call message("WARNING",
     $           "e.g. DM initialization from atomic data")
         if (harrisfun) call die("Harris functional run needs " //
     $                              "'force-aux-cell T'")

      else if (folding) then
         if ( .not. use_aux_cell ) then
               call message("INFO","Gamma-point calculation " //
     $                      "with multiply-connected orbital pairs")
               call message("INFO",
     $              "Folding of H and S implicitly performed")
               call check_cohp()
         else
            write(6,"(a,/,a)") "Non Gamma-point calculation " //
     $           "with multiply-connected orbital pairs " //
     $           "in auxiliary supercell.",
     $           "Possible internal error. " //
     $           "Use 'debug-folding T' to debug."
            call die("Inadequate auxiliary supercell")
         endif
      endif
!
      call globalize_sum(maxnh,nnz)
      if (cml_p) then
         call cmlStartPropertyList(mainXML,title='Orbital info')
         call cmlAddProperty(xf=mainXML, value=no_u,
     $        title='Number of orbitals in unit cell',
     $        dictref='siesta:no_u', units="cmlUnits:countable")
         call cmlAddProperty(xf=mainXML, value=nnz,
     $        title='Number of non-zeros',
     $        dictref='siesta:nnz', units="cmlUnits:countable")
         call cmlEndPropertyList(mainXML)
      endif
      !
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
         ! We still need to re-create Julian Gale's
         ! indexing for O(N) in parallel.
         print "(a5,i3,a20,3i8)",
     $         "Node: ", Node, "no_u, no_l, maxnh: ", no_u, no_l, maxnh
         call setup_ordern_indexes(no_l, no_u, Nodes)
      endif

      ! I would like to skip this alloc/move/dealloc/attach
      ! by allowing sparsity to have pointer targets.
      ! However, this poses a problem with intel compilers,
      ! as it apparently errors out when de-allocating a target pointer
      write(oname,"(a,i0)") "sparsity for geom step ", istep
      call newSparsity(sparse_pattern,no_l,no_u,maxnh,
     &     numh,listhptr,listh, name = oname)
      deallocate(numh,listhptr,listh)
      call attach(sparse_pattern, 
     &     n_col = numh, list_ptr = listhptr, list_col = listh )

      ! In case the user requests to create the connectivity graph
      if ( write_GRAPHVIZ > 0 ) then
         ! first create the unit-cell sparsity pattern
         call crtSparsity_SC(sparse_pattern, g_Sp, UC=.true.)
         ! next move to global sparsity pattern
         call Sp_to_Spglobal(block_dist, g_Sp, g_Sp)
         if ( IONode ) then
            if ( write_GRAPHVIZ /= 2 )
     &           call sp2graphviz(trim(slabel)//'.ORB.gv', g_Sp)
            ! Convert to atomic 
            if ( write_GRAPHVIZ /= 1 ) then
               call SpOrb_to_SpAtom(single_dist,g_Sp,na_u,lasto,g_Sp)
               call sp2graphviz(trim(slabel)//'.ATOM.gv', g_Sp)
            end if
         end if
         call delete(g_Sp)
      end if

             
      ! Copy over xijo array (we can first do it here... :( )
      call newdData2D(tmp_2D,xijo,'xijo')
      deallocate(xijo)
      write(oname,"(a,i0)") "xijo at geom step ", istep
      call newdSpData2D(sparse_pattern,tmp_2D,block_dist,xij_2D, 
     &     name=oname)
      call delete(tmp_2D) ! decrement container...
      xijo => val(xij_2D)

!     Convert the xijo array into a super cell offset array
!     isc_off (located in siesta_geom).
!     This array is primarily used in TranSiesta but may easily
!     be used elsewhere to figure out orbital/atomic placements
!     in the sparsity pattern.
      if ( .not. use_aux_cell ) then
         ! Here we create the super-cell offsets
         call re_alloc(isc_off,1,3,1,1)
         isc_off(:,:) = 0
      else
         call xij_offset(ucell,nsc, na_u,xa,lasto, 
     &        xij_2D, isc_off, Bcast=.true.)
      end if

      
      ! When the user requests to only do an analyzation, we can call
      ! appropriate routines and quit
      if ( TS_Analyze ) then

         ! Force the creation of the full sparsity pattern
         call ts_sparse_init(slabel,IsVolt, N_Elec, Elecs, 
     &        ucell, nsc, na_u, xa, lasto, block_dist, sparse_pattern, 
     &        .not. use_aux_cell, isc_off)

         ! create the tri-diagonal matrix
         call ts_tri_analyze( block_dist, sparse_pattern , N_Elec,
     &        Elecs, ucell, na_u, lasto, nsc, isc_off,
     &        BTD_method )

         ! Print-out timers
         call timer('TS-analyze',3)

         ! Bye also waits for all processors
         call bye('transiesta analyzation performed')
         
      end if


      write(oname,"(a,i0)") "EDM at geom step ", istep
      call newdSpData2D(sparse_pattern,spin%EDM,block_dist,EDM_2D,
     &     name=oname)
      !if (ionode) call print_type(EDM_2D)
      Escf => val(EDM_2D)

      call re_alloc(Dold,1,maxnh,1,spin%DM,name='Dold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      call re_alloc(Hold,1,maxnh,1,spin%H,name='Hold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      if ( converge_EDM ) then
         call re_alloc(Eold,1,maxnh,1,spin%EDM,name='Eold',
     .     routine='sparseMat',copy=.false.,shrink=.true.)
      end if

!     Allocate/reallocate storage associated with Hamiltonian/Overlap matrix
      write(oname,"(a,i0)") "H at geom step ", istep
      call newdSpData2D(sparse_pattern,spin%H,block_dist,H_2D,
     &                  name=oname)
      !if (ionode) call print_type(H_2D)
      H => val(H_2D)

      write(oname,"(a,i0)") "H_vkb at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,H_vkb_1D,name=oname)
      !if (ionode) call print_type(H_vkb_1D)

      write(oname,"(a,i0)") "H_kin at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,H_kin_1D,name=oname)
      !if (ionode) call print_type(H_kin_1D)

      if ( switch_dftu ) then
         write(oname,"(a,i0)") "H_dftu at geom step ", istep
         call newdSpData2D(sparse_pattern,spin%spinor,
     &        block_dist,H_dftu_2D,name=oname)
         
         ! Initialize to 0, LDA+U may re-calculate
         !   this matrix sporadically doing the SCF.
         ! Hence initialization MUST be performed upon
         ! re-allocation.
         call init_val(H_dftu_2D)
         if ( inicoor /= istep ) then
            ! Force initialization of the LDA+U
            ! when changing geometry
            ! For the first geometry this is controlled
            ! by the user via an fdf-key
            dftu_init = .true.
         end if
      end if
      
      if ( spin%SO ) then
         write(oname,"(a,i0)") "H_so at geom step ", istep
         call newdSpData2D(sparse_pattern,spin%H - 2,
     &        block_dist,H_so_2D,name=oname)
      end if

      write(oname,"(a,i0)") "S at geom step ", istep
      call newdSpData1D(sparse_pattern,block_dist,S_1D,name=oname)
      if (ionode) call print_type(S_1D)
      S => val(S_1D)

!>    Before proceeding we need to "fix" a few things before we can
!>    successfully read the new DM.
!>    1. Cases where the atomic displacements yields a new sparse
!>    pattern it is vital to remove the elements that does not exist.
!>    Such cases are often encountered because atoms move in/out of
!>    neighbouring atoms orbital ranges. Everytime two orbitals meet,
!>    new sparse elements are added, and everytime two orbitals flee
!>    sparse elements are removed.
!>    2. If the supercell has changed size it is necessary to
!>    remove/translate the old/new supercells such that they are
!>    conforming with the new sparse pattern.
!>    For details regarding the sparsity pattern, see sparse_matrices.F90
      do i = 1, n_items(DM_history)
        pair => get_pointer(DM_history,i)
        call secondp(pair,tmp_Sp2D)
        if ( auxchanged ) then
!         Transfer the old mixing matrix to the new supercell pattern
          call correct_supercell_SpD(nsc_old, tmp_Sp2D, nsc)
        end if
!       Only retain the new orbital interactions
        call restruct_dSpData2D(tmp_Sp2D, sparse_pattern, tmp_Sp2D,
     &      show_warning = .false.)
      end do

      
!     Find overlap matrix 
      call overlap( na_u, na_s, no_s, scell, xa, indxua, rmaxo, maxnh,
     &              lasto, iphorb, isa, numh, listhptr, listh, S )

#ifdef NCDF_4
!     At this point the sparsity pattern, overlap matrix and some
!     other details.
!     Before continuing we will create the CDF file output
!
!     Initialize the NC file
      if ( write_cdf ) then
        if ( inicoor == istep ) then
          call cdf_init_file(trim(slabel)//'.nc', is_FC=(idyn==6))
#ifdef MPI
          call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

!         Save the basis set (only once)
          call cdf_save_basis(trim(slabel)//'.nc')
#ifdef MPI
          call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

          d_sav = d_sav//('xa'.kv.1)//('cell'.kv.1)
          d_sav = d_sav//('sp'.kv.1)//('S'.kv.1)//('xij'.kv.1)
          d_sav = d_sav//('isc_off'.kv.1)//('nsc'.kv.1)
          call cdf_save_state(trim(slabel)//'.nc',d_sav)
          call delete(d_sav)
          
        end if
      end if
#endif

!     
!     Here we could also read a Hamiltonian, either to proceed to
!     the analysis section (with nscf=0) or to start a mix-H scf cycle.
!
!     Initialize density matrix
!     The resizing of Dscf is done inside new_dm
      call new_DM(auxchanged, DM_history, DM_2D, EDM_2D)
      Dscf => val(DM_2D)
      Escf => val(EDM_2D)
      if (spin%H > 1) call print_spin(dummy_qspin)

!     The number of old supercells is used in new_DM. Hence we may first
!     update the oldnsc here.
      if ( auxchanged ) nsc_old(:) = nsc(:)

      ! Initialize energy-density matrix to zero for first call to overfsm
      ! Only part of Escf is updated in TS, so if it is put as zero here
      ! a continuation run gives bad forces.
      if ( .not. TSrun ) then
         call normalize_DM( first= .true. )
         Escf(:,:) = 0.0_dp
      end if

#ifdef TEST_IO
      ! We test the io-performance here
      call time_io(spin%H,nsc,H_2D)
#endif

#ifdef SIESTA__FLOOK
      call dict_repopulate_sparse_matrices()
#endif
      
!     Write out the ORB_INDX file (if requested)
      if ( save_ORB_INDX ) then
        if ( IOnode .and. istep == inicoor ) then
          call write_orb_indx( na_u, na_s, no_u, no_s, isa, xa,
     .        iaorb, iphorb, indxuo, nsc, ucell )
        end if
      end if
      
!     If onlyS, Save overlap matrix and exit
      if ( onlyS ) then
         fname = fname_TSHS(slabel, onlyS = .true. )
         ! We include H as S, well-knowing that we only write one of
         ! them, there is no need to allocate space for no reason!
         call ts_write_tshs(fname, 
     &        .true., .not. use_aux_cell, ts_Gamma_file,
     &        ucell, nsc, isc_off, na_u, no_s, spin%H,
     &        ts_kscell_file, ts_kdispl_file,
     &        xa, lasto, 
     &        H_2D, S_1D, indxuo, 
     &        dummyEf, dummyQtot, Temp,0,0)
         call bye( 'Save overlap matrix and exit' ) ! Exit siesta
         
      endif 

#ifdef CDF
      if (writedm_cdf) then
         call setup_dm_netcdf_file( maxnh, no_l, spin%H,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh)
      endif
      if (writedm_cdf_history) then
         call setup_dm_netcdf_file( maxnh, no_l, spin%H,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              istep)
      endif
      if (writedmhs_cdf) then
         call setup_dmhs_netcdf_file( maxnh, no_l, spin%H,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              s)
      endif
      if(writedmhs_cdf_history) then
         call setup_dmhs_netcdf_file( maxnh, no_l, spin%H,
     &                              no_s, indxuo,
     &                              numh, listhptr, listh,
     &                              s,
     &                              istep)
      endif
#endif
      call timer( 'state_init', 2 )

      END subroutine state_init

      subroutine check_cohp()
      use siesta_options, only: write_coop
      use sys,            only: message
      
      if (write_coop) then
         call message("WARNING","There are multiply-connected "//
     $                          "orbitals.")
         call message("WARNING","Your COOP/COHP analysis might " //
     $                          "be affected by folding.")
         call message("WARNING",'Use "force-aux-cell T "' //
     $                          'or k-point sampling')
      endif
      end subroutine check_cohp

      END module m_state_init
