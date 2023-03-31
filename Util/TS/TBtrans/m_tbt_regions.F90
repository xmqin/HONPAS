! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This module contains
! the different regions used in tbtrans

! Coded by Nick Papior Andersen
! 2014

module m_tbt_regions

  use precision, only: dp
  use m_region
  use class_OrbitalDistribution
  use class_Sparsity
  ! To re-use as much from transiesta
  use m_ts_method, only : r_aBuf, r_oBuf
  use m_ts_method, only : r_aDev => r_aC, r_oDev => r_oC

  implicit none

  private
  save

  ! The full tbtrans sparsity region in the unit-cell equivalent
  ! region, this is used to accomodate the Hamiltonian and overlap
  ! matrices. In principle this could be made obsolete.
  type(Sparsity), public :: sp_uc ! TBT-GLOBAL (UC)

  ! The SC sparsity pattern in the device region
  type(Sparsity), public :: sp_dev_sc

  ! the different regions that connects to the equivalent
  ! electrodes.
  ! I.e. it is the regions that goes from the electrode
  ! and down to the central region (without any central
  ! region overlap)
  type(tRgn), allocatable, target, public :: r_aEl(:), r_oEl(:)
  type(tRgn), allocatable, target, public :: r_oElpD(:)

  ! the device region (the calculated GF region)
  ! Note that these arrays are pivoted indices
  public :: r_aDev, r_oDev

  ! The buffer region, just for completeness
  public :: r_aBuf, r_oBuf

  public :: tbt_init_regions
  public :: tbt_region_options
  public :: tbt_print_regions

contains

  subroutine tbt_init_regions(N_Elec, Elecs, cell, na_u, xa, lasto, &
       dit, sp, &
       nsc, isc_off)

    use fdf
    use fdf_extra
    use parallel, only : IONode, Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World, MPI_Barrier
#endif
    use files, only: slabel

    use m_char, only: lcase
    use m_pivot
    use m_pivot_methods, only : sp2graphviz

    use geom_helper, only : iaorb
    use intrinsic_missing, only : SPC_PROJ, VNORM, VEC_PROJ
    use create_Sparsity_SC
    use create_Sparsity_Union

    use m_ts_electype
    use m_ts_method, only : atom_type, TYP_DEVICE, TYP_BUFFER

    use m_ts_sparse, only : ts_Sparsity_Global
    use m_ts_pivot, only : ts_pivot

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    
    ! Number of electrodes
    integer, intent(in) :: N_Elec
    ! electrodes
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! The device region unit-cell
    real(dp), intent(in) :: cell(3,3)
    ! Last orbital of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! The distribution for the sparsity pattern
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The supercell information
    integer, intent(in) :: nsc, isc_off(3,nsc)

    integer :: iEl, jEl

    ! A temporary sparsity pattern
    type(Sparsity) :: sp_tmp

    ! the different regions that becomes the electrodes
    type(tRgn) :: r_aEl_alone(N_Elec), r_oEl_alone(N_Elec)

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=128) :: g, csort
    integer :: i, ia1, ia2, no_u
    integer :: init_nz
    type(tRgn) :: r_tmp, r_tmp2, r_tmp3, r_Els, priority
    real(dp) :: tmp
    
    no_u = lasto(na_u)

    call timer('init-region+sp', 1)

    ! Create the sparsity pattern and remove the buffer atoms...
    if ( r_oBuf%n > 0 ) then
       sp_tmp = sp
       call Sp_remove_region(dit,sp_tmp,r_oBuf,sp)
       call delete(sp_tmp)
    end if

#ifdef TRANSIESTA_DEBUG
    open(file='NO_BUF_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Create all the "alone" electrode regions
    do iEl = 1 , N_Elec

       ! Create electrode region
       ia1 = Elecs(iEl)%idx_a
       ia2 = ia1 - 1 + TotUsedAtoms(Elecs(iEl))
       call rgn_range(r_aEl_alone(iEl), ia1, ia2)

       ia1 = Elecs(iEl)%idx_o
       ia2 = ia1 - 1 + TotUsedOrbs(Elecs(iEl))
       call rgn_range(r_oEl_alone(iEl), ia1, ia2)
       
       ! Check that we have a legal region
       if ( rgn_overlaps(r_aEl_alone(iEl),r_aBuf) ) then
          write(*,*)'Overlapping electrode: '//trim(Elecs(iEl)%name)
          call die('Buffer region overlaps with an electrode &
               &please correct your input!')
       end if

    end do

    ! Delete to be ready to populate the device
    ! Read in electrode down-folding regions
    call rgn_delete(r_aDev, r_Els)

    ! Read in device region via the block or list
    if ( fdf_islist('TBT.Atoms.Device') ) then
       
       ! Query size of list
       i = -1
       call rgn_init(r_aDev, 1)
       call fdf_list('TBT.Atoms.Device', i, r_aDev%r)
       call rgn_init(r_aDev, i)
       call fdf_list('TBT.Atoms.Device', r_aDev%n, r_aDev%r)
       
    else if ( fdf_block('TBT.Atoms.Device',bfdf) ) then

       ! read by line and set them to be buffer atoms
       do while ( fdf_bline(bfdf,pline) ) 

          ! empty line
          if ( fdf_bnnames(pline) == 0 ) cycle
       
          g = fdf_bnames(pline,1)

          ! Actual device atoms...
          if ( leqi(g,'atom') .or. leqi(g,'position') ) then
             call fdf_brange(pline, r_tmp, 1, na_u)
             if ( r_tmp%n == 0 ) &
                  call die('Could not read in any atoms &
                  &in line of TBT.Atoms.Device')
             call rgn_union(r_aDev,r_tmp,r_aDev)
             
          end if

          ! Atoms NOT in the device region...
          if ( leqi(g,'not-atom') .or. leqi(g,'not-position') .or. &
               leqi(g,'-atom') .or. leqi(g,'-position') ) then
             call fdf_brange(pline, r_tmp, 1, na_u)
             if ( r_tmp%n == 0 ) &
                  call die('Could not read in any atoms &
                  &in line of TBT.Atoms.Device')
             call rgn_union(r_Els,r_tmp,r_Els)
             
          end if
          
       end do
       call rgn_delete(r_tmp)       

    end if

    if ( r_aDev%n == 0 ) then
       ! populate the device region with all atoms
       call rgn_range(r_aDev, 1, na_u)
    end if
    ! Just to speed up the following stuff
    call rgn_sort(r_aDev)
    
    ! remove the buffer and electrode atoms
    if ( r_aBuf%n > 0 ) then
       call rgn_complement(r_aBuf,r_aDev,r_aDev)
    end if
    do iEl = 1 , N_Elec
       call rgn_complement(r_aEl_alone(iEl),r_aDev,r_aDev)
    end do

    ! remove the except region from the device region
    call rgn_complement(r_Els,r_aDev,r_aDev)
    ! Clean atomic downfolding region
    call rgn_delete(r_Els)

    if ( r_aDev%n == 0 ) then
       call die('Zero atoms are in the device region...?')
    end if

    ! Create device region
    call rgn_Atom2Orb(r_aDev,na_u,lasto,r_oDev)

    if ( IONode ) then
       write(*,'(/,a)')'tbt: Analyzing electrode sparsity &
            &pattern and electrode pivot-tables'
    end if

    ! In case the user wants "a correct DOS"
    ! in this region, we extend it
    if ( fdf_get('TBT.Atoms.Device.Connect',.false.) ) then

      ! TBTrans will truncate connections at electrode interfaces.
      call rgn_sp_connect(r_oDev, dit, sp, r_tmp)
      if ( r_tmp%n == 0 ) &
          call die('No orbitals connect to the specified device &
          &region. This is not allowed.')

      ! Convert connecting region to atoms (this also
      ! folds supercell orbitals to the correct atoms)
      call rgn_Orb2Atom(r_tmp, na_u, lasto, r_tmp2)
      call rgn_delete(r_tmp)

      ! Remove buffer atoms (in case the electrode is too small)
      if ( r_aBuf%n > 0 ) then
        call rgn_complement(r_aBuf, r_tmp2, r_tmp2)
      end if

      ! Retain number of atoms before removing electrodes
      i = r_tmp2%n

      ! Remove all electrodes from the region
      do iEl = 1 , N_Elec
        call rgn_complement(r_aEl_alone(iEl), r_tmp2, r_tmp2)
      end do

      ! Append
      call rgn_append(r_aDev, r_tmp2, r_aDev)

      ! If we have removed some atoms from the electrodes, it means
      ! there are connections to the electrodes from the device region.
      if ( IONode .and. r_tmp2%n < i ) then
        ! It only connects to electrodes
        write(*,'(a)')'tbt: Device regions &
            &connects directly with electrodes'
        write(*,'(a)')'tbt: If the overlap is large this might &
            &produce spurious effects in DOS calculations'
      end if

      call rgn_delete(r_tmp2)

      ! In its current state we force the entire atoms
      ! to be in the orbital connection scheme (even though
      ! some orbitals might not connect...)
      call rgn_Atom2Orb(r_aDev, na_u, lasto, r_oDev)

    end if

    ! Makes searching a little faster
    call rgn_sort(r_aDev)
    call rgn_sort(r_oDev)

#ifdef TRANSIESTA_DEBUG
    open(file='FULL_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Allocate the different regions
    allocate(r_aEl(N_Elec),r_oEl(N_Elec))
    allocate(r_oElpD(N_Elec))

    do iEl = 1 , N_Elec

       ! Remove the connections that cross the boundary
       ! starting from this electrode
       call rgn_sp_connect(r_oEl_alone(iEl), dit, sp, r_tmp)
       call rgn_union(r_oEl_alone(iEl), r_tmp, r_tmp2)
       ! o_inD is used in the subsequent pivoting routines
       call rgn_copy(r_oEl_alone(iEl),Elecs(iEl)%o_inD)

       select case ( Elecs(iEl)%t_dir )
       case ( 4 ) ! B-C
         i = Elecs(iEl)%pvt(2)
         call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case ( 5 ) ! A-C
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case ( 6 ) ! A-B
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(2)
       case ( 7 ) ! A-B-C
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(2)
         call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case default
         i = Elecs(iEl)%pvt(Elecs(iEl)%t_dir)
       end select
       
       ! Remove connections from this electrode across the boundary...
       call Sp_remove_crossterms(dit,sp,nsc,isc_off, i, sp, r = r_tmp2)

       ! Check that the device region does not overlap
       if ( rgn_overlaps(r_aEl_alone(iEl),r_aDev) ) then
          write(*,*)'Overlapping electrode: '//trim(Elecs(iEl)%name)
          call die('Device region overlaps with an electrode &
               &please correct your input!')
       end if

    end do

    ! Retrieve initial number of non-zero elements
    ! Note that this number of non-zero elements is _after_
    ! we have removed those terms that connect across
    ! the cell boundaries.
    init_nz = nnzs(sp)
    
    do iEl = 1 , N_Elec - 1

       ! in order to get the correct connections
       ! i.e. without connections across the device region
       ! we need to remove the connections to the other 
       ! electrodes
       ! Step 1. build a unified region of all the following
       !         electrodes!
       call rgn_delete(r_tmp)
       do jEl = iEl + 1 , N_Elec
          call rgn_append(r_tmp, r_oEl_alone(jEl), r_tmp)
       end do
       ! Speeds up stuff
       call rgn_sort(r_tmp)

       ! First we update the sparsity pattern to remove any connections
       ! between the electrode and the other ones
       ! It will NOT remove connections between the central region and
       ! the other electrodes!
       sp_tmp = sp
       call Sp_remove_region2region(dit,sp_tmp,r_oEl_alone(iEl),r_tmp,sp)
       call delete(sp_tmp)

    end do

    ! Issue a warning if there are any removed elements
    ! from the above loop
    ! That would mean that there are connections across the device region
    ! which would mean a leak in the transmission.
    i = nnzs(sp)
    if ( i < init_nz .and. IONode ) then
       write(*,'(/a,i0,a/)')'*** WARNING! Removed ',init_nz - i, ' elements &
            &which connect electrodes across the device region!'
    end if
       

#ifdef TRANSIESTA_DEBUG
    open(file='NO_ELECTRODE_CONNECTIONS_SP',unit=1400,form='formatted')
    call sp_to_file(1400,sp)
    close(1400)
#endif

    ! Create the temporary unit-cell sparsity pattern
    call crtSparsity_SC(sp, sp_tmp, UC = .true. )

    call timer('init-region+sp', 2)

    ! Create the electrode down-folding regions.
    ! Note that sorting according to the more advanced methods
    ! is not directly applicable as the methods involve non-stringent
    ! ending elements.
    call timer('pivot-elec', 1)

    ! Collect buffer and device orbitals to one list
    call rgn_append(r_oBuf, r_oDev, r_tmp)
    call rgn_sort(r_tmp)

    do iEl = 1 , N_Elec

       if ( mod(iEl-1,Nodes) /= Node ) cycle

       ! Create pivoting region (except buffer+device)
       call rgn_range(r_oEl(iEl), 1, no_u)
       call rgn_complement(r_tmp, r_oEl(iEl), r_oEl(iEl))

       ! Sort according to the connectivity of the electrode
       ! This will also reduce the pivoting table (r_oEl) to
       ! _only_ have orbitals from the electrode up to the
       ! device region.
       ! First get the actual sub region that connects
       ! from the electrode to the device
       g = 'atom+'//trim(Elecs(iEl)%name)
       csort = fdf_get('TBT.BTD.Pivot.Elecs',trim(g))
       csort = fdf_get('TBT.BTD.Pivot.Elec.'//&
            trim(Elecs(iEl)%name),trim(csort))

       ! If the electrode is in the pivoting scheme we
       ! are for sure doing a connectivity graph.
       if ( index(csort, trim(g)) > 0 ) then
          ! the requested sorting algorithm
          ! is the connectivity graph, then we can do
          ! everything in one go!
          call ts_pivot(dit, sp_tmp, &
               1, Elecs(iEl:iEl), &
               cell, na_u, xa, lasto, &
               r_oEl(iEl), csort, extend = .false.)
       else
          ! First find the connectivity graph to limit the
          ! electrode region
          ! re-create g
          g = 'atom+'//trim(Elecs(iEl)%name)
          call ts_pivot(dit, sp_tmp, &
               1, Elecs(iEl:iEl), &
               cell, na_u, xa, lasto, &
               r_oEl(iEl), g, extend = .false.)
          ! now r_oEl(iEl) contains the connetivity graph
          ! from electrode iEl (without crossing the device region)
          call ts_pivot(dit, sp_tmp, &
               1, Elecs(iEl:iEl), &
               cell, na_u, xa, lasto, &
               r_oEl(iEl), csort, extend = .false.)
       end if

       ! Print out the pivoting scheme that was used for this electrode
       if ( IONode ) then
          write(*,'(4a)')'tbt: BTD pivoting scheme for electrode (', &
               trim(Elecs(iEl)%name),'): ', trim(csort)
       end if

       if ( .not. leqi(g, csort) ) then
          ! if the majority of the electrode pivoting
          ! elements are in the upper half, then reverse
          ! Note that this will (can) only happen if the
          ! method is not the connectivity graph from it-self
          ! hence we have the above if-clause
          
          tmp = sum(rgn_pivot(r_oEl(iEl),r_oEl_alone(iEl)%r))
          ! average position in pivoting array
          tmp = tmp / r_oEl_alone(iEl)%n
          if ( tmp > 0.5_dp * r_oEl(iEl)%n ) then
             ! the actual electrode position
             ! is primarily placed closer to the device,
             ! hence we reverse the list to get a
             ! better sorting
             call rgn_reverse(r_oEl(iEl))
          end if

       end if

       ! This aligns the atoms in the same way the orbitals 
       ! introduce the atoms.
       call rgn_Orb2Atom(r_oEl(iEl), na_u, lasto , r_aEl(iEl))
       call rgn_sort(r_aEl(iEl))

       ! Create the region that connects the electrode-followed
       ! region to the central region
       call rgn_sp_connect(r_oEl(iEl), dit, sp_tmp, Elecs(iEl)%o_inD)
       ! To streamline the memory layout in the tri-diagonal
       ! blocks we might as well sort this
       call rgn_sort(Elecs(iEl)%o_inD)
       ! Append the found region that is connecting out to the device region
       call rgn_append(r_oEl(iEl), Elecs(iEl)%o_inD, r_oElpD(iEl))

       ! We now know how many orbitals that we are down-folding the 
       ! electrode self-energy to
       if ( Elecs(iEl)%o_inD%n == 0 ) then
          ! Print-out for debugging purposes
          call rgn_print(r_oDev)
          call rgn_print(r_oEl(iEl))
          call die('The electrode down-folding region is 0 in the device &
               &region. Please expand your device region.')
       end if

    end do

    ! Possibly deleting it
    call rgn_delete(r_tmp2, r_Els)
    call rgn_copy(r_oDev, r_tmp)
    call rgn_sort(r_tmp) ! copy and sort to faster find pivoting

    do iEl = 1 , N_Elec

#ifdef MPI
       i = mod(iEl-1,Nodes)
       ! Bcast the regions
       call rgn_MPI_Bcast(r_aEl(iEl),i)
       call rgn_MPI_Bcast(r_oEl(iEl),i)
       call rgn_MPI_Bcast(Elecs(iEl)%o_inD,i)
       call rgn_MPI_Bcast(r_oElpD(iEl),i)
#endif

       ! Check that the region does not overlap with any previous 
       ! electrode region...
       if ( rgn_overlaps(r_tmp2, r_aEl(iEl)) ) then
         if ( Node == 0 ) then
           do jEl = 1 , N_Elec
             ! We are dying anyway, so might as well sort to make it easier
             ! to debug
             call rgn_sort(r_aEl(jEl))
             call rgn_print(r_aEl(jEl))
           end do
         end if
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World, jEl)
#endif
         call die('Electrode regions connect across the device region, &
             &please increase your device region!')
       end if
       if ( iEl < N_Elec ) then
         call rgn_append(r_tmp2, r_aEl(iEl), r_tmp2)
         call rgn_sort(r_tmp2)
       end if

       ! Set the names
       r_aEl(iEl)%name = '[A]-'//trim(Elecs(iEl)%name)//' folding region'
       r_oEl(iEl)%name = '[O]-'//trim(Elecs(iEl)%name)//' folding region'
       r_oElpD(iEl)%name = '[O]-'//trim(Elecs(iEl)%name)//' folding El + D'
       Elecs(iEl)%o_inD%name = '[O]-'//trim(Elecs(iEl)%name)//' in D'

       ! Prepare the inDpvt array (will be filled in tri_init
       ! as we sort each block individually)
       call rgn_copy(Elecs(iEl)%o_inD, Elecs(iEl)%inDpvt)
       
       ! Check that the electrode down-folded self-energy is fully contained
       ia1 = minval(rgn_pivot(r_tmp, Elecs(iEl)%o_inD%r))
       if ( ia1 <= 0 ) then
          call die('A downfolded region is not existing. Programming error')
       end if

       ! Collect all electrode down-fold regions into one
       ! The downfolded regions *must* not overlap.
       ! So using union may interfere.
       call rgn_append(r_Els, r_oEl(iEl), r_Els)

       ! If the user requests GRAPHVIZ output
       if ( fdf_get('TBT.BTD.Pivot.Graphviz',.false.) .and. IONode ) then
          csort = trim(slabel) // '.TBT.' // trim(Elecs(iEl)%name) // '.gv'
          call sp2graphviz(csort, sp_tmp, pvt=r_oEl(iEl))
       end if

       ! Enlarge the sparse pattern by adding all electrode self-energy terms in
       ! the central region.
       ! Otherwise we do not have a DENSE part of the self-energies in the central
       ! region. This is _very_ important.
       call crtSparsity_Union(dit, sp_tmp, Elecs(iEl)%o_inD, sp_tmp)

    end do

    call timer('pivot-elec', 2)

    call timer('pivot-device', 1)

    if ( IONode ) then
       write(*,'(/a)')'tbt: Analyzing device sparsity pattern and &
            &pivot-table'
    end if
    
    ! We sort the device region based on the
    ! first orbital that the electrode connects the most too.
    ! This seems like a good choice as we know 
    ! it will be on the boundary between one electrode
    ! and the device region

    ! Re-create the device region from a sorting algorithm.
    ! However, in cases where the user is sure that
    ! the device region is correctly sorted 
    ! then we need not re-sort it.
    ! This seems like the best choice when looking
    ! at TB models where number of connections is the same

    ! Prepare the check-regions
    ! Default the device region pivoting to start from the electrode
    ! which has the largest overlap in the device region.
    ! This seems like the option that will provide the best BTD
    ! due to the larger initial base.
    iEl = 1
    do i = 2, N_Elec
       if ( rgn_size(Elecs(iEl)%o_inD) < rgn_size(Elecs(i)%o_inD) ) then
          iEl = i
       end if
    end do
    
    if ( .not. fdf_get('TBT.Analyze', .false.) ) then

       ! Only perform the pivoting *if* we do not analyze the
       ! sparsity pattern
       csort = 'atom+'//trim(Elecs(iEl)%name)
       csort = fdf_get('TS.BTD.Pivot',trim(csort))
       csort = fdf_get('TBT.BTD.Pivot',trim(csort))
       csort = fdf_get('TBT.BTD.Pivot.Device',trim(csort))
       call ts_pivot(dit, sp_tmp, &
            N_Elec, Elecs, &
            cell, na_u, xa, lasto, &
            r_oDev, csort)

       ! Print out what we found
       if ( IONode ) then
          write(*,'(a)')'tbt: BTD pivoting scheme in device: '//trim(csort)
       end if

    end if

    call timer('pivot-device', 2)

    ! Check that there is no overlap with the other regions
    call rgn_sort(r_Els)
    if ( rgn_overlaps(r_Els, r_oDev) ) then
       call rgn_print(r_oDev)
       call rgn_print(r_Els)
       print *,'Overlapping device, down-folding region(s)...'
       call die('tbt_regions: Error in programming, electrode down')
    end if

    if ( rgn_overlaps(r_oBuf, r_oDev) ) then
       call rgn_print(r_oDev)
       call rgn_print(r_oBuf)
       print *,'Overlapping device and buffer region...'
       call die('tbt_regions: Error in programming, buffer')
    end if

    ! Ensure that the number of device orbitals + electrode downfolding
    ! + buffer orbitals equal the full system
    if ( r_Els%n + r_oDev%n + r_oBuf%n /= no_u ) then
       r_Els%name = 'Electrodes'
       r_oDev%name = 'Device'
       call rgn_print(r_oBuf)
       call rgn_print(r_oDev)
       call rgn_print(r_Els)
       call die('tbt_regions: Error in programming, total')
    end if

    call rgn_Orb2Atom(r_oDev,na_u,lasto,r_aDev)
    r_oDev%name = '[O]-device'
    r_aDev%name = '[A]-device'
    
    ! If the user requests GRAPHVIZ output
    if ( fdf_get('TBT.BTD.Pivot.Graphviz',.false.) .and. IONode ) then
       csort = trim(slabel) // '.TBT.gv'
       call sp2graphviz(csort, sp_tmp, pvt=r_oDev)
    end if

    ! Before we proceed we should create the Hamiltonian
    ! sparsity pattern used for setting up the Green function
    call ts_Sparsity_Global(dit, sp, N_Elec, Elecs, sp_uc)
    
    ! Do a final check that all regions are correctly setup
    ! We know that the sum of each segment has to be the 
    ! total number of orbitals in the region.
    i = r_oDev%n
    do iEl = 1 , N_Elec
       i = i + r_oEl(iEl)%n
    end do
    i = i + r_oBuf%n
    if ( i /= no_u ) then
       if ( Node == 0 ) then
          write(*,'(a,i0)')'Buffer orbitals: ',r_oBuf%n
          write(*,'(a,i0)')'Device orbitals: ',r_oDev%n
          do iEl = 1 , N_Elec
             write(*,'(a,i0)')trim(Elecs(iEl)%name)//' orbitals: ',r_oEl(iEl)%n
          end do
          if ( i > no_u ) then
             ! find the overlapping orbitals
             call rgn_union(r_oBuf,r_oDev,r_tmp)
             do iEl = 1 , N_Elec
                call rgn_union(r_tmp,r_oEl(iEl),r_tmp)
             end do
             r_tmp%name = 'Double counted orbitals'
             call rgn_sort(r_tmp)
             call rgn_print(r_tmp)
          else
             ! find the missing orbitals
             call rgn_range(r_tmp,1,no_u)
             if ( r_oBuf%n > 0 ) then
                call rgn_complement(r_oBuf,r_tmp,r_tmp)
             end if
             call rgn_complement(r_oDev,r_tmp,r_tmp)
             do iEl = 1 , N_Elec
                call rgn_complement(r_oEl(iEl),r_tmp,r_tmp)
             end do
             r_tmp%name = 'Missing orbitals'
             call rgn_sort(r_tmp)
             call rgn_print(r_tmp)
          end if
          write(*,'(a,2(tr1,i0))')'Total number of orbitals vs. counted:',no_u,i

          write(*,'(/,a)')'Missing/Excess orbitals can happen if your device &
               &region is ill-formatted.'
          write(*,'(a)')'Suppose you create a device region which disconnects &
               &certain non-device region orbitals from the electrode regions.'
          write(*,'(a)')'Then this will occur, please ensure that you have &
               &defined your device region such that the above does not occur.'
       end if
#ifdef MPI
       call MPI_Barrier(MPI_Comm_World,i)
#endif
       call die('Something went wrong when asserting the &
            &total number of orbitals. Have you requested &
            &something not applicable?')
    end if

    ! Clean-up
    call rgn_delete(r_tmp,r_tmp2,r_tmp3,r_Els,priority)
    call delete(sp_tmp)

    if ( IONode ) then
       write(*,'(a)')'tbt: Done analyzing electrode and device sparsity pattern and pivot-tables'
    end if

    ! Clean up un-used memory
    do iEl = 1, N_Elec
      call rgn_delete(r_aEl_alone(iEl))
      call rgn_delete(r_oEl_alone(iEl))
    end do
    
  contains

    function sort_contain(str,name) result(contain)
      use m_char, only : lcase
      character(len=*), intent(in) :: str, name
      logical :: contain

      character(len=len(str)) :: lstr
      character(len=len_trim(name)) :: lname

      integer :: i

      contain = .false.

      lstr = lcase(str)
      lname = lcase(trim(name))

      ! check whether it is in this stuff
      i = index(lstr,lname)
      if ( i > 1 ) then
         contain = scan(lstr(i-1:i-1),'+ ') == 1
         i = i + len_trim(lname)
         if ( i <= len(str) ) then
            contain = contain .and. scan(lstr(i:i),'+ ') == 1
         end if
      else if ( i == 1 ) then
         i = i + len_trim(lname)
         if ( i <= len(str) ) then
            contain = scan(lstr(i:i),'+ ') == 1
         else
            ! it was found and the string is too short
            contain = .true.
         end if
      end if
      
    end function sort_contain

  end subroutine tbt_init_regions


  subroutine tbt_region_options( sp, save_DATA )
    use dictionary
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    use m_sparsity_handling, only : Sp_retain_region, Sp_sort, Sp_union
#ifdef NCDF_4
    use m_tbt_delta, only: read_delta_Sp
    use m_tbt_dH, only: use_dH, dH
#endif

    type(Sparsity), intent(inout) :: sp
    type(dictionary_t), intent(in) :: save_DATA
#ifdef NCDF_4
    type(OrbitalDistribution) :: fdit

    integer :: no_u
    type(Sparsity) :: sp_dH
#endif

    ! Make sure to initialize the device region
    ! sparsity pattern
    call delete(sp_dev_sc)
#ifdef NCDF_4
    if ( ('orb-current' .in. save_DATA) .or. &
         ('proj-orb-current' .in. save_DATA) .or. &
         ('DM-Gf' .in. save_DATA) .or. ('DM-A' .in. save_DATA) .or. &
         ('COOP-Gf' .in. save_DATA) .or. ('COHP-Gf' .in. save_DATA) .or. &
         ('COOP-A' .in. save_DATA) .or. ('COHP-A' .in. save_DATA) ) then

       call attach(sp,nrows_g=no_u)
#ifdef MPI
       call newDistribution(no_u,MPI_Comm_Self,fdit,name='TBT-fake dist')
#else
       call newDistribution(no_u,-1           ,fdit,name='TBT-fake dist')
#endif
       call Sp_retain_region(fdit,sp,r_oDev,sp_dev_sc)
       ! Note that the delta-Sigma is not necessary because
       ! the self-energy does not add to bond-currents, etc.
       if ( use_dH ) then
         call read_delta_Sp(dH,no_u,sp_dH)
         call Sp_retain_region(fdit,sp_dH,r_oDev,sp_dH)
         call Sp_union(fdit,sp_dev_sc,sp_dH,sp_dev_sc)
         call delete(sp_dH)
       end if
       call Sp_sort(sp_dev_sc)
       call delete(fdit)

    end if
#endif

  end subroutine tbt_region_options

  subroutine tbt_print_regions(na_u, lasto, N_Elec, Elecs)

    use parallel, only : Node
    use m_verbosity, only : verbosity
    use m_ts_electype

    integer, intent(in) :: na_u
    integer, intent(in) :: lasto(0:na_u)
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer :: i
    type(tRgn) :: r

    if ( Node /= 0 ) return

    if ( verbosity < 3 ) return

    ! Print out the buffer regions
    if ( r_aBuf%n > 0 ) then
       call local_print(r_aBuf, .false. )
       call local_print(r_oBuf, .true. )
    end if

    ! Print out the device region
    write(*,'(a,i0)')'tbt: # of device region orbitals: ',r_oDev%n

    call local_print(r_aDev, .false. )
    call local_print(r_oDev, .true. )

    ! Print out all the electrodes + their projection region
    do i = 1 , N_Elec
      write(*,*) ! new-line
      write(*,'(3a,i0)')'tbt: # of ',trim(Elecs(i)%name), &
          ' downfolding orbitals: ',r_oElpD(i)%n
      write(*,'(3a,i0)')'tbt: # of ',trim(Elecs(i)%name), &
          ' device orbitals: ',Elecs(i)%o_inD%n
      call local_print(r_aEl(i), .false.)
      call local_print(r_oEl(i), .true.)
      if ( verbosity > 3 ) then
        ! Create the atom equivalent regions
        call rgn_Orb2Atom(Elecs(i)%o_inD, na_u, lasto, r)
        call rgn_sort(r)
        r%name = '[A]-'//trim(Elecs(i)%name)//' folding in D'
        call local_print(r, .false.)
      end if
      if ( verbosity > 7 ) then
        call rgn_intersection(r_oElpD(i),r_oDev,r)
        r%name = '[O]-'//trim(Elecs(i)%name)//' folding in D'
        call local_print(r, .true.)
      end if
    end do

    ! Clean-up
    call rgn_delete(r)

  contains
    
    subroutine local_print(r_in,is_orb)
      type(tRgn), intent(in) :: r_in
      logical, intent(in) :: is_orb
      type(tRgn) :: r
      integer :: seq, mid, high

      if ( is_orb ) then
         seq = 10
         mid = 7
         high = 9
      else
         seq = 12
         mid = 4
         high = 6
      end if

      if ( verbosity > high ) then
         ! print-unsorted
         call rgn_print(r_in, seq_max = seq )
      else if ( verbosity > mid ) then
         call rgn_copy(r_in,r)
         call rgn_sort(r)
         call rgn_print(r, seq_max = seq )
         call rgn_delete(r)
      end if

    end subroutine local_print

  end subroutine tbt_print_regions

end module m_tbt_regions
