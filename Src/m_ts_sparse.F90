!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! A module that supplements the reduced memory TranSiesta version.
! It greatly reduces the memory requirement of transiesta, as well
! as making the code more clearer in intent.
! This module is probably the most *hard* to understand part.

module m_ts_sparse

  use class_Sparsity
  use class_iSpData1D

  use precision, only : dp

  implicit none

  ! The full Transiesta sparsity region in the unit-cell equivalent
  ! region, this is used to accomodate the Hamiltonian and overlap
  ! matrices. In principle this could be made obsolete. However,
  ! that task seems more cumbersome than worthy of notice (for
  ! the moment).
  ! This is a *SORTED* sparse matrix
  type(Sparsity), save :: ts_sp_uc ! TS-GLOBAL (UC)

  ! We will save the "update region sparsity"
  ! Note that this is NOT the same as the sparsity pattern
  ! provided by SIESTA!
  ! After using this sparsity pattern, we do not need the listud[g] arrays
  ! as this sparsity pattern is a reflection of the mask used by listud[g]
  ! Lastly the usage of a MASKED sparsity pattern, will reduce the clutter
  ! between MPI and non-MPI codes as it will be the same array in both 
  ! circumstances.
  ! This is a *SORTED* sparse matrix
  type(Sparsity), save :: tsup_sp_uc ! TS-update-GLOBAL (UC)

  ! We will save the local "update region sparsity"
  ! Note that this is NOT the same as the sparsity pattern
  ! provided by SIESTA!
  ! After using this sparsity pattern, we do not need the listud[g] arrays
  ! as this sparsity pattern is a reflection of the mask used by listud[g]
  ! This reflects the local sparsity pattern updated elements.
  ! This is a *SORTED* sparse matrix
  type(Sparsity), save :: ltsup_sp_sc ! TS-update-local (SC)

  ! This is an index array which points from ltsup_sp_sc to the local siesta
  ! sparsity pattern.
  ! TODO: check how much speed we gain from not searching in the sparsity pattern
  type(iSpData1D), save :: ltsup_sc_pnt

  ! The offsets for the supercells
  real(dp), pointer, save :: sc_off(:,:) => null()

contains

! This routine setups what-ever is needed to do the
! memory reduced TranSiesta code.
! This means collecting information about which region needs
! update, etc.
  subroutine ts_sparse_init(slabel, &
       IsVolt, N_Elec, Elecs, &
       ucell, nsc, na_u,xa,lasto, dit,sparse_pattern, Gamma, &
       isc_off)
    
    use class_OrbitalDistribution

    use alloc
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif 
    use parallel, only: IONode

    use create_Sparsity_SC
    use m_ts_electype
    use m_ts_method
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
    use parallel, only: Node
#endif
    use m_region
    use m_sparsity_handling

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: slabel
    logical, intent(in) :: IsVolt ! bias calculation
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! Unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of super-cells in each direction
    integer, intent(in) :: nsc(3)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! The distribution for the sparsity-pattern
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sparse_pattern
    ! Whether we have xij or not
    logical, intent(in) :: Gamma
    ! supercell indices
    integer, intent(in) :: isc_off(3,product(nsc))

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution) :: fdit
    ! We create a temporary sparsity pattern which removes
    ! all cross-connections across the electrode transport direction.
    type(Sparsity) :: tmp_sp
    ! Temporary arrays for knowing the electrode size
    integer :: no_u_TS, i

    ! Number of orbitals in TranSiesta
    no_u_TS = nrows_g(sparse_pattern) - no_Buf

    ! Do a crude check of the sizes
    if ( no_u_TS <= sum(TotUsedOrbs(Elecs)) ) then
      call die("The contact region size is &
          &smaller than the electrode size. Please correct.")
    end if

    ! Create the ts-offsets
    if ( Gamma ) then
      ! Initialize the sc_off array
      call re_alloc(sc_off,1,3,1,1)
      sc_off(:,:) = 0._dp
    else
      i = product(nsc)
      call re_alloc(sc_off,1,3,1,i)
      sc_off(:,:) = matmul(ucell,isc_off)
    end if
    
    ! Removes all buffer atoms, cross-terms between electrodes
    ! and electrode cross-boundaries.
    call ts_Sp_calculation(dit,sparse_pattern,N_Elec,Elecs, &
        ucell, nsc, isc_off, tmp_sp)

    if ( IONode .and. Gamma ) then
      write(*,'(a,a)')'transiesta: ',&
          'We cannot assure cross-boundary connections in Gamma calculations.'
      write(*,'(a,a)')'transiesta: ',&
          'Ensure that the electrode says: Principal cell is perfect'
    end if

    ! Create the global transiesta H(k), S(k) sparsity pattern
    call ts_Sparsity_Global(dit, tmp_sp, N_Elec, Elecs, ts_sp_uc)
    
    if ( IONode ) then
      write(*,'(/,a)') 'transiesta: created H and S sparsity pattern:'
      call print_type(ts_sp_uc)
    end if
    
#ifdef TRANSIESTA_DEBUG
    if ( IONode ) write(*,*) 'Created TS-Global HS (100)'
    call sp_to_file(100+Node,ts_sp_uc)
#endif

    if ( IsVolt ) then 

      ! Create the update region (a direct subset of the local sparsity pattern)
      ! Hence it is still a local sparsity pattern.
      call ts_Sparsity_Update(dit, tmp_sp, N_Elec, Elecs, ltsup_sp_sc)

      if ( IONode ) then
        write(*,'(/,a)') 'transiesta: created local update sparsity pattern:'
        call print_type(ltsup_sp_sc)
      end if
      
#ifdef TRANSIESTA_DEBUG
      if ( IONode ) write(*,*)'Created TS-local UP (300)'
      call sp_to_file(300+Node,ltsup_sp_sc)
#endif

      ! Create the pointer from the local transiesta update sparsity 
      ! to the local siesta sparsity
      call ts_Sparsity_Subset_pointer(dit,sparse_pattern,ltsup_sp_sc, &
          ltsup_sc_pnt)

    end if

    ! In order to ensure that the electrodes are in the
    ! tri-diagonal sparsity pattern, we can easily create
    ! the full sparsity pattern with the electrodes included
    ! and then recreate the tri-diagonal sparsity pattern
    ! This is probably the crudest way of doing it.
#ifdef MPI
    call newDistribution(nrows_g(ts_sp_uc),MPI_Comm_Self,fdit, &
        name='TranSiesta UC distribution')
#else    
    call newDistribution(nrows_g(ts_sp_uc),-1,fdit, &
        name='TranSiesta UC distribution')
#endif

    ! The update sparsity pattern can be simplied to the H,S sparsity
    ! pattern, if all electrodes have certain options to be the same.
    if ( ( all(Elecs(:)%Bulk) .and. all(Elecs(:)%DM_update == 1) ) .or. &
        ( all(.not. Elecs(:)%Bulk) .and. all(Elecs(:)%DM_update == 2) ) ) then
      ! The sparsity patterns are the same, i.e. the
      ! update and Hamiltonian sparse matrices are the same.
      ! We can re-use the sparsity pattern.
      
      tsup_sp_uc = ts_sp_uc
      if ( IONode ) then
        write(*,'(/a)') 'transiesta: update sparsity pattern same as H and S.'
      end if

    else

      ! In this case we cannot be sure that the update sparsity pattern
      ! is a direct sub-set of the Hamiltonian.
      ! Hence, we need to create the correct one.
      ! In this case is is simply a global sparsity of the 
      !   'ltsup_sp_sc', in case there is an applied bias
      if ( IsVolt ) then

        ! Convert to unit-cell sparsity pattern
        call crtSparsity_SC(ltsup_sp_sc, tmp_sp, UC=.TRUE.)

      else

        ! Create update region (in local sparsity pattern)
        call ts_Sparsity_Update(dit, tmp_sp, N_Elec, Elecs, tsup_sp_uc)
        ! Convert to unit-cell sparsity pattern
        call crtSparsity_SC(tsup_sp_uc, tmp_sp, UC=.TRUE.)

      end if

      ! Convert to the global one
      call Sp_to_Spglobal(dit, tmp_sp, tsup_sp_uc)

      if ( IONode ) then
        write(*,'(/,a)') 'transiesta: created update sparsity pattern:'
        call print_type(tsup_sp_uc)
      end if

    end if

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created TS-Global update (200)'
    call print_type(tsup_sp_uc)
    call sp_to_file(200+Node,tsup_sp_uc)
#endif

    ! Clean-up
    call delete(tmp_sp)

  end subroutine ts_sparse_init


  ! This routine returns a sparsity pattern which only contains
  ! orbitals that are in the device + electrode region.
  ! thus we have removed the following orbitals:
  !   - buffer
  !   - cross terms between:
  !     * electrode-electrode
  !     * electrode-cross boundaries, say if the electrode
  !       basis range is extending across the other side electrode
  !       this ensures no "ghost-transmission"
  subroutine ts_Sp_calculation(dit,s_sp,N_Elec,Elecs, &
       ucell, nsc, isc_off, &
       ts_sp)

    use parallel, only: IONode
    
    use class_OrbitalDistribution
    use class_Sparsity
    use m_region

    use m_ts_electype

    use m_sparsity_handling
    use m_ts_method

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: s_sp ! the local sparse pattern
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))
    ! the global sparse pattern for the calculation region (in SC format)
    type(Sparsity), intent(inout) :: ts_sp

    type(tRgn) :: r_oE(N_Elec), r_tmp1, r_tmp2
    integer :: iEl, i

    integer :: init_nz, old_nz

    if ( r_oBuf%n > 0 ) then
      ! Remove buffer atoms...
      call Sp_remove_region(dit,s_sp,r_oBuf,ts_sp)
    else
      ts_sp = s_sp
    end if

    ! Get initial number of nnz
    init_nz = nnzs(ts_sp)

    ! Remove all electrode to other side connections
    ! this only has effect when we cross 
    do iEl = 1 , N_Elec

       ! Create electrode region
       i = Elecs(iEl)%idx_o
       call rgn_range(r_oE(iEl), i, i+TotUsedOrbs(Elecs(iEl))-1)

       ! Remove the connections that cross the boundary
       ! starting from this electrode
       call rgn_sp_connect(r_oE(iEl), dit, ts_sp, r_tmp1)
       call rgn_union(r_oE(iEl), r_tmp1, r_tmp2)

       ! Calculate the transport direction in the device cell.
       ! We expect there to be only one and thus find the transport
       ! direction in the big cell
       old_nz = nnzs(ts_sp)
 
       select case ( Elecs(iEl)%t_dir )
       case ( 4 ) ! B-C
         i = Elecs(iEl)%pvt(2)
         call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case ( 5 ) ! A-C
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case ( 6 ) ! A-B
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(2)
       case ( 7 ) ! A-B-C
         i = Elecs(iEl)%pvt(1)
         call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(2)
         call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)
         i = Elecs(iEl)%pvt(3)
       case default
         i = Elecs(iEl)%pvt(Elecs(iEl)%t_dir)
       end select

       ! Remove connections from this electrode across the boundary...
       call Sp_remove_crossterms(dit,ts_sp,product(nsc),isc_off, i, ts_sp, r = r_tmp2)

       ! Update init_nz for the valid removed elements
       init_nz = init_nz - (old_nz - nnzs(ts_sp))

       ! We also be sure to remove all direct connections
       if ( iEl > 1 ) then
          call rgn_delete(r_tmp1,r_tmp2)
          do i = 1 , iEl - 1
             call rgn_copy(r_tmp2,r_tmp1)
             call rgn_union(r_oE(i),r_tmp1,r_tmp2)
          end do
          call Sp_remove_region2region(dit,ts_sp,r_oE(iEl),r_tmp2,ts_sp)
       end if

       call rgn_delete(r_tmp1,r_tmp2)
       
    end do

    i = nnzs(ts_sp)
    if ( i < init_nz .and. IONode ) then
       write(*,'(/a,i0,a/)')'*** WARNING! Removed ',init_nz - i, ' elements &
           &which connect electrodes across the device region!'
    end if
    
    do iEl = 1 , N_Elec
      call rgn_delete(r_oE(iEl))
    end do

  end subroutine ts_Sp_calculation


  ! Returns the sparsity pattern of the Hamiltonian and the overlap
  ! matrix that is used in the Green function calculation.
  ! It requires the passed sparsity pattern to have removed
  ! any cross-terms.
  ! If they are not already removed a more simpler algorithm will
  ! remove those terms via direct interaction.
  subroutine ts_Sparsity_Global(dit,s_sp, &
       N_Elec, Elecs, ts_sp)

#ifdef TRANSIESTA_DEBUG
    use parallel, only : IONode
#endif
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use intrinsic_missing, only : SORT_QUICK
    use create_Sparsity_SC
    use m_ts_electype
    use m_ts_method
    use m_sparsity_handling

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
    use parallel, only: Node
#endif

! **********************
! * INPUT variables    *
! **********************
    ! The SIESTA distribution of the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! Sparsity patterns of SIESTA (local)
    type(Sparsity), intent(inout) :: s_sp
    ! All the electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the returned update region.    
    type(Sparsity), intent(inout) :: ts_sp

! **********************
! * LOCAL variables    *
! **********************
    ! We need temporary sparsity patterns which will be deleted
    ! We need this to generate the Transiesta sparsity
    type(Sparsity) :: sp_global, sp_uc

    ! to globalize from the local sparsity pattern (SIESTA)
    ! and afterwards used as the looping mechanisms
    ! to create the mask for creating the UC transiesta pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_l, no_u, uc_n_nzs, n_nzsg

    ! search logical to determine the update region...
    logical, allocatable :: l_HS(:)

    ! Loop-counters
    integer :: io, jo, ind
    integer :: ict, jct
    logical :: UseBulk

    ! Initialize
    call delete(ts_sp)

    call attach(s_sp,nrows=no_l,nrows_g=no_u)

    ! Create the (local) SIESTA-UC sparsity...
#ifdef MPI
    call crtSparsity_SC(s_sp, sp_global, UC=.TRUE.)
    uc_n_nzs = nnzs(sp_global)

    ! point to the local (SIESTA-UC) sparsity pattern arrays
    call Sp_to_Spglobal(dit,sp_global,sp_uc)

    ! Delete the local UC sparsity pattern
    call delete(sp_global)

#else
    call crtSparsity_SC(s_sp, sp_uc, UC=.TRUE.)
    uc_n_nzs = nnzs(sp_uc)
#endif

#ifdef TRANSIESTA_DEBUG
    if(IONode)write(*,*)'Created UC SIESTA sparsity (400)'
    call print_type(sp_uc)
    call sp_to_file(400+Node, sp_uc)
#endif

    ! Now we have the globalized SIESTA Unit-cell pattern, make it 
    ! only to the Transiesta sparsity pattern

    ! Immediately point the global arrays to their respective parts
    call attach(sp_uc,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nnzs=n_nzsg)

    ! allocate space for the MASK to create the TranSiesta GLOBAL region
    allocate(l_HS(n_nzsg))
    call memory('A','L',n_nzsg,'transiesta')

    ! Initialize
    l_HS(:) = .false.

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
!$OMP parallel do default(shared), &
!$OMP&private(io,jo,ind,ict,jct,UseBulk)
    do io = 1 , no_u

       ict = orb_type(io)
       if ( ict /= TYP_BUFFER ) then

       ! The index in the pointer array is retrieved
       do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)
          
          ! The unit-cell column index (remember we are looping a UC SP)
          jo = l_col(ind)

          ! If we are in the buffer region, cycle (l_HS(ind) =.false. already)
          ! note, that we have already *checked* ic
          jct = orb_type(jo)
          if ( jct == TYP_BUFFER ) cycle

          UseBulk = .false.
          ! Specify bulk according to electrode
          ! In order to allow to have the update sparsity pattern
          ! as a subset, we require that the Hamiltonian also
          ! has the electrode interconnects
          if ( ict > TYP_DEVICE ) UseBulk = Elecs(ict)%Bulk
          if ( jct > TYP_DEVICE ) UseBulk = UseBulk .or. Elecs(jct)%Bulk

          if ( UseBulk ) then
             ! here we create the Hamiltonian matrix on these criterias:
             !  1) no electrode-electrode connections
             !  2) no only-electrode connections
             !  3) add cross-terms

             l_HS(ind) = ict == TYP_DEVICE .or. jct == TYP_DEVICE

          else
             
             ! If not bulk we only want to save the Hamiltonian elements
             ! for the same electrode, otherwise everything is needed
             if ( ict > TYP_DEVICE .and. jct > TYP_DEVICE ) then
                l_HS(ind) = ict == jct
             else
                l_HS(ind) = .true.
             end if
             
          end if

       end do

       end if

    end do
!$OMP end parallel do

    ! We now have a MASK of the actual needed TranSiesta sparsity pattern
    ! We create the TranSiesta sparsity pattern
    call crtSparsity_SC(sp_uc,ts_sp,MASK=l_HS)

    ! Ensure that it is sorted
    ! Sorting the sparsity pattern greatly speeds up the Hamiltonian
    ! and overlap matrix searches
    call attach(ts_sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
!$OMP parallel do default(shared), private(io)
    do io = 1 , no_u
       if ( l_ncol(io) /= 0 ) then
          call sort_quick(l_ncol(io), l_col(l_ptr(io)+1:))
       end if
    end do
!$OMP end parallel do

    ! clean-up
    deallocate(l_HS)
    call memory('D','L',n_nzsg,'transiesta')
    
    ! Furthermore, we dont need the SIESTA UC sparsity global...
    call delete(sp_uc)

  end subroutine ts_Sparsity_Global


  ! Returns the global sparsity pattern of the density density
  ! matrix that is to be updated.
  ! This may be a sub-set of the Hamiltonian and overlap matrices
  ! used for creating the Green function, or it could be different,
  ! depending on the bulkiness and update settings for transiesta.
  subroutine ts_Sparsity_Update(dit,s_sp, N_Elec, Elecs, &
       tsup_sp)

    use geom_helper, only : UCORB
    use create_Sparsity_SC
    use class_OrbitalDistribution
    use m_ts_method
    use m_ts_electype
! **********************
! * INPUT variables    *
! **********************
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: s_sp
    ! the electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(Sparsity), intent(inout) :: tsup_sp

! **********************
! * LOCAL variables    *
! **********************
    ! to globalize from the local sparsity pattern (SIESTA)
    ! and afterwards used as the looping mechanisms
    ! to create the mask for creating the UC transiesta pattern
    integer, pointer :: l_ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    ! Also used in non-MPI (to reduce dublicate code)
    integer :: no_l, no_u, n_nzs

    ! search logical to determine the update region...
    logical, allocatable :: lup_DM(:)
    logical :: DM_bulk, DM_cross

    ! Loop-counters
    integer :: lio, io, ind
    integer :: ict, jct
    
    ! Logical for determining the region
    logical :: i_in_C, j_in_C

    ! Initialize the tsup
    call delete(tsup_sp)

    call attach(s_sp,nrows=no_l,nrows_g=no_u,nnzs=n_nzs)

    ! allocate space for the MASK to create the TranSiesta UPDATE region
    allocate(lup_DM(n_nzs))
    call memory('A','L',n_nzs,'transiesta')

    ! Initialize
    lup_DM(:) = .false.

    call attach(s_sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    ! We do not need to check the buffer regions...
    ! We know they will do NOTHING! :)
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ict,jct,ind,DM_bulk,DM_cross), &
!$OMP&private(i_in_C,j_in_C)
    do lio = 1 , no_l

       ! Shift out of the buffer region
       io = index_local_to_global(dit,lio)

       ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
       ict = orb_type(io)
       if ( ict /= TYP_BUFFER ) then

       ! Loop the index of the pointer array
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! If we are in the buffer region, cycle (lup_DM(ind) =.false. already)
          ! note, that we have already *checked* io
          jct = orb_type(l_col(ind))
          if ( jct == TYP_BUFFER ) cycle

          ! Remark,
          !   DM_update == 0 (none)
          !   DM_update == 1 (cross-terms)
          !   DM_update == 2 (all)
          DM_cross = .true.

          if ( ict > TYP_DEVICE ) then
            ! Assign that the density matrix, should not be updated
            DM_bulk = Elecs(ict)%DM_update < 2
            ! If DM_cross is different from 0, the entire
            ! electrode region is updated
            DM_cross = Elecs(ict)%DM_update /= 0
          end if
          if ( jct > TYP_DEVICE ) then
            DM_bulk = Elecs(jct)%DM_update < 2
            DM_cross = Elecs(jct)%DM_update /= 0
          end if

          ! We check whether it is electrode-connections. 
          ! If, so, they are not used in transiesta:
          if      ( ict > TYP_DEVICE .and. jct > TYP_DEVICE ) then

             ! Remove connections between electrodes
             ! but maintain same electrode updates if not bulk dm
             if ( ict == jct .and. .not. DM_bulk ) then
                lup_DM(ind) = .true.
             end if

          else

             ! Note that this if-statement will never be reached
             ! if both ict and jct are electrodes.
             ! We will only be here if either of them is an electrode.
             
             i_in_C = ict == TYP_DEVICE
             j_in_C = jct == TYP_DEVICE
             
             if ( DM_cross ) then
                ! the user has requested to also update cross-terms
                lup_DM(ind) = i_in_C .or. j_in_C
             else
                ! the user has requested to ONLY update the central region
                lup_DM(ind) = i_in_C .and. j_in_C
             end if

          end if

       end do

       end if

    end do
!$OMP end parallel do

    ! We now have a MASK of the actual needed TranSiesta sparsity pattern
    ! We create the TranSiesta sparsity pattern
    call crtSparsity_SC(s_sp,tsup_sp,MASK=lup_DM)

    call memory('D','L',n_nzs,'transiesta')
    deallocate(lup_DM)

  end subroutine ts_Sparsity_Update

  subroutine ts_Sparsity_Subset_pointer(dit,sp,sub_sp,ipnt)
    use class_OrbitalDistribution

! **********************
! * INPUT variables    *
! **********************
    ! The distribution pattern for everything
    type(OrbitalDistribution), intent(inout) :: dit
    ! The sparsity pattern we wish to point to
    type(Sparsity), intent(inout) :: sp
    ! The sparsity pattern we wish to point from
    type(Sparsity), intent(inout) :: sub_sp
    ! The pointer index
    type(iSpData1D), intent(inout) :: ipnt

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: sub_ncol(:), sub_ptr(:), sub_col(:)
    integer, pointer :: pnt(:)

    integer :: no_l, io, j, sub_ind, ind

    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l)
    call attach(sub_sp,n_col=sub_ncol,list_ptr=sub_ptr,list_col=sub_col, &
         nrows=io)
    if ( io /= no_l ) call die('Could not do index matching due to &
         &inconsistent sparsity patterns')

    ! Clear and create
    call delete(ipnt)
    call newiSpData1D(sub_sp,dit,ipnt,name='TS pointer')
    ! Point to the array
    pnt => val(ipnt)
    ! Initialize to check that we can locate all indices
    pnt(:) = 0

    ! Loop in the subset sparsity pattern
!$OMP parallel do default(shared), &
!$OMP&private(io,j,sub_ind,ind)
    do io = 1 , no_l

       ! Loop number of entries in the row...
       do j = 1 , sub_ncol(io)

          ! The index in the pointer array is retrieved
          sub_ind = sub_ptr(io) + j

          ! Loop in the super-set sparsity pattern
          super_idx: do ind = l_ptr(io) + 1 , l_ptr(io) + l_ncol(io)

             ! If we have the same column index it must be
             ! the same entry they represent
             if ( sub_col(sub_ind) == l_col(ind) ) then
                pnt(sub_ind) = ind
                exit super_idx
             end if

          end do super_idx
       end do
    end do
!$OMP end parallel do

    if ( any(pnt(:) == 0) ) then
       call die('An index could not be located in the super-set &
            &sparsity pattern. Are you surely having the correct &
            &sparsity?')
    end if

    if ( any(pnt(:) > nnzs(sp)) ) then
       call die('An index could not be located in the super-set &
            &sparsity pattern. Are you surely having the correct &
            &sparsity?')
    end if

  end subroutine ts_Sparsity_Subset_pointer

end module m_ts_sparse
