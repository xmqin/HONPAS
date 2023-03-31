! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This particular solution method relies on solving the GF
! with the tri-diagonalization routine.
! This will leverage memory usage and also the execution time.

module m_tbt_tri_init

  use precision, only : dp, i8b
  use m_region

  implicit none

  public :: tbt_tri_init
  public :: tbt_tri_print_opti

  type(tRgn), save, allocatable, target :: ElTri(:)
  type(tRgn), save :: DevTri

  public :: ElTri, DevTri
  public :: fold_elements

  private
  
contains

  subroutine tbt_tri_init_elec( dit , sp )

    use parallel, only : Node, Nodes
    use class_OrbitalDistribution
    use class_Sparsity
    use create_Sparsity_Union

    use m_ts_tri_common, only: ts_pivot_tri_sort_El
    use m_ts_rgn2trimat
    use m_ts_electype
    use m_ts_method, only: TS_BTD_A_COLUMN, TS_BTD_A_PROPAGATION

#ifdef MPI
    use mpi_siesta
#endif
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    use m_tbt_options, only : N_Elec, Elecs, BTD_method
    use m_tbt_regions

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp

    type(Sparsity) :: tmpSp1, tmpSp2
    type(tRgn) :: tmp_roEl
    integer :: i, iEl

    call timer('tri-init-elec',1)

    ! This works as creating a new sparsity deletes the previous
    ! and as it is referenced several times it will not be actually
    ! deleted...
    allocate(ElTri(N_Elec))
    do i = 1 + Node , N_Elec , Nodes

       ! Retain region
       call Sp_retain_region(dit,sp,r_oElpD(i),tmpSp2)

       ! Add the self-energy of the electrode (in its original position)
       call rgn_range(tmp_roEl, Elecs(i)%idx_o, Elecs(i)%idx_o + TotUsedOrbs(Elecs(i)) - 1)
       call crtSparsity_Union(dit,tmpSp2, tmp_roEl,tmpSp1)
       call delete(tmpSp2)

#ifdef TRANSIESTA_DEBUG
       open(file='ELEC_'//trim(Elecs(i)%name)//'_SP',unit=1400,form='formatted')
       call sp_to_file(1400,tmpSp1)
       close(1400)
#endif

       ! Create tri-diagonal parts for this electrode
       ! IF parts == 0 will create new partition
       call ts_rgn2TriMat(1, Elecs(i:i), .false., &
            dit, tmpSp1, r_oElpD(i), ElTri(i)%n, ElTri(i)%r, &
            BTD_method, last_eq = Elecs(i)%o_inD%n , par = .false. )
       call delete(tmpSp1)

    end do

    ! Clean memory
    call rgn_delete(tmp_roEl)

    ! The i'th processor has the following electrodes
    do iEl = 1 , N_Elec
                 
#ifdef MPI
       ! The node having this electrode is
       i = mod(iEl-1,Nodes)

       ! B-cast the tri-diagonal matrix from the
       ! processor
       call rgn_MPI_Bcast(ElTri(iEl),i)
#endif

       ! Set the name 
       ElTri(iEl)%name = '[TRI] '//trim(Elecs(iEl)%name)

       ! Sort the tri-diagonal blocks
       call ts_pivot_tri_sort_El(nrows_g(sp), r_oElpD(iEl), 1, Elecs(iEl:iEl), ElTri(iEl))
       
    end do

    call timer('tri-init-elec',2)
    
  end subroutine tbt_tri_init_elec

  subroutine tbt_tri_init( dit , sp , cell, na_u, xa, lasto, proj )

    use sys, only : bye
    use fdf, only: fdf_get
    use parallel, only : IONode
#ifdef MPI
    use mpi_siesta, only: MPI_Bcast, MPI_Comm_World
#endif
    use class_OrbitalDistribution
    use class_Sparsity

    use create_Sparsity_SC, only: crtSparsity_SC
    use create_Sparsity_Union, only: crtSparsity_Union

    use m_ts_rgn2trimat
    use m_ts_tri_common, only : ts_pivot_tri_sort_El
    use m_ts_tri_common, only : nnzs_tri_i8b
    use m_ts_electype
#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    use m_sparsity_handling
    use m_tbt_options, only : N_Elec, Elecs, BTD_method
    use m_tbt_regions

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)
    ! An array of additional projection regions
    ! which determines the projection of a molecule
    ! onto seperate regions
    type(tRgn), intent(in), optional :: proj(:)

    type(Sparsity) :: tmpSp1, tmpSp2
    type(tRgn) :: r_tmp
    integer :: i, n, iEl
    integer(i8b) :: els
    integer, allocatable :: rpvt(:)

    call timer('tri-init',1)

    if ( IONode ) then
       write(*,'(/,a)')'tbt: Creating electrode tri-diagonal matrix blocks'
    end if
    
    ! Copy over sparsity pattern (in UC format)
    call crtSparsity_SC(sp, tmpSp1, UC = .true.)
    
    do iEl = 1 , N_Elec

       ! Add the self-energy of the electrode in the projected position
       ! of the "device" region.
       call crtSparsity_Union(dit,tmpSp1,Elecs(iEl)%o_inD,tmpSp1)

    end do

    ! We have now already added the projected position of the
    ! self-energies. Create the tri-diagonal matrices
    ! for the electrode down-folding regions
    call tbt_tri_init_elec( dit , tmpSp1 )

    if ( IONode ) then
       write(*,'(a)')'tbt: Creating device tri-diagonal matrix blocks'
    end if

    if ( present(proj) ) then
       do i = 1 , size(proj)

          ! Add the self-energy of the electrode in the projected position
          ! of the "device" region.
          call crtSparsity_Union(dit,tmpSp1,proj(i),tmpSp1)
          
       end do
    end if

    ! Create the device region sparsity pattern by removing everything
    ! else....
    call Sp_retain_region(dit,tmpSp1,r_oDev,tmpSp2)
    call delete(tmpSp1)

    if ( fdf_get('TBT.Analyze', .false.) ) then
       
       call tbt_tri_analyze(dit, tmpSp2, N_Elec, Elecs, &
           cell, na_u, xa, lasto, r_oDev, BTD_method)

       call timer('TBT-analyze', 3)
       
       call bye('Stopping TBtrans on purpose after analyzation step...')
       
    end if

#ifdef TRANSIESTA_DEBUG
    open(file='DEV_FULL_SP',unit=1400,form='formatted')
    call sp_to_file(1400,tmpSp2)
    close(1400)
#endif

    call rgn_delete(DevTri)

    
    ! Create tri-diagonal parts for this one...
    call ts_rgn2TriMat(N_Elec, Elecs, .true., &
       dit, tmpSp2, r_oDev, DevTri%n, DevTri%r, &
       BTD_method, last_eq = 0, par = .true. )
    call delete(tmpSp2) ! clean up

    i = nrows_g(sp)
    ! Sort the tri-diagonal blocks
    call ts_pivot_tri_sort_El(i, r_oDev, N_Elec, Elecs, DevTri)

    ! Create back-pivoting table to easily find indices in the device
    ! pivoting table.
    allocate(rpvt(i))
    do i = 1, r_oDev%n
      rpvt(r_oDev%r(i)) = i
    end do

    ! The down-folded region can "at-will" be sorted
    ! in the same manner it is seen in the device region.
    ! We enforce this as it increases the chances of consecutive 
    ! memory layout.
    do iEl = 1 , N_Elec

      ! Associate so we put pivoting table directly where it should be
      call rgn_assoc(r_tmp, Elecs(iEl)%inDpvt)
      r_tmp%n = 0

      ! Put electrode in device orbitals in the region
      do i = 1, Elecs(iEl)%o_inD%n
        if ( .not. rgn_push(r_tmp, rpvt(Elecs(iEl)%o_inD%r(i))) ) &
            call die('tri_init: error on pushing values')
      end do

      ! Sort according to device input (this sorts inDpvt)
      call rgn_sort(r_tmp)
      ! It is clear that the inDpvt array is a sorted array
      Elecs(iEl)%inDpvt%sorted = .true.

      if ( r_tmp%n /= Elecs(iEl)%inDpvt%n ) then
        call die('tri_init: error on inDpvt update')
      end if

      ! Copy back the device in sorted device region order
      do i = 1, Elecs(iEl)%o_inD%n
        Elecs(iEl)%o_inD%r(i) = r_oDev%r(r_tmp%r(i))
      end do

      ! Copy this information to the El + D
      i = r_oElpD(iEl)%n
      n = Elecs(iEl)%o_inD%n
      r_oElpD(iEl)%r(i-n+1:i) = Elecs(iEl)%o_inD%r(1:n)

    end do
    call rgn_nullify(r_tmp)

    ! Clean-up memory
    deallocate(rpvt)

    DevTri%name = '[TRI] device region'

    if ( IONode ) then
       
       ! Print out stuff
       call rgn_print(DevTri, seq_max = 8 , repeat = .true.)
       ! Print out memory estimate
       els = nnzs_tri_i8b(DevTri%n,DevTri%r)
       ! check if there are overflows
       if ( els > huge(1) ) then
         write(*,'(a,i0)') 'Elements: ', els
         write(*,'(a,i0)') 'Max: ', huge(1)
         call die('tbt: Memory consumption is too large, try &
             &another pivoting scheme.')
       end if
       write(*,'(a,i0)') 'tbt: Matrix elements in BTD: ', els

       write(*,'(/,a)') 'tbt: Electrodes tri-diagonal matrices'
       do i = 1 , N_Elec
          call rgn_print(ElTri(i), seq_max = 8 , repeat = .true.)
       end do
       
    end if

    call timer('tri-init',2)

  end subroutine tbt_tri_init

  subroutine tbt_tri_print_opti(na_u,lasto,r_oDev,N_Elec)
    use parallel, only : IONode
    use m_region
    use m_verbosity, only : verbosity
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(tRgn), intent(in) :: r_oDev
    integer, intent(in) :: N_Elec

    integer :: cum_sum, li, i, off
    type(tRgn) :: ro, ra

    if ( .not. IONode ) return
    if ( verbosity < 2 ) return
    if ( N_Elec > 2 ) return
    if ( DevTri%n <= 2 ) return

    ! Currently we don't do this, it is rarely used, I guess
    ! and it would be better for people to determine it externally.
    return

    ! In case we have more than two tri-mat regions we can advice
    ! the user to a minimal tri-mat matrix

    ! In fact we should do this analysis for a sparsity pattern
    ! without any self-energies added, this would correctly
    ! get the size of the couplings and can shrink the 
    ! transmission region based on the mininum connections
    ! the self-energies will follow.

    ! The first thing is to find the two parts where the electrodes are
    ! "living".
    ! Then finally we can look at regions in between the two
    ! parts and suggest the minimal one.

    ! 1. r_oDev is the pivoted array
    !    So we simply need to find the first/last index of both the
    !    electrodes.
    

    ! We take the minimal region in the middle
    li = 0
    cum_sum = huge(1)
    do i = 2 , DevTri%n - 1
       if ( DevTri%r(i) < cum_sum ) then
          li = i
          cum_sum = DevTri%r(i)
       end if
    end do

    if ( DevTri%r(1) + DevTri%r(2) < cum_sum ) then
       li = 0
       cum_sum = DevTri%r(1) + DevTri%r(2)
       call rgn_list(ro,cum_sum,r_oDev%r)
    end if

    if ( DevTri%r(DevTri%n-1) + DevTri%r(DevTri%n) < cum_sum ) then
       li = 0
       cum_sum = DevTri%r(DevTri%n-1) + DevTri%r(DevTri%n)
       call rgn_list(ro,cum_sum,r_oDev%r(r_oDev%n-cum_sum+1:))
    end if

    if ( li > 0 ) then
       off = 1
       do i = 1 , li - 1
          off = off + DevTri%r(i)
       end do
       call rgn_list(ro,cum_sum,r_oDev%r(off:))
    end if

    call rgn_Orb2Atom(ro,na_u,lasto,ra)
    call rgn_delete(ro)

    ! Sort transmission region atoms
    call rgn_sort(ra)

    write(*,*) ''
    write(*,'(a)') 'tbt: Suggested atoms for fastest transmission calculation:'

    ra%name = '[A]-Fast transmission'
    call rgn_print(ra, seq_max = 12)
    write(*,*) ''

    ! Clean-up
    call rgn_delete(ra)

  end subroutine tbt_tri_print_opti

  function fold_elements(N_tri,tri) result(elem)
    integer, intent(in) :: N_tri, tri(N_tri)
    integer :: elem, i, tmp

    elem = 0
    tmp = 0
    do i = 1 , N_tri - 1
       tmp = tri(i)**2
       tmp = tmp + tri(i)*( tri(i) + 2 * tri(i+1) )
       tmp = tmp + tri(i+1) ** 2
       elem = max(elem,tmp)
    end do

  end function fold_elements


  subroutine tbt_tri_analyze( dit, sp, N_Elec, Elecs, &
       cell, na_u, xa, lasto, r_pvt, method)
    
    use parallel, only : IONode
    use fdf, only : fdf_get, leqi
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use class_OrbitalDistribution
    use class_Sparsity
    
    use m_sparsity_handling

    use m_pivot

    use m_ts_pivot, only: crt_el_priority
    use m_ts_electype
    use m_ts_sparse, only : ts_sp_calculation

    use m_ts_tri_common
    use m_ts_rgn2trimat

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp ! the local sparse pattern
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)
    type(tRgn), intent(inout) :: r_pvt

    ! The method used to partition the BTD format
    integer, intent(in) :: method
                                 
    type(Sparsity) :: tmpSp1, tmpSp2

    integer :: iEl, no

    integer :: n, n_nzs
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    character(len=*), parameter :: fmt = '(/,''TBT.BTD.Pivot.Device '',a,''+'',a)'
    character(len=64) :: fmethod
    character(len=4) :: corb
    
    ! Regions used for sorting the device region
    type(tRgn) :: r_tmp, start, r_El, full, priority, r_apvt
    integer :: orb_atom
    logical :: one_orb

    ! Capture the min memory pivoting scheme
    character(len=64) :: min_mem_method
    real(dp) :: min_mem

    call timer('TBT-analyze', 1)

    ! Write out all pivoting etc. analysis steps
    if ( IONode ) write(*,'(/,a)') 'tbt: BTD analysis'

    ! Copy over the sparse matrix to tmpSp1
    tmpSp1 = sp

    min_mem = huge(1._dp)
    min_mem_method = 'TOO LARGE'

    call rgn_Orb2Atom(r_pvt, na_u, lasto, r_apvt)

    ! Attach the sparsity pattern of the orbitals
    ! later (tmpSp2) may be atom
    call attach(tmpSp1, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nrows_g = no , nnzs = n_nzs )
    ! Check whether we are dealing with a 1-orbital system
    one_orb = no == na_u

    fmethod = 'orb+none'
    if ( IONode ) write(*,fmt) 'orb','none'
    call tri(r_pvt)

    orb_atom_switch: do orb_atom = 1 , 2
    ! The user can skip the orbital analysis if it takes too long
    if ( orb_atom == 1 ) then
       corb = 'orb'
          
       ! We default to only looking at the atomic sparsity
       ! pattern. This is *much* faster and does provide
       ! a very near optimal sparse pattern. 
       ! The user can select to do both.
       if ( one_orb ) then
          ! Same as 'atom'
          corb = 'atom'
        else
          if ( leqi(fdf_get('TBT.BTD.Analyze','atom'),'atom') ) cycle orb_atom_switch
       end if

       call rgn_copy(r_pvt,full)

       tmpSp2 = tmpSp1

    else
       if ( one_orb ) exit
       corb = 'atom'

       ! Convert the sparsity pattern to the atom
       call SpOrb_to_SpAtom(dit,tmpSp1,na_u,lasto,tmpSp2)
       ! *** the distribution will always
       !     be bigger than for the atoms, hence we need
       !     not re-construct it ***

       ! Reduce the searching place of atoms
       call rgn_copy(r_apvt, full)

    end if

    ! Sort the device region
    call rgn_sort(full)

    n = nrows_g(tmpSp2)

    ! Create priority list
    call rgn_init(priority,n)
    call crt_El_priority(N_Elec,Elecs,priority, &
         na_u,lasto,is_orb = orb_atom == 1 )

    ! Attach the sparsity pattern of the orbitals
    call attach(tmpSp1, n_col = ncol, list_ptr = l_ptr, &
         list_col = l_col , nrows_g = no , nnzs = n_nzs )


    ! *** Start analyzing sparsity pattern
    do iEl = 1, N_Elec
       
       if ( orb_atom == 1 ) then
          call rgn_copy(Elecs(iEl)%o_inD, start)
       else
          ! transfer to atom
          call rgn_Orb2Atom(Elecs(iEl)%o_inD,na_u,lasto,start)
       end if

       fmethod = trim(corb)//'+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),trim(Elecs(iEl)%name)
       call sp_pvt(n,tmpSp2,r_tmp, PVT_CONNECT, sub = full, start = start)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+rev-'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'rev-'//trim(Elecs(iEl)%name)
       call rgn_reverse(r_tmp)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+CM+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'CM+'//trim(Elecs(iEl)%name)
       call sp_pvt(n,tmpSp2,r_tmp, PVT_CUTHILL_MCKEE, sub = full, start = start)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+rev-CM+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'rev-CM+'//trim(Elecs(iEl)%name)
       call rgn_reverse(r_tmp)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+CM+priority+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'CM+priority+'//trim(Elecs(iEl)%name)
       call sp_pvt(n,tmpSp2,r_tmp, PVT_CUTHILL_MCKEE, sub = full, start = start, &
            priority = priority%r)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+rev-CM+priority+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'rev-CM+priority+'//trim(Elecs(iEl)%name)
       call rgn_reverse(r_tmp)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+PCG+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'PCG+'//trim(Elecs(iEl)%name)
       call sp_pvt(n,tmpSp2,r_tmp, PVT_PCG, sub = full, start = start)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+rev-PCG+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'rev-PCG+'//trim(Elecs(iEl)%name)
       call rgn_reverse(r_tmp)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+PCG+priority+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'PCG+priority+'//trim(Elecs(iEl)%name)
       call sp_pvt(n,tmpSp2,r_tmp, PVT_PCG, sub = full, start = start, priority = priority%r)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

       fmethod = trim(corb)//'+rev-PCG+priority+'//trim(Elecs(iEl)%name)
       if ( IONode ) write(*,fmt) trim(corb),'rev-PCG+priority+'//trim(Elecs(iEl)%name)
       call rgn_reverse(r_tmp)
       if ( orb_atom == 1 ) then
          call tri(r_tmp)
       else
          call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
          call tri(r_El)
       end if

    end do

    call rgn_delete(start)


    fmethod = trim(corb)//'+GPS'
    if ( IONode ) write(*,fmt) trim(corb),'GPS'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GPS'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GPS'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+GPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'GPS+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GPS, sub = full, priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GPS+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

#if defined(TS_PVT_GGPS) || defined(TBT_PVT_GGPS)
    fmethod = trim(corb)//'+GGPS'
    if ( IONode ) write(*,fmt) trim(corb),'GGPS'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GGPS, sub = full)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GGPS'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GGPS'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+GGPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'GGPS+priority'
    call sp_pvt(n,tmpSp2,r_tmp, PVT_GGPS, sub = full, priority = priority%r)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if

    fmethod = trim(corb)//'+rev-GGPS+priority'
    if ( IONode ) write(*,fmt) trim(corb),'rev-GGPS+priority'
    call rgn_reverse(r_tmp)
    if ( orb_atom == 1 ) then
       call tri(r_tmp)
    else
       call rgn_atom2orb(r_tmp,na_u,lasto,r_El)
       call tri(r_El)
    end if
#endif

    end do orb_atom_switch

    call rgn_delete(r_tmp,r_El,full,priority)

    call delete(tmpSp1) ! clean up
    call delete(tmpSp2)

    if ( IONode ) then
       write(*,*) ! new-line
       write(*,*) ! new-line
       write(*,'(a)') ' **********'
       write(*,'(a)') ' *  NOTE  *'
       write(*,'(a)') ' **********'
       if ( trim(min_mem_method) == 'TOO LARGE' ) then
         write(*,'(a)') ' All pivoting methods requires more elements than can be allocated'
         write(*,'(a)') ' Therefore you cannot run your simulation using TBtrans'
         write(*,'(a)') ' If you could reduce the device region [TBT.Atoms.Device] you'
         write(*,'(a)') ' may be able to run this system.'
       else
         write(*,'(a)') ' This minimum memory pivoting scheme may not necessarily be the'
         write(*,'(a)') ' best performing algorithm!'
         write(*,'(a,/)') ' You should analyze the pivoting schemes!'
         write(*,'(a)') ' Minimum memory required pivoting scheme:'
         write(*,'(a,a)') '  TBT.BTD.Pivot.Device ', trim(min_mem_method)
         write(*,'(a,en11.3,a)') '  Memory: ', min_mem, ' GB'
       end if
       write(*,*) ! new-line
    end if
    
    call timer('TBT-analyze', 2)

  contains

    ! Print out all relevant information for this
    ! pivoting scheme
    subroutine tri(r_pvt)
      use m_pivot_methods, only : bandwidth, profile
      type(tRgn), intent(in) :: r_pvt

      type(tRgn) :: cur, cTri

      integer :: bw
      ! Possibly very large numbers
      integer(i8b) :: prof, els
      logical :: is_suitable
      real(dp) :: total

      call rgn_copy(r_pvt, cur)

      ! Create a new tri-diagonal matrix, do it in parallel
      call ts_rgn2TriMat(N_Elec, Elecs, .true., &
           dit, tmpSp1, cur, ctri%n, ctri%r, &
           method, 0, par = .true. )
      
      ! Sort the pivoting table for the electrodes
      ! such that we reduce the Gf.Gamma.Gf
      ! However, this also makes it easier to
      ! insert the self-energy as they become consecutive
      ! in index, all-in-all, win-win!
      call ts_pivot_tri_sort_El(nrows_g(tmpSp1), cur, N_Elec, Elecs, ctri)

      bw   = bandwidth(no,n_nzs,ncol,l_ptr,l_col,cur)
      prof = profile(no,n_nzs,ncol,l_ptr,l_col,cur)
      if ( IONode ) then
         write(*,'(tr3,a,t23,i10,/,tr3,a,t13,i20)') &
              'Bandwidth: ',bw,'Profile: ',prof
      end if

      ! Calculate size of the tri-diagonal matrix
      els = nnzs_tri_i8b(ctri%n,ctri%r)
      ! check if there are overflows
      if ( els > huge(1) ) then
        write(*,'(tr3,a,i0,'' / '',i0)')'*** Number of elements exceeds integer limits [elements / max] ', &
            els, huge(1)
        write(*,'(tr3,a)')'*** Will not be able to use this pivoting scheme!'
        is_suitable = .false.
      else
        is_suitable = .true.
      end if

      if ( IONode ) then
        call rgn_print(ctri, name = 'BTD partitions' , &
            seq_max = 10 , indent = 3 , repeat = .true. )

        write(*,'(tr3,a,i0,'' / '',f10.3)') &
            'BTD matrix block size [max] / [average]: ', &
            maxval(ctri%r), sum(real(ctri%r)) / ctri%n

        total = real(r_pvt%n, dp) ** 2
        write(*,'(tr3,a,f9.5,'' %'')') &
            'BTD matrix elements in % of full matrix: ', &
            real(els,dp)/total * 100._dp
      end if

      total = size2gb(els) * 2
      if ( IONode ) then
        write(*,'(tr3,a,t39,en11.3,a)') 'Rough estimation of MEMORY: ', &
            total,' GB'
      end if
      if ( total < min_mem .and. is_suitable ) then
        min_mem = total
        min_mem_method = fmethod
      end if

      call rgn_delete(ctri, cur)
      
    end subroutine tri

    function size2gb(i) result(b)
      integer(i8b) :: i
      real(dp) :: b
      b = real(i, dp) * 16._dp / 1024._dp ** 3
    end function size2gb

  end subroutine tbt_tri_analyze
  
end module m_tbt_tri_init
