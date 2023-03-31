! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication
module m_tbt_kregions
  
  use precision, only : dp
  use m_region

  implicit none

  private
  save 

  type :: kRegion
     ! The atoms in this periodic region
     type(tRgn) :: atm
     ! The current k-point
     integer :: ik
     ! A list of assigned k-points and weights
     real(dp), pointer :: kpt(:,:), wkpt(:)
     ! Point to the k-region which is the offset
     ! basis
     type(kRegion), pointer :: off => null()
  end type kRegion
  public :: kRegion

  ! A list of different k-regions
  ! This is a list of 0:n
  ! where 0 contains the 'Gamma' region
  ! and all subsequent n regions contain
  ! different regions which will get different k
  type(kRegion), allocatable, target :: r_k(:)
  public :: r_k
  integer, public :: n_k = 0

  ! Routines
  public :: tbt_init_kRegions, tbt_print_kRegions
  public :: kregion_k, kregion_step
  public :: calc_GS_k

contains

  subroutine tbt_init_kRegions(r_aBuf,N_Elec, Elecs, cell, &
       dit,sp, na_u, xa, lasto, nsc, isc_off)
    
    use fdf
    use fdf_extra, only : fdf_bnext
    use parallel, only : IONode
    use precision, only : dp

    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_electype

    use m_tbt_kpoint, only : kpoint, kweight, read_kgrid

    ! The buffer region
    type(tRgn), intent(in) :: r_aBuf
    ! Number of electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The unit cell
    real(dp), intent(in) :: cell(3,3)
    ! the orbital distribution for the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Atomic orbital configuration
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! Atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! the supercell information
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))

    ! ** local variables
    logical :: TRS
    type(tRgn) :: r1, r2
    integer :: i, il, iEl
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g
    real(dp) :: k(3)

    n_k = 0 

    ! If the user does not request k->k', do not read in
    if ( .not. fdf_get('TBT.Region.k',.false.) ) return

     ! Warn the user about possible mis-use of regions
    if ( IONode ) then
       write(*,'(/,a)') 'tbt: WARNING'
       write(*,'(a)') 'tbt: Using k-regions is not recommended &
            &unless you really know what you are doing.'
       write(*,'(a)') 'tbt: WARNING'
    end if

    ! First change original k-points to reciprocal units
    n_k = size(kweight)
    do i = 1 , n_k
       k(:) = kpoint(:,i)
       call kpoint_convert(cell,k,kpoint(:,i),1)
    end do
    n_k = 0

    if ( .not. fdf_block('TBT.Atoms.k',bfdf) ) &
         call die('tbt: k-regions, could not read &
         &different periodic regions. See %block TBT.Atoms.k')
    call crt_read_k(bfdf,na_u)

    ! Get the atoms which connects to a super-cell to check
    ! whether they are all in a k-region.
    call get_connect_k(N_Elec, Elecs, &
         dit, sp, na_u, xa, lasto, nsc, isc_off, r1)

    ! Assert that the regions all have at least one atom
    do il = 1 , n_k
       ! Sort to fasten the search, and stream-line
       ! their usage
       call rgn_sort(r_k(il)%atm)
       if ( r_k(il)%atm%n == 0 ) then
          call die('Region: '//trim(r_k(il)%atm%name)//' &
               &is empty, this is not allowed.')
       end if
       i = r1%n
       do while ( i > 0 ) 
          if ( in_rgn(r_k(il)%atm,r1%r(i)) ) then
             if ( 0 > rgn_pop(r1,i) ) iEl = 0
          end if
          i = i - 1
       end do
    end do

    ! If there are still elements in r1 we definitely have
    ! an error
    if ( r1%n /= 0 ) then
       r1%name = 'Periodic atoms'
       call rgn_print(r1)
       call die('Could not assert all periodic atoms to be &
            &k-region.')
    end if
    call rgn_delete(r1)

    ! Assert that each electrode "only" exists in one
    ! periodicity region (we do not allow the user
    ! doing weird partially non-periodic electrodes)
    do iEl = 1 , N_Elec
       ! Create a region of this electrode atoms
       i = Elecs(iEl)%idx_a
       call rgn_range(r1,i,i-1 + TotUsedAtoms(Elecs(iEl)) )
       do il = 1 , n_k
          ! Check that the entire electrode exists in
          ! any one of the regions.
          call rgn_union(r_k(il)%atm,r1,r2)
          ! if there is no overlap
          if ( r2%n == r_k(il)%atm%n + r1%n ) cycle
          ! if there is complete overlap
          if ( r2%n == r_k(il)%atm%n ) exit
          ! This means that r2%n > 0 /= r1%n
          ! this is partially periodic electrode and not
          ! allowed.
          if ( IONode ) then
          do i = 1 , n_k
             call rgn_print(r_k(i)%atm)
          end do
          write(*,*)'tbt: Periodicity regions are ill-formatted'
          write(*,*)'Each electrode MUST only be present in any one &
               &periodicity region at any time'
          write(*,*)'Dividing the electrode into two periodicity &
               &regions will create erroneous self-energy terms.'
          write(*,*)'Please assert input.'
          end if
          call die('Periodicity regions ill-formatted, check output')
       end do
    end do

    ! At this point we have all different k-regions
    ! We can now create the Gamma-region
    ! When the user creates k-regions, the user
    ! will intrinsically force the Gamma-point on atoms
    ! not specified.
    ! This may come as a surprise but is trivially accomplished
    ! while allowing greater flexibility. 
    ! If the user does not want this, simply add all atoms
    ! to list
    call rgn_copy(r_k(1)%atm,r1)
    do il = 2 , n_k
       call rgn_append(r1,r_k(il)%atm,r1)
    end do
    call rgn_range(r_k(0)%atm,1,na_u)
    call rgn_complement(r1,r_k(0)%atm,r_k(0)%atm)
    ! Probably we should remove the buffer atoms
    r_k(0)%atm%name = 'default'
    call rgn_delete(r1,r2)

    ! Sort all regions (makes the setup of the Hamiltonian faster)
    do il = 0 , n_k
       if ( r_aBuf%n > 0 ) then
          ! Remove all buffer atoms
          call rgn_complement(r_aBuf,r_k(il)%atm,r_k(il)%atm)
       end if
       call rgn_sort(r_k(il)%atm)
    end do
    
    ! The 0 region contains the kgrid used as "backend"
    r_k(0)%kpt  => kpoint(:,:)
    r_k(0)%wkpt => kweight(:)

    TRS = .not. fdf_get('SpinSpiral',.false.)
    TRS = fdf_get('TBT.Symmetry.TimeReversal',TRS)

    ! First we initialize all k-regions to
    ! be the regular grid
    do il = 1 , n_k
       
       g = 'TBT.Atoms.k.'//trim(r_k(il)%atm%name)
       if ( fdf_block(trim(g),bfdf) ) then

          ! First we figure out if this is an
          ! offset, if it is, then TRS cannot be
          ! applied. Else TRS can be applied as any
          ! normal standard. This of course
          ! assumes that TRS applies when 
          ! momentum losses occur, which might not be
          ! the case.
          ! This is totally the users responsibility ;)
          do while ( fdf_bnext(bfdf,pline) ) 
             if ( fdf_bnnames(pline) == 0 ) cycle
             g = fdf_bnames(pline,1)
             ! If we have offset in it,
             if ( leqi(g,'offset') .or. leqi(g,'off') ) then
                g = fdf_bnames(pline,2)
                if ( leqi(g,'default') ) then
                   ! The offset is the default k-grid
                   r_k(il)%off => r_k(0)
                else
                   ! find the name in the list
                   do i = 1 , n_k
                      if ( leqi(g,r_k(i)%atm%name) ) then
                         r_k(il)%off => r_k(i)
                         exit
                      end if
                   end do
                   if ( .not. associated(r_k(il)%off) ) then
                      call die('Could not find k-region '//trim(g)// &
                           ' please ensure that the region with the given &
                           &name exists.')
                   end if
                end if
             end if
          end do

          g = 'TBT.Atoms.k.'//trim(r_k(il)%atm%name)

          ! The k-specification exists
          ! Read in the information
          if ( associated(r_k(il)%off) ) then
             ! The k-points using offsets cannot be ensured
             ! TRS reduced. Hence we force it to be false
             ! if using offset.
             call read_kgrid(trim(g), &
                  .false.,cell,r_k(il)%kpt,r_k(il)%wkpt, &
                  is_b = .true. )
          else
             call read_kgrid(trim(g), &
                  TRS,cell,r_k(il)%kpt,r_k(il)%wkpt, &
                  is_b = .true. )
          end if

       else

          ! We default a non-mentioned region
          ! to have the same grid as the default setup
          r_k(il)%off => r_k(0)
          nullify(r_k(il)%kpt) ; allocate(r_k(il)%kpt(3,1))
          r_k(il)%kpt = 0._dp
          nullify(r_k(il)%wkpt) ; allocate(r_k(il)%wkpt(1))
          r_k(il)%wkpt = 1._dp

       end if

    end do

  end subroutine tbt_init_kRegions

  subroutine get_connect_k(N_Elec, Elecs, &
       dit,sp, na_u, xa, lasto, nsc, isc_off, r_per)
    
    use fdf
    use precision, only : dp

    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_electype

    ! Number of electrodes
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the orbital distribution for the sparsity pattern
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Atomic orbital configuration
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! Atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    ! the supercell information
    integer, intent(in) :: nsc(3), isc_off(3,product(nsc))
    ! The atoms which are periodic
    type(tRgn), intent(inout) :: r_per

    ! ** local variables
    integer :: no_l, no_u
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: na
    integer, allocatable :: a_list(:)
    logical :: periodic
    integer :: io, ia

    call attach(sp,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    if ( no_u /= no_l ) call die('Error in fold creation')

    ! Gather a list of all atoms which only have Gamma connections
    na = 0
    allocate(a_list(na_u))
    do ia = 1 , na_u

       periodic = .false. ! in-case l_ncol == 0
       do io = lasto(ia-1) + 1 , lasto(ia)
          if ( l_ncol(io) == 0 ) cycle
          periodic = any(l_col(l_ptr(io)+1:l_ptr(io)+l_ncol(io)) > no_u)
          if ( periodic ) exit
       end do
       if ( periodic ) then
          na = na + 1
          a_list(na) = ia
       end if

    end do

    call rgn_list(r_per,na,a_list)

    deallocate(a_list)

    ! Sort it to search better
    call rgn_sort(r_per)

  end subroutine get_connect_k

  subroutine crt_read_k(bfdf,na_u)
    
    use fdf
    use fdf_extra
    use precision, only : dp

    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_electype

    ! INPUT
    type(block_fdf), intent(inout) :: bfdf
    integer, intent(in) :: na_u

    ! ** local variables
    type(tRgn) :: r1
    integer :: i, il, ic
    type(parsed_line), pointer :: pline => null()
    character(len=50) :: g
    character(len=50), allocatable :: rlist(:)
    logical :: found

    call fdf_brewind(bfdf)

    ! We assume that the user will never create more than
    ! the initial number of different k-regions plus 50
    do while ( fdf_bnext(bfdf,pline) ) 
       il = il + 1
    end do
    allocate(rlist(il))
    call fdf_brewind(bfdf)
    
    ! first count number of differently named regions
    n_k = 0
    do while ( fdf_bnext(bfdf,pline) ) 

       found = .false.
       if ( n_k > 0 ) then
          g = fdf_bnames(pline,1)
          do i = 1 , n_k
             if ( leqi(g,rlist(i)) ) then
                found = .true.
                exit
             end if
          end do
       end if
       if ( .not. found ) then
          n_k = n_k + 1
          if ( n_k > size(rlist) ) then
             call die('Number of differently specified k-regions exceeds &
                  &limit. Change 50 in m_tbt_k.F90')
          end if
          rlist(n_k) = g
       end if

    end do

    ! Clean-up
    deallocate(rlist)
    call fdf_brewind(bfdf)

    allocate(r_k(0:n_k))
    
    il = 0
    do while ( fdf_bnext(bfdf,pline) ) 
       
       g = fdf_bnames(pline,1)
       if ( leqi(trim(g),'default') ) then
          call die('You MUST not name a k-region default, &
               &it is used internally. Please change the name &
               &of the region.')
       end if
       
       ! Check if the name already has been read (then
       ! we accumulate the atoms)
       found = .false.
       ic = il + 1
       if ( il > 0 ) then
          do i = 1 , il
             if ( leqi(g,r_k(i)%atm%name) ) then
                ic = i
                exit
             end if
          end do
       end if
       if ( ic == il + 1 ) then
          ! we have a new name
          il = ic
          
          ! We can read in a range
          call fdf_brange(pline,r1, 1, na_u)
          if ( r1%n == 0 ) &
               call die('Could not read in anything in k-regions')
          call rgn_union(r_k(il)%atm,r1,r_k(il)%atm)
          r_k(il)%atm%name = trim(g)
          
       else
          
          call fdf_brange(pline,r1, 1, na_u)
          if ( r1%n == 0 ) &
               call die('Could not read in anything in k-regions')
          call rgn_union(r_k(ic)%atm,r1,r_k(ic)%atm)
          r_k(ic)%atm%name = trim(g)
          
       end if
       
    end do

    call rgn_delete(r1)

  end subroutine crt_read_k

    ! Subroutine for reading in both the left and right next energy point
  subroutine calc_GS_k(ispin, cE, &
       N_Elec, Elecs, k_idx, &
       nzwork, zwork)

#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum, MPI_Integer
#endif

    use m_ts_electype
    use m_ts_cctype
    use m_ts_electrode, only: calc_next_GS_Elec
      
    integer, intent(in) :: ispin
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    ! Index of the k-region that belongs to the electrode
    integer, intent(in) :: k_idx(N_Elec)
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: zwork(nzwork)

    integer :: i
    type(ts_c_idx) :: c
    real(dp) :: bkpt(3)

    c = cE

    ! TODO Move reading of the energy points
    ! directly into the subroutines which need them
    ! In this way we can save both GAA, Sigma AND Gamma arrays!!!!
    ! However, this will probably come at the expense 
    ! of doing the same "repetition" expansion twice, we can live with
    ! that!
    do i = 1 , N_Elec

       if ( Elecs(i)%out_of_core ) then
          call die('Using the k-region option requires the &
               &self-energies to be calculated when used.')
       end if

       ! If the index for the contour is negative
       ! It means that we are dealing with a Fermi
       ! charge correction
       ! If it is different from 1 we do not have an equilibrium contour
       ! In this case the energy is the eta value of the electrode
       if ( Elecs(i)%Eta > 0._dp ) then
#ifdef TBT_PHONON
         c%e = cmplx(real(cE%e,dp)**2,Elecs(i)%Eta, dp)
#else
         c%e = cmplx(real(cE%e,dp),Elecs(i)%Eta, dp)
#endif
       else
#ifdef TBT_PHONON
         c%e = cmplx(real(cE%e,dp)**2,aimag(cE%e)**2, dp)
#else
         c%e = cE%e
#endif
       end if

       ! Get k-point
       call kregion_k(k_idx(i),bkpt)

       ! This routine will automatically check
       ! (and SET) the k-point for the electrode.
       ! This is necessary for the expansion to work.
       call calc_next_GS_Elec(Elecs(i),ispin,bkpt,c%e, &
            nzwork, zwork)
    end do

  end subroutine calc_GS_k

  ! For il in [1:n_k] it returns the k-point
  ! for that k-region.
  ! For il == 0 it returns the Gamma point
  ! for il < 0 it returns the current k-point
  ! for the DEFAULT region.
  subroutine kregion_k(il,k,w)
    integer, intent(in) :: il
    real(dp), intent(out) :: k(3)
    real(dp), intent(out), optional :: w

    type(kRegion), pointer :: rK

    if ( il == 0 ) then
       k(:) = 0._dp
       if ( present(w) ) &
            w = 1._dp
       return
    else if ( il < 0 ) then
       k(:) = r_k(0)%kpt(:,r_k(0)%ik)
       if ( present(w) ) &
            w = r_k(0)%wkpt(r_k(0)%ik)
       return
    end if

    k(:) = r_k(il)%kpt(:,r_k(il)%ik)
    if ( present(w) ) &
         w = r_k(il)%wkpt(r_k(il)%ik)
    rK => r_k(il)%off
    do while ( associated(rK) )
       k(:) = k(:) + rK%kpt(:,rK%ik)
       if ( present(w) ) &
            w = w * rK%wkpt(rK%ik)
       rK => rK%off
    end do

  end subroutine kregion_k

  subroutine kregion_step()
    integer :: i
    ! step the k-region
    do i = n_k , 0 , -1
       if ( r_k(i)%ik < size(r_k(i)%wkpt) ) then
          r_k(i)%ik = r_k(i)%ik + 1
          exit
       else
          r_k(i)%ik = 1
       end if
    end do
  end subroutine kregion_step

  subroutine tbt_print_kRegions( cell )

    use parallel, only : IONode
    use m_verbosity, only : verbosity

    use m_tbt_kpoint, only : tbt_iokp

    real(dp), intent(in) :: cell(3,3)

    integer :: il, nperm
    character(len=30) :: g

    if ( .not. IONode ) return
    if ( n_k == 0 ) return

    ! Print out information regarding the k regions

    if ( verbosity > 3 ) then 
       write(*,'(/,a)') 'tbt: k-regions print out, in the permutation order'
       write(*,'(/,a)') 'tbt: Gamma region'
       write(g,'(a,i0)') 'nkpt = ',1
       call rgn_print(r_k(0)%atm, name = g , seq_max = 12 )
       
       write(*,'(/,a)') 'tbt: k-regions'
       nperm = size(r_k(0)%wkpt)
       do il = 1 , n_k
          write(g,'(a,i0)') 'nkpt = ',size(r_k(il)%wkpt)
          if ( associated(r_k(il)%off) ) then
             g = trim(g)//', offsets: '//trim(r_k(il)%off%atm%name)
          end if
          call rgn_print(r_k(il)%atm, name = g , seq_max = 12 )
          nperm = nperm * size(r_k(il)%wkpt)
       end do
    else
       write(*,'(/,a)') 'tbt: k-regions permutation order'
       nperm = size(r_k(0)%wkpt)
       do il = 1 , n_k
          nperm = nperm * size(r_k(il)%wkpt)
          write(g,'(a,i0)') ', nkpt = ',size(r_k(il)%wkpt)
          if ( associated(r_k(il)%off) ) then
             g = trim(g)//', offsets: '//trim(r_k(il)%off%atm%name)
          end if
          write(*,'(3a)')'tbt: ',trim(r_k(il)%atm%name), trim(g)
       end do
    end if
       
    write(*,'(/,2(a,i0))')'tbt: k-point &
         &permutations [default/total]: ',size(r_k(0)%wkpt),' / ',nperm

    ! Save the k-points
    nperm = size(r_k(0)%wkpt)
    call tbt_iokp(nperm,r_k(0)%kpt,r_k(0)%wkpt, &
         fend = 'TBT.DEFAULT.BKP' )
    
    do il = 1 , n_k
       nperm = size(r_k(il)%wkpt)
       call tbt_iokp(nperm,r_k(il)%kpt,r_k(il)%wkpt, &
            fend = 'TBT.'//trim(r_k(il)%atm%name)//'.BKP' )
    end do

  end subroutine tbt_print_kRegions
      
end module m_tbt_kregions
