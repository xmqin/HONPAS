! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!
! This code segment has been fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_sparse

  use precision, only : dp
  use geom_helper, only : ucorb, iaorb
  use parallel, only : Node, Nodes

  implicit none

  private

  public :: nsc_to_offset
  public :: offset2idx

  public :: calc_nsc

  interface list_col_correct
     module procedure list_col_correct_sp
  end interface list_col_correct
  public :: list_col_correct

  interface xij_offset
     module procedure xij_offset_sp
  end interface xij_offset
  public :: xij_offset

  interface offset_xij
     module procedure offset_xij_sp
  end interface offset_xij
  public :: offset_xij

  integer, parameter, public :: TRANSFER_ALL = -999999

contains

  subroutine nsc_to_offset(nsc,sc_off)
    ! Number of supercells in each direction
    integer, intent(in) :: nsc(3)
    ! index of offsets
    integer, pointer :: sc_off(:,:)

    ! ** local variables
    integer :: n_s
    integer :: ia, ib, ic, is

    nullify(sc_off)

    n_s = product(nsc)
    allocate(sc_off(3,n_s))

    ! Create offsets
    do ic = -nsc(3)/2 , nsc(3)/2
    do ib = -nsc(2)/2 , nsc(2)/2
    do ia = -nsc(1)/2 , nsc(1)/2
       is = offset2idx(nsc,(/ia,ib,ic/))
       sc_off(1,is) = ia
       sc_off(2,is) = ib
       sc_off(3,is) = ic
    end do
    end do
    end do
    
  end subroutine nsc_to_offset

  ! From lnsc we calculate the supercell index
  function offset2idx(nsc,tm) result(idx)
    integer, intent(in) :: nsc(3), tm(3)
    integer :: idx
    integer :: ia, ib, ic

    ! We need to decipher all the stuff
    idx = 1
    if ( all(tm == 0) ) then
       ! We still have the basic unit cell in
       ! the first index!
       return
    end if

    do ic = -nsc(3)/2 , nsc(3)/2
    do ib = -nsc(2)/2 , nsc(2)/2
    do ia = -nsc(1)/2 , nsc(1)/2
       ! correct for the misalignment of the unit-cells...
       if ( all((/ic,ib,ia/) == 0 ) ) idx = idx - 1
       idx = idx + 1
       if ( tm(1) == ia .and. &
            tm(2) == ib .and. &
            tm(3) == ic ) return
    end do
    end do
    end do
    idx = 0

  end function offset2idx

  subroutine calc_nsc(cell,na_u,xa,lasto,xij_2D,nsc,Bcast)
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: cell(3,3) ! The unit cell of system
    integer, intent(in)  :: na_u ! Unit cell atoms
    real(dp), intent(in) :: xa(3,na_u) ! Atomic coordinates (needed for RemZConnection & RemUCellDistances)
    integer, intent(in)  :: lasto(0:na_u) ! last orbital number of equivalent atom
    type(dSpData2D), intent(inout) :: xij_2D ! differences with unitcell, differences with unitcell
    logical, intent(in), optional :: Bcast
! ***********************
! * OUTPUT variables    *
! ***********************
    integer, intent(out) :: nsc(3)

! ***********************
! * LOCAL variables     *
! ***********************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: xij(:,:)

    logical :: lBcast
    real(dp) :: recell(3,3)
    real(dp) :: xijo(3), xc
    integer :: ia, ja
    integer :: no_l, no_u, lio, io, jo, ind
#ifdef MPI
    integer :: MPIerror, tnsc(3)
#endif

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    dit => dist(xij_2D)
    sp  => spar(xij_2D)
    xij => val (xij_2D)

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    ! Initialize the transfer cell to:
    nsc(:) = 0

    ! Prepare the cell to calculate the index of the atom
    call reclat(cell,recell,0) ! Without 2*Pi
    
    do lio = 1 , no_l
#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(dit,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)
          jo = ucorb(l_col(ind),no_u)
          ja = iaorb(jo,lasto)
          xijo(:) = xij(:,ind) - ( xa(:,ja)-xa(:,ia) )

          ! Loop over directions
          ! recell is already without 2*Pi
          xc = sum(xijo(:) * recell(:,1))
          nsc(1) = max(nsc(1),abs(nint(xc)))
          xc = sum(xijo(:) * recell(:,2))
          nsc(2) = max(nsc(2),abs(nint(xc)))
          xc = sum(xijo(:) * recell(:,3))
          nsc(3) = max(nsc(3),abs(nint(xc)))
       end do
    end do
    
#ifdef MPI
    if ( lBcast ) then
       tnsc = nsc
       call MPI_AllReduce(tnsc,nsc,3,MPI_Integer, &
            MPI_Max,MPI_Comm_World,MPIerror)
    end if
#endif
    
    ! Calculate the actual number of super-cells
    nsc(:) = 2 * nsc(:) + 1
    
  end subroutine calc_nsc


  ! Create an index array containing the unit-cell expansions
  ! This means we can do this:
  !   xij(:,ind) = ucell(:,1) * tm(1) + ucell(:,2) * tm(2) + ucell(:,3) * tm(3) + xa(:,iaorb(jo))-xa(:,iaorb(io))
  ! to create the full xij array...
  subroutine list_col_correct_sp(cell, nsc, na_u, xa, lasto, &
       xij_2D, Bcast)

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: cell(3,3) ! The unit cell of system
    ! Number of supercells in each direction ** MUST be corrected using nsc_to_offests
    integer, intent(in)  :: nsc(3) ! Number of supercells in each direction
    integer, intent(in)  :: na_u ! Unit cell atoms
    integer, intent(in)  :: lasto(0:na_u) ! Last orbital of atom
    real(dp), intent(in) :: xa(3,na_u) ! atom positions
    type(dSpData2D), intent(inout) :: xij_2D
    logical, intent(in), optional :: Bcast
! ***********************
! * LOCAL variables     *
! ***********************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: xij(:,:)

    logical :: lBcast
    integer :: no_l, no_u, lio, io, jo, ind, ia, ja, is, tm(3)
    real(dp) :: xijo(3), rcell(3,3)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    dit => dist(xij_2D)
    sp  => spar(xij_2D)
    xij => val (xij_2D)

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)
    
    ! Prepare the cell to calculate the index of the atom
    call reclat(cell,rcell,0) ! Without 2*Pi

    !err = 0._dp
    do lio = 1 , no_u
#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(dit,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)
       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          jo = ucorb(l_col(ind),no_u)
          ja = iaorb(jo,lasto)

          xijo(:) = xij(:,ind) - ( xa(:,ja) - xa(:,ia) )
          tm(:) = nint( matmul(xijo,rcell) )

          ! get supercell index
          is = offset2idx(nsc,tm) - 1

          if ( is < 0 ) then
             ! Index not found
             call die('Error in index, possible shifting')
          end if

          ! Correct index
          l_col(ind) = is * no_u + jo

          ! The error on this conversion
          ! tends to be on the order of 1e-14
          ! xijo = ucell(:,1) * tm(1) &
          !      + ucell(:,2) * tm(2) &
          !      + ucell(:,3) * tm(3) &
          !      + xa(:,ja) - xa(:,ia)
          ! err = max(maxval(abs(xijo - xij(:,ind))),err)

       end do
    end do

  end subroutine list_col_correct_sp

  ! Routine used for calculation the supercell offsets using
  ! the xij_2D array.
  ! Basically this routine runs through all xijo elements
  ! and calculates the cell distance:
  !    R(col(ind) % no_u) = xij(ind) - (xa(j) - xa(i))
  ! This should result in vectors R(...) which are integer
  ! multiples of the lattice vectors:
  !   cell X = R
  ! returns X as integers.
  ! This routine will quit/exit if two supercell indices
  !   col(ind) % no_u does not yield the same R vector.
  subroutine xij_offset_sp(cell,nsc,na_u,xa,lasto, &
       xij_2D,isc_off,Bcast)
    
    use intrinsic_missing, only : VNORM
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use alloc
#ifdef MPI
    use mpi_siesta
#endif

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: cell(3,3)
    ! Number of super-cells in each direction
    integer, intent(in) :: nsc(3)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    ! vectors from i->J
    type(dSpData2D), intent(inout) :: xij_2D
    logical, intent(in), optional :: Bcast

! **********************
! * OUTPUT variables   *
! **********************
    integer, pointer :: isc_off(:,:)

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: xij(:,:)

    logical :: lBcast
#ifdef MPI
    integer, allocatable :: ioff_2(:,:)
    integer :: MPIerror, iNode
#endif
    integer :: n_s, tm(3), is
    integer :: no_l, no_u
    integer :: lio, io, ind, ia, ja
    integer :: i1, i2, i3, isc(3), s1, s2, s3
    real(dp) :: xijo(3), rcell(3,3)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    dit => dist(xij_2D)
    sp  => spar(xij_2D)
    xij => val(xij_2D)

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    ! Number of super cells
    n_s = product(nsc)
    ! Calculate the offsets (instead of using xij)
    call re_alloc(isc_off,1,3,1,n_s)

    ! Initialize the integer offset
    isc_off(:,:) = 0
    
    ! Prepare the cell to calculate the index of the atom
    call reclat(cell,rcell,0) ! Without 2*Pi

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(dit,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from one)
          is = (l_col(ind) - 1)/no_u + 1

          xijo(:) = xij(:,ind) - ( xa(:,ja) - xa(:,ia) )

          tm(:) = nint( matmul(xijo,rcell) )

          if ( any(tm(:) /= isc_off(:,is)) .and. any(isc_off(:,is)/=0) ) then
             write(*,'(2(a10,2(tr1,i6)))')'r,C',io,l_col(ind),'ia,ja',ia,ja
             write(*,'(a,i3,2(tr2,3(tr1,i3)))') 'is, tm, old_tm: ',is, tm, isc_off(:,is)
             write(*,'(a,4(tr1,f10.5))') 'xij - ucell*tm, |V|: ',xijo(:) - &
                  cell(:,1) * tm(1) - cell(:,2) * tm(2) - cell(:,3) * tm(3), vnorm(xijo)
             write(*,'(2(a10,3(tr1,f10.5)))') 'xijo: ',xijo(:), &
                  'ucell: ',cell(:,1) * tm(1) + cell(:,2) * tm(2) + cell(:,3) * tm(3)
             call die('Error on transfer matrix indices...')
          else
             isc_off(:,is) = tm(:)
          end if

       end do
    end do

#ifdef MPI
    if ( lBcast ) then
    ! Reduce all sc_off
    allocate(ioff_2(3,n_s))

    do iNode = 0 , Nodes - 1
       if ( iNode == Node ) ioff_2(:,:) = isc_off(:,:)
       call MPI_Bcast(ioff_2(1,1),3*n_s, &
            MPI_Integer, iNode, MPI_Comm_World, MPIerror)

       ! Search for differing sc_off
       do is = 1 , n_s
          if ( all(isc_off(:,is) == 0) ) then
             ! it hasn't been found locally, just copy
             isc_off(:,is) = ioff_2(:,is)
          else if ( all(ioff_2(:,is) == 0) ) then
             ! do nothing, this is not set on other node
          else if ( any(isc_off(:,is) /= ioff_2(:,is)) ) then
             print '(2(tr1,a,3(tr1,i2)))','Local:',isc_off(:,is),'Found:',ioff_2(:,is)
             ! We do not have the same supercell indices
             call die('Supercell reduction is erroneous. &
                  &Please consult the developers.')
          end if
          
       end do

    end do

    deallocate(ioff_2)
    end if
#endif

    ! Correct isc_off for missing elements, otherwise
    ! we could have an ambiguity.
    ! Here we force the unit-cell to be the first index
    do is = 2 , n_s
       if ( all(isc_off(:,is) == 0) ) then
          lBcast = .false.
          ! We need to find a non-existing one
          ! It appears siesta sets it up "backwards"
          ! In a rather obscure way...
          do i3 = 0 , nsc(3)/2
          do s3 = 1 , -1 , -2
             isc(3) = i3 * s3
             do i2 = 0, nsc(2)/2
             do s2 = 1 , -1 , -2
                isc(2) = i2 * s2
                do i1 = 0, nsc(1)/2
                sc_1: do s1 = 1 , -1 , -2
                   isc(1) = i1 * s1
                   if ( all(isc == 0) ) cycle
                   do io = 1 , n_s
                      if ( all(isc == isc_off(:,io)) ) cycle sc_1
                   end do
                   isc_off(:,is) = isc
                   lBcast = .true.
                   exit
                end do sc_1
                if ( lBcast ) exit
                end do
             if ( lBcast ) exit
             end do
             if ( lBcast ) exit
             end do
          if ( lBcast ) exit
          end do
          if ( lBcast ) exit
          end do
       end if
    end do
    
  end subroutine xij_offset_sp

  subroutine offset_xij_sp(cell,n_s,isc_off, na_u,xa,lasto, &
       xij_2D,Bcast)
    
    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use alloc

! **********************
! * INPUT variables    *
! **********************
    ! Unit cell
    real(dp), intent(in) :: cell(3,3)
    ! Number of super-cells
    integer, intent(in) :: n_s
    ! Integer supercells
    integer, intent(in) :: isc_off(3,0:n_s-1)
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Last orbital of the equivalent unit-cell atom
    real(dp), intent(in) :: xa(3,na_u)
    ! Last orbital of the equivalent unit-cell atom
    integer, intent(in) :: lasto(0:na_u)
    logical, intent(in), optional :: Bcast
! **********************
! * OUTPUT variables   *
! **********************
    ! vectors from i->J
    type(dSpData2D), intent(inout) :: xij_2D

! **********************
! * LOCAL variables    *
! **********************
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    real(dp), pointer :: xij(:,:)

    logical :: lBcast
    integer :: tm(3), is
    integer :: no_l, no_u
    integer :: lio, io, ind, ia, ja
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast

    dit => dist(xij_2D)
    sp  => spar(xij_2D)
    xij => val(xij_2D)

    ! Attach to the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u)

    do lio = 1 , no_l

       if ( ncol(lio) == 0 ) cycle

#ifdef MPI
       if ( lBcast ) then
          ! Transfer to global index
          io = index_local_to_global(dit,lio,Node)
       else
          io = lio
       end if
#else
       io = lio
#endif
       ia = iaorb(io,lasto)

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + ncol(lio)

          ja = iaorb(ucorb(l_col(ind),no_u),lasto)

          ! the supercell index (counting from zero)
          is = (l_col(ind) - 1)/no_u

          tm(:) = isc_off(:,is)
          xij(:,ind) = tm(1)*cell(:,1)+tm(2)*cell(:,2)+tm(3)*cell(:,3)
          xij(:,ind) = xij(:,ind) + xa(:,ja) - xa(:,ia)

       end do
    end do

  end subroutine offset_xij_sp

end module m_sparse
  
