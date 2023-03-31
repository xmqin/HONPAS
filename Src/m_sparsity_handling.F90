! Module to handle sparsity formats
! We supply routines for nullifying specific atoms
! in a sparsity format. 
! We also provide routines for removal of those entries in the equivalent
! data arrays.

! Fully created by Nick Papior Andersen, 2014

module m_sparsity_handling

  use precision, only : dp
  use class_OrbitalDistribution
  use class_Sparsity
  use class_dData1D
  use class_dSpData1D
  use class_dData2D
  use class_dSpData2D
  use geom_helper, only : iaorb, ucorb
  use m_region
  use intrinsic_missing, only: SORT_QUICK

  implicit none

  private

  public :: SpOrb_to_SpAtom
  public :: Sp_remove_crossterms
  public :: Sp_union
  public :: Sp_sort
  public :: Sp_remove_region
  public :: Sp_retain_region
  public :: Sp_remove_region2region
  public :: Sp_to_Spglobal
  interface SpData_to_Sp
    module procedure dSpData1D_to_Sp
    module procedure dSpData2D_to_Sp
  end interface SpData_to_Sp
  public :: SpData_to_Sp

  interface SpData_interp
    module procedure dSpData2D_interp
  end interface SpData_interp
  public :: SpData_interp

contains

  ! Creates a sparsity pattern in the atomic subspace.
  ! I.e. we first reduce the sparsity pattern to the Gamma
  ! sparsity pattern and create the sparsity pattern for the 
  ! atoms (it will most likely be relatively dense, but...)
  subroutine SpOrb_to_SpAtom(dit,in,na_u,lasto,out)
    ! the sparsity pattern may be distributed
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced to a sparsity pattern for
    ! atomic connections
    type(Sparsity), intent(inout) :: in
    ! number of orbitals per atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    ! The out-put sparsity is _never_ distributed
    ! downfolding to atoms is a huge decrease, and 
    ! it makes no sense to distribue it...
    type(Sparsity), intent(inout) :: out

    ! The new sparsity pattern
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, allocatable :: num(:), listptr(:), list(:)
    logical, allocatable :: la_c(:)
    integer :: ia, ja, io, lio, ind
    integer :: no_l, no_u, n_nzs
    integer :: ic

    call attach(in,nrows=no_l,nrows_g=no_u, nnzs=n_nzs, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    if ( no_u /= no_l ) call die('Error in conversion SpOrb2SpAt')

    ! Instead of doing these loops 2 times we simply allocate
    ! as though the sparsity pattern was an atomic sparsity pattern.
    ! Then we only fill in new entries as needed.
    ! This has a slight overhead of memory but should be negligeble.
    allocate(la_c(na_u))
    la_c(:) = .true.
    allocate(num(na_u))
    allocate(list(n_nzs))
    n_nzs = 0
    do ia = 1 , na_u
      ! we figure out the number of atom
      ! connections for each atom

      ! Initialize count
      ic = 0
      ! Initialize to specify where to accept entries
      do io = lasto(ia-1) + 1 , lasto(ia)

        ! Check the local orbital
        lio = index_global_to_local(dit,io)
        if ( lio <= 0 ) cycle
        if ( l_ncol(lio) == 0 ) cycle

        do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          ! Get connecting atom
          ja = orb_to_atom(l_col(ind), na_u, lasto, no_u)
          if ( la_c(ja) ) then
            la_c(ja) = .false.
            n_nzs = n_nzs + 1
            list(n_nzs) = ja
            ic = ic + 1
          end if

        end do

      end do
      
      num(ia) = ic

      ! Reset true values
      do ja = 0, ic - 1
        la_c(list(n_nzs-ja)) = .true.
      end do

    end do

    deallocate(la_c)

    ! Create listptr
    allocate(listptr(na_u))
    listptr(1) = 0
    do ia = 2 , na_u
      listptr(ia) = listptr(ia-1) + num(ia-1)
    end do

    ! Create new sparsity pattern and copy over
    call newSparsity(out,na_u,na_u,n_nzs,num,listptr,list, &
        name='Atomic ('//trim(name(in))//')', &
        ncols=na_u,ncols_g=na_u)

    ! Clean up
    deallocate(num,listptr,list)

  contains

    function orb_to_atom(orb, na_u, lasto, no_u) result(atom)
      use geom_helper, only : iaorb
      integer, intent(in) :: orb, na_u, lasto(0:na_u), no_u
      integer :: atom

      atom = (orb-1)/no_u
      atom = iaorb(orb, lasto) + atom * na_u
      
    end function orb_to_atom

  end subroutine SpOrb_to_SpAtom

  subroutine Sp_remove_crossterms(dit,in,nsc,isc_off,dir,out,r)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The supercell connections (in integers)
    integer, intent(in) :: nsc, isc_off(3,0:nsc-1)
    ! The direction in the unit-cell we wish to remove connections to
    ! All OUT OF unit-cell connections will be removed!
    integer, intent(in) :: dir
    ! The output sparsity pattern
    type(Sparsity), intent(inout) :: out
    ! The region that we restrict our limitation to
    type(tRgn), intent(in), optional :: r
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, n_nzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, indx, is

    logical, allocatable :: log_r(:)

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    allocate(log_r(no_u))
    if ( present(r) ) then
      log_r(:) = .false.
      call rgn_2logical(r, log_r)
    else
      log_r(:) = .true.
    end if

    allocate(num(no_l))
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ind,is)
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0

       if ( l_ncol(lio) /= 0 ) then

       io = index_local_to_global(dit,lio)
       if ( .not. log_r(io) ) then
          ! we are not asked to remove these cross-terms
          num(lio) = l_ncol(lio)

       else

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          if ( log_r(ucorb(l_col(ind),no_u)) ) then

             is = (l_col(ind)-1)/no_u
             if ( isc_off(dir,is) /= 0 ) cycle

          end if
          
          ! The orbital exists on the atom
          num(lio) = num(lio) + 1
          
       end do

       end if
       end if
       
    end do
!$OMP end parallel do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    n_nzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(n_nzs))
    indx = 0
    do lio = 1 , no_l

       ! This will be zero if l_ncol(lio) is... :)
       if ( num(lio) == 0 ) cycle

       io = index_local_to_global(dit,lio)
       if ( .not. log_r(io) ) then
          ! we are not asked to remove these cross-terms
          ind = l_ptr(lio)
          list(indx+1:indx+l_ncol(lio)) = l_col(ind+1:ind+l_ncol(lio))
          indx = indx + l_ncol(lio)
          cycle
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          if ( log_r(UCORB(l_col(ind),no_u)) ) then
             
             is = (l_col(ind)-1)/no_u
             if ( isc_off(dir,is) /= 0 ) cycle

          end if

          indx = indx + 1

          list(indx) = l_col(ind)
          
       end do
       
    end do
    if ( n_nzs /= indx ) then
       call die('Sp_remove_crossterms: could not ensure sparsity pattern')
    end if

    deallocate(log_r)

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,n_nzs,num,listptr,list, &
         name='T '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))
    
    ! Clean up
    deallocate(num,listptr,list)

  end subroutine Sp_remove_crossterms

  ! Create a unified sparsity pattern
  subroutine Sp_union(dit,sp1,sp2,out)
    ! The distribution the two sparsity patterns lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be unified
    type(Sparsity), intent(inout) :: sp1, sp2
    ! sparsity pattern to be unified too
    type(Sparsity), intent(inout) :: out
    
    integer, pointer :: ncol1(:), ptr1(:), col1(:)
    integer, pointer :: ncol2(:), ptr2(:), col2(:)
    integer, allocatable :: ncol(:), ptr(:)
    integer :: no1, no2, no_u
    integer :: nc, ncg
    type(tRgn) :: r1, r2, rcol
    integer :: io

    if ( .not. initialized(sp1) ) then
       if ( .not. initialized(sp2) ) then
          call delete(out)
       else
          ! Point out to sp2
          out = sp2
       end if
       
       return
    else if ( .not. initialized(sp2) ) then
       if ( .not. initialized(sp1) ) then
          call delete(out)
       else
          ! Point out to sp1
          out = sp1
       end if

       return
    end if

    ! Attach the two sparsity patterns
    call attach(sp1,n_col=ncol1,list_ptr=ptr1,list_col=col1, &
         nrows=no1,nrows_g=no_u,ncols=nc,ncols_g=ncg)
    call attach(sp2,n_col=ncol2,list_ptr=ptr2,list_col=col2, &
         nrows=no2)
    
    if ( no1 /= no2 ) then
       call die('Sp_union: Different distributions for unifying not allowed.')
    end if

    ! Allocate count of orbitals
    allocate(ncol(no1),ptr(no1))

    ! Make room for *everything*
    ! Yes, this may be much too large, but it shouldn't matter too much
    call rgn_init(rcol, size(col1) + size(col2))
    rcol%n = 0

    ! Count the new sparsity pattern
    ! Loop over all the entries 
    do io = 1 , no1

       ! Create sorted list
       call rgn_list(r1,ncol1(io),col1(ptr1(io)+1:ptr1(io)+ncol1(io)))
       call rgn_sort(r1)
       call rgn_list(r2,ncol2(io),col2(ptr2(io)+1:ptr2(io)+ncol2(io)))
       call rgn_union(r1, r2, r2)
       ! Sort column indices, because it will most likely
       ! speed up everything
       call rgn_sort(r2)

       ncol(io) = r2%n
       ptr(io) = rcol%n

       ! Simultaneously we build col
       ! We do this after assigning ptr as we then 
       ! count that easily
       if ( .not. rgn_push(rcol, r2) ) &
            call die('Sp_Union: push-error')

    end do

    call rgn_delete(r1,r2)

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no1,no_u,rcol%n,ncol,ptr,rcol%r, &
         name=trim(name(sp1))//' U '//trim(name(sp2)), &
         ncols=nc,ncols_g=ncg)

    ! Clean up
    deallocate(ncol,ptr)
    
    call rgn_delete(rcol)

  end subroutine Sp_union

  ! Create a unified sparsity pattern
  subroutine Sp_sort(sp)
    ! sparsity pattern to be sorted
    type(Sparsity), intent(inout) :: sp
    
    integer, pointer :: ncol(:), ptr(:), col(:)
    integer :: io, no

    if ( .not. initialized(sp) ) then
       return
    end if

    ! Attach the sparsity pattern
    call attach(sp,n_col=ncol,list_ptr=ptr,list_col=col, &
         nrows=no)
    
    ! Count the new sparsity pattern
    ! Loop over all the entries 
    do io = 1 , no

       ! Create sorted list
       call sort_quick(ncol(io), col(ptr(io)+1:))

    end do

  end subroutine Sp_sort

  subroutine Sp_remove_region(dit,in,rr,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The region we wish to remove
    type(tRgn), intent(in) :: rr
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, n_nzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, jo, indx, n_rr
    logical, allocatable :: log_rr(:)

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)

    allocate(log_rr(no_u))
    log_rr(:) = .false.
    n_rr = rr%n
    call rgn_2logical(rr, log_rr)

    allocate(num(no_l))
!$OMP parallel do default(shared), private(lio,io,ind,jo)
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0

       if ( l_ncol(lio) /= 0 ) then

       io = index_local_to_global(dit,lio)

       if ( .not. log_rr(io) ) then

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          jo = ucorb(l_col(ind),no_u)
          if ( .not. log_rr(jo) ) then
       
            ! The orbital exists on the atom
            num(lio) = num(lio) + 1
          end if
       end do

       end if
       end if
       
    end do
!$OMP end parallel do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    n_nzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(n_nzs))
    indx = 0
    do lio = 1 , no_l

       io = index_local_to_global(dit,lio)
       if ( .not. log_rr(io) ) then

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

          jo = ucorb(l_col(ind),no_u)
          if ( .not. log_rr(jo) ) then

            indx = indx + 1
            list(indx) = l_col(ind)
          end if
          
       end do

       end if
       
    end do
    if ( n_nzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    deallocate(log_rr)
    
    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,n_nzs,num,listptr,list, &
         name='T '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))

    ! Clean up
    deallocate(num,listptr,list)
    
  end subroutine Sp_remove_region

  subroutine Sp_retain_region(dit,in,rr,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The region to retain
    type(tRgn), intent(in) :: rr
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    type(tRgn) :: full, rem_r
    integer :: no_u

    call attach(in,nrows_g=no_u)

    ! instead of creating two different
    ! algorithms we create the complement region
    ! and then remove that!
    ! HOW genius! :)

    call rgn_range(full,1,no_u)
    
    ! create complement of rem_r
    call rgn_complement(rr,full,rem_r)
    ! clean-up
    call rgn_delete(full)

    call Sp_remove_region(dit,in,rem_r,out)

    call rgn_delete(rem_r)
    
  end subroutine Sp_retain_region

  subroutine Sp_remove_region2region(dit,in,r1,r2,out)
    ! The distribution this sparsity pattern lives in
    type(OrbitalDistribution), intent(in) :: dit
    ! sparsity pattern to be reduced
    type(Sparsity), intent(inout) :: in
    ! The cross-terms "from" and "to", or "to" and "from"
    type(tRgn), intent(in) :: r1, r2
    ! The region of orbitals that will be removed
    type(Sparsity), intent(inout) :: out
    
    ! ** local variables
    ! tssp variables
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, n_nzs
    integer, allocatable :: num(:), listptr(:), list(:)
    integer :: lio, io, ind, jo, indx, ridx
    logical, allocatable :: log_r(:,:)

    ! set the connections to zero for the supplied atoms
    ! We do NOT reduce the sparsity row-size, ONLY the number of
    ! entries.

    call attach(in,nrows=no_l,nrows_g=no_u, &
        n_col=l_ncol,list_ptr=l_ptr,list_col=l_col)
    allocate(log_r(no_u, 2))
    log_r(:,:) = .false.
    do io = 1, r1%n
      log_r(r1%r(io), 1) = .true.
    end do
    do io = 1, r2%n
      log_r(r2%r(io), 2) = .true.
    end do

    allocate(num(no_l))
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ridx,ind,jo)
    do lio = 1 , no_l

       ! Initialize sparsity to 0 entries
       num(lio) = 0
       
       if ( l_ncol(lio) /= 0 ) then

         io = index_local_to_global(dit,lio)

         if ( log_r(io, 1) ) then
           ! io is in region-1, the checked region
           ! for jo will then be in region two
           ridx = 2
         else if ( log_r(io, 2) ) then
           ridx = 1
         else
           ridx = 0
           num(lio) = l_ncol(lio)
         end if

         if ( ridx /= 0 ) then

           do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             jo = ucorb(l_col(ind),no_u)
             ! If ridx is 2, then io is in region 1
             ! Else, if ridx is 1, then io is in region 2
             if ( log_r(jo, ridx) ) cycle

             ! The orbital exists on the atom
             num(lio) = num(lio) + 1

           end do

         end if
       end if

    end do
!$OMP end parallel do

    ! Create listptr
    allocate(listptr(no_l))
    listptr(1) = 0
    do lio = 2 , no_l
       listptr(lio) = listptr(lio-1) + num(lio-1)
    end do

    n_nzs = listptr(no_l) + num(no_l)

    ! Create actual sparsity pattern 
    allocate(list(n_nzs))
    indx = 0
    do lio = 1 , no_l

       io = index_local_to_global(dit,lio)
       if ( log_r(io, 1) ) then
          ridx = 2
       else if ( log_r(io, 2) ) then
          ridx = 1
       else
          ridx = 0
       end if

       do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

         if ( ridx /= 0 ) then
           ! See above for explanation
           jo = ucorb(l_col(ind),no_u)
           if ( log_r(jo, ridx) ) cycle
         end if

          indx = indx + 1
          list(indx) = l_col(ind)
          
       end do
       
    end do
    if ( n_nzs /= indx ) then
       call die('Could not ensure sparsity pattern')
    end if

    deallocate(log_r)

    ! Create new sparsity pattern and copy over
    call newSparsity(out,no_l,no_u,n_nzs,num,listptr,list, &
         name='T '//trim(name(in)), &
         ncols=ncols(in),ncols_g=ncols_g(in))

    ! Clean up
    deallocate(num,listptr,list)
    
  end subroutine Sp_remove_region2region

  subroutine Sp_to_Spglobal(dit,in,out)

#ifdef MPI
    use mpi_siesta
#endif

    type(OrbitalDistribution), intent(in) :: dit
    ! The sparsity that needs to be globalized
    type(Sparsity), intent(inout) :: in
    ! The globalized sparsity pattern
    type(Sparsity), intent(inout) :: out

    integer :: lio, io, no_l, no_u, n_nzs, n_nzsg, Bn, ind, lind
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, allocatable :: ncol(:), n_ptr(:), col(:)
    integer :: Node, comm
#ifdef MPI
    integer :: MPIerror
#endif

    if ( dist_nodes(dit) == 1 ) then
       out = in
       return
    end if

#ifndef MPI
    call die('Error in this, you are using non-MPI &
         &yet have a distribution... Weird')
#else

    Node = dist_node(dit)
    comm = dist_comm(dit)
    
    call attach(in ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_l,nrows_g=no_u,nnzs=n_nzs)

    ! Grab number of non-zero elements
    call MPI_AllReduce(n_nzs,n_nzsg, 1, MPI_Integer, MPI_SUM, &
         comm,MPIerror)

    allocate(ncol(no_u),n_ptr(no_u),col(n_nzsg))

    n_ptr(1) = 0
    ind = 0
    do io = 1 , no_u

       Bn = node_handling_element(dit,io)
       if ( Node == Bn ) then
          lio = index_global_to_local(dit,io)
          ncol(io) = l_ncol(lio)
       end if

       ! First grab number of non-zero elements
       call MPI_Bcast(ncol(io),1,MPI_Integer, &
            Bn, comm,MPIerror)

       ! Then b-cast the number of elements
       if ( Node == Bn .and. ncol(io) > 0 ) then
          lind = l_ptr(lio)
          col(ind+1:ind+ncol(io)) = l_col(lind+1:lind+l_ncol(lio))
       end if

       if ( ncol(io) > 0 ) then
          call MPI_Bcast(col(ind+1),ncol(io),MPI_Integer, &
               Bn, comm, MPIerror)
       end if

       ! Update index-pointer
       ind = ind + ncol(io)

       if ( io > 1 ) then
          n_ptr(io) = n_ptr(io-1) + ncol(io-1)
       end if

    end do

    if ( ind /= n_nzsg ) then
       call die('Error in sparsity pattern globalization')
    end if

    call newSparsity(out,no_u,no_u,n_nzsg,ncol,n_ptr,col, &
         name='G ('//name(in)//')')
    
#endif

  end subroutine Sp_to_Spglobal

  subroutine dSpData1D_to_Sp(A,out,B)

    ! the sparse data that needs to be transferred
    ! to the new sparsity format (out)
    type(dSpData1D), intent(inout) :: A,B
    ! The sparsity it needs to be transferred to
    type(Sparsity), intent(inout) :: out

    type(dData1D) :: dat
    real(dp), pointer :: inA(:), outA(:)

    ! get lists of sparsity patterns
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: n_ncol(:), n_ptr(:), n_col(:)
    integer :: lio, io, nr, lnr, lind, rind, nind

    ! get value and distribution
    inA => val(A)
    dit => dist(A)
    sp  => spar(A)

    ! Allocate the new data array
    io = nnzs(out)
    call newdData1D(dat,io,trim(name(A))//' reduced')
    outA => val(dat)

    call attach(sp ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(out,n_col=n_ncol,list_ptr=n_ptr,list_col=n_col, &
         nrows=lio,nrows_g=io)
    if ( lnr /= lio .or. nr /= io ) then
       call die('Error in reduced sparsity pattern')
    end if
    
!$OMP parallel do default(shared), &
!$OMP&private(lio,lind,rind,nind)
    do lio = 1 , lnr

       nind = 0
       
       if ( l_ncol(lio) /= 0 ) then
       if ( n_ncol(lio) /= 0 ) then
       
       do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
          
          find: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
             if ( l_col(lind) == n_col(rind) ) then
                outA(rind) = inA(lind)
                nind = nind + 1
                exit find
             end if
          end do find

       end do

       ! Check that copying goes as planned
       if ( nind /= n_ncol(lio) ) then
          call die('Error in sparsity copying')
       end if

       end if
       end if
       
    end do
!$OMP end parallel do

    ! copy over data
    call newdSpData1D(out,dat,dit,B,trim(name(A))//' reduced')

    ! Delete the data (we have copied it...)
    call delete(dat)
    
  end subroutine dSpData1D_to_Sp

  subroutine dSpData2D_to_Sp(A,out,B)

    ! the sparse data that needs to be transferred
    ! to the new sparsity format (out)
    type(dSpData2D), intent(inout) :: A, B
    ! The sparsity it needs to be transferred to
    type(Sparsity), intent(inout) :: out

    type(dData2D) :: dat
    real(dp), pointer :: inA(:,:), outA(:,:)

    ! get lists of sparsity patterns
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:) , l_ptr(:) , l_col(:)
    integer, pointer :: n_ncol(:), n_ptr(:), n_col(:)
    integer :: lio, io, nr, lnr, lind, rind, nind
    integer :: sp_dim
    integer :: d

    ! get value and distribution
    inA => val(A)
    dit => dist(A)
    sp  => spar(A)

    ! Create the data array for the reduced size
    sp_dim = spar_dim(A)
    io = nnzs(out)
    if ( sp_dim == 1 ) then
       d = size(inA,dim=2)
       call newdData2D(dat,io,d,trim(name(A))//' reduced')
    else
       d = size(inA,dim=1)
       call newdData2D(dat,d,io,trim(name(A))//' reduced')
    end if
    outA => val(dat)

    call attach(sp ,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=lnr,nrows_g=nr)
    call attach(out,n_col=n_ncol,list_ptr=n_ptr,list_col=n_col, &
         nrows=lio,nrows_g=io)
    if ( lnr /= lio .or. nr /= io ) then
       call die('Error in reduced sparsity pattern')
    end if
    
!$OMP parallel do default(shared), &
!$OMP&private(lio,lind,rind,nind)
    do lio = 1 , lnr

       nind = 0
       
       if ( l_ncol(lio) /= 0 ) then
       if ( n_ncol(lio) /= 0 ) then
       
       select case ( sp_dim )
       case ( 1 )
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             find1: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
                if ( l_col(lind) == n_col(rind) ) then
                   outA(rind,:) = inA(lind,:)
                   nind = nind + 1
                   exit find1
                end if
             end do find1
          end do
       case ( 2 )
          do lind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)
             
             find2: do rind = n_ptr(lio) + 1 , n_ptr(lio) + n_ncol(lio)
                if ( l_col(lind) == n_col(rind) ) then
                   outA(:,rind) = inA(:,lind)
                   nind = nind + 1
                   exit find2
                end if
             end do find2
          end do
       end select

       ! Check that copying goes as planned
       if ( nind /= n_ncol(lio) ) then
          call die('Error in sparsity copying')
       end if

       end if
       end if
       
    end do
!$OMP end parallel do

    ! copy over data
    call newdSpData2D(out,dat,dit,B,trim(name(A))//' reduced', &
        sparsity_dim=sp_dim)

    ! Delete the data (we have copied it...)
    call delete(dat)
    
  end subroutine dSpData2D_to_Sp

  subroutine dSpData2D_interp(N,A_2Ds,x,x0)
    use m_interpolate
    integer, intent(in) :: N
    type(dSpData2D), intent(inout) :: A_2Ds(N)
    ! The real values that we wish to interpolate from
    real(dp), intent(in) :: x(N), x0 ! and to

    ! Container to do array assignments
    type :: A2D
       real(dp), pointer :: v(:,:)
    end type A2D
    type(A2D) :: array(N)

    integer :: io, i, j, n_nzs, sp_dim, dim2
    real(dp) :: y(N)

    ! prep the array assignments...
    do i = 1 , N
       array(i)%v => val(A_2Ds(i))
    end do

    ! Get the dimensionality and the sparsity dimension
    sp_dim = spar_dim(A_2Ds(1))
    if ( sp_dim == 1 ) then
      dim2 = size(array(1)%v, dim=2)
    else
      dim2 = size(array(1)%v, dim=1)
    end if
    n_nzs = nnzs(A_2Ds(1))

!$OMP parallel default(shared), private(j,io,i,y)

    select case ( sp_dim )
    case ( 1 )
      do j = 1 , dim2
!$OMP do
        do io = 1 , n_nzs
          do i = 1 , N
            y(i) = array(i)%v(io,j)
          end do
          call interp_spline(N,x,y,x0, array(1)%v(io,j))
        end do
!$OMP end do nowait
      end do
    case ( 2 )
!$OMP do
      do io = 1 , n_nzs
        do j = 1 , dim2
          do i = 1 , N
            y(i) = array(i)%v(j,io)
          end do
          call interp_spline(N,x,y,x0, array(1)%v(j,io))
        end do
      end do
!$OMP end do nowait
    end select
!$OMP end parallel

  end subroutine dSpData2D_interp

end module m_sparsity_handling
