! Module for easy reading of sparsity patterns and 
! data distributed in sparsity patterns.

! Ideally all IO routines should utilize these routines

! Fully implemented by Nick Papior Andersen

! SIESTA io-module
module m_io_s

  use class_OrbitalDistribution
  use class_Sparsity
  use precision, only : sp, dp
  use parallel, only : Node
#ifdef MPI
  use mpi_siesta, only : MPI_Bcast, MPI_AllReduce, MPI_Sum, MPI_Max
  use mpi_siesta, only : MPI_Comm_World, MPI_Comm_Self
  use mpi_siesta, only : MPI_Send, MPI_Recv
  use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
  use mpi_siesta, only : MPI_Success, MPI_Status_Size
  use mpi_siesta, only : MPI_REQUEST_NULL
#endif

  implicit none 

  private

  public :: Node_Sp_gncol
  public :: io_read_Sp
  public :: io_read_d1D, io_read_d2D
  public :: io_write_Sp
  public :: io_write_d1D, io_write_d2D

  ! The counting functions
  public :: count_blocks
  public :: count_consecutive, count_consecutive_sum
  public :: max_consecutive, max_consecutive_sum

contains

  ! Returns a consecutive number of contributions
  ! starting from the specified index
  function count_blocks(dit,n_t) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t
    integer :: n
    ! Local variables
    integer :: cur_node, i

    n = 1
    cur_node = node_handling_element(dit,1)
    do i = 2 , n_t
       ! if the idx is not present, just return
       if ( cur_node /= node_handling_element(dit,i) ) then
          n = n + 1
          cur_node = node_handling_element(dit,i)
       end if
    end do

  end function count_blocks

  ! Returns a consecutive number of contributions
  ! starting from the specified index
  function count_consecutive(dit,n_t,idx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, idx
    integer :: n
    ! Local variables
    integer :: cur_node, i

    n = 1
    cur_node = node_handling_element(dit,idx)
    do i = idx + 1 , n_t
       ! if the idx is not present, just return
       if ( cur_node /= node_handling_element(dit,i) ) return
       n = n + 1
    end do

  end function count_consecutive

  function count_consecutive_sum(dit,n_t,nidx,idx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, nidx(n_t), idx
    integer :: n

    n = count_consecutive(dit,n_t,idx)
    n = sum(nidx(idx:idx-1+n))

  end function count_consecutive_sum

  function max_consecutive(dit,n_t) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t
    integer :: n
    ! Local variables
    integer :: i, cur_n

    i = 1
    n = 0
    do while ( i <= n_t )
       
       ! Count number of consecutive numbers
       cur_n = count_consecutive(dit,n_t,i)
       if ( cur_n > n ) then
          n = cur_n
       end if

       ! Step counter
       i = i + cur_n

    end do
    
  end function max_consecutive

  function max_consecutive_sum(dit,n_t,nidx) result(n)
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: n_t, nidx(n_t)
    integer :: n
    ! Local variables
    integer :: i, cur_n, tmp

    i = 1
    n = 0
    do while ( i <= n_t )
       
       ! Count number of consecutive numbers
       cur_n = count_consecutive(dit,n_t,i)
       tmp = sum(nidx(i:i-1+cur_n))
       if ( tmp > n ) then
          n = tmp
       end if

       ! Step counter
       i = i + cur_n

    end do
    
  end function max_consecutive_sum

  subroutine Node_Sp_gncol(lNode,sp,dit,no,gcol)
    integer, intent(in) :: lNode
    type(Sparsity), intent(inout) :: sp
    type(OrbitalDistribution), intent(in) :: dit
    integer, intent(in) :: no
    integer, intent(out) :: gcol(no)

    integer, pointer :: ncol(:)
    integer :: gio, n, io, nb, ib
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif
    
    ! grab ncol
    call attach(sp,n_col=ncol)
    gcol(1) = 1 ! assure it is not performed again

#ifndef MPI
    gcol = ncol
#else

    nb = count_blocks(dit,no)
    allocate(ibuf(nb))
    ibuf(:) = MPI_REQUEST_NULL
          
    gio = 1
    ib = 0
    do while ( gio <= no )
       ib = ib + 1
          
       BNode = node_handling_element(dit,gio)

       ! Get number of consecutive orbitals
       n = count_consecutive(dit,no,gio)
       
       if ( Node == BNode ) then
          io = index_global_to_local(dit,gio,Node)
          if ( Node == lNode ) then
             gcol(gio:gio-1+n) = ncol(io:io-1+n)
          else
             call MPI_ISSend( ncol(io) , n, MPI_Integer, &
                  lNode, gio, MPI_Comm_World, ibuf(ib), MPIerror)
          end if
       else if ( Node == lNode ) then
          call MPI_IRecv( gcol(gio) , n, MPI_Integer, &
               BNode, gio, MPI_Comm_World, ibuf(ib), MPIerror )
       end if
       gio = gio + n
    end do

    do ib = 1 , nb
       if( ibuf(ib) /= MPI_REQUEST_NULL ) &
            call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
    end do
    deallocate(ibuf)
#endif
    
  end subroutine Node_Sp_gncol
    

  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine io_read_Sp(iu, no, sp, tag, dit, Bcast, gncol)

    ! File handle
    integer, intent(in) :: iu
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(no)

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, ind, nl, n, i, nb, ib
    logical :: ldit, lBcast, lIO
    integer, pointer :: lncol(:) => null()
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: MPIerror, MPIstatus(MPI_STATUS_SIZE), BNode
#endif

#ifdef MPI
    ldit = present(dit)
#else
    ldit = .false.
#endif
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    ! If one of them is provided, it will only be Node == 0
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    if ( lIO ) then

       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
       end if

       ! First read in number of non-zero 
       ! entries per orbital
       read(iu) lncol

    end if

    nl = no

    ! If a distribution is present, then do something
#ifdef MPI
    if ( ldit ) then

       ! First count number of local entries
       nl = 0
       do gio = 1 , no
          BNode = node_handling_element(dit,gio)
          if ( BNode == Node ) nl = nl + 1
       end do

       ! allocate room for the number of columns in
       ! each row
       allocate(ncol(nl))

       nb = count_blocks(dit,no)

       ! allocate all requests
       allocate(ibuf(nb))
       ibuf(:) = MPI_REQUEST_NULL

       ! Distribute it
       gio = 1
       ib = 0
       do while ( gio <= no ) 
          ib = ib + 1

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          ! belong to the same node...
          n = count_consecutive(dit,no,gio)

          if ( BNode == Node ) then

             io = index_global_to_local(dit,gio,Node)

             if ( Node == 0 ) then
                ncol(io:io-1+n) = lncol(gio:gio-1+n)
             else
                call MPI_IRecv( ncol(io) , n, MPI_Integer, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
             end if

          else if ( Node == 0 ) then

             call MPI_ISSend( lncol(gio) , n, MPI_Integer, &
                  BNode, ib, MPI_Comm_World, ibuf(ib), MPIerror)

          end if

          gio = gio + n

       end do

       do ib = 1 , nb
          if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
               call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
       end do
       deallocate(ibuf)

    else if ( lBcast ) then

       ! Everything should be b-casted
       if ( Node == 0 ) then
          ncol => lncol
       else
          allocate(ncol(nl))
       end if
       
       ! Bcast everything
       call MPI_Bcast(ncol(1),nl,MPI_Integer, &
            0,MPI_Comm_World,MPIError)
       
    else if ( lIO ) then
       ncol => lncol

    end if
#else
    ! Point to the buffer
    ncol => lncol
#endif

    ! Allocate pointer
    allocate(l_ptr(nl))
    
    l_ptr(1) = 0
    do io = 2 , nl
       l_ptr(io) = l_ptr(io-1) + ncol(io-1)
    end do
    
    ! Number of local non-zero elements
    ! (also works for any bcast methods)
    n_nzs = l_ptr(nl) + ncol(nl)

    ! Allocate space
    allocate(l_col(n_nzs))

#ifdef MPI
    if ( ldit ) then
       
       ! We have a distributed read
       if ( Node == 0 ) then
          max_n = max_consecutive_sum(dit,no,lncol)
          allocate(ibuf(max_n))
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

       ! Read in columns
       ind = 0
       gio = 1
       ib = 0
       do while ( gio <= no ) 
          ib = ib + 1

          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          ! belong to the same node...
          n = count_consecutive(dit,no,gio)

          if ( BNode == Node ) then

             ! Get the local orbital
             io = index_global_to_local(dit,gio,Node)

             if ( Node == 0 ) then

                do i = io , io - 1 + n
                   read(iu) l_col(ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do

             else

                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( l_col(ind+1) , i, MPI_Integer, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
                
             end if

          else if ( Node == 0 ) then

             i = 0
             do io = gio , gio - 1 + n
                read(iu) ibuf(i+1:i+lncol(io))
                i = i + lncol(io)
             end do

             call MPI_Send( ibuf(1) , i, MPI_Integer, &
                  BNode, ib, MPI_Comm_World, MPIerror)
             
          end if

          gio = gio + n

       end do

       if ( Node == 0 ) then
          if ( .not. present(gncol) ) deallocate(lncol)
       else
          do ib = 1 , nb
             if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
       end if
       deallocate(ibuf)
       
    else if ( lBcast ) then

       if ( Node == 0 ) then
          
          ind = 0
          do gio = 1 , no
             read(iu) l_col(ind+1:ind+ncol(gio))
             ind = ind + ncol(gio)
          end do

       end if

       ! Bcast
       call MPI_Bcast(l_col(1),n_nzs,MPI_Integer, &
            0,MPI_Comm_World,MPIError)

    else if ( lIO ) then
#endif

       ind = 0
       do io = 1 , no
          read(iu) l_col(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do

#ifdef MPI       
    end if
#endif

    ! Create the sparsity pattern
    call newSparsity(sp,nl,no, &
         n_nzs, ncol, l_ptr, l_col, trim(tag))

    ! de-allocate
    deallocate(l_ptr,l_col)
    if ( ldit ) deallocate(ncol)
    if ( lBcast .and. Node /= 0 ) deallocate(ncol)
    if ( lIO .and. .not. present(gncol) ) deallocate(lncol)

  end subroutine io_read_Sp

  ! Writes a sparsity pattern at the
  ! current position in the file (iu)
  ! If dist is supplied it will write a distributed sparsity pattern
  subroutine io_write_Sp(iu, sp, dit, gncol)

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! distribution
    type(OrbitalDistribution), intent(in), optional :: dit
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    integer, pointer :: lncol(:) => null()
    integer, pointer :: ncol(:), l_col(:) => null()

    integer :: lno, no, io, max_n, ind, n, i, nb, ib
    logical :: ldit
#ifdef MPI
    integer, allocatable :: ibuf(:)
    integer :: gio
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    ! Get the sparsity sizes
    call attach(sp,n_col=ncol, list_col=l_col, nrows=lno,nrows_g=no)

    ldit = present(dit)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then

#ifdef MPI
       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
          lncol(1) = -1
       end if
       if ( lncol(1) < 0 ) then
          call Node_Sp_gncol(0,sp,dit,no,lncol)
       end if

#else
       call die('Error in code, non-full contained sp')
#endif

    else
       lncol => ncol
    end if
    
    if ( Node == 0 ) then
       
       write(iu) lncol

    end if

#ifdef MPI
    ! Write the list_col array
    if ( ldit ) then

       nb = count_blocks(dit,no)

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
          ! Retrive the maximum number of non-zero
          ! elements in each row
          max_n = max_consecutive_sum(dit,no,lncol)
          allocate(ibuf(max_n))
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

       ! Loop size
       ind = 0
       gio = 1 
       ib = 0
       do while ( gio <= no )
          ib = ib + 1
          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
#ifdef TEST_IO
                i = sum(ncol(io:io-1+n))
                write(iu) l_col(ind+1:ind+i)
                ind = ind + i
#else
                do i = io , io - 1 + n
                   write(iu) l_col(ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do
#endif
             else
                i = sum(ncol(io:io-1+n))
                call MPI_ISSend( l_col(ind+1) , i, MPI_Integer, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror)
                ind = ind + i
             end if
          else if ( Node == 0 ) then
             call MPI_Recv( ibuf(1) , max_n, MPI_Integer, &
                  BNode, ib, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: io_write_Sp')
#ifdef TEST_IO
             i = sum(lncol(gio:gio-1+n))
             write(iu) ibuf(1:i)
#else
             i = 0
             do io = gio , gio - 1 + n
                write(iu) ibuf(i+1:i+lncol(io))
                i = i + lncol(io)
             end do
#endif

          end if
          gio = gio + n
       end do
       
       if ( .not. present(gncol) ) deallocate(lncol)
       if ( Node /= 0 ) then
          do ib = 1 , nb
             if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
       end if
       deallocate(ibuf)
    else

       ind = 0
       do io = 1 , no
          write(iu) l_col(ind+1:ind+lncol(io))
          ind = ind + lncol(io)
       end do
       
    end if

#else
    
    ind = 0
    do io = 1 , no
       write(iu) l_col(ind+1:ind+lncol(io))
       ind = ind + lncol(io)
    end do

#endif

  end subroutine io_write_Sp

  subroutine io_read_d1D(iu, sp, dSp1D, tag, dit, Bcast, gncol)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, n_nzs, n, i, nb, ib
    logical :: ldit, lBcast, lIO
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: max_n, gio
    integer :: MPIerror, BNode, MPIstatus(MPI_STATUS_SIZE)
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    ! If one of them is provided, it will only be Node == 0
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
       call newdSpData1D(sp,dit,dSp1D,name=trim(tag))
       
#ifdef MPI
       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
          lncol(1) = -1
       end if
       if ( lncol(1) < 0 ) then
          call Node_Sp_gncol(0,sp,dit,no,lncol)
       end if
#else
       call die('Error in distribution, io_read_d1D')
#endif

    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData1D(sp,fdit,dSp1D,name=trim(tag))
       ! Clean up the distribution again
       call delete(fdit)
    end if

    a => val(dSp1D)

    if ( ldit ) then
#ifdef MPI

       nb = count_blocks(dit,no)

       ! Allocate the maximum number of entries
       if ( Node == 0 ) then
          max_n = max_consecutive_sum(dit,no,lncol)
          allocate(buf(max_n))
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

       ! Loop size
       ind = 0
       gio = 1 
       ib = 0
       do while ( gio <= no )
          ib = ib + 1
          BNode = node_handling_element(dit,gio)

          ! Get number of consecutive orbitals
          ! belong to the same node...
          n = count_consecutive(dit,no,gio)

          if ( BNode == Node ) then

             ! Get the local orbital
             io = index_global_to_local(dit,gio,Node)

             if ( Node == 0 ) then

                do i = io , io - 1 + n
                   read(iu) a(ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do

             else

                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( a(ind+1) , i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
                
             end if
             
          else if ( Node == 0 ) then

             i = 0
             do io = gio , gio - 1 + n
                read(iu) buf(i+1:i+lncol(io))
                i = i + lncol(io)
             end do

             call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIerror)
             
          end if

          gio = gio + n

       end do

       if ( .not. present(gncol) ) deallocate(lncol)
       if ( Node == 0 ) then
          deallocate(buf)
       else
          do ib = 1 , nb
             if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
          deallocate(ibuf)
       end if
#else
       call die('Error in distribution for, io_read_d1D')
#endif
    else if ( lIO ) then

       ind = 0
       do io = 1 , no
          read(iu) a(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do

    end if

#ifdef MPI
    if ( lBcast ) then

       call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)

    end if
#endif

  end subroutine io_read_d1D

  subroutine io_write_d1D(iu, dSp1D, gncol)

    use class_dSpData1D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData1D), intent(inout) :: dSp1D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:) => null()
    integer, pointer :: lncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, n, i, nb, ib
    logical :: ldit
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    call attach(sp,nrows=lno,nrows_g=no,n_col=ncol)

    ! If they are different we should 
    ! use the distribution setting
    ldit = lno /= no

    ! Retrieve data
    a => val(dSp1D)

    if ( ldit ) then

#ifdef MPI
       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
          lncol(1) = -1
       end if
       if ( lncol(1) < 0 ) then
          call Node_Sp_gncol(0,sp,dit,no,lncol)
       end if
#else
       call die('Error in distribution, io_write_d1D')
#endif

#ifdef MPI

       nb = count_blocks(dit,no)

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
          max_n = max_consecutive_sum(dit,no,lncol)
          allocate(buf(max_n))
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

       ! Write the data...
       ind = 0
       gio = 1
       ib = 0 
       do while ( gio <= no )
          ib = ib + 1
          BNode = node_handling_element(dit,gio)
          
          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
#ifdef TEST_IO
                i = sum(ncol(io:io-1+n))
                write(iu) a(ind+1:ind+i)
                ind = ind + i
#else
                do i = io , io - 1 + n
                   write(iu) a(ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do
#endif
             else
                i = sum(ncol(io:io-1+n))
                call MPI_ISSend( a(ind+1) , i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror)
                ind = ind + i
             end if
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: io_write_d1D')
#ifdef TEST_IO
             i = sum(lncol(gio:gio-1+n))
             write(iu) buf(1:i)
#else
             i = 0
             do io = gio , gio - 1 + n
                write(iu) buf(i+1:i+lncol(io))
                i = i + lncol(io)
             end do
#endif
             
          end if
          gio = gio + n
       end do
          
       if ( .not. present(gncol) ) deallocate(lncol)
       if ( Node == 0 ) then
          deallocate(buf)
       else 
          ! Wait for the last one to not send
          ! two messages with the same tag...
          do ib = 1 , nb
             if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
          deallocate(ibuf)
       end if
#endif
    else

       ind = 0
       do io = 1 , no
          write(iu) a(ind+1:ind+ncol(io))
          ind = ind + ncol(io)
       end do
       
    end if
    
  end subroutine io_write_d1D

  subroutine io_read_d2D(iu, sp, dSp2D, dim2, tag, &
       sparsity_dim, dit, Bcast, gncol)
    
    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! Data array with sparsity pattern
    type(dSpData2D), intent(inout) :: dSp2D
    ! The non-sparse dimension
    integer, intent(in) :: dim2
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! This denotes the sparsity dimension (either 1 or 2) 1=default
    integer, intent(in), optional :: sparsity_dim
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(in), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(OrbitalDistribution) :: fdit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lncol(:) => null(), ncol(:)

    integer :: io, lno, no, s, ind, n_nzs, n, i, nb, ib
    integer :: sp_dim
    logical :: ldit, lBcast, lIO
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: MPIerror, BNode, MPIstatus(MPI_STATUS_SIZE)
#endif

    ldit = present(dit)
    lBcast = .false.
    if ( present(Bcast) ) lBcast = Bcast
    ! If one of them is provided, it will only be Node == 0
    lIO = .not. (lBcast .or. ldit)
    if ( .not. lIO ) lIO = (Node == 0)

    sp_dim = 1
    if ( present(sparsity_dim) ) sp_dim = sparsity_dim

    call attach(sp,nrows=lno, nrows_g=no, n_col=ncol,nnzs=n_nzs)
    if ( ldit ) ldit = lno /= no

    if ( ldit ) then
       call newdSpData2D(sp,dim2,dit,dSp2D,name=trim(tag), &
            sparsity_dim=sp_dim)

#ifdef MPI
       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
          lncol(1) = -1
       end if
       if ( lncol(1) < 0 ) then
          call Node_Sp_gncol(0,sp,dit,no,lncol)
       end if
#else
       call die('Error in distribution, io_read_d2D')
#endif

    else
       ! Create the Fake distribution
#ifdef MPI
       call newDistribution(no,MPI_Comm_Self,fdit,name='Fake dist')
#else
       call newDistribution(no,-1           ,fdit,name='Fake dist')
#endif
       call newdSpData2D(sp,dim2,fdit,dSp2D,name=trim(tag), &
            sparsity_dim=sp_dim)
       call delete(fdit)
    end if

    a => val(dSp2D)

    if ( ldit ) then

#ifdef MPI

       nb = count_blocks(dit,no)

       ! Allocate maximum number of entries
       if ( Node == 0 ) then
          max_n = max_consecutive_sum(dit,no,lncol)
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

    if ( sp_dim == 2 ) then ! collapsed IO

       if ( Node == 0 ) then
          allocate(buf(max_n*dim2))
       end if

       ! Read sparse blocks and distribute
       ind = 0
       gio = 1 
       ib = 0
       do while ( gio <= no )
          ib = ib + 1
          BNode = node_handling_element(dit,gio)
          ! Get number of consecutive orbitals
          ! belong to the same node...
          n = count_consecutive(dit,no,gio)
          if ( BNode == Node ) then
             ! Get the local orbital
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                do i = io , io - 1 + n
                   read(iu) a(1:dim2,ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do
             else
                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( a(1,ind+1) ,dim2*i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
             end if
          else if ( Node == 0 ) then
             i = 0
             do io = gio , gio - 1 + n
                read(iu) buf(i+1:i+dim2*lncol(io))
                i = i + dim2*lncol(io)
             end do
             call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIerror)
          end if
          gio = gio + n
       end do

    else ! non-collapsed IO

       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       do s = 1 , dim2
          ind = 0
          gio = 1 
          ib = 0
          do while ( gio <= no )
             ib = ib + 1
             BNode = node_handling_element(dit,gio)
             ! Get number of consecutive orbitals
             ! belong to the same node...
             n = count_consecutive(dit,no,gio)
             if ( BNode == Node ) then
                ! Get the local orbital
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   do i = io , io - 1 + n
                      read(iu) a(ind+1:ind+ncol(i),s)
                      ind = ind + ncol(i)
                   end do
                else
                   ! count the number of received entities
                   i = sum(ncol(io:io-1+n))
                   call MPI_IRecv( a(ind+1,s) ,i, MPI_Double_Precision, &
                        0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                   ind = ind + i
                end if
             else if ( Node == 0 ) then
                i = 0
                do io = gio , gio - 1 + n
                   read(iu) buf(i+1:i+lncol(io))
                   i = i + lncol(io)
                end do
                call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                     BNode, ib, MPI_Comm_World, MPIerror)
             end if
             gio = gio + n
          end do

          if ( Node /= 0 .and. s < dim2 ) then
             do ib = 1 , nb
                if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                     call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
                ibuf(ib) = MPI_REQUEST_NULL
             end do
          end if

       end do

    end if
    
       if ( Node == 0 ) then
          deallocate(buf)
       else
          do ib = 1 , nb
             if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
          deallocate(ibuf)
       end if

       if ( .not. present(gncol) ) deallocate(lncol)

#else
       call die('Error in distribution for, io_read_d2D')
#endif
    else if ( lIO ) then

       if ( sp_dim == 2 ) then ! collapsed IO
          ind = 0
          do io = 1 , no
             read(iu) a(1:dim2,ind+1:ind+ncol(io))
             ind = ind + ncol(io)
          end do
       else ! non-collapsed IO
          do s = 1 , dim2 
             ind = 0
             do io = 1 , no
                read(iu) a(ind+1:ind+ncol(io),s)
                ind = ind + ncol(io)
             end do
          end do
       end if

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

       call MPI_Bcast(a(1,1),dim2*n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)

    end if
#endif

  end subroutine io_read_d2D

  subroutine io_write_d2D(iu, dSp2D, gncol)

    use class_dSpData2D

    ! File handle
    integer, intent(in) :: iu
    ! Data array with sparsity pattern (and an attached distribution)
    type(dSpData2D), intent(inout) :: dSp2D
    integer, intent(inout), target, optional :: gncol(:)

    ! Local variables for reading the values
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: a(:,:) => null()
    integer, pointer :: lncol(:) => null(), ncol(:)

    integer :: io, lno, no, ind, nb, ib
    integer :: s, dim2, sp_dim, n_nzs, n, i
    logical :: ldit
#ifdef MPI
    real(dp), allocatable :: buf(:)
    integer, allocatable :: ibuf(:)
    integer :: gio, max_n
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    call attach(sp,nrows=lno,nrows_g=no, n_col=ncol,nnzs=n_nzs)
    
    ldit = lno /= no
    
    ! Retrieve data
    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

    if ( ldit ) then

#ifdef MPI
       if ( present(gncol) ) then
          lncol => gncol
       else
          allocate(lncol(no))
          lncol(1) = -1
       end if
       if ( lncol(1) < 0 ) then
          call Node_Sp_gncol(0,sp,dit,no,lncol)
       end if

#else
       call die('Error in distribution, io_write_d2D')
#endif
       
#ifdef MPI

       nb = count_blocks(dit,no)

       if ( Node == 0 ) then
          max_n = max_consecutive_sum(dit,no,lncol)
       else
          allocate(ibuf(nb))
          ibuf(:) = MPI_REQUEST_NULL
       end if

       ! Loop size
    if ( sp_dim == 1 ) then

       ! The ionode now has the maximum retrieved array
       if ( Node == 0 ) then
         allocate(buf(max_n))
       end if

       do s = 1 , dim2

          ! Write the data...
          ind = 0
          gio = 1
          ib = 0
          do while ( gio <= no )
             ib = ib + 1
             BNode = node_handling_element(dit,gio)
          
             ! Get number of consecutive orbitals
             n = count_consecutive(dit,no,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
#ifdef TEST_IO
                   i = sum(ncol(io:io-1+n))
                   write(iu) a(ind+1:ind+i,s)
                   ind = ind + i
#else
                   do i = io , io - 1 + n
                      write(iu) a(ind+1:ind+ncol(i),s)
                      ind = ind + ncol(i)
                   end do
#endif
                else
                   i = sum(ncol(io:io-1+n))
                   call MPI_ISSend( a(ind+1,s) , i, MPI_Double_Precision, &
                        0, ib, MPI_Comm_World, ibuf(ib), MPIerror)
                   ind = ind + i
                end if
             else if ( Node == 0 ) then
                call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                     BNode, ib, MPI_Comm_World, MPIstatus, MPIerror )
                if ( MPIerror /= MPI_Success ) &
                     call die('Error in code (1): io_write_d2D')
#ifdef TEST_IO
                i = sum(lncol(gio:gio-1+n))
                write(iu) buf(1:i)
#else
                i = 0
                do io = gio , gio - 1 + n
                   write(iu) buf(i+1:i+lncol(io))
                   i = i + lncol(io)
                end do
#endif
             
             end if
             gio = gio + n
          end do
          if ( Node /= 0 .and. s < dim2 ) then
             do ib = 1 , nb
                if ( ibuf(ib) /= MPI_REQUEST_NULL ) &
                     call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
                ibuf(ib) = MPI_REQUEST_NULL
             end do
          end if
       end do

    else

       if ( Node == 0 ) then
          allocate(buf(max_n*dim2))
       end if
       
       ind = 0
       gio = 1
       ib = 0
       do while ( gio <= no )
          ib = ib + 1
          BNode = node_handling_element(dit,gio)
          
          ! Get number of consecutive orbitals
          n = count_consecutive(dit,no,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
#ifdef TEST_IO
                i = sum(ncol(io:io-1+n))
                write(iu) a(1:dim2,ind+1:ind+i)
                ind = ind + i
#else
                do i = io , io - 1 + n
                   write(iu) a(1:dim2,ind+1:ind+ncol(i))
                   ind = ind + ncol(i)
                end do
#endif
             else
                i = sum(ncol(io:io-1+n))
                call MPI_ISSend( a(1,ind+1) , dim2*i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror)
                ind = ind + i
             end if
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code (2): io_write_d2D')
#ifdef TEST_IO
             i = sum(lncol(gio:gio-1+n))
             write(iu) buf(1:i)
#else
             i = 0
             do io = gio , gio - 1 + n
                write(iu) buf(i+1:i+dim2*lncol(io))
                i = i + lncol(io)
             end do
#endif
             
          end if
          gio = gio + n
       end do

    end if
    
       if ( .not. present(gncol) ) deallocate(lncol)
       if ( Node == 0 ) then
          deallocate(buf)
       else
          do ib = 1 , nb
             if( ibuf(ib) /= MPI_REQUEST_NULL ) &
                  call MPI_Wait(ibuf(ib),MPIstatus,MPIerror)
          end do
          deallocate(ibuf)
       end if
#else
       call die('Error in io_write_d2D')
#endif
    else
       
       if ( sp_dim == 1 ) then
          do s = 1 , dim2
             ind = 0
             do io = 1 , no
                write(iu) a(ind+1:ind+ncol(io),s)
                ind = ind + ncol(io)
             end do
          end do
       else
          ind = 0
          do io = 1 , no
             write(iu) a(1:dim2,ind+1:ind+ncol(io))
             ind = ind + ncol(io)
          end do
       end if
       
    end if
    
  end subroutine io_write_d2D
  
end module m_io_s
  
