! Module for easy writing of all information
! in a single netCDF file.

! We rely on the NetCDF4 library and can
! utilise parallel IO or standard IO

module m_ncdf_io

  use precision, only : sp, dp, grid_p
  use parallel, only : Node, Nodes

  use class_OrbitalDistribution
  use class_Sparsity
  use m_io_s, only: max_consecutive_sum, Node_Sp_gncol
  use m_io_s, only: count_consecutive, count_blocks

#ifdef NCDF_4
  use netcdf_ncdf, ncdf_parallel => parallel
#endif
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  private

#ifdef NCDF_4
  public :: cdf_init_mesh
  public :: cdf_r_grid
  public :: cdf_w_grid

  public :: cdf_w_sp
  public :: cdf_w_d1D, cdf_w_d2D

  public :: cdf_r_sp
  public :: cdf_r_d1D, cdf_r_d2D

#endif

  ! Create the type to hold the mesh distribution data
  type :: tMeshDist
     private
     ! number of mesh divisions in each axis
     integer :: m(3)
     ! mesh box bounds of each node in each direction
            ! box(1,iAxis,iNode)=lower bounds
            ! box(2,iAxis,iNode)=upper bounds
     integer, allocatable :: box(:,:,:)
  end type tMeshDist
  type(tMeshDist), save :: distr

contains

#ifdef NCDF_4

  subroutine cdf_init_mesh(mesh,nsm)
    
    use parallel, ProcY => ProcessorY

    ! mesh divisions (fine points), number of fine points per big points
    integer, intent(in) :: mesh(3), nsm

    ! Local quantities
    integer :: nm(3), nmesh, ntot
    integer :: ProcZ, blocY, blocZ, nremY, nremZ
    integer :: dimX, dimY, dimZ
    integer :: PP, iniY, iniZ, PY, PZ

    if ( .not. allocated(distr%box) ) then
       allocate(distr%box(2,3,0:Nodes-1))
    end if
    
    nm(1:3) = mesh(1:3) / nsm
    nmesh = product(mesh)
    distr%m(1:3) = mesh(1:3)
    
    ProcZ = Nodes/ProcY
    
    blocY = nm(2)/ProcY
    nremY = nm(2) - blocY*ProcY
    blocZ = nm(3)/ProcZ
    nremZ = nm(3) - blocZ*ProcZ

    dimX = nm(1) * nsm
    
    ntot = 0
    
    PP   = 0
    iniY = 1
    do PY = 1, ProcY
       
       dimY = blocY
       if ( PY <= nremY ) dimY = dimY + 1  ! Add extra points starting from the first nodes
       dimY = dimY * nsm                 ! For fine points
       
        iniZ = 1
        do PZ = 1, ProcZ
           dimZ = blocZ
           if ( PZ <= nremZ ) dimZ = dimZ + 1
           dimZ = dimZ*nsm                 ! For fine points
           
           distr%box(1,1,PP) = 1
           distr%box(2,1,PP) = dimX
           distr%box(1,2,PP) = iniY
           distr%box(2,2,PP) = iniY + dimY - 1
           distr%box(1,3,PP) = iniZ
           distr%box(2,3,PP) = iniZ + dimZ - 1

           ntot = ntot + dimX * dimY * dimZ
           
           iniZ = iniZ + dimZ
           PP   = PP + 1
           
        end do

        iniY = iniY + dimY

     end do
     
     if (ntot /= nmesh) then
        if (Node == 0) then
           write(6,*) "Nominal npt: ", nmesh, " /= assigned npt:", ntot
        end if
        call die()
     end if
     
  end subroutine cdf_init_mesh

  subroutine cdf_w_Sp(ncdf,dit,sp)
  
    type(hNCDF), intent(inout) :: ncdf
    type(OrbitalDistribution), intent(inout) :: dit
    type(Sparsity), intent(inout) :: sp

    ! Local variables
    integer, pointer :: ncol(:), l_col(:)
    integer, allocatable :: ibuf(:), gncol(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n, j
#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif
    integer :: n_nnzs, g_nzs
    integer :: nb, ib

    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

#ifdef MPI
    if ( no_l /= no_u ) then
       call MPI_Reduce(n_nzs,g_nzs,1,MPI_Integer, &
            MPI_Sum, 0, MPI_Comm_World, MPIerror)
    else
       g_nzs = n_nzs
    end if
#else
    g_nzs = nnzs(sp)
#endif

    ! Read dimension nnzs
    call ncdf_inq_dim(ncdf,'nnzs',len=n_nnzs)
    if ( Node == 0 .and. n_nnzs /= g_nzs ) then
       call die('Number of non-zero elements is not equivalent.')
    end if
    
    if ( parallel_io(ncdf) ) then

       ! we write it using MPI
       call ncdf_par_access(ncdf,name='n_col',access=NF90_INDEPENDENT)

       allocate(ibuf(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,ibuf)
#ifdef MPI
       ! Globalize n_col
       call MPI_Bcast(ibuf,no_u,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif
       if ( Node == 0 ) then
          call ncdf_put_var(ncdf,'n_col',ibuf)
       end if

       call ncdf_par_access(ncdf,name='list_col',access=NF90_COLLECTIVE)

       ! Determine the number of blocks
       nb = count_blocks(dit, no_u)
       ind = 1
       io = 1
       do ib = 1 , nb
         if ( io <= no_l ) then
         
           ! We loop on each segment until no more
           ! segment exists
           ! Get the current inset point
           gio = index_local_to_global(dit,io)

           if ( gio > 1 ) then
             gind = sum(ibuf(1:gio-1)) + 1
           else
             gind = 1
           end if

           ! count number of orbitals in this block
           j = count_consecutive(dit,no_u,gio)
          
           ! Figure out how many this corresponds to in the 
           ! list_col array
           max_n = sum(ncol(io:io+j-1))
          
           call ncdf_put_var(ncdf,'list_col',l_col(ind:ind+max_n-1), &
               start=(/gind/), count=(/max_n/) )

           ! Update sparse index
           ind = ind + max_n

           ! Update row index
           io = io + j

         else
           ! Fake access
           call ncdf_put_var(ncdf,'list_col',l_col(1:1), &
               start=(/1/), count=(/0/) )
         end if

       end do

       deallocate(ibuf)

    else

#ifdef MPI
       allocate(gncol(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,gncol)
       call ncdf_put_var(ncdf,'n_col',gncol)

       if ( Node == 0 ) then
          max_n = maxval(gncol)
          allocate(ibuf(max_n))
       end if

       ! Write list_col

       ! Loop size
       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)

          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 .and. ncol(io) > 0 ) then
                call ncdf_put_var(ncdf,'list_col',l_col(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else if ( ncol(io) > 0 ) then
                call MPI_Send( l_col(ind+1) , ncol(io), MPI_Integer, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 .and. gncol(gio) > 0 ) then
             call MPI_Recv( ibuf(1) , max_n, MPI_Integer, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: cdf_w_Sp')
             call MPI_Get_Count(MPIstatus, MPI_Integer, io, MPIerror)
             call ncdf_put_var(ncdf,'list_col',ibuf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

       if ( Node == 0 ) then
          deallocate(ibuf)
       end if
       deallocate(gncol)

#else
       call ncdf_put_var(ncdf,'n_col',ncol)
       call ncdf_put_var(ncdf,'list_col',l_col)
#endif

    end if

  end subroutine cdf_w_Sp

  ! Reads in a sparsity pattern at the
  ! current position in the file (iu)
  ! The sparsity pattern "sp" will be returned
  ! as populated.
  ! If dist is supplied it will distribute
  ! the sparsity pattern as supplied (this implies Bcast = .true.)
  ! Else if Bcast is true it will b-cast the sparsity 
  ! pattern fully.
  subroutine cdf_r_Sp(ncdf, no, sp, tag, dit, Bcast, gncol)

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Number of orbitals readed
    integer, intent(in) :: no
    ! Sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The tag of the sparsity pattern
    character(len=*), intent(in) :: tag
    ! distribution if needed to be b-cast in a non-global
    ! fashion
    type(OrbitalDistribution), intent(inout), optional :: dit
    ! Bcast the values?
    logical, intent(in), optional :: Bcast
    integer, intent(inout), target, optional :: gncol(no)

    ! Local variables for reading the values
    integer, pointer :: ncol(:) => null()
    integer, pointer :: l_ptr(:) => null()
    integer, pointer :: l_col(:) => null()

    integer :: io, n_nzs, gind, ind, nl, n, i, nb, ib
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
       call ncdf_get_var(ncdf,'n_col',lncol)

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
       gind = 1
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

                i = sum(ncol(io:io-1+n))
                call ncdf_get_var(ncdf,'list_col', l_col(ind+1:ind+i), &
                     start = (/gind/) , count = (/i/) )
                ind  =  ind + i
                gind = gind + i

             else

                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( l_col(ind+1) , i, MPI_Integer, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
                
             end if

          else if ( Node == 0 ) then

             i = sum(lncol(gio:gio-1+n))
             call ncdf_get_var(ncdf,'list_col', ibuf(1:i), &
                  start = (/gind/) , count = (/i/) )

             call MPI_Send( ibuf(1) , i, MPI_Integer, &
                  BNode, ib, MPI_Comm_World, MPIerror)

             gind = gind + i
             
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

          ind = sum(ncol(1:no))
          call ncdf_get_var(ncdf,'list_col',l_col(1:ind), &
               count = (/ind/) )
          if ( ind /= n_nzs ) then
             call die('Error in reading sparsity pattern, &
                  &size not equivalent.')
          end if
       end if

       ! Bcast
       call MPI_Bcast(l_col(1),n_nzs,MPI_Integer, &
            0,MPI_Comm_World,MPIError)

    else
#endif

       ind = sum(ncol(1:no))
       call ncdf_get_var(ncdf,'list_col',l_col(1:ind), &
            count = (/ind/) )

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

  end subroutine cdf_r_Sp

  subroutine cdf_w_d1D(ncdf,vname,dSp1D)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData1D), intent(inout) :: dSp1D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n, j
    integer, allocatable :: gcol(:)
    real(dp), pointer :: a(:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif
    integer :: nb, ib

    dit => dist(dSp1D)
    sp => spar(dSp1D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp1D)

    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       allocate(gcol(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,gcol)
#ifdef MPI
       ! Globalize n_col
       call MPI_Bcast(gcol,no_u,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

       ! Determine the number of blocks
       nb = count_blocks(dit, no_u)
       ind = 1
       io = 1
       do ib = 1 , nb
         if ( io <= no_l ) then

           ! We loop on each segment until no more
           ! segment exists
           ! Get the current inset point
           gio = index_local_to_global(dit,io)

           if ( gio > 1 ) then
             gind = sum(gcol(1:gio-1)) + 1
           else
             gind = 1
           end if

           ! count number of orbitals in this block
           j = count_consecutive(dit,no_u,gio)
          
           ! Figure out how many this corresponds to in the 
           ! list_col array
           max_n = sum(ncol(io:io+j-1))
          
           call ncdf_put_var(ncdf,vname,a(ind:ind+max_n-1), &
               start=(/gind/), count=(/max_n/) )

           ! Update sparse index
           ind = ind + max_n

           ! Update row index
           io = io + j

         else
           ! Fake access
           call ncdf_put_var(ncdf,vname,a(1:1), &
               start=(/1/), count=(/0/) )
         end if
          
       end do

       deallocate(gcol)

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,vname,a(ind+1:ind+ncol(io)), &
                     count=(/ncol(io)/),start=(/gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(ind+1) , ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code: cdf_w_d1D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             call ncdf_put_var(ncdf,vname,buf(1:io), &
                  count=(/io/),start=(/gind/))
             gind = gind + io
          end if
       end do

       if ( Node == 0 ) then
         deallocate(buf)
       end if
#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d1D

  subroutine cdf_r_d1D(ncdf,vname,sp, dSp1D, tag, dit, Bcast, gncol)

    use class_dSpData1D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(Sparsity), intent(inout) :: sp
    type(dSpData1D), intent(inout) :: dSp1D
    ! The tag of the data format
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

    integer :: io, lno, no, gind, ind, n_nzs, n, i, nb, ib
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
       call die('Error in distribution, cdf_r_d1D')
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
       gind = 1
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

                ! Count number of elements
                i = sum(ncol(io:io-1+n))
                call ncdf_get_var(ncdf,vname, &
                     a(ind+1:ind+i),start=(/gind/),count=(/i/))
                ind = ind + i
                gind = gind + i

             else

                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( a(ind+1) , i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
                
             end if
             
          else if ( Node == 0 ) then
             
             i = sum(lncol(gio:gio-1+n))
             call ncdf_get_var(ncdf,vname, &
                  buf(1:i),start=(/gind/),count=(/i/))

             call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIerror)

             gind = gind + i
             
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
       call die('Error in distribution for, cdf_r_d1D')
#endif
    else if ( lIO ) then

       i = sum(ncol(1:no))
       call ncdf_get_var(ncdf,vname,a(:), count = (/i/) )

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

       call MPI_Bcast(a(1),n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)

    end if
#endif
    
  end subroutine cdf_r_d1D

  subroutine cdf_w_d2D(ncdf,vname,dSp2D)

    use class_dSpData2D

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: vname
    type(dSpData2D), intent(inout) :: dSp2D

    ! Local variables
    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    integer, pointer :: ncol(:), l_col(:)
    integer :: no_l, no_u, n_nzs, gio, io, ind, gind, max_n, j
    integer :: is, id2, dim2, sp_dim
    integer, allocatable :: gcol(:)
    real(dp), pointer :: a(:,:)
    real(dp), allocatable :: buf(:)

#ifdef MPI
    integer :: BNode, MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif
    integer :: nb, ib

    dit => dist(dSp2D)
    sp => spar(dSp2D)
    ! Write the sparsity to the file...
    call attach(sp,nrows=no_l, nrows_g=no_u, &
         n_col=ncol,list_col=l_col,nnzs=n_nzs)

    a => val(dSp2D)
    sp_dim = spar_dim(dSp2D)
    if ( sp_dim == 1 ) then
      dim2 = size(a, dim=2)
    else
      dim2 = size(a, dim=1)
    end if

    if ( parallel_io(ncdf) ) then

       call ncdf_par_access(ncdf,name=vname,access=NF90_COLLECTIVE)

       allocate(gcol(no_u))
       call Node_Sp_gncol(0,sp,dit,no_u,gcol)
#ifdef MPI
       ! Globalize n_col
       call MPI_Bcast(gcol,no_u,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

       ! Determine the number of blocks
       nb = count_blocks(dit, no_u)
       ind = 1
       io = 1
       do ib = 1 , nb
         if ( io <= no_l ) then

           ! We loop on each segment until no more
           ! segment exists
           ! Get the current inset point
           gio = index_local_to_global(dit,io)

           if ( gio > 1 ) then
             gind = sum(gcol(1:gio-1)) + 1
           else
             gind = 1
           end if

           ! count number of orbitals in this block
           j = count_consecutive(dit,no_u,gio)

           ! Figure out how many this corresponds to in the 
           ! list_col array
           max_n = sum(ncol(io:io+j-1))
          
           if ( sp_dim == 1 ) then
             do is = 1 , dim2
               call ncdf_put_var(ncdf,vname,a(ind:ind+max_n-1,is), &
                   start=(/gind,is/), count=(/max_n/))
             end do
           else
             call ncdf_put_var(ncdf,vname,a(:,ind:ind+max_n-1), &
                 start=(/1,gind/), count=(/dim2,max_n/) )
           end if

           ! Update sparse index
           ind = ind + max_n

           ! Update row index
           io = io + j

         else
           ! Fake access
           if ( sp_dim == 1 ) then
             do is = 1 , dim2
               call ncdf_put_var(ncdf,vname,a(1:1,is), &
                   start=(/1,is/), count=(/0/))
             end do
           else
             call ncdf_put_var(ncdf,vname,a(:,1:1), &
                 start=(/1,1/), count=(/dim2,0/) )
           end if
         end if

       end do

       deallocate(gcol)

    else

#ifdef MPI
       io = maxval(ncol)
       call MPI_Reduce(io,max_n,1,MPI_Integer, &
            MPI_Max, 0, MPI_Comm_World, MPIerror)
       max_n = max_n * dim2
       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

    if ( sp_dim == 1 ) then

       do id2 = 1 , dim2
          ind = 0
          gind = 1
          do gio = 1 , no_u
             BNode = node_handling_element(dit,gio)
             
             if ( Node == BNode ) then
                io = index_global_to_local(dit,gio,Node)
                if ( Node == 0 ) then
                   call ncdf_put_var(ncdf,trim(vname),a(ind+1:ind+ncol(io),id2), &
                        count=(/ncol(io)/),start=(/gind,id2/))
                   gind = gind + ncol(io)
                else
                   call MPI_Send( a(ind+1,id2) , ncol(io), MPI_Double_Precision, &
                        0, gio, MPI_Comm_World, MPIerror)
                   
                end if
                ind = ind + ncol(io)
             else if ( Node == 0 ) then
                call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                     BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
                if ( MPIerror /= MPI_Success ) &
                     call die('Error in code (1): cdf_w_d2D')
                call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
                call ncdf_put_var(ncdf,trim(vname),buf(1:io), &
                     count=(/io/),start=(/gind,id2/))
                gind = gind + io
             end if
          end do ! gio
       end do ! id2

    else

       ind = 0
       gind = 1
       do gio = 1 , no_u
          BNode = node_handling_element(dit,gio)
          
          if ( Node == BNode ) then
             io = index_global_to_local(dit,gio,Node)
             if ( Node == 0 ) then
                call ncdf_put_var(ncdf,trim(vname),a(1:dim2,ind+1:ind+ncol(io)), &
                     count=(/dim2,ncol(io)/),start=(/1,gind/))
                gind = gind + ncol(io)
             else
                call MPI_Send( a(1,ind+1) , dim2*ncol(io), MPI_Double_Precision, &
                     0, gio, MPI_Comm_World, MPIerror)
             end if
             ind = ind + ncol(io)
          else if ( Node == 0 ) then
             call MPI_Recv( buf(1) , max_n, MPI_Double_Precision, &
                  BNode, gio, MPI_Comm_World, MPIstatus, MPIerror )
             if ( MPIerror /= MPI_Success ) &
                  call die('Error in code (2): cdf_w_d2D')
             call MPI_Get_Count(MPIstatus, MPI_Double_Precision, io, MPIerror)
             io = io / dim2
             call ncdf_put_var(ncdf,trim(vname),reshape(buf(1:io*dim2),(/dim2,io/)), &
                  count=(/dim2,io/),start=(/1,gind/))
             gind = gind + io
          end if
       end do

    end if

       if ( Node == 0 ) then
         deallocate(buf)
       end if
#else
       call ncdf_put_var(ncdf,trim(vname),a)
#endif

    end if
    
  end subroutine cdf_w_d2D

  subroutine cdf_r_d2D(ncdf,vname, sp, dSp2D, dim2, tag, &
       sparsity_dim, dit, Bcast, gncol)

    use class_dSpData2D

    ! File handle
    type(hNCDF), intent(inout) :: ncdf
    ! Variable name
    character(len=*), intent(in) :: vname
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

    integer :: io, lno, no, s, gind, ind, n_nzs, n, i, nb, ib
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
       call die('Error in distribution, cdf_r_d2D')
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
       gind = 1
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
                i = sum(ncol(io:io-1+n))
                call ncdf_get_var(ncdf,vname,a(1:dim2,ind+1:ind+i), &
                     start = (/1,gind/) , count = (/dim2,i/) )
                ind = ind + i
                gind = gind + i
             else
                ! count the number of received entities
                i = sum(ncol(io:io-1+n))
                call MPI_IRecv( a(1,ind+1) ,dim2*i, MPI_Double_Precision, &
                     0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                ind = ind + i
             end if
          else if ( Node == 0 ) then
             i = sum(lncol(gio:gio-1+n))
             call ncdf_get_var(ncdf,vname,buf(1:i), &
                  start = (/1,gind/) , count = (/dim2,i/) )
             call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                  BNode, ib, MPI_Comm_World, MPIerror)
             gind = gind + i
          end if
          gio = gio + n
       end do

    else ! non-collapsed IO

       if ( Node == 0 ) then
          allocate(buf(max_n))
       end if

       do s = 1 , dim2
          gind = 1
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
                   i = sum(ncol(io:io-1+n))
                   call ncdf_get_var(ncdf,vname,a(ind+1:ind+i,s), &
                        start = (/gind,s/) , count = (/i/) )
                   ind = ind + i
                   gind = gind + i
                else
                   ! count the number of received entities
                   i = sum(ncol(io:io-1+n))
                   call MPI_IRecv( a(ind+1,s) ,i, MPI_Double_Precision, &
                        0, ib, MPI_Comm_World, ibuf(ib), MPIerror )
                   ind = ind + i
                end if
             else if ( Node == 0 ) then
                i = sum(lncol(gio:gio-1+n))
                call ncdf_get_var(ncdf,vname,buf(1:i), &
                     start = (/gind,s/) , count = (/i/) )
                call MPI_Send( buf(1) , i, MPI_Double_Precision, &
                     BNode, ib, MPI_Comm_World, MPIerror)
                gind = gind + i
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
       call die('Error in distribution for, cdf_r_d2D')
#endif
    else if ( lIO ) then

       i = sum(ncol(1:no))
       if ( sp_dim == 2 ) then ! collapsed IO
          call ncdf_get_var(ncdf,vname,a(:,1:i), count = (/dim2,i/) )
       else ! non-collapsed IO
          call ncdf_get_var(ncdf,vname,a(1:i,:), &
               start = (/1,1/) , count = (/i,dim2/) )
       end if

    end if

    ! If a distribution is present, then do something
#ifdef MPI
    if ( lBcast ) then

       call MPI_Bcast(a(1,1),dim2*n_nzs,MPI_Double_Precision, &
            0, MPI_Comm_World, MPIError)

    end if
#endif

  end subroutine cdf_r_d2D
  
  subroutine cdf_w_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(in) :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), allocatable :: gb(:)
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

       ! Ensure collective writing
       call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

       lb(:)  = distr%box(1,:,Node)
       nel(:) = distr%box(2,:,Node) - lb(:) + 1
       if ( .not. all(nel == nmeshl) ) then
          call die('cdf_r_grid: cannot assert the grid-size from the &
               &stored grid and the denoted size')
       end if

       if ( present(idx) ) then
          call ncdf_put_var(ncdf,name,grid, &
               start=(/lb(1),lb(2),lb(3),idx/), &
               count=(/nel(1),nel(2),nel(3),1/) )          
       else
          call ncdf_put_var(ncdf,name,grid, start=lb, count=nel )
       end if
   
    else

       mnpt = 0
       do iN = 0 , Nodes - 1
          nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
          mnpt = max(mnpt,product(nel))
       end do
    
       ! The main node can safely write the data...
       if ( Node == 0 ) then
          
          allocate(gb(mnpt))
          
          ! First save it's own data
          lb(:) = distr%box(1,:,0) 
          nel(:) = distr%box(2,:,0) - lb(:) + 1
          
          if ( present(idx) ) then
             call ncdf_put_var(ncdf,name,grid, &
                  start=(/lb(1),lb(2),lb(3),idx/), &
                  count=(/nel(1),nel(2),nel(3),1/) )
          else
             call ncdf_put_var(ncdf,name,grid, &
                  start=lb, count=nel )
          end if

          ! Loop on the remaining nodes
          do iN = 1 , Nodes - 1

             ! we retrieve data from the iN'th node.
             call MPI_Recv(gb,mnpt,MPI_grid_real,iN,iN, &
                  MPI_Comm_World,MPIstat, MPIerror)
             
             ! Just make sure we only pass the correct size
             call MPI_Get_Count(MPIstat, MPI_Grid_Real, inpt, MPIerror)
             
             lb(:) = distr%box(1,:,iN) 
             nel(:) = distr%box(2,:,iN) - lb(:) + 1
             if ( inpt /= product(nel) ) then
                call die('Error when receiving the distributed &
                     &grid data for writing to the NetCDF file.')
             end if
             
             if ( present(idx) ) then
                call ncdf_put_var(ncdf,name,gb(1:inpt), &
                     start=(/lb(1),lb(2),lb(3),idx/), &
                     count=(/nel(1),nel(2),nel(3),1/) )
             else
                call ncdf_put_var(ncdf,name,gb(1:inpt), &
                     start=lb, count=nel )
             end if
          
          end do
          
          deallocate(gb)
       else
          mnpt = product(nmeshl)
          call MPI_Send(grid,mnpt,MPI_grid_real,0, &
               Node,MPI_Comm_World,MPIerror)
       end if
       
    end if

#else
    if ( present(idx) ) then
       call ncdf_put_var(ncdf,name,grid,start=(/1,1,1,idx/), &
            count=nmeshl)
    else
       call ncdf_put_var(ncdf,name,grid, &
            count=nmeshl)
    end if
#endif
  
  end subroutine cdf_w_grid

  subroutine cdf_r_grid(ncdf,name,nmeshl,grid,idx)

    type(hNCDF), intent(inout) :: ncdf
    character(len=*), intent(in) :: name
    integer, intent(in) :: nmeshl(3)
    real(grid_p), intent(inout), target :: grid(:)
    integer, intent(in), optional :: idx

#ifdef MPI
    integer :: MPIstat(MPI_STATUS_SIZE)
    integer :: MPIerror, mnpt
    integer :: lb(3), nel(3), iN, inpt
    real(grid_p), pointer :: gb(:) => null()
#endif

#ifdef MPI

    if ( parallel_io(ncdf) ) then

       ! Ensure collective writing
       call ncdf_par_access(ncdf,name=name,access=NF90_COLLECTIVE)

       lb(:)  = distr%box(1,:,Node)
       nel(:) = distr%box(2,:,Node) - lb(:) + 1
       if ( .not. all(nel == nmeshl) ) then
          call die('cdf_r_grid: cannot assert the grid-size from the &
               &stored grid and the denoted size')
       end if

       if ( present(idx) ) then
          call ncdf_get_var(ncdf,name,grid, &
               start=(/lb(1),lb(2),lb(3),idx/), &
               count=(/nel(1),nel(2),nel(3),1/) )          
       else
          call ncdf_get_var(ncdf,name,grid, start=lb, count=nel )
       end if
   
    else

       mnpt = 0
       do iN = 0 , Nodes - 1
          nel(:) = distr%box(2,:,iN) - distr%box(1,:,iN) + 1
          mnpt = max(mnpt,product(nel))
       end do
       
       ! The main node can safely read the data...
       if ( Node == 0 ) then
          
          if ( mnpt > product(nmeshl) ) then
             ! allocate, the IO node have a smaller
             ! space than the others
             allocate(gb(mnpt))
          else
             gb => grid
          end if
                    
          ! Loop on the remaining nodes (we may possibly
          ! reuse the current grid array, and hence
          ! we first read the other nodes...)
          do iN = 1 , Nodes - 1
     
             lb(:) = distr%box(1,:,iN) 
             nel(:) = distr%box(2,:,iN) - lb(:) + 1

             ! get total number of points
             inpt = product(nel)
             
             if ( present(idx) ) then
                call ncdf_get_var(ncdf,name,gb(1:inpt), &
                     start=(/lb(1),lb(2),lb(3),idx/), &
                     count=(/nel(1),nel(2),nel(3),1/) )
             else
                call ncdf_get_var(ncdf,name,gb(1:inpt), &
                     start=lb, count=nel )
             end if

             ! we send data to the iN'th node.
             call MPI_Send(gb,inpt,MPI_grid_real,iN,iN, &
                  MPI_Comm_World,MPIerror)
             
          end do

          if ( mnpt > product(nmeshl) ) then
             deallocate(gb)
          end if

          ! Retrieve it's own data
          lb(:) = distr%box(1,:,0) 
          nel(:) = distr%box(2,:,0) - lb(:) + 1

          if ( present(idx) ) then
             call ncdf_get_var(ncdf,name,grid, &
                  start=(/lb(1),lb(2),lb(3),idx/), &
                  count=(/nel(1),nel(2),nel(3),1/) )
          else
             call ncdf_get_var(ncdf,name,grid, &
                  start=lb, count=nel )
          end if
          
       else
          mnpt = product(nmeshl)
          call MPI_Recv(grid,mnpt,MPI_grid_real,0, &
               Node,MPI_Comm_World,MPIstat,MPIerror)
       end if
       
    end if

#else
    if ( present(idx) ) then
       call ncdf_get_var(ncdf,name,grid,start=(/1,1,1,idx/), &
            count=nmeshl)
    else
       call ncdf_get_var(ncdf,name,grid,start=(/1,1,1/), &
            count=nmeshl)
    end if
#endif
  
  end subroutine cdf_r_grid

#else
  subroutine dummy()
    !pass
  end subroutine dummy
#endif

end module m_ncdf_io
