!
module class_Distribution
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  character(len=*), parameter :: mod_name="class_Distribution.F90"

#ifndef MPI
  ! To allow compilation in serial mode
  integer, parameter, private :: MPI_GROUP_NULL = -huge(1)
  integer, parameter, private :: MPI_COMM_NULL = -huge(1)
  integer, parameter, private :: MPI_UNDEFINED = -huge(1)
#endif

  !---------------------------------------------
  integer, parameter, public          :: TYPE_NULL = -1
  integer, parameter, public          :: TYPE_BLOCK_CYCLIC = 1
  integer, parameter, public          :: TYPE_PEXSI = 2
  !---------------------------------------------

  type Distribution_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null Dist"
     !------------------------
     integer       :: dist_type
     !------------------------------------------
     integer  :: ref_comm  = MPI_COMM_NULL
     integer  :: group     = MPI_GROUP_NULL      ! MPI group
     integer  :: node      = MPI_UNDEFINED       ! MPI rank in group  (my_proc)
     integer, allocatable :: ranks_in_ref_comm(:)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   ! 
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
  end type Distribution_

  type Distribution
     type(Distribution_), pointer :: data => null()
  end type Distribution

  type(Distribution_), pointer :: obj => null()

  public :: newDistribution, distType
  public :: ref_comm, get_ranks_in_ref_comm
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newDistribution_
  end interface newDistribution

  interface distType
     module procedure dist_type
  end interface distType

  interface ref_comm
     module procedure ref_comm_
  end interface ref_comm

  interface get_ranks_in_ref_comm
     module procedure get_ranks_in_ref_comm_
  end interface get_ranks_in_ref_comm

  interface num_local_elements
     module procedure num_local_elements_
  end interface num_local_elements
  interface index_local_to_global
     module procedure index_local_to_global_
  end interface index_local_to_global
  interface index_global_to_local
     module procedure index_global_to_local_
  end interface index_global_to_local
  interface node_handling_element
     module procedure node_handling_element_
  end interface node_handling_element

  !====================================    
#define TYPE_NAME Distribution
#include "basic_type.inc"
  !====================================    

  subroutine delete_Data(spdata)
    type(Distribution_) :: spdata
    integer :: ierr

#ifdef MPI    
    if (spdata%group /= MPI_GROUP_NULL) then
       call MPI_Group_free(spdata%group,ierr)
    endif
#endif    
    if (allocated(spdata%ranks_in_ref_comm)) then
       deallocate(spdata%ranks_in_ref_comm)
    endif

  end subroutine delete_Data

  function dist_type(this) result(dtype)
    type(Distribution)  :: this
    integer :: dtype

    obj => this%data

    dtype = obj%dist_type

  end function dist_type

  !
  subroutine newDistribution_(this,Ref_Comm,Ranks_in_Ref_Comm, &
       dist_type,Blocksize,name)
    !........................................
    ! Constructor
    !........................................
    type (Distribution), intent(inout) :: this
    integer, intent(in)                       :: ref_comm
    integer, intent(in)                       :: ranks_in_ref_comm(:)
    integer, intent(in)                       :: dist_type
    integer, intent(in)                       :: Blocksize
    character(len=*), intent(in), optional    :: name

    integer :: error, gsize, ref_group

    call init(this)

    obj => this%data

    obj%dist_type = dist_type
    obj%ref_comm = ref_comm
    obj%blocksize = Blocksize

    ! Do we need to allocate, or a simple assignment would suffice? (F2003 only)
    gsize = size(ranks_in_ref_comm)
    if (allocated(obj%ranks_in_ref_comm)) deallocate(obj%ranks_in_ref_comm)
    allocate(obj%ranks_in_ref_comm(gsize))
    obj%ranks_in_ref_comm(:)  = ranks_in_ref_comm(:)

    ! Define a group internally only for the purposes of
    ! getting the ranks consistently

#ifdef MPI
    call MPI_Comm_Group(ref_comm,ref_group,error)
    call MPI_Group_Incl(ref_group,gsize,obj%ranks_in_ref_comm,obj%group,error)
    call MPI_Group_Rank( obj%Group, obj%node, error )
    call MPI_Group_Size( obj%Group, obj%nodes, error )
    call MPI_Group_free(ref_group, error)
#else
    obj%node = 0
    obj%nodes = 1
#endif
!    print "(i0,a,i0,1x,4i2)", obj%dist_type,  &
!                " node, ranks in ref: ", obj%node, obj%ranks_in_ref_comm(:)
    if (present(name)) then
       obj%name = trim(name)
    else
       obj%name = "(Distribution from BlockSize and Ranks)"
    endif
    call tag_new_object(this)

  end subroutine newDistribution_

  !-----------------------------------------------------------
  subroutine get_ranks_in_ref_comm_(this,ranks)
    type(Distribution), intent(in)            :: this
    integer, allocatable, intent(out) :: ranks(:)

    integer :: rsize

    obj => this%data
    rsize = size(obj%ranks_in_ref_comm)
    allocate(ranks(rsize))
    ranks(:) = obj%ranks_in_ref_comm(:)

  end subroutine get_ranks_in_ref_comm_
  !-----------------------------------------------------------
  function ref_comm_(this) result(comm)
    type(Distribution), intent(in)  :: this
    integer                 :: comm

    obj => this%data
    comm =  obj%ref_comm

  end function ref_comm_
  !-----------------------------------------------------------
  function num_local_elements_(this,nels,Node) result(nl)
    type(Distribution), intent(in)  :: this
    integer, intent(in)                    :: nels
    integer, intent(in)                    :: Node
    integer                                :: nl

#ifdef MPI
    integer, external :: numroc
#endif
    integer :: remainder

    obj => this%data

    select case(obj%dist_type)
    case (TYPE_BLOCK_CYCLIC)
       if (Node >= obj%Nodes) then
          nl = 0
       else
#ifdef MPI          
          nl = numroc(nels,obj%blocksize,Node,  &
               obj%isrcproc,obj%nodes)
#else
          nl = nels
#endif          
       endif
    case (TYPE_PEXSI) 
       remainder = nels - obj%blocksize * obj%nodes
       if (Node == (obj%Nodes - 1)) then
          nl = obj%blocksize + remainder
       else if (Node >= obj%Nodes) then
          nl = 0
       else 
          nl = obj%blocksize
       endif

    case default
       nl = -huge(1)   ! Other signal?
    end select

  end function num_local_elements_

  function index_local_to_global_(this,il,Node) result(ig)
    type(Distribution), intent(in)  :: this
    integer, intent(in)                    :: il
    integer, intent(in)                    :: Node
    integer                                :: ig
    
#ifdef MPI
    integer, external :: indxl2g
#endif

    obj => this%data

    select case(obj%dist_type)
    case (TYPE_BLOCK_CYCLIC)
       if (Node >= obj%Nodes) then
          ig = 0
       else
#ifdef MPI          
          ig = indxl2g(il,obj%blocksize,Node, &
               obj%isrcproc,obj%nodes)
#else
          ig = il
#endif          
       endif
    case (TYPE_PEXSI) 
       if (Node >= obj%Nodes) then
          ig = 0
       else
          ig = obj%blocksize*Node + il
       endif
    case default
       ig = 0   ! Other signal?
    end select

  end function index_local_to_global_

  function index_global_to_local_(this,ig,Node) result(il)
    type(Distribution), intent(in)  :: this
    integer, intent(in)                    :: ig
    integer, intent(in), optional          :: Node
    integer                                :: il

    integer :: owner, myrank, error
#ifdef MPI
    integer, external :: indxg2l
#endif

    obj => this%data

    select case(obj%dist_type)
    case (TYPE_BLOCK_CYCLIC)
       if (present(Node)) then
          if (Node >= obj%nodes) then
             il = 0
          else
#ifdef MPI
             il = indxg2l(ig,obj%blocksize,Node,  &
                  obj%isrcproc,obj%nodes)
#else
             il = ig
#endif             
          endif
       else
          ! Assume that we only want a non-zero value if the orb
          ! is local to this node
          owner = node_handling_element_(this,ig)
          if (owner == obj%node) then
#ifdef MPI
             il = indxg2l(ig,obj%blocksize,myrank,  &
                  obj%isrcproc,obj%nodes)
#else
             il = ig
#endif             
          else
             il = 0
          endif
       endif
    case (TYPE_PEXSI) 
       if (present(Node)) then
          if (Node >= obj%nodes) then
             il = 0
          else
             il = ig - obj%blocksize*Node 
          endif
       else
          ! Assume that we only want a non-zero value if the orb
          ! is local to this node
          owner = node_handling_element_(this,ig)
          if (owner == obj%node) then
             il = ig - obj%blocksize * obj%node
          else
             il = 0
          endif
       endif
    case default
       il = 0   ! Other signal?
    end select

  end function index_global_to_local_

  function node_handling_element_(this,ig) result(proc)
    type(Distribution), intent(in)  :: this
    integer, intent(in)                    :: ig
    integer                                :: proc

    integer, parameter :: dummy_Node = 0
#ifdef MPI
    integer, external  :: indxg2p
#endif
    
    obj => this%data

    select case(obj%dist_type)
    case (TYPE_BLOCK_CYCLIC)
#ifdef MPI
       proc = indxg2p(ig,obj%blocksize,dummy_Node,  &
            obj%isrcproc,obj%nodes)
#else
       proc = 0
#endif       

    case (TYPE_PEXSI) 
       ! Assume bs=2, norbs=13, nodes=5
       ! Then, the distribution is 2, 2, 2, 2, 5 for procs 0,1,2,3,4
       ! For ig=13, proc=6 according to the naive calculation (first line)
       ! We have to correct this to assign this orb to proc 4.
       ! Same for ig=10: proc=5
       ! In this case the load balancing is quite bad.

       proc = (ig-1) / obj%blocksize
       if (proc > obj%Nodes-1) proc = obj%Nodes - 1

    case default
       proc = MPI_UNDEFINED   ! Other signal?
    end select

  end function node_handling_element_

end module class_Distribution
