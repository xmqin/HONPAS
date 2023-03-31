!

  module class_Dist
#ifdef MPI
    use mpi
#endif

    use class_BlockCyclicDist
    use class_PEXSIDist

  implicit none

  character(len=*), parameter :: mod_name=__FILE__

!---------------------------------------------
  type, private :: distKind
     private
     integer :: i
  end type distKind

  type(distKind), parameter, public          :: TYPE_NULL = distKind(-1)
  type(distKind), parameter, public          :: TYPE_BLOCK_CYCLIC = distKind(1)
  type(distKind), parameter, public          :: TYPE_PEXSI = distKind(2)

  public operator(==)
  interface operator(==)
     module procedure isequal_
  end interface
!---------------------------------------------

  type Dist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null Dist"
     !------------------------
     integer  :: group = -1      ! MPI group
     integer  :: node = -1       ! MPI rank in group  (my_proc)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     type(distKind)           :: dist_type
     type(BlockCyclicDist) :: bdist
     type(PEXSIDist)       :: pdist
     !------------------------------------------
  end type Dist_

  type Dist
     type(Dist_), pointer :: data => null()
  end type Dist

  public :: newDistribution, distType, group
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newDistribution_
  end interface
  interface distType
     module procedure dist_type
  end interface
  interface group
     module procedure group_
  end interface
  interface num_local_elements
     module procedure num_local_elements_
  end interface
  interface index_local_to_global
     module procedure index_local_to_global_
  end interface
  interface index_global_to_local
     module procedure index_global_to_local_
  end interface
  interface node_handling_element
     module procedure node_handling_element_
  end interface

!====================================    
#define TYPE_NAME Dist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(Dist_) :: spdata
        if (spdata%dist_type == TYPE_BLOCK_CYCLIC) then
           call delete(spdata%bdist)
        else if (spdata%dist_type == TYPE_PEXSI) then
           call delete(spdata%pdist)
        else if (spdata%dist_type == TYPE_NULL) then
           ! do nothing
        else
           ! do nothing
        end if

     end subroutine delete_Data

  function dist_type(this) result(dtype)
    type(Dist)  :: this
    type(distKind) :: dtype
  
    dtype = this%data%dist_type

  end function dist_type

  function isequal_ (a,b) result(iseq)
    type(distKind), intent(in) :: a, b
    logical                 :: iseq
    
    iseq = (a%i == b%i)
  end function isequal_
!
  subroutine newDistribution_(Blocksize,Group,this,dist_type,name)
     !........................................
     ! Constructor
     !........................................
     type (Dist), intent(inout) :: this
     integer, intent(in)                       :: Blocksize
     integer, intent(in)                       :: Group
     type (distKind), intent(in)               :: dist_type
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)

     if (dist_type == TYPE_BLOCK_CYCLIC) then
        call newDistribution(Blocksize,Group,this%data%bdist,name)
     else if (dist_type == TYPE_PEXSI) then
        call newDistribution(Blocksize,Group,this%data%pdist,name)
     else if (dist_type == TYPE_NULL) then
        ! do nothing
     else
        ! do nothing
     end if

     this%data%dist_type = dist_type

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(Abstract Distribution from BS and Group)"
     endif
     call tag_new_object(this)

   end subroutine newDistribution_

!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        nl =  num_local_elements(this%data%bdist,nels,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        nl =  num_local_elements(this%data%pdist,nels,Node)
     else if (dist_type(this) == TYPE_NULL) then
        nl = 0
     else
        nl = 0
     end if

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        ig = index_local_to_global(this%data%bdist,il,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        ig = index_local_to_global(this%data%pdist,il,Node)
     else if (dist_type(this) == TYPE_NULL) then
        ig = 0
     else
        ig = 0
     end if

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        il = index_global_to_local(this%data%bdist,ig,Node)
     else if (dist_type(this) == TYPE_PEXSI) then
        il = index_global_to_local(this%data%pdist,ig,Node)
     else if (dist_type(this) == TYPE_NULL) then
        il = 0
     else
        il = 0
     end if
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(Dist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        proc = node_handling_element(this%data%bdist,ig)
     else if (dist_type(this) == TYPE_PEXSI) then
        proc = node_handling_element(this%data%pdist,ig)
     else if (dist_type(this) == TYPE_NULL) then
        proc = MPI_UNDEFINED
     else
        proc = MPI_UNDEFINED
     end if

   end function node_handling_element_

!-----------------------------------------------------------
   function group_(this) result(g)
     type(Dist), intent(in)  :: this
     integer                 :: g

     if (dist_type(this) == TYPE_BLOCK_CYCLIC) then
        g =  this%data%bdist%data%group
     else if (dist_type(this) == TYPE_PEXSI) then
        g =  this%data%pdist%data%group
     else if (dist_type(this) == TYPE_NULL) then
        g = MPI_GROUP_NULL
     else
        g = MPI_GROUP_NULL
     end if

   end function group_


end module class_Dist
