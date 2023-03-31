!
  module class_PEXSIDist
#ifdef MPI  
    use mpi
#endif
  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  type PEXSIDist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null PEXSIDist"
     !------------------------
     integer  :: group = -1      ! MPI group
     integer  :: node = -1       ! MPI rank in group  (my_proc)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
     !
  end type PEXSIDist_

  type PEXSIDist
     type(PEXSIDist_), pointer :: data => null()
  end type PEXSIDist

  public :: newDistribution
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newPEXSIDistribution
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
#define TYPE_NAME PEXSIDist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(PEXSIDist_) :: spdata
      ! do nothing
     end subroutine delete_Data

  subroutine newPEXSIDistribution(Blocksize,Group,this,name)
     !........................................
     ! Constructor
     !........................................
     type (PEXSIDist), intent(inout) :: this
     integer, intent(in)                       :: Blocksize
     integer, intent(in)                       :: Group
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)

     this%data%blocksize = Blocksize
     this%data%group      = Group

#ifdef MPI
     if (Group == MPI_GROUP_NULL) then
        this%data%node = MPI_PROC_NULL
     else
        call MPI_Group_Rank( Group, this%data%node, error )
        call MPI_Group_Size( Group, this%data%nodes, error )
     endif
#else
     this%data%node = 0
     this%data%nodes = 1
#endif
     this%data%node_io = 0

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(PEXSI Dist from BlockSize and Group)"
     endif
     call tag_new_object(this)

   end subroutine newPEXSIDistribution

!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     integer :: remainder

     remainder = nels - this%data%blocksize * this%data%nodes
     if (Node == (this%data%Nodes - 1)) then
        nl = this%data%blocksize + remainder
     else if (Node >= this%data%Nodes) then
        nl = 0
     else 
        nl = this%data%blocksize
     endif

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     if (Node >= this%data%Nodes) then
        ig = 0
     else
        ig = this%data%blocksize*Node + il
     endif

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: owner, myrank, error

     if (present(Node)) then
        if (Node >= this%data%nodes) then
           il = 0
        else
           il = ig - this%data%blocksize*Node 
        endif
     else
        ! Assume that we only want a non-zero value if the orb
        ! is local to this node
        owner = node_handling_element_(this,ig)
        if (owner == this%data%node) then
           il = ig - this%data%blocksize * this%data%node
        else
           il = 0
        endif
     endif
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(PEXSIDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     ! Assume bs=2, norbs=13, nodes=5
     ! Then, the distribution is 2, 2, 2, 2, 5 for procs 0,1,2,3,4
     ! For ig=13, proc=6 according to the naive calculation (first line)
     ! We have to correct this to assign this orb to proc 4.
     ! Same for ig=10: proc=5
     ! In this case the load balancing is quite bad.

     proc = (ig-1) / this%data%blocksize
     if (proc > this%data%Nodes-1) proc = this%data%Nodes - 1

   end function node_handling_element_


end module class_PEXSIDist
