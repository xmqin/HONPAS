!
  module class_BlockCyclicDist
#ifdef MPI
    use mpi
#endif
  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  type BlockCyclicDist_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null BlockCyclicDist"
     !------------------------
     integer  :: group = -1      ! MPI group
     integer  :: node = -1       ! MPI rank in group  (my_proc)
     integer  :: nodes = 0       ! MPI size of group  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   ! 
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
  end type BlockCyclicDist_

  type BlockCyclicDist
     type(BlockCyclicDist_), pointer :: data => null()
  end type BlockCyclicDist

  public :: newDistribution
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local

  interface newDistribution
     module procedure newBlockCyclicDistribution
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
#define TYPE_NAME BlockCyclicDist
#include "basic_type.inc"
!====================================    

     subroutine delete_Data(spdata)
      type(BlockCyclicDist_) :: spdata
      ! do nothing
     end subroutine delete_Data

  subroutine newBlockCyclicDistribution(Blocksize,Group,this,name)
     !........................................
     ! Constructor
     !........................................
     type (BlockCyclicDist), intent(inout) :: this
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
        this%data%name = "(Distribution from BlockSize and Group)"
     endif
     call tag_new_object(this)

   end subroutine newBlockCyclicDistribution

!-----------------------------------------------------------
   function num_local_elements_(this,nels,Node) result(nl)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in)                    :: Node
     integer                                :: nl

     if (Node >= this%data%Nodes) then
        nl = 0
     else
        nl = numroc(nels,this%data%blocksize,Node,  &
                 this%data%isrcproc,this%data%nodes)
     endif

   end function num_local_elements_

   function index_local_to_global_(this,il,Node) result(ig)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in)                    :: Node
     integer                                :: ig

     if (Node >= this%data%Nodes) then
        ig = 0
     else
        ig = indxl2g(il,this%data%blocksize,Node, &
                  this%data%isrcproc,this%data%nodes)
     endif

   end function index_local_to_global_

   function index_global_to_local_(this,ig,Node) result(il)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: owner, myrank, error

     if (present(Node)) then
        if (Node >= this%data%nodes) then
           il = 0
        else
           il = indxg2l(ig,this%data%blocksize,Node,  &
                     this%data%isrcproc,this%data%nodes)
        endif
     else
        ! Assume that we only want a non-zero value if the orb
        ! is local to this node
        owner = node_handling_element_(this,ig)
        if (owner == this%data%node) then
           il = indxg2l(ig,this%data%blocksize,myrank,  &
                        this%data%isrcproc,this%data%nodes)
        else
           il = 0
        endif
     endif
 
   end function index_global_to_local_

   function node_handling_element_(this,ig) result(proc)
     type(BlockCyclicDist), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     integer :: dummy_Node = 0

     proc = indxg2p(ig,this%data%blocksize,dummy_Node,  &
                    this%data%isrcproc,this%data%nodes)


   end function node_handling_element_

!========================================================================
       INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
       INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  NUMROC computes the NUMber of Rows Or Columns of a distributed
!  matrix owned by the process indicated by IPROC.
!
!  Arguments
!  =========
!
!  N         (global input) INTEGER
!            The number of rows/columns in distributed matrix.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local input) INTEGER
!            The coordinate of the process whose local array row or
!            column is to be determined.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row or column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Local Scalars ..
       INTEGER              EXTRABLKS, MYDIST, NBLOCKS
!     ..
!     .. Intrinsic Functions ..
       INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
!     Figure PROC's distance from source process
!
       MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
!
!     Figure the total number of whole NB blocks N is split up into
!
       NBLOCKS = N / NB
!
!     Figure the minimum number of rows/cols a process can have
!
       NUMROC = (NBLOCKS/NPROCS) * NB
!
!     See if there are any extra blocks
!
       EXTRABLKS = MOD( NBLOCKS, NPROCS )
!
!     If I have an extra block
!
       IF( MYDIST.LT.EXTRABLKS ) THEN
           NUMROC = NUMROC + NB
!
!         If I have last block, it may be a partial block
!
       ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
           NUMROC = NUMROC + MOD( N, NB )
       END IF
!
       RETURN
!
!     End of NUMROC
!
       END function numroc
!==========================================================================
      INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXLOC, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXL2G computes the global index of a distributed matrix entry
!  pointed to by the local index INDXLOC of the process indicated by
!  IPROC.
!
!  Arguments
!  =========
!
!  INDXLOC   (global input) INTEGER
!            The local index of the distributed matrix entry.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local input) INTEGER
!            The coordinate of the process whose local array row or
!            column is to be determined.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
      INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) +     &
                MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1
!
      RETURN
!
!     End of INDXL2G
!
      END function indxl2g

!===========================================================
      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXG2P computes the process coordinate which posseses the entry of a
!  distributed matrix specified by a global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the element.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )
!
      RETURN
!
!     End of INDXG2P
!
      END function indxg2p

!================================================================

      INTEGER FUNCTION INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXG2L computes the local index of a distributed matrix entry
!  pointed to by the global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the distributed matrix entry.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1
!
      RETURN
!
!     End of INDXG2L
! 
      END function indxg2l

end module class_BlockCyclicDist
