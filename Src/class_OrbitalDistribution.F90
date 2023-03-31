! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This should be used to specify the orbital distribution, in the spirit
! of the Scalapack descriptors
! Closer to Siesta, it should mimic the "order-N arrays" of J. Gale
!
!      subroutine GetNodeOrbs( NOrb, Node, Nodes, NOrbNode)
!            calculates the number of orbitals stored on the (local)? Node
!            Could be in principle that of any node.
!            NOTE: NOrb is not always no_u:   nocc in savepsi... (ksv)
!            NOrbNode = nOrbPerNode(Node+1)
!      subroutine GlobalToLocalOrb( GOrb, Node, Nodes, LOrb)
!            ditto, returning zero if Gorb is not handled in local node
!            LOrb = nG2L(GOrb)
!      subroutine LocalToGlobalOrb( LOrb, Node, Nodes, GOrb)
!            easy in principle... but Julian has the option to have the
!            info from any node...
!            GOrb = nL2G(LOrb,Node+1)
!           It is used in iodm_netcdf, but it could be emulated by having
!           each node send also the info about the global counterparts.
!      subroutine WhichNodeOrb( GOrb, Nodes, Node)
!            returns the Node number that handles GOrb
!            Node = nNode(GOrb)

  module class_OrbitalDistribution
  
  implicit none

  character(len=*), parameter :: mod_name="class_OrbitalDistribution"

  type OrbitalDistribution_
     integer              :: refCount = 0
     character(len=36)    :: id
     !------------------------
     character(len=256)   :: name = "null OrbitalDistribution"
     !------------------------
     integer  :: comm = -1       ! MPI communicator
     integer  :: node = -1       ! MPI rank in comm  (my_proc)
     integer  :: nodes = 0       ! MPI size of comm  (nprocs)
     integer  :: node_io = -1    ! Node capable of IO
     !------------------------
     integer  :: blocksize = 0   ! Flag to determine whether the dist is block-cyclic
     !                   
     integer  :: isrcproc = 0    ! Processor holding the first element (unused)
     !------------------------------------------
     !
     !                   If the distribution is not block-cyclic,
     !                   we need explicit information
     !
     !
     ! Just for sanity: number of global elements (see function num_local_elements)
     ! This will have to be set at the time of filling in the arrays
     ! for non-block-cyclic distributions
     integer                         :: n = -1
     ! Number of elements held by each proc
     integer, pointer                :: nroc_proc(:)  => null() 
     ! Global index from Local index (in my_proc)   
     integer, pointer                :: nl2g(:)       => null() 
     ! Local index (in my_proc) from global index  (zero if not held by my_proc) 
     integer, pointer                :: ng2l(:)       => null()  
     ! Processor number from global index (base 0)
     integer, pointer                :: ng2p(:)       => null()    
     !------------------------------------------
  end type OrbitalDistribution_

  type OrbitalDistribution
     type(OrbitalDistribution_), pointer :: data => null()
  end type OrbitalDistribution

  public :: newDistribution
  public :: dist_comm, dist_node, dist_nodes
  public :: num_local_elements, node_handling_element
  public :: index_local_to_global, index_global_to_local
  public :: print_type

  interface newDistribution
     module procedure newBlockCyclicDistribution
  end interface

  interface print_type
     module procedure printOrbitalDistribution
  end interface

  interface dist_comm
     module procedure comm_
  end interface
  interface dist_node
     module procedure dist_node_
  end interface
  interface dist_nodes
     module procedure dist_nodes_
  end interface

!==================================== 
#define TYPE_NAME OrbitalDistribution
#include "basic_type.inc"
!==================================== 

  subroutine delete_Data(spdata)
    type(OrbitalDistribution_) :: spdata
    if ( associated(spdata%nroc_proc) ) &
         deallocate( spdata%nroc_proc)
    if ( associated(spdata%nl2g) ) &
         deallocate( spdata%nl2g)
    if ( associated(spdata%ng2l) ) &
         deallocate( spdata%ng2l)
    if ( associated(spdata%ng2p) ) &
         deallocate( spdata%ng2p)
  end subroutine delete_Data

  subroutine newBlockCyclicDistribution(Blocksize,Comm,this,name)
     !........................................
     ! Constructor
     !........................................
     type (OrbitalDistribution), intent(inout) :: this
     integer, intent(in)                       :: Blocksize
     integer, intent(in)                       :: Comm
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)

     this%data%blocksize = Blocksize
     this%data%comm      = Comm

#ifdef MPI
     call MPI_Comm_Rank( Comm, this%data%node, error )
     call MPI_Comm_Size( Comm, this%data%nodes, error )
#else
     this%data%node = 0
     this%data%nodes = 1
#endif
     this%data%node_io = 0

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(Distribution from BlockSize and Comm)"
     endif
     call tag_new_object(this)

   end subroutine newBlockCyclicDistribution

!-----------------------------------------------------------
   function num_local_elements(this,nels,Node) result(nl)
#ifdef MPI
     use mpi_siesta, only: MPI_COMM_Self
#endif
     type(OrbitalDistribution), intent(in)  :: this
     integer, intent(in)                    :: nels
     integer, intent(in), optional          :: Node
     integer                                :: nl

     integer :: lNode
     integer :: Remainder, MinPerNode, RemainderBlocks
     
     if ( present(Node) ) then
        lNode = Node
     else
        lNode = this%data%Node
     end if

     if (this%data%blocksize == 0) then
        ! Not a block-cyclic distribution
        if (nels /= this%data%n) then
           call die("Cannot figure out no_l if nels/=n")
           nl = -1
        endif
        if (.not. associated(this%data%nroc_proc)) then
           call die("Dist arrays not setup")
           nl = -1
        endif
        nl = this%data%nroc_proc(lNode)
#ifdef MPI
      else if ( this%data%Comm == MPI_Comm_Self ) then ! globalized local distribution

         ! No matter the node, it is the same orbital
         ! The interface for creating such an orbital distribution
         ! is 
         ! call newDistribution(nrows_g(sp),MPI_Comm_Self,dit)
         ! TODO Talk with Alberto about adding: "lproc = this%data%nodes - MOD(proc,this%data%nodes)"
         if ( nels /= this%data%blocksize ) &
              call die('Contact Nick Papior Andersen, nickpapior@gmail.com nel')
         nl = this%data%blocksize
#endif
      else  ! block-cyclic distribution

          ! Use Julian's code for now, although this is 
          ! really Scalapack's numroc  (copied at the end)
          ! nl = numroc(nels,this%data%blocksize,Node,this%data%isrcproc,this%data%nodes)

          MinPerNode = nels / (this%data%nodes*this%data%blocksize)
          Remainder = nels - MinPerNode * this%data%nodes*this%data%blocksize
          RemainderBlocks = Remainder / this%data%blocksize
          Remainder = Remainder - RemainderBlocks * this%data%blocksize
          nl = MinPerNode*this%data%blocksize
          if (lNode.lt.RemainderBlocks) nl = nl + this%data%blocksize
          if (lNode.eq.RemainderBlocks) nl = nl + Remainder

      endif
      end function num_local_elements

      function dist_node_(this) result(node)
        type(OrbitalDistribution), intent(in) :: this
        integer :: node
        node = this%data%node
      end function dist_node_
      function dist_nodes_(this) result(nodes)
        type(OrbitalDistribution), intent(in) :: this
        integer :: nodes
        nodes = this%data%nodes
      end function dist_nodes_

   function index_local_to_global(this,il,Node) result(ig)
#ifdef MPI
     use mpi_siesta, only: MPI_COMM_Self
#endif
     type(OrbitalDistribution), intent(in)  :: this
     integer, intent(in)                    :: il
     integer, intent(in), optional          :: Node
     integer                                :: ig

     integer :: lNode
     integer :: LBlock, LEle

     lNode = this%data%node
     if ( present(Node) ) lNode = Node

     if (this%data%blocksize == 0) then
        ! Not a block-cyclic distribution
        if (lNode /= this%data%node) then
           call die("Cannot figure out ig if Node/=my_proc")
           ig = -1
        endif
        if (.not. associated(this%data%nl2g)) then
           call die("Dist arrays not setup")
           ig = -1
        endif
        ig = this%data%nl2g(il)
#ifdef MPI
      else if ( this%data%Comm == MPI_Comm_Self ) then ! globalized local distribution

         ! No matter the node, it is the same orbital
         ! The interface for creating such an orbital distribution
         ! is 
         ! call newDistribution(nrows_g(sp),MPI_Comm_Self,dit)
         if ( il > this%data%blocksize ) &
              call die('Contact Nick Papior Andersen, nickpapior@gmail.com l2g')

         ig = il
#endif
      else  ! block-cyclic distribution

          ! Use Julian's code for now, although this is 
          ! really Scalapack's indxl2g  (copied at the end)
          ! ig = indxl2g(il,this%data%blocksize,Node,this%data%isrcproc,this%data%nodes)

          !  Find local block number
          LBlock = ((il -1)/this%data%blocksize)
          !  Substract local base line to find element number within the block
          LEle = il - LBlock*this%data%blocksize
          !  Calculate global index
          ig = (LBlock*this%data%nodes + lNode)*this%data%blocksize + LEle

      endif
      end function index_local_to_global


   function index_global_to_local(this,ig,Node) result(il)
#ifdef MPI
     use mpi_siesta, only: MPI_COMM_Self
#endif
     type(OrbitalDistribution), intent(in)  :: this
     integer, intent(in)                    :: ig
     ! In case Node is not supplied we expect it to request its
     ! own index
     integer, intent(in), optional          :: Node
     integer                                :: il

     integer :: lNode
     integer :: LBlock, LEle, GBlock, OrbCheck

     lNode = this%data%node
     if ( present(Node) ) lNode = Node

     if (this%data%blocksize == 0) then
        ! Not a block-cyclic distribution
        if (lNode /= this%data%node) then
           call die("Cannot figure out il if Node/=my_proc")
           il = -1
        endif
        if (.not. associated(this%data%ng2l)) then
           call die("Dist arrays not setup")
           il = -1
        endif
        il = this%data%ng2l(ig)
        ! Alternatively: Use an extended ng2p, in which "local" elements are negative.
        ! if (ng2p(ig)) < 0 then il = -that, else 0
        ! Example:  0 0 1 1 2 2 -1 -2 4 4 5 5 0 0 1 1 2 2 -3 -4 5 5 ..... (nb=2, np=6)
#ifdef MPI
      else if ( this%data%Comm == MPI_Comm_Self ) then ! globalized local distribution

         ! No matter the node, it is the same orbital
         ! The interface for creating such an orbital distribution
         ! is 
         ! call newDistribution(nrows_g(sp),MPI_Comm_Self,dit)
         if ( ig > this%data%blocksize ) &
              call die('Contact Nick Papior Andersen, nickpapior@gmail.com g2l')

         il = ig
#endif
      else  ! block-cyclic distribution

          ! Use Julian's code for now, although this is 
          ! really a combination of Scalapack's indxg2p and indxg2l  (copied at the end)
          ! owner = indxg2p(ig,this%data%blocksize,Node,this%data%isrcproc,this%data%nodes)
          ! if (owner == Node) then
          !   il = indxg2l(ig,this%data%blocksize,Node,this%data%isrcproc,this%data%nodes)
          ! else
          !   il = 0
          ! endif 

          !  Find global block number
          GBlock = ((ig -1)/this%data%blocksize)
          !  Substract global base line to find element number within the block
          LEle = ig - GBlock*this%data%blocksize
          !  Find the block number on the local node
          LBlock = ((GBlock - lNode)/this%data%nodes)
          !  Generate the local orbital pointer based on the local block number
          il = LEle + LBlock*this%data%blocksize
          !  Check that this is consistent - if it is not then this
          !  local orbital is not on this node and so we return 0
          !  to indicate this.
          OrbCheck = (LBlock*this%data%nodes + lNode)*this%data%blocksize + LEle
          if (OrbCheck.ne.ig) il = 0
      endif

      end function index_global_to_local

   function node_handling_element(this,ig) result(proc)
#ifdef MPI
     use mpi_siesta, only: MPI_COMM_Self
#endif
     type(OrbitalDistribution), intent(in)  :: this
     integer, intent(in)                    :: ig
     integer                                :: proc

     integer :: GBlock

     if (this%data%blocksize == 0) then
        ! Not a block-cyclic distribution
        if (.not. associated(this%data%ng2p)) then
           call die("Dist arrays not setup")
           proc = -1
        endif
        proc = this%data%ng2p(ig)
        ! Alternatively: Use an extended ng2p, in which "local" elements are negative.
        ! if (ng2p(ig)) < 0 then il = -that, else 0
        ! Example:  0 0 1 1 2 2 -1 -2 4 4 5 5 0 0 1 1 2 2 -3 -4 5 5 ..... (nb=2, np=6)
#ifdef MPI
      else if ( this%data%Comm == MPI_Comm_Self ) then ! globalized local distribution

         ! No matter the node, it is the same orbital
         ! The interface for creating such an orbital distribution
         ! is 
         ! call newDistribution(nrows_g(sp),MPI_Comm_Self,dit)
         if ( ig > this%data%blocksize ) &
              call die('Contact Nick Papior Andersen, nickpapior@gmail.com nhe')

         proc = this%data%node
#endif
      else  ! block-cyclic distribution

          ! Use Julian's code for now, although this is 
          ! really Scalapack's indxg2p (copied at the end)
          ! proc = indxg2p(ig,this%data%blocksize,dummy_Node,  &
          !                this%data%isrcproc,this%data%nodes)

          !  Find global block number
          GBlock = ((ig -1)/this%data%blocksize)
          !  Find the Node number that has this block
          proc = mod(GBlock,this%data%nodes)

      endif

      end function node_handling_element

   ! Returns the communicator assigned to this
   ! distribution
   function comm_(this) result(comm)
     type(OrbitalDistribution), intent(in) :: this
     integer :: comm
     comm = this%data%comm
   end function comm_
   
      subroutine printOrbitalDistribution(this)
        type(OrbitalDistribution), intent(in) :: this

        if (.not. initialized(this) ) then
           print "(a)", "Orbital distribution Not Associated"
           RETURN
        endif

        print "(a,/,t4,5(a,i0),a)", &
             "  <orb-dist:"//trim(this%data%name), &
             " comm=", this%data%comm, &
             " node/nodes=", this%data%node,' / ',this%data%nodes, &
             " blocksize=",this%data%blocksize, &
             ", refcount: ", refcount(this), ">"
      end subroutine printOrbitalDistribution

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

end module class_OrbitalDistribution
