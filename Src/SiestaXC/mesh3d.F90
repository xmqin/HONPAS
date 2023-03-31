!!@LICENSE
!
!******************************************************************************
! MODULE mesh3D
! Defines and handles different parallel distributions of a 3D mesh of points.
! Written by J.M.Soler. Jan-July.2008
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! addMeshData   : Adds the data in a box array to the equiv. unit-cell points
! associateMeshTask : Associates a communication task to mesh distr.
! copyMeshData  : Copies data in a unit-cell array to another mesh box array
! fftMeshDistr  : Sets four mesh distributions to perform FFTs 
! freeMeshDistr : Frees a distribution index (they are limited to maxDistr)
! freeMeshTask  : Frees a task index (they are limited to maxTasks)
! myMeshBox     : Returns mesh box of my processor node
! nodeMeshBox   : Returns mesh box of any given processor node
! redistributeMeshData : Copies data between different mesh distributions
! sameMeshDistr : Finds if two mesh distributions are equal
! setMeshDistr  : Defines a distribution of mesh points over parallel nodes
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use alloc,      only: de_alloc        ! De-allocation routine
! use alloc,      only: re_alloc        ! Re-allocation routine
! use sys,        only: die             ! Termination routine
! use sorting,    only: ordix           ! Finds the sorting index of a vector
! use parallel,  only: parallel_init    ! Sets node and nodes variables
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
! use precision,  only: gp=>grid_p      ! Grid-data real precision type
! use parallel,   only: myNode=>node    ! This process node
! use parallel,   only: totNodes=>nodes ! Total number of processor nodes
!
! Copy and add for limiting temporary arrays
! use m_array, only: array_copy, array_add
					
!   USED MPI procedures, interfaces and types:
! use mpi_siesta, only: MPI_AllGather
! use mpi_siesta, only: MPI_Send
! use mpi_siesta, only: MPI_Recv
! use mpi_siesta, only: MPI_COMM_WORLD
! use mpi_siesta, only: MPI_STATUS_SIZE
! use mpi_siesta, only: MPI_Integer
! use mpi_siesta, only: MPI_grid_real   ! Grid-data real precision type
!
!   EXTERNAL procedures used:
! none
!
!   Private fixed module parameters (see/set values in code):
! integer maxDistr      : Max. mesh distributions allocated at any given time
! integer maxTasks      : Max. communication tasks allocated at any given time
! integer maxDistrID    : Max. IDs assigned to the same distrib.
! integer maxTaskID     : Max. IDs assigned to the same task
! integer maxParts      : Max. parts of a mesh box, in different unit cells
! integer maxDistrTasks : Max. commun. tasks associated to a mesh distribution
!
!******************************************************************************
!
! SUBROUTINE addMeshData( nMesh, srcBox, srcData, dstDistr, dstData, task )
!
! Adds (reduces) the values defined in a 3D mesh box to the equivalent mesh 
!   points in a periodic unit cell. The input box may be entirely inside, or
!   partially or totally outside the the unit cell, and it is folded to add 
!   all its values to the corresponding unit-cell points.
! In parallel execution, the routine handles the interprocessor transfers
!   needed to send the box data to the nodes that will store them.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)    : Global mesh size in each axis
! integer srcBox(2,3) : Source mesh box: box(1,iAxis)=lower box limits
!                                        box(2,iAxis)=upper box limits
!                       The box may be partially or totally outside the first
!                       unit cell (FUC). It is folded to the FUC in order to
!                       assign the values of srcData to those of dstData.
! real(gp) srcData(0:,0:,0:,:) : Source data values to be copied. 
!                       The first three upper array bounds must be equal or 
!                       larger than srcBox(2,iAxis)-src(1,iAxis). Real type gp 
!                       is identical to grid_p, defined in module precision.
! integer dstDistr    : Distrib. ID of destination data, or zero for serial
!                       execution
!----------------------- INPUT and OUTPUT -------------------------------------
! real(gp) dstData(0:,0:,0:,:) : Destination data array. Its upper array bounds
!                       must be enough to contain the node box in the dstDistr
!                       parallel distrib or nMesh(iAxis)-1 in serial execution
!                       It is NOT initialized: its input values are incremented
!                       with the source array values.
!-------------------- OPTIONAL INPUT and OUTPUT -------------------------------
! integer task : ID of communication task. It allows to save information on an
!                optimized interprocessor communication pattern, required to
!                exchange data between source and destination mesh boxes.
!                This info is used in future calls with the same srcBox and
!                dstDistr. task value must be saved between calls. See also
!                associateMeshTask.
!----------------------------- UNITS ------------------------------------------
! Units of srcData and dstData are arbitrary but they must be equal
!----------------------------- USAGE ------------------------------------------
!   Typical usage for parallel execution is:
! use precision,    only: gp=>grid_p
! use mesh3D, only: addMeshData, associateMeshTask, myMeshBox, setMeshDistr
! integer:: box1(2,3), box2(2,3), myBox(2,3), nData, nMesh(3)
! integer,save:: myDistr=-1, task1=-1, task2=-1
! real(gp),allocatable:: myData(:,:,:,:), data1(:,:,:,:), data2(:,:,:,:)
! ... Find/define nMesh
! call setMeshDistr( myDistr, nMesh )
! call myMeshBox( nMesh, myDistr, myBox )
! allocate( myData(myBox(1,1):myBox(2,1), &
!                  myBox(1,2):myBox(2,2), &
!                  myBox(1,3):myBox(2,3),1:nData) )
! myData = 0
! ... Find/define box1 and box2
! allocate( data1(box1(1,1):box1(2,1), &
!                 box1(1,2):box1(2,2), &
!                 box1(1,3):box1(2,3),1:nData) )
! allocate( data2(box2(1,1):box2(2,1), &
!                 box2(1,2):box2(2,2), &
!                 box2(1,3):box2(2,3),1:nData) )
! ... Find data1 and data2
! call associateMeshTask( task1, myDistr )
! call associateMeshTask( task2, myDistr )
! call addMeshData( nMesh, box1, data1, myDistr, myData, task1 )
! call addMeshData( nMesh, box2, data2, myDistr, myData, task2 )
!   This usage is also valid in serial execution, although it could be
! simplified in this case (e.g. the call to setMeshDistr could be replaced by
! myDistr=0).
!   The optimization of the communication task is done at the end of the first
! call with given srcBox and dstDistr. Therefore, the use of the task argument
! is not advised unless the same srcBox and dstDistr are used repeatedly.
!
! In using some mesh3D routines, it is essential to understand that, in
!   fortran90, assumed-shape array arguments (like srcData and dstData in 
!   addMeshData) are passed with their shape (i.e. their dimension in each
!   index) but NOT with the lower or upper bounds of these indexes. Thus,
!   element a(6,1) of an array dimensioned as a(5:10,2) in the calling program,
!   and as a(0:,0:) in the called subroutine, will be accessed as a(1,0) in the
!   latter. This is why arrays dimensioned as srcData(0:,0:,0:,:) may be (and
!   are expected to be) dimensioned in the calling program as 
!   srcData(box(1,1):box(2,1),box(1,2):box(2,2),box(1,3):box(2,3),nData)
!---------------------------- BEHAVIOUR ---------------------------------------
! If dstDistr==0, it assumes serial execution, i.e. that dstData contains the
!   destination data values at all the mesh points in the unit-cell.
! It stops with an error message in any of these cases:
! - size(srcData,1:3) < srcBox(2,1:3)-srcBox(1,1:3)+1
! - size(dstData,1:3) < dstBox(2,1:3)-dstBox(1,1:3)+1, with dstBox as 
!   returned by call myMeshBox(nMesh,dstDistr,dstBox)
! - size(srcData,4) /= size(dstData,4)
! If task is present, it is automatically associated to dstDistr (see
!   associateMeshDistr). If dstDistr changes between calls (see freeMeshDistr),
!   task info is automatically erased.
! If task is present and it contains a predefined value (returned by a previous
!   call to addMeshData or copyMeshData), but srcBox has changed, it stops with
!   an error message (since the communication task has also changed). To avoid
!   this, reset task to a negative value. If srcBox derives from a defined
!   mesh distribution, and task had been associated to it, task info will be 
!   automatically deleted when the distrib changes, and no error will occur.
!--------------------------- ALGORITHMS ---------------------------------------
! The same as in copyMeshData, except that the destination boxes, rather than
!   the source boxes, are previously obtained from dstDistr.
!
!******************************************************************************
!
! SUBROUTINE associateMeshTask( taskID, distrID1, distrID2 )
!
! Associates a communication task to one or two parallel mesh distributions
!------------------------------ INPUT -----------------------------------------
! integer distrID1 : ID of first distribution to be associated with task
!-------------------------- OPTIONAL INPUT ------------------------------------
! integer distrID2 : ID of second distrribution to be associated
!----------------------- INPUT and OUTPUT -------------------------------------
! integer taskID   : ID of task to be associated with distribution(s)
!----------------------------- USAGE ------------------------------------------
! use mesh3D, only: associateMeshTask, setMeshDistr
! integer:: nMesh(3)
! integer,save:: distr1=-1, distr2=-1, task1=-1, task2=-1, task12=-1, task21=-1
! real(gp),allocatable:: workload(:,:,:)
! ... Find nMesh. Allocate and find workload
! call setMeshDistr( distr1, nMesh )
! call setMeshDistr( distr2, nMesh, workload )
! call associateMeshTask( task1, distr1 )
! call associateMeshTask( task2, distr2 )
! call associateMeshTask( task12, distr1, distr2 )
! call associateMeshTask( task21, distr1, distr2 )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution (totNodes==1) it does nothing
! If taskID is already associated with distrID1 (and distrID2, if present),
!   nothing is done. Thus, it is safe (and not inefficient) to make the same
!   call repeatedly.
! If the taskID value has not been defined nor associated with any distribution,
!   it is initialized before associating it, and it returns changed.
! If taskID is defined, but it is associated with two distributions different
!   from distrID1 (or from distrID1 and distrID2, if the latter is present),
!   the taskID is erased (freed) and reinitialized.
! It stops with an error message in the following cases:
!   - If distrID1 (or distrID2 if present) has not been defined or if it has
!     been erased (freed).
!   - If the available maxDistrTasks slots, for tasks associated to distrID1
!     (or distrID2, if present), is overflooded.
!--------------------------- ALGORITHMS ---------------------------------------
! The taskID is an identifier for a parallel communication "task", defined as
!   a set of source and destination mesh boxes that belong to each processor.
! The same task can be identified by different IDs (for example by different
!   calling programs, but the information is stored only once. This is ensured
!   by comparing the boxes of a new taskID with those of all previously defined
!   tasks. If they coincide, the taskID is simply assigned to that task.
! When a taskID is erased (freed), the task itself is not erased, unless all 
!   its IDs have been erased. There may be up to maxTasks possible simultaneous 
!   tasks, each with up to maxTaskID IDs.
! A task may be associated with up to two mesh distributions. This association 
!   state is stored in both the task and the distribution(s) info. The 
!   association means that the task's source or destination boxes are linked
!   to the distribution. This does not necessarily mean that the boxes are the 
!   node's boxes in the distribution. For example, they might be "wings" to
!   the left or to the right of those boses. The key point is that the task
!   boxes will change when the distribution changes.
! The purpose of the task-distribution association is to erase automatically
!   the task when the distribution is erased (freed). This is convenient
!   (and safer) because it allows not to have to free the task explicitly.
!
!******************************************************************************
!
! SUBROUTINE copyMeshData( nMesh, srcDistr, srcData, dstBox, dstData, task )
!
! Copies a box of values defined in a 3D mesh of a periodic unit cell.
! The box limits may be entirely inside, or partially or totally outside 
!   the unit cell, which is periodically repeated to fill the box.
! In parallel execution, the routine handles the interprocessor transfers
!   needed to bring the requested box data from the nodes that store them.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)    : Global mesh size in each axis
! integer srcDistr    : Distrib. ID of source data, or zero for serial exec.
! real(gp) srcData(0:,0:,0:,:) : Source data values to be copied. Its upper
!                       array bounds must be enough to contain the node box in
!                       the srcDistr parallel distrib, or nMesh(iAxis)-1 in
!                       serial execution. Real type gp is identical to grid_p, 
!                       defined in module precision.
! integer dstBox(2,3) : Destination mesh box: box(1,iAxis)=lower box limits
!                                             box(2,iAxis)=upper box limits
!                       The box may be partially or totally outside the first
!                       unit cell (FUC). It is folded to the FUC in order to
!                       assign the values of srcData to those of dstData.
!----------------------------- OUTPUT -----------------------------------------
! real(gp) dstData(0:,0:,0:,:) : Destination data array (different from srcData)
!                       The first three upper array bounds must be equal or 
!                       larger than dstBox(2,iAxis)-dst(1,iAxis). 
!-------------------- OPTIONAL INPUT and OUTPUT -------------------------------
! integer task : ID of communication task. It allows to save information on an
!                optimized interprocessor communication pattern, required to
!                exchange data between source and destination mesh boxes.
!                This info is used in future calls with the same srcBox and
!                dstDistr. task value must be saved between calls. See also
!                associateMeshTask.
!----------------------------- UNITS ------------------------------------------
! Units of srcData and dstData are arbitrary but they must be equal
!----------------------------- USAGE ------------------------------------------
!   A typical usage for data redistribution in parallel execution is:
! use precision,    only: gp=>grid_p
! use mesh3D, only: associateMeshTask, copyMeshData, myMeshBox, setMeshDistr
! integer:: box1(2,3), box2(2,3), nData, nMesh(3)
! integer,save:: distr1=-1, distr2=-1, task=-1
! real(gp),allocatable:: data1(:,:,:,:), data2(:,:,:,:)
! ... Find/define nMesh, box1 and box2 (so that they are two partitions of mesh)
! call setMeshDistr( distr1, nMesh, box1 )
! call setMeshDistr( distr2, nMesh, box2 )
! allocate( data1(box1(1,1):box1(2,1), &
!                 box1(1,2):box1(2,2), &
!                 box1(1,3):box1(2,3),1:nData) )
! allocate( data2(box2(1,1):box2(2,1), &
!                 box2(1,2):box2(2,2), &
!                 box2(1,3):box2(2,3),1:nData) )
! ... Find data1 in box1
! call associateMeshTask( task, distr1, distr2 )
! call copyMeshData( nMesh, distr1, data1, box2, data2, task )
!   This usage is still valid in serial execution, although it makes little 
! sense in this case, since both mesh distributions would be identical.
!   The optimization of the communication task is done at the end of the first
! call with given srcDistr and dstBox. Therefore, the use of the task argument
! is not advised unless the same srcDistr and dstBox are used repeatedly.
!---------------------------- BEHAVIOUR ---------------------------------------
! If srcDistr==0, it assumes serial execution, i.e. that srcData contains the
!   source data values at all the mesh points in the unit-cell.
! It stops with an error message in any of these cases:
! - size(srcData,1:3) < srcBox(2,1:3)-srcBox(1,1:3)+1, with srcBox as 
!   returned by call myMeshBox(nMesh,srcDistr,srcBox)
! - size(dstData,1:3) < dstBox(2,1:3)-dstBox(1,1:3)+1
! - size(srcData,4) /= size(dstData,4)
! If task is present, it is automatically associated to srcDistr (see
!   associateMeshDistr). If srcDistr changes between calls (see freeMeshDistr),
!   task info is automatically erased.
! If task is present and it contains a predefined value (returned by a previous
!   call to addMeshData or copyMeshData), but dstBox has changed, it stops with
!   an error message (since the communication task has also changed). To avoid
!   this, reset task to a negative value. If dstBox derives from a defined
!   mesh distribution (as in the example above), and task had been associated
!   to it, task info will be automatically deleted when the distrib changes,
!   and no error will occur.
!--------------------------- ALGORITHMS ---------------------------------------
! - The source boxes are first obtained from srcDistr.
! - All the nodes' source and destination boxes are gathered, redistributed, 
!   and stored in each other node. These boxes are also identified by a taskID,
!   which is returned, for retrieval in future calls with the same task 
!   (defined by all the source and destination boxes themselves).
! - A universal (dependent only on totNodes) ordered sequence of all-to-all
!   inter-node communication transfers is found. It is designed to minimize
!   communication collisions. It is based on timing together transfers
!   between processors separated by a given index span. 
! - In the first call, the true, nonempty transfers for a given task are 
!   selected. Then, the order of the true transfers is (will be) reoptimized 
!   based on a graph-theory algorithm. This optimized sequence is then stored 
!   as part of the task information.
! - For each transfer (true or temptative):
!   - The source box (my node's in sending transfers, or that of the sending 
!     node, in receiving transfers) and the destination box (mine when 
!     receiving, or the receiving node's when sending) are each divided in
!     nParts, each in a different periodic unit cell. These part boxes are 
!     then shifted to the first unit cell (FUC).
!   - For each pair of the partial source and destination boxes, the
!     intersection between their FUC images is found. If the intersection is
!     empty (only happens if task info is not available), nothing is done.
!   - When sending, the corresponding source data values are added to a buffer,
!     that is finally sent with all parts together. When receiving, the buffer
!     is splitted in the corresponding intersection parts and added to the
!     destination data array.
!
!******************************************************************************
!
! SUBROUTINE fftMeshDistr( nMesh, fftDistr, axisDistr )
!
! Creates and returns a homogeneaous 3D mesh distribution fftDistr, and three
! 2D distributions axisDistr, whose boxes extend wholly over one of the three
! cell axes. These axisDistr allow fully local 1D FFTs in one of the axes.
! The fftDistr and axisDistr are designed to make already optimal (without any
! collisions) the default sequence of inter-node communications, returned by
! all2allTransferOrder. This is beacause all simultaneous communications occur
! between nodes separated by the same span (difference between node indexes).
! Thus, no task arguments are required to optimize such communications.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3) : Mesh divisions in each axis (including siesta "subpoints")
!------------------------- INPUT and OUTPUT -----------------------------------
! integer fftDistr : ID of parallel 3D mesh distribution to be used in 3D FFTs
!--------------------- OPTIONAL INPUT and OUTPUT ------------------------------
! integer axisDistr(3) : 2D distributions for the 1D FFTs required in 3D FFTs
!----------------------------- USAGE ------------------------------------------
!   Typical usage to find my box of k points:
! use mesh3D, only: fftMeshDistr, myMeshBox
! integer:: kBox(2,3), nMesh(3)
! integer,save:: fftDistr=-1
! ...Find nMesh
! call fftMeshDistr( nMesh, fftDistr )
! call myMeshBox( nMesh, fftDistr, kBox )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution (totNodes==1) it simply returns fftDistr=0.
! If the input value of fftDistr (and axisDistr, if present) are already 
!   defined for that nMesh, it returns without changing them. In this case, it
!   does NOT check that the distributions are suited for FFTs.
! Since it calls setMeshDistr, some of the behaviour of that routine is
!   inherited by fftMeshDistr.
!--------------------------- ALGORITHMS ---------------------------------------
! First, it uses setMeshDistr, with default input parameters, to create an 
!   optimized uniform 3D distr fftDistr. 
! Then, if axisDistr is present: for the x axis, it divides each row of node
!   boxes along that axis (say nNodeX boxes) uniformly in the other two axes.
!   Then it assigns one of these new long boxes (that extend along the whole x
!   axis) to each of the nNodeX processors, that owned the primitive 3D boxes.
!   This new boxes define axisDistr(1). The same is done for the y and z axes
!   (stricly, the axes are the cell axes, not x, y, and z).
!
!******************************************************************************
!
! SUBROUTINE freeMeshDistr( distrID )
!
! Erases a mesh distribution ID and frees the distribution slot if 
! no other IDs are assigned to it
!
!------------------------------ INPUT -----------------------------------------
! integer distrID  ! ID of an allocated distribution
!----------------------------- USAGE ------------------------------------------
! use mesh3D, only: free MeshDistr, setMeshDistr
! integer:: nMesh(3)
! integer,save:: myDistr=-1
! ...Find nMesh
! call setMeshDistr( nMesh, myDistr )
! ...Use myDistr
! call freeMeshDistr( myDistr )
!---------------------------- BEHAVIOUR ---------------------------------------
! If distrID is not defined, it returns without doing anything
! If the distribution has only that distrID, its information is erased.
!   Otherwise, only that particular ID is removed from it. 
!--------------------------- ALGORITHMS ---------------------------------------
! See setMeshDistr
!
!******************************************************************************
!
! SUBROUTINE myMeshBox( nMesh, distrID, box )

! Finds the mesh box of the local processor in a parallel mesh distribution.
! Equivalent to nodeMeshBox with node=myNode
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3) : Mesh divisions in each axis
! integer distrID  : Mesh-distribution ID
!----------------------------- OUTPUT -----------------------------------------
! integer box(2,3) : Mesh box: box(1,:)=lower bounds, box(2,:)=upper bounds
!----------------------------- USAGE ------------------------------------------
! Same as nodeMeshBox but without need of using myNode
!---------------------------- BEHAVIOUR ---------------------------------------
! Same as nodeMeshBox
!--------------------------- ALGORITHMS ---------------------------------------
! It simply calls nodeMeshBox with node=myNode
!
!******************************************************************************
!
! SUBROUTINE nodeMeshBox( nMesh, distrID, node, box )
!
! Finds the mesh box belonging to any node in a parallel mesh distribution,
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3) : Mesh divisions in each axis
! integer distrID  : Mesh-distribution ID
!----------------------------- OUTPUT -----------------------------------------
! integer box(2,3) : Mesh box: box(1,:)=lower bounds, box(2,:)=upper bounds
!----------------------------- USAGE ------------------------------------------
!   Typical usage to find my box of mesh points:
! use mesh3D, only: nodeMeshBox
! use parallel,     only: myNode=>node
! integer:: myBox(2,3), nMesh(3)
! integer,save:: myDistr=-1
! ...Find nMesh and myDistr
! call nodeMeshBox( nMesh, myDistr, myNode, myBox )
!---------------------------- BEHAVIOUR ---------------------------------------
! If iDistr=0, it returns box(1,:)=0 and  box(2,:)=nMesh(:)-1 
!   Note: although nMesh is available from the distr. data, it is also included
!   in the routine arguments, to handle directly the case iDistr=0.
! If node stores no mesh points, its (empty) box returns with box(1,:)>box(2,:)
!--------------------------- ALGORITHMS ---------------------------------------
! The corresponding box is simply retrived from the ditrib. data stored.
!
!******************************************************************************
!
! SUBROUTINE redistributeMeshData( srcDistr, srcData, dstDistr, dstData, task )
!
! Copies data values on the mesh with a new mesh distribution among nodes
!------------------------------ INPUT -----------------------------------------
! integer  srcDistr          : Distrib. ID of srcData
! real(gp) srcData(:,:,:,:)  : Source data. Must be declared target or pointer.
! integer  dstDistr          : Distrib. ID of dstData
!------------------------- OUTPUT POINTER -------------------------------------
! real(gp) dstData(:,:,:,:)  : Output destination array 
!                            : (it must be different from srcData)
!--------------------- OPTIONAL INPUT and OUTPUT ------------------------------
! integer, task              : ID of communication task
!----------------------------- UNITS ------------------------------------------
! Units of srcData and dstData are arbitrary but they must be equal
!----------------------------- USAGE ------------------------------------------
!   Typical usage to 'translate' between two distributions:
! use precision,    only: gp=>grid_p
! use mesh3D, only: associateMeshTask, redistributeMeshData
! integer:: nData, newBox(2,3), nMesh(3)
! integer,save:: newDistr=-1, old2new=-1, oldDistr=-1
! real(gp),pointer:: newData(:,:,:,:), oldData(:,:,:,:)
! ...Find nMesh, nData, oldDistr. Allocate and find oldData
! call setMeshBox( newDistr, nMesh )
! call myMeshBox( nMesh, newDistr, newBox )
! allocate( newData( newBox(1,1):newBox(2,1), &
!                    newBox(1,2):newBox(2,2), &
!                    newBox(1,3):newBox(2,3), nData )
! call associateMeshTask( old2new, oldDistr, newDistr )
! call redistributeMeshData( oldDistr, oldData, newDistr, newData, old2new )
!
!   In serial execution, srcData and dstData will be the same physical array
! upon return. This allows to write a common serial-parallel code, without
! increasing the required memory in the serial case. But this must be handled
! with much care, to avoid erasing or overwritting dstData inadvertently.
! Therefore, the use of copyMeshData is prefferable in most other situations.
!---------------------------- BEHAVIOUR ---------------------------------------
! If srcDistr==0 and dstDistr==0 (serial execution) it simply makes dstData 
!   point to srcData. This allows to save memory, but it must be handled with 
!   care in calling program.
! If the two distributions are equal (sameMeshDistr(srcDistr,dstDistr)==.true.)
!   it simply (re)allocates dstData and copies dstData=srcData. Othewise,
!   it calls copyMeshData after the (re)allocation of dstData.
!--------------------------- ALGORITHMS ---------------------------------------
! First, the dstDistr mesh box is obtained, then dstData is (re)allocated and
!   finally, copyMeshData is used to obtain dstData from srcData
!
!******************************************************************************
!
! logical FUNCTION sameMeshDistr( ID1, ID2 )
!
! Compares two distribution IDs
!------------------------------ INPUT -----------------------------------------
! integer ID1, ID2  ! Mesh distribution IDs to be compared
!----------------------------- USAGE ------------------------------------------
! use mesh3D, only: sameMeshDistr
! integer,save:: oldDistr=-1, newDistr=-1
! ... Find oldDistr and newDistr
! if (.not.sameMeshDistr(oldDistr,newDistr)) do something
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution (ID1=ID2=0), it returns .true.
! If ID1 and ID2 are equal, but not defined distributions, it returns .false.
!--------------------------- ALGORITHMS ---------------------------------------
! In parallel execution, it looks for the distributions that own the given IDs,
!   then it checks whether all their mesh boxes are the same.
!
!******************************************************************************
!
! SUBROUTINE setMeshDistr( distrID, nMesh, box, firstNode, nNodes, &
!                          nNodesX, nNodesY, nNodesZ, nBlock, &
!                          wlDistr, workload )
!
! Defines a parallel distribution of mesh points over processor nodes.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)  : Mesh divisions in each axis (including siesta "subpoints")
!------------------------- OPTIONAL INPUT -------------------------------------
! integer box(2,3)  : Mesh box of my processor node: 
!                     box(1,iAxis)=lower box limits, in range (0:nMesh(iAxis)-1)
!                     box(2,iAxis)=upper box limits, in same range
! integer firstNode : First node in the mesh distr.
! integer nNodes    : Total nodes in the mesh distr
! integer nNodesX   : Nodes in the X (first) axis.
!                     Must be present if present(nNodesYZ)
! integer nNodesY   : Nodes in the Y (second) axis.
!                     Must be present if present(nNodesZ)
! integer nNodesZ   : Nodes in the Z (third) axis
! integer nBlock    : Size of blocks of mesh points, in each axis, which are
!                     not splitted in different nodes. It must be a factor of
!                     all of nMesh(1:3). If box is also present, nBlock must be
!                     a factor of all box(1,1:3) and box(2,1:3)+1.
!                     nBlock corresponds to nsm (lateral size of "superpoints")
!                     in siesta mesh terminology.
! integer wlDistr   : Distr. index of workload array
! real(gp) workload(0:,0:,0:) ! Approx. relative workload of mesh points. 
!                     Must be nonnegative at all points and have nonzero sum.
!-------------------------- INPUT and OUTPUT ----------------------------------
! integer distrID : ID assigned to the mesh distrib.
!----------------------------- USAGE ------------------------------------------
!    Arguments box, firstNode, nNodes, nNodesXYZ, and nBlock are provided to
! force compatibility with user distributions, but they should be avoided
! otherwise, since they difficult the optimal distribution of mesh points.
!
!    Typical usage to create a uniform 3D mesh distribution:
! use mesh3D, only: setMeshDistr
! integer:: nMesh(3)
! integer,save:: myDistr=-1
! ... Find nMesh
! call setMeshDistr( myDistr, nMesh )
!
!    Typical usage to distribute evenly the workload of each node:
! use mesh3D, only: setMeshDistr
! integer:: nMesh(3)
! integer,save:: newDistr=-1, oldDistr=-1
! real(gp),allocatable:: workload(:,:,:)
! ... Find nMesh
! call setMeshDistr( oldDistr, nMesh )
! call myMeshBox( oldDistr, nMesh, box )
! allocate( workload(box(1,1):box(2,1), &
!                    box(1,2):box(2,2), &
!                    box(1,3):box(2,3)) )
! ... Find approximate CPU workload associated to each point of box
! call setMeshDistr( newDistr, nMesh, wlDistr=oldDistr, workload=workload )
! 
!   The mesh distribution of J.D.Gale, used up to siesta-2.5, is set by:
! use parallel, only: ProcessorY
! use mesh,     only: nsm
! use mesh3D,   only: setMeshDistr
! integer:: nMesh(3)
! integer,save:: JDGdistr=-1
! call setMeshDistr( JDGdistr, nMesh, nNodesX=1, nNodesY=ProcessorY, 
!                    nBlock=nsm )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution (totNodes==1) it simply returns distrID=0, irrespective
!   of all arguments. Parameter totNodes is obtained from module parallel.
! If the input distribution ID is still valid (i.e. consistent with the other
!   input arguments), the same value is returned. If it points to an existing
!   distribution that is no longer consistent, the old distribution ID is 
!   freed before returning with a new distrID. This makes it convenient to 
!   make succesive calls with the same distrID but different other arguments
!   (e.g. different workloads in different iterations).
! New IDs are never repeated, even if they identify the same distribution.
! If box is present, all other optional arguments are ignored. The different
!   node mesh boxes should be a nonoverlapping partition of all the mesh 
!   points, but this is NOT checked.
! If firstNode is not present, firstNode=0 is assumed.
! If nNodes is not present all nodes are used (nNodes=totNodes).
! It stops with an error message in the following cases:
!   - nNodes is present and it is smaller than 1 or larger than totNodes.
!   - nNodesZ is present but either nNodesX or nNodesY are not present.
!   - nNodesY is present but nNodesX is not present.
!   - nBlock is present and it is not a factor of one of nMesh(1:3).
!   - box and nBlock are present, and nBlock is not a factor of any of 
!     box(1,1:3) or box(2,1:3)+1.
!   - workload is present but wlDistr is not present.
!   - workload is present and negative at any mesh point, or zero at all 
!     points of the unit cell mesh.
! If nNodesXYZ are not present, the distribution of nodes over each axis is
!   optimized, in the sense of leading to node mesh boxes as cubic as possible. 
!   If only nNodesX is present, it is optimized over the Y and Z axes.
! If nBlock is not present, nBlock=1 is assumed.
! If workload is present, its distribution over processors is optimized, in the
!   sense of load balance. If it is not present, the even distribution of mesh 
!   points is optimized (equivalent to using the same workload for all points).
!--------------------------- ALGORITHMS ---------------------------------------
! For uniform distributions (workload not present):
! - nNodes is first fatorized in its prime factors.
! - All possible distributions of (products of) factors over axes are tried,
!   and that with lowest dispersion of nMesh(axis)/factor(axis) is selected.
!   This is done only over axes not constrained by nNodesXYZ, if present, and
!   respecting that axis factors must be multiples of nBlock, if present.
! - nMesh(axis) is divided in nNodes(axis) boxes. If there is a rest, the first
!   rest(axis) of the nNodes(axis) are given an extra point. Again, this is 
!   done with blocks of nBlock points (if present) rather than single points.
!
! For nonuniform distributions (workload present):
! - nNodes is first fatorized in its prime factors.
! - Beginning with a box of the entire cell size, each box is divided in 
!   nParts=factor (in order of decreasing factors). The division axis is that
!   along which the spatial dispersion of the box workload is maximum. The 
!   division surfaces are chosen to split the workload as evenly as possible.
!   If the division surface is in a void region (with zero workload), it is
!   placed at the midpoint between the nonzero workload regions. 
! - The divided boxes are assigned to a contiguous set of 
!   boxNodes=boxNodes/nParts. This process is repeated until boxNodes=1.
!
! A stored distribution is defined by nMesh and its node boxes. There may be up
!   to maxDistr distributions defined simultaneously. Each distribution may be
!   identified simultaneously by up to maxDistrID identifiers (for example by 
!   different calling routines), but its information is stored only once. 
!   This is ensured by comparing the boxes of a new distrID with those of all
!   previously defined distributions. If they coincide, the distrID is simply
!   added to the existing distribution.
! Each distrID is never repeated. When a distrID is erased (freed), the
!   distribution itself is not erased, unless all its IDs have been erased.
!   This prevents that a calling routine may erase a distribution that is
!   still being used by other routines.
!
!******************************************************************************

MODULE mesh3D

! Used module procedures
  use alloc,     only: de_alloc        ! De-allocation routine
  use alloc,     only: re_alloc        ! Re-allocation routine
  use sys,       only: die             ! Termination routine
  use sorting,   only: ordix           ! Sorting routine
  use parallel,  only: parallel_init   ! Sets node and nodes variables

  ! Copy and add for limiting temporary arrays
  use m_array,   only: array_add, array_copy

! Used module parameters
  use precision, only: dp              ! Real double precision type
  use precision, only: gp=>grid_p      ! Grid-data real precision type
  use parallel,  only: myNode=>node    ! This process node
  use parallel,  only: totNodes=>nodes ! Total number of processor nodes
#ifdef DEBUG_XC
  use debugXC,   only: udebug          ! Output file unit for debug info
#endif /* DEBUG_XC */

! Used MPI procedures and types
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

! Public procedures
PUBLIC:: &
  addMeshData,   &! Adds the data in a box array to the equiv. unit-cell points
  associateMeshTask, &! Associates a communication task to mesh distr.
  copyMeshData,  &! Copies data in a unit-cell array to another mesh box array
  fftMeshDistr,  &! Sets four mesh distributions to perform FFTs 
  freeMeshDistr, &! Frees a distribution index (they are limited to maxDistr)
  freeMeshTask,  &! Frees a task index (they are limited to maxTasks)
  myMeshBox,     &! Returns mesh box of my processor node
  nodeMeshBox,   &! Returns mesh box of any given processor node
  redistributeMeshData, &! Copies data between different mesh distributions
  sameMeshDistr, &! Finds if two mesh distributions are equal
  setMeshDistr    ! Defines a distribution of mesh points over parallel nodes

! Public types, parameters, or variables:
! none

PRIVATE ! Nothing is declared public beyond this point

  ! Private module parameters
  character(len=*),parameter:: moduleName = 'mesh3D '
  integer,parameter:: maxDistr   = 20 ! Max. number of mesh distributions
                                      !   allocated at any given time
  integer,parameter:: maxTasks  = 100 ! Max. number of communication tasks
                                      !   allocated at any given time
  integer,parameter:: maxDistrID = 20 ! Max. IDs assigned to the same distrib.
  integer,parameter:: maxTaskID  = 10 ! Max. IDs assigned to the same task
  integer,parameter:: maxParts  = 125 ! Max. parts of a mesh box
  integer,parameter:: maxDistrTasks = 50 ! Max. communication tasks associated
                                         !   to a parallel mesh distribution

  ! Private type to hold mesh distribution data
  type distrType
    private
    logical:: defined=.false.   ! Has this distribution been defined?
    integer:: ID(maxDistrID)=-1 ! ID numbers assigned to the distribution
    integer:: nNodes=-1         ! Number of nodes participating in the distrib.
    integer:: firstNode=-1      ! First of the consecutive nNodes in the distr.
    integer:: nMesh(3)=-1       ! Number of total mesh divisions in each axis
    integer:: task(maxDistrTasks)=-1 ! Communic. tasks associated to distr.
    integer,pointer:: box(:,:,:)=>null() ! Mesh box bounds of each node:
                                         ! box(1,iAxis,iNode)=lower bounds
                                         ! box(2,iAxis,iNode)=upper bounds
  end type distrType

  ! Private type to hold a mesh communication task
  type taskType
    private
    logical:: defined=.false.    ! Has this task been defined?
    logical:: optimized=.false.  ! Has the transfer order been optimized?
    logical:: associated=.false. ! Has this task been associated to any distr?
    integer:: ID(maxTaskID)=-1   ! ID numbers assigned to the distribution
    integer:: distr(2)=-1        ! Distr. indexes to which task is associated
    integer:: nTrsf=0            ! Number of needed communication transfers
    integer,pointer:: srcBox(:,:,:)=>null() ! Source mesh box of each node.
    integer,pointer:: dstBox(:,:,:)=>null() ! Destination mesh box of each node.
    integer,pointer:: trsfNode(:) ! Sequence of to/from transfer nodes
    integer,pointer:: trsfDir(:) ! Sequence of to/from transfer directions
  end type taskType

  ! Private array that will hold the different distributions data
  type(distrType),target,save:: storedMeshDistr(maxDistr)

  ! Private array that will hold the different communication tasks
  type(taskType),target,save:: storedMeshTask(maxTasks)

  ! Private counters
  integer,save:: nDistrID = 0  ! Number of assigned distribution IDs
  integer,save:: nTaskID = 0  ! Number of assigned mesh communication task IDs

#ifdef MPI
  integer:: MPIerror, MPIstatus(MPI_STATUS_SIZE)
#endif

CONTAINS

!-------------------------------------------------------------------------------

subroutine addMeshData( nMesh, srcBox, srcData, dstDistr, dstData, task )

! Adds (reduces) the values defined in a 3D mesh box to the equivalent mesh 
! points in a periodic unit cell. 

  implicit none
  integer, intent(in) :: nMesh(3)          ! Global mesh size in each axis
  integer, intent(in) :: srcBox(2,3)       ! Source mesh box:
                                           ! box(1,iAxis)=lower box limits
                                           ! box(2,iAxis)=upper box limits
                                           ! Box may be larger than the
                                           ! (0:nMesh-1) unit-cell range
  real(gp),intent(in) :: srcData(0:,0:,0:,:) ! Source data values to be copied.
                                           ! The first three upper array bounds 
                                           ! must be equal or larger than
                                           ! srcBox(2,iAxis)-src(1,iAxis)
  integer, intent(in) :: dstDistr          ! Distrib. ID of destination data
                                           ! or zero for serial execution
  real(gp),intent(inout):: dstData(0:,0:,0:,:) ! Destination data array
                                           ! Its upper array bounds must be
                                           ! enough to contain the node box in
                                           ! the dstDistr parallel distrib.
                                           ! or nMesh(iAxis)-1 in serial exe
  integer,intent(inout),optional:: task    ! ID of communication task

  integer :: dstBox(2,3)

! Associate task to distribution
  if (present(task)) call associateMeshTask( task, dstDistr )

! Find box limits of destination data in this node
  call myMeshBox( nMesh, dstDistr, dstBox )

! Add srcData to dstData
  call reduceData( nMesh, srcBox, srcData, dstBox, dstData, taskID=task )

end subroutine addMeshData

!-------------------------------------------------------------------------------

subroutine all2allTransferOrder( nNodes, node, nTrsf, trsfNode, trsfDir )

! Finds an ordered sequence of matching all-to-all send-receive transfers
! that avoids collisions and that approximately minimizes total transfer time
! Algorithm:
! At any given 'time' (index of the ordered sequence), half of the nodes send
! right-wards (to higher nodes) and the other half receives from leftwards,
! with a fixed internode span, modulo nNodes2 (defined below). The sending/
! receiving nodes belong to odd/even blocks (or viceversa) of consecutive nodes,
! whose blockSize is the power-of-two factor of the span (i.e. span is the
! product of the blocSize times an odd number). nNodes2 is the smallest 
! multiple of 2*blockSize equal to, or larger than nNodes. This ensures opposite
! parities of the blocks of sendNode and recvNode=modulo(sendNode+span,nNodes2).

  implicit none
  integer,intent(in) :: nNodes   ! Number of parallel processor nodes. The
                                 ! number of transfers per node is 2*nNodes-1
  integer,intent(in) :: node     ! Local processor node, in range (0:nNodes-1)
  integer,intent(out):: nTrsf    ! Number of required transfers
  integer,intent(out):: trsfNode(0:2*(nNodes-1)) ! To/from transfer node.
                                 ! First 'transfer' is local: trsfNode(0)=node
  integer,intent(out):: trsfDir(0:2*(nNodes-1))  ! Transfer direction:
                                 ! trsfDir=+1 => from node to trsfNode
                                 ! trsfDir=-1 => from trsfNode to node
                                 ! trsfDir=0  => local 'transfer' within node

  integer:: blockSize, iTrsf, nNodes2, nodePar, oddFac, pow2, span, trsfPar

  ! Set num. of transfers needed (a send and a receive to/from each other node)
  nTrsf = 2*(nNodes-1)

  ! Special local transfer
  trsfNode(0) = node
  trsfDir(0)  = 0
  if (nNodes==1) return

  ! Loop on block sizes (powers of 2). iTrsf is the transfer-order index.
  iTrsf = 1
  do pow2 = 0,int(log(real(nNodes))/log(2.))+1
    blockSize = 2**pow2
    ! Find smallest multiple of 2*blockSize equal to, or larger than nNodes
    nNodes2 = ((nNodes-1)/(2*blockSize)+1) * (2*blockSize)
    ! Find parity of node's block
    nodePar = mod(node/blockSize,2)
    ! Loop on odd factors of span. Maximum span is nNodes2-1
    ! Fake nodes are 'added' between nNodes and nNodes2-1
    do oddFac = 1,nNodes2/blockSize,2
      span = blockSize * oddFac
      ! Loop on transfer 'parity' (direction)
      do trsfPar = 0,1
        ! Find transfer sense, and the node sending/receiving to/from node
        ! All transfers are rightwards, with transfer node modulo mNodes2
        if (nodePar==trsfPar) then ! node sends
          trsfNode(iTrsf) = modulo( node+span, nNodes2 )
          trsfDir(iTrsf) = +1
        else ! node receives
          trsfNode(iTrsf) = modulo( node-span, nNodes2 )
          trsfDir(iTrsf) = -1
        end if ! (nodePar==trsfPar)
        ! If trsfNode is within range (not a fake node), accept this transfer
        if (trsfNode(iTrsf)<nNodes) iTrsf = iTrsf + 1
        ! Finish when all transfers have been found
        if (iTrsf>2*(nNodes-1)) return
      end do ! trsfPar
    end do ! oddFact
  end do ! pow2

end subroutine all2allTransferOrder

!-------------------------------------------------------------------------------

subroutine associateMeshTask( taskID, distrID1, distrID2 )

! Associates a communication task to one or two parallel mesh distributions

  implicit none
  integer,         intent(inout):: taskID   ! ID of task to be associated
  integer,         intent(in)   :: distrID1 ! ID of first distr, to be assoc.
  integer,optional,intent(in)   :: distrID2 ! ID of second distr. to be assoc.

  character(len=*),parameter:: myName = 'associateMeshTask '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  logical:: found
  integer:: iDistr1, iDistr2, it, iTask
  type(taskType), pointer :: task
  type(distrType), pointer :: distr1, distr2

! In serial execution, do nothing
  if (totNodes<2) return

! Find distribution(s)
  iDistr1 = indexDistr( distrID1 )
  if (iDistr1<0) call die(errHead//'distrID1 not defined')
  distr1 => storedMeshDistr(iDistr1)  ! Just a shorter name
  if (present(distrID2)) then
    iDistr2 = indexDistr( distrID2 )
    if (iDistr2<0) call die(errHead//'distrID2 not defined')
    distr2 => storedMeshDistr(iDistr2)
  end if ! (present(distrID2))

! Find task
  found = .false.
  iTask = indexTask( taskID )
  if (iTask>0) then  ! Task exists already
    task => storedMeshTask(iTask)  ! Just a shorter name
    found = .true.
    if (.not.any(iDistr1==task%distr)) found = .false.
    if (present(distrID2)) then
      if (.not.any(iDistr2==task%distr)) found = .false.
    end if
    if (.not.found) & ! Free existing task ID, since it has new associations
      call freeMeshTask( taskID )
  end if ! (iTask>0)

! Create new task, if necessary
  if (found) then
    return  ! Since task is alredy associated to same distribution(s)
  else      ! Task did not exist or has been eliminated => create it
    call initTask( taskID )
    iTask = indexTask( taskID )
    task => storedMeshTask(iTask)
  end if ! (found)

! Store association to first distribution
  task%distr(1) = iDistr1  ! Associate distribution to task
  if (.not.any(distr1%task==iTask)) then  ! Associate task to distribution
    found = .false.
    do it = 1,maxDistrTasks
      if (distr1%task(it)<0) then  ! Empty slot
        distr1%task(it) = iTask
        found = .true.
        exit ! it loop
      end if
    end do ! it
    if (.not.found) call die(errHead//'parameter maxDistrTasks too small')
  end if ! (.not.any(distr1%task==iTask))

! Store association to second distribution, if different from first one
  if (present(distrID2) .and. iDistr1/=iDistr2) then
    task%distr(2) = iDistr2
    if (.not.any(distr2%task==iTask)) then
      found = .false.
      do it = 1,maxDistrTasks
        if (distr2%task(it)<0) then  ! Empty slot
          distr2%task(it) = iTask
          found = .true.
          exit ! it loop
        end if
      end do ! it
      if (.not.found) call die(errHead//'parameter maxDistrTasks too small')
    end if ! (.not.any(distr2%task==iTask))
  end if ! (present(distrID2))

! Mark task as already associated
  task%associated = .true.

#ifdef DEBUG_XC
!  if (present(distrID2)) then
!    write(udebug,'(a,2i4,2(2x,2i4))') &
!      myName//'taskID,iTask,distrID,iDistr=', &
!      taskID, iTask, distrID1, iDistr1, distrID2, iDistr2
!    write(udebug,'(a,i3,2x,25i3,/,(36x,25i3))') myName//'iDistr,tasks=', &
!      iDistr1, pack(distr1%task,distr1%task>0)
!    write(udebug,'(a,i3,2x,25i3,/,(36x,25i3))') myName//'iDistr,tasks=', &
!      iDistr2, pack(distr2%task,distr2%task>0)
!  else
!    write(udebug,'(a,2i4,2x,2i4)') &
!      myName//'taskID,iTask,distrID,iDistr=', &
!      taskID, iTask, distrID1, iDistr1
!    write(udebug,'(a,i3,2x,25i3,/,(36x,25i3))') myName//'iDistr,tasks=', &
!      iDistr1, pack(distr1%task,distr1%task>0)
!  end if
#endif /* DEBUG_XC */

end subroutine associateMeshTask

!-------------------------------------------------------------------------------

subroutine commonBox( nMesh, aBox, bBox, aComBox, bComBox, sizeSum, nParts )

! Finds the common intersection(s) of two mesh boxes, not necessarily contained
! in first unit cell (FUC). Each box is folded in parts contained in FUC.
! Then the intersection is obtained between each part or the folded aBox and
! each part of the folded bBox and returned in nParts intersection parts.

  implicit none
  integer,intent(in) :: nMesh(3)       ! Mesh divisions in each axis
  integer,intent(in) :: aBox(2,3)      ! First box, not necessarily within FUC
  integer,intent(in) :: bBox(2,3)      ! Second box, not necessarily within FUC
  integer,intent(out):: aComBox(:,:,:) ! Partial intersection boxes, each one
                                       ! contained in a single periodic cell, 
                                       ! relative to the aBox origin
  integer,intent(out):: bComBox(:,:,:) ! Same partial intersection boxes, but 
                                       ! relative to the bBox origin
  integer,intent(out):: sizeSum(0:)    ! Sum of the sizes of the partial boxes
  integer,intent(out):: nParts         ! Number of partial boxes

  character(len=*),parameter:: myName = moduleName//'commonBox '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: aPart, aParts, aPartBox(2,3,maxParts), aPartBoxFUC(2,3,maxParts), &
            bPart, bParts, bPartBox(2,3,maxParts), bPartBoxFUC(2,3,maxParts), &
            comBox(2,3), comSize, mxParts

! Find maximum number of parts allowed by size of output arrays
  mxParts = min( size(aComBox,3), size(bComBox,3), size(sizeSum)-1 )

! Divide aBox and bBox in parts contained in a single unit cell, as well as
! the equivalent part boxes in first unit cell (FUC)
  call unitCellParts( nMesh, aBox, aPartBox, aPartBoxFUC, aParts )
  call unitCellParts( nMesh, bBox, bPartBox, bPartBoxFUC, bParts )

! Find intersection boxes common to aBox and bBox
#ifdef DEBUG_XC
!  write(udebug,'(a,3(2i4,2x))') myName//'       aBox=', aBox
!  write(udebug,'(a,3(2i4,2x))') myName//'       bBox=', bBox
#endif /* DEBUG_XC */
  nParts = 0
  sizeSum(0) = 0
  do bPart = 1,bParts
  do aPart = 1,aParts
    comBox(1,:) = max( aPartBoxFUC(1,:,aPart), bPartBoxFUC(1,:,bPart) )
    comBox(2,:) = min( aPartBoxFUC(2,:,aPart), bPartBoxFUC(2,:,bPart) )
    if (all(comBox(1,:)<=comBox(2,:))) then ! Nonempty intersection box
      nParts = nParts + 1
      if (nParts > mxParts) &
        call die(errHead//'size of aComBox, bComBox, or sizeSum too small')
      ! Sum size of common box
      comSize = product(comBox(2,:)-comBox(1,:)+1)
      sizeSum(nParts) = sizeSum(nParts-1) + comSize
      ! Find common box in the unfolded unit cell of aPartBox
      aComBox(:,:,nParts) = comBox(:,:) + &
                            aPartBox(:,:,aPart) - aPartBoxFUC(:,:,aPart)
      ! Find common box in aBox origin. Index=1 in srcBox in both lines!!!
      aComBox(1,:,nParts) = aComBox(1,:,nParts) - aBox(1,:)
      aComBox(2,:,nParts) = aComBox(2,:,nParts) - aBox(1,:)
      ! Find common box in bBox origin
      bComBox(:,:,nParts) = comBox(:,:) + &
                            bPartBox(:,:,bPart) - bPartBoxFUC(:,:,bPart)
      bComBox(1,:,nParts) = bComBox(1,:,nParts) - bBox(1,:)
      bComBox(2,:,nParts) = bComBox(2,:,nParts) - bBox(1,:)
#ifdef DEBUG_XC
!      write(udebug,'(a,3(2i4,2x))') myName//'   aPartBox=', &
!                                                   aPartBox(:,:,aPart)
!      write(udebug,'(a,3(2i4,2x))') myName//'     comBox=', &
!        comBox(:,:) - aPartBox(:,:,aPart) + aPartBoxFUC(:,:,aPart)
!      write(udebug,'(a,3(2i4,2x))') myName//'    aComBox=', &
!                                                    aComBox(:,:,nParts)
!      write(udebug,'(a,3(2i4,2x))') myName//'    bComBox=', &
!                                                    bComBox(:,:,nParts)
#endif /* DEBUG_XC */
    end if
  end do ! aPart
  end do ! bPart

end subroutine commonBox

!-------------------------------------------------------------------------------

subroutine copyMeshData( nMesh, srcDistr, srcData, dstBox, dstData, task )

! Copies a box of values defined in a 3D mesh of a periodic unit cell.

  implicit none
  integer, intent(in) :: nMesh(3)          ! Global mesh size in each axis
  integer, intent(in) :: srcDistr          ! Distrib. ID of source data
                                           ! or zero for serial execution
  real(gp),intent(in) :: srcData(0:,0:,0:,:) ! Source data values to be copied
                                           ! Its upper array bounds must be
                                           ! enough to contain the node box in
                                           ! the srcDistr parallel distrib.
                                           ! or nMesh(iAxis)-1 in serial exe
  integer, intent(in) :: dstBox(2,3)       ! Destination mesh box:
                                           ! box(1,iAxis)=lower box limits
                                           ! box(2,iAxis)=upper box limits
                                           ! Box may be larger than the
                                           ! (0:nMesh-1) unit-cell range
  real(gp),intent(out):: dstData(0:,0:,0:,:) ! Destination data values.
                                           ! The upper array bounds must be 
                                           ! dstBox(2,iAxis)-dst(1,iAxis)
                                           ! srcData and dstData must be 
                                           ! different arrays
  integer,intent(inout),optional:: task    ! ID of communication task

  integer :: srcBox(2,3)

! Associate task to distribution
  if (present(task)) call associateMeshTask( task, srcDistr )

! Find box limits of source data in this node
  call myMeshBox( nMesh, srcDistr, srcBox )

! Check that srcData and dstData are not identical already
!  Avoid this in parallel, since different nodes might do different things
!  if (all(srcBox==dstBox) .and. all(shape(srcData)==shape(dstData))) then
!    if (all(srcData==dstData)) return
!  end if

! Copy data
  dstData = 0
  call reduceData( nMesh, srcBox, srcData, dstBox, dstData, taskID=task )

end subroutine copyMeshData

!-------------------------------------------------------------------------------

subroutine divideBox1D( box, nParts, partBox, blockSize, workload )

! Divides a 1D box according to workload (or uniformly if workload is absent)

  implicit none
  integer,          intent(in) :: box(2)
  integer,          intent(in) :: nParts
!
! Avoid the creation of a temporary array by passing descriptor info
!!!  integer,          intent(out):: partBox(2,nParts)
  integer,          intent(out):: partBox(:,:)  
!
  integer, optional,intent(in) :: blockSize
  real(gp),optional,intent(in) :: workload(0:)

  character(len=*),parameter:: myName = moduleName//'divideBox1D '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  logical :: nextPart
  integer :: blckSize, boxSize, i, iPart, largerSize, last, &
             nBlocks, nLargerParts, partSize
  real(dp):: partWkld, wlSum

! Trap a trivial case
  if (nParts==1) then
    partBox(:,1) = box
    return
  end if

! Find box size
  boxSize = box(2) - box(1) + 1

! Set block size
  if (present(blockSize)) then
    if (mod(boxSize,blockSize)/=0) &
      call die(errHead//'boxSize not a multiple of blockSize')
    blckSize = blockSize
  else
    blckSize = 1
  end if

! Choose division mode
  if (present(workload)) then ! Divide box according to workload

    ! Check size of workload array
    if (size(workload)/=boxSize) call die(errHead//'size(workload)/=boxSize')

    ! Find intended workload of each part
    wlSum = sum(workload)
    partWkld = wlSum / nParts

    ! Loop on parts
    partBox(1,1) = box(1)
    do iPart = 1,nParts-1
      last = 0  ! Last point with nonzero workload (or zero if none yet)
      do i = 0,boxSize-1
        if ( sum(workload(0:i)) > iPart*partWkld ) then  ! Limit surpassed
          if (i==0 .or. last < i-1) then   ! Place limit at center of vacuum
            partBox(2,iPart) = box(1) + (last+i)/2
          else if ( i==boxSize-1 .or. &
                    abs(sum(workload(0:i-1))-iPart*partWkld) &
                  < abs(sum(workload(0:i))  -iPart*partWkld) ) then
            partBox(2,iPart) = box(1) + i-1
          else
            partBox(2,iPart) = box(1) + i
          end if
          partBox(1,iPart+1) = partBox(2,iPart) + 1
          ! Avoid dividing blocks: parts must begin in a multiple of blockSize
          partBox(1,iPart+1) = partBox(1,iPart+1) &
                             - mod( partBox(1,iPart+1), blckSize )
          ! Prevent that iPart and iPart+1 overlap
          partBox(1,iPart+1) = max( partBox(1,iPart+1), &
                                    partBox(1,iPart)+blckSize )
          partBox(2,iPart) = partBox(1,iPart+1) - 1
          ! Go to next part
          nextPart = .true.
        else ! Limit not surpassed
          nextPart = .false.
        end if ! ( sum(workload(axis,0:i)) > iPart*partWkld )
        if (workload(i)>0._gp) last = i
        if (nextPart) exit ! i loop
      end do ! i
    end do ! iPart
    partBox(2,nParts) = box(2)

  else ! (.not.present(workload)) => Divide box uniformly

    nBlocks = boxSize / blckSize              ! Number of blocks of points
    partSize = nBlocks / nParts               ! Blocks per part
    nLargerParts = nBlocks - partSize*nParts  ! Number of larger parts
    partSize = partSize * blckSize            ! Points per part
    largerSize = partSize + blckSize ! One more block for larger parts
    partBox(1,1) = box(1)
    do iPart = 1,nLargerParts        ! largerParts of size partSize+blockSize
      if (iPart>1) partBox(1,iPart) = partBox(2,iPart-1) + 1
      partBox(2,iPart) = partBox(1,iPart) + largerSize - 1
    end do
    do iPart = nLargerParts+1,nParts  ! remaining parts of size partSize
      if (iPart>1) partBox(1,iPart) = partBox(2,iPart-1) + 1
      partBox(2,iPart) = partBox(1,iPart) + partSize - 1
    end do

  end if ! (present(workload))

end subroutine divideBox1D

!-------------------------------------------------------------------------------

subroutine divideBox3D( nMesh, wlDistr, workload, box, &
                        nParts, blockSize, axis, partBox )

! Divides a mesh box according to a given workload and axis. If axis<1 on
! input, it is chosen to maximize the spatial dispersion of workload

  implicit none
  integer, intent(in)   :: nMesh(3)  ! Global mesh size in each axis
  integer, intent(in)   :: wlDistr   ! Mesh distribution ID of workload array
  real(gp),intent(in)   :: workload(0:,0:,0:) ! Work load of each mesh point
  integer, intent(in)   :: box(2,3)  ! Mesh box to be divided
  integer, intent(in)   :: nParts    ! Number of parts of the divided box
  integer, intent(in)   :: blockSize ! Size of blocks of mesh points in each 
                                     ! axis that cannot be split among nodes
  integer, intent(inout):: axis      ! Axis along which division is made
  integer, intent(out)  :: partBox(2,3,nParts) ! Divided part boxes

  character(len=*),parameter:: myName = moduleName//'divideBox3D '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer :: boxShape(3), i, iPart, ix, maxShape
  real(dp):: maxDisp, wlDisp, wlSum, wlSum1, wlSum2
!  logical :: nextPart
!  integer :: largerSize, last, nBlocks, nLargerParts, partSize
!  real(dp):: partWkld
  real(gp),allocatable:: prjWkld(:,:)

! Initialize part boxes as total box
  forall (iPart=1:nParts) partBox(:,:,iPart) = box(:,:)

! Trap a trivial case
  if (nParts<2) return

! Check that nMesh and box are consistent with blockSize
  if (any(mod(nMesh,blockSize)/=0)) &
    call die(errHead//'nMesh and blockSize are not consistent')
  if (any(mod(box(1,:),blockSize)/=0) .or. &
      any(mod(box(2,:)+1,blockSize)/=0))   &
    call die(errHead//'box and blockSize are not consistent')

! Check that workload is nonnegative
  if (any(workload<0._gp)) call die(errHead//'negative workload')

! Allocate array for the projection of workload along each axis
  boxShape(:) = box(2,:) - box(1,:) + 1
  allocate( prjWkld(3,0:maxval(boxShape)-1) )

! Find projection of workload along each axis
  prjWkld = 0
  call projectMeshData( nMesh, wlDistr, workload, box, prjWkld )

#ifdef DEBUG_XC
!  write(udebug,*) myName//'projected workload:'
!  write(udebug,'(i4,3f15.6)') (i,prjWkld(:,i),i=0,maxval(boxShape)-1)
#endif /* DEBUG_XC */

! Find total workload in box
  wlSum = sum(prjWkld) / 3  ! Since the sum of each projection must be equal

! Select axis for division
  if (axis<1 .or. axis>3) then ! Otherwise, accept input axis 
    if (wlSum==0._dp) then  ! Choose longest axis
      maxShape = maxval(boxShape)
      do ix = 3,1,-1  ! If equivalent, take largest ix (z better than x)
        if (boxShape(ix) == maxShape) then
          axis = ix
          exit ! ix loop
        end if
      end do
    else ! Choose axis with maximum spatial dispersion of workload
      maxDisp = 0
      do ix = 3,1,-1
        wlSum1 = 0
        wlSum2 = 0
        do i = 0,boxShape(ix)-1
          wlSum1 = wlSum1 + prjWkld(ix,i) * i
          wlSum2 = wlSum2 + prjWkld(ix,i) * i**2
        end do
        wlDisp = wlSum2/wlSum - (wlSum1/wlSum)**2
        if (wlDisp>maxDisp) then
          axis = ix
          maxDisp = wlDisp
        end if
      end do ! ix
    end if ! (wlSum==0._gp)
  end if ! (axis<1 .or. axis>3)

! Divide box along axis
  if (wlSum==0._dp) then  ! Divide box uniformly
    call divideBox1D( box(:,axis), nParts, partBox(:,axis,:), blockSize )
  else                    ! Divide projected workload uniformly
    call divideBox1D( box(:,axis), nParts, partBox(:,axis,:), blockSize, &
                      prjWkld(axis,0:boxShape(axis)-1) )
  end if ! (wlSum==0._dp)

#ifdef DEBUG_XC
! Check that all box limits are consistent
  if ( any(box(1,axis) > partBox(1,axis,:)) .or. &
       any(partBox(2,axis,1:nParts-1) >= partBox(1,axis,2:nParts)) .or. &
       any(partBox(1,axis,:) > partBox(2,axis,:)) .or. &
       any(partBox(2,axis,:) > box(2,axis)) ) then
    write(udebug,'(a,3(2i6,2x))') errHead//'box=', box
    write(udebug,'(a,/,(3(2i6,2x)))') errHead//'partBox=', partBox
    call die(errHead//'inconsistent partBox limits')
  end if
#endif /* DEBUG_XC */

  deallocate( prjWkld )

end subroutine divideBox3D

!-------------------------------------------------------------------------------

subroutine fftMeshDistr( nMesh, fftDistr, axisDistr )

! Creates and returns a homogeneaous 3D mesh distribution fftDistr, and three
! 2D distributions axisDistr, whose boxes extend wholly over one of the three
! cell axes. These axisDistr allow fully local 1D FFTs in one of the axes.
! The fftDistr and axisDistr are designed to make already optimal (without any
! collisions) the default sequence of inter-node communications, returned by
! all2allTransferOrder. This is beacause all simultaneous communications occur
! between nodes separated by the same span (difference between node indexes).
! Thus, no task arguments are required to optimize such communications.

  implicit none
  integer,         intent(in)   :: nMesh(3)     ! Mesh divisions of each axis
  integer,         intent(inout):: fftDistr     ! 3D mesh distr. used in 3D FFT
  integer,optional,intent(inout):: axisDistr(3) ! 2D distr. used for the 1D FFTs

  character(len=*),parameter:: myName = 'fftMeshDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: axis, axis1, axis2, axis3, axisNodes(3), box0(2,3), boxNodes(3), &
            iDistr, iNode, iNode1, iNode2, iNode3, jNode2, jNode3, &
            node0, nodeSpan(3), oldDistr, rowMesh(3), rowNodes
  integer,allocatable:: subBox(:,:,:)
  logical:: found
  type(distrType),pointer:: distr

! Trap the serial case (totNodes available from module parallel)
  if (totNodes==1) then
    fftDistr = 0
    if (present(axisDistr)) axisDistr(:) = 0
    return
  end if

! Check if the input FFT distribution IDs are already defined as such
  iDistr = indexDistr( fftDistr )
  if (iDistr>0 .and. iDistr<=maxDistr) then
    distr => storedMeshDistr(iDistr)
    if (distr%defined .and. all(distr%nMesh==nMesh)) then
      found = .true.  ! But now check the axis distributions
      if (present(axisDistr)) then
        do axis = 1,3
          iDistr = indexDistr( axisDistr(axis) )
          if (iDistr>0 .and. iDistr<=maxDistr) then
            distr => storedMeshDistr(iDistr)
            if (.not.distr%defined .or. any(distr%nMesh/=nMesh)) &
              found = .false.
          else
            found = .false.
          end if
        end do ! axis
      end if ! (present(axisDistr))
      if (found) return ! Since the input values are still valid
    end if ! (distr%defined .and. ...)
  end if ! (iDistr>0)

! Find optimal distribution of nodes on axes
  call optimizeNodeDistr( nMesh, totNodes, axisNodes )

! Create homogeneous 3D distribution
#ifdef DEBUG_XC
!  write(udebug,*) myName//'calling setMeshDistr with fftDistr'
#endif /* DEBUG_XC */
  call setMeshDistr( fftDistr, nMesh, nNodesX=axisNodes(1), &
                     nNodesY=axisNodes(2), nNodesZ=axisNodes(3) )

! Return already if axis distributions are not required
  if (.not.present(axisDistr)) return

! Find span between nodes along each axis
!  nodeSpan(1) = 1
!  nodeSpan(2) = axisNodes(1)
!  nodeSpan(3) = axisNodes(1) * axisNodes(2)
! Third axis is innermost for node distribution (see setMeshDistr)
  nodeSpan(3) = 1
  nodeSpan(2) = axisNodes(3)
  nodeSpan(1) = axisNodes(3) * axisNodes(2)

#ifdef DEBUG_XC
!    write(udebug,'(a,3i4)') myName//'axisNodes=', axisNodes
!    write(udebug,'(a,3i4)') myName//'nodeSpan=', nodeSpan
#endif /* DEBUG_XC */

! Allocate a small temporary array
  allocate( subBox(2,0:maxval(axisNodes)-1,3) )

! Loop on the three cell axes
  do axis1 = 1,3

    ! (Re)initialize the axis distribution for the 1D FFT transforms
    oldDistr = axisDistr(axis1)
    call freeMeshDistr( axisDistr(axis1) )
    call initDistr( axisDistr(axis1), nMesh, 0, totNodes )

    ! Find the two other axes
    axis2 = modulo(axis1,3) + 1
    axis3 = modulo(axis2,3) + 1

    ! Point towards the apropriate axis distribution
    iDistr = indexDistr( axisDistr(axis1) )
    distr => storedMeshDistr(iDistr)

    ! Loop on the node boxes over the two other axes
    do iNode3 = 0,axisNodes(axis3)-1
    do iNode2 = 0,axisNodes(axis2)-1

      ! Find the first node in the row of nodes along axis1
      node0 = nodeSpan(axis2)*iNode2 + nodeSpan(axis3)*iNode3

      ! Find the box limits along axis2 and axis3
      call nodeMeshBox( nMesh, fftDistr, node0, box0 )

      ! Find the total box of a row of nodes along axis1
      rowNodes = axisNodes(axis1)
      rowMesh(1) = nMesh(axis1)
      rowMesh(2) = box0(2,axis2) - box0(1,axis2) + 1
      rowMesh(3) = box0(2,axis3) - box0(1,axis3) + 1

      ! Find an optimal (re)distribution of rowNodes along axis2 and/or axis3
      boxNodes(1) = 1
      call optimizeNodeDistr( rowMesh(2:3), rowNodes, boxNodes(2:3) )

      ! Divide the row box along axis2 and/or axis3
      subBox(1,0,1) = 0
      subBox(2,0,1) = nMesh(axis1) - 1
      call divideBox1D( box0(:,axis2), boxNodes(2), subBox(:,0:boxNodes(2)-1,2))
      call divideBox1D( box0(:,axis3), boxNodes(3), subBox(:,0:boxNodes(3)-1,3))

      ! Set the new boxes
      iNode = node0
      do jNode3 = 0,boxNodes(3)-1
      do jNode2 = 0,boxNodes(2)-1
        distr%box(1,axis1,iNode) = 0
        distr%box(2,axis1,iNode) = nMesh(axis1)-1
        distr%box(:,axis2,iNode) = subBox(:,jNode2,2)
        distr%box(:,axis3,iNode) = subBox(:,jNode3,3)
        iNode = iNode + nodeSpan(axis1)
      end do ! jNode2
      end do ! jNode3

#ifdef DEBUG_XC
!      write(udebug,'(a,2i4,3(2x,2i4))') &
!        myName//'axis1,node0,box0=', axis1, node0, box0
!      write(udebug,'(a,3i4)') myName//'boxNodes=', boxNodes
!      do jNode2 = 0,boxNodes(2)-1
!        write(udebug,'(a,i4,2x,3i4)') &
!          myName//'axis2,subBox=', axis2, subBox(:,jNode2,2)
!      end do
!      do jNode3 = 0,boxNodes(3)-1
!        write(udebug,'(a,i4,2x,3i4)') &
!          myName//'axis3,subBox=', axis3, subBox(:,jNode3,3)
!      end do
#endif /* DEBUG_XC */

    end do ! iNode2
    end do ! iNode3

    ! Check if this distribution was already defined
    call reduceDistr( axisDistr(axis1) )

#ifdef DEBUG_XC
    iDistr = indexDistr( axisDistr(axis1) )
    distr => storedMeshDistr(iDistr)
!    write(udebug,'(a,3i4)') myName//'axis1,axis2,axis3=', axis1, axis2, axis3
    write(udebug,'(a,6x,3i4,3(2x,2i4))') &
      myName//'    old/newDistrID,iDistr,myBox=', &
      oldDistr, axisDistr(axis1), iDistr, distr%box(:,:,myNode)
#endif /* DEBUG_XC */

  end do ! axis1

! Deallocate temporary array
  deallocate( subBox )

end subroutine fftMeshDistr

!-------------------------------------------------------------------------------

subroutine freeMeshDistr( distrID )

! Erases a mesh distribution ID and frees the distribution slot if 
! no other IDs are assigned to it

  implicit none
  integer,intent(in):: distrID  ! ID of an allocated distribution

  character(len=*),parameter:: myName = 'freeMeshDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: iDistr, iID, iNode, it, iTask, taskID
  logical:: found
  type(distrType),pointer :: distr
  type(taskType),pointer :: task

! Find internal distribution index
  iDistr = indexDistr( distrID )

! Check that distribution exists
  if (iDistr<1 .or. iDistr>maxDistr) return
  distr => storedMeshDistr(iDistr) ! Just a shorter name

! Erase ID from the distribution
  do iID = 1,maxDistrID
    if (distr%ID(iID)==distrID) then
      distr%ID(iID) = -1
      exit ! iID loop
    end if
  end do ! iID

! Free distribution if no other IDs are assigned to it
  if (all(distr%ID<0)) then
    ! First free all distribution tasks
#ifdef DEBUG_XC
!    if (any(distr%task>0)) &
!      write(udebug,'(a,i3,2x,25i3,/,(32x,25i3))') &
!        myName//'iDistr,tasks=', iDistr, pack(distr%task,distr%task>0)
#endif /* DEBUG_XC */
    do it = 1,maxDistrTasks
      iTask = distr%task(it)
      if (iTask>0) then
        task => storedMeshTask(iTask)
        found = .false.
        do iID = 1,maxTaskID  ! Find a valid ID of task, to call freeMeshTask
          taskID = task%ID(iID)
          if (taskID>0) then
            found = .true.
            call freeMeshTask( taskID )
            exit ! iID loop
          end if ! (taskID>0)
        end do ! iID
        if (.not.found) call die(myName//'ERROR: no valid task ID found')
      end if ! (iTask>0)
    end do ! it
    ! Finally free distribution itself
    deallocate( distr%box )
    distr%defined = .false.
  end if

end subroutine freeMeshDistr

!-------------------------------------------------------------------------------

subroutine freeMeshTask( taskID )

! Erases a mesh task ID and frees the task slot if 
! no other IDs are assigned to it

  implicit none
  integer,intent(in):: taskID  ! ID of an allocated task

  character(len=*),parameter:: myName = 'freeMeshTask '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: distrID, i, iDistr, it, iTask, iID
  logical:: found
  type(distrType),pointer :: distr
  type(taskType),pointer :: task

! Find internal task index
  iTask = indexTask( taskID )

! Check that task exists
  if (iTask<1 .or. iTask>maxTasks) return
  task => storedMeshTask(iTask) ! Just a shorter name

#ifdef DEBUG_XC
!  write(udebug,'(a,2i4,2x,2i4)') &
!    myName//'taskID,iTask,task%distr=', taskID, iTask, task%distr
#endif /* DEBUG_XC */

! Erase ID from the task ID list
  do iID = 1,maxTaskID
    if (task%ID(iID)==taskID) then
      task%ID(iID) = -1
      exit ! iID loop
    end if
  end do ! iID

! Free task if no other IDs are assigned to it
  if (all(task%ID<0)) then
    ! First eliminate task from all its associated distributions
    do i = 1,2
      iDistr = task%distr(i)
      if (iDistr>0) then
        distr => storedMeshDistr(iDistr)
        found = .false.
        do it = 1,maxDistrTasks
          if (distr%task(it)==iTask) then
            distr%task(it) = -1
            found = .true.
          end if
        end do ! it
        if (.not.found) call die(errHead//'task-distr association not found')
      end if ! (iDistr>0)
    end do ! i
    ! Finally free task itself
    deallocate( task%trsfDir, task%trsfNode, task%dstBox, task%srcBox )
    task%associated = .false.
    task%optimized = .false.
    task%defined = .false.
    task%distr = -1
    task%ID = -1
  end if

end subroutine freeMeshTask

!-------------------------------------------------------------------------------

subroutine gatherBoxes( box, boxes )

! Gathers the mesh boxes of all processors

  implicit none
  integer,intent(in) :: box(2,3)                 ! My node's mesh box
  integer,intent(out):: boxes(2,3,0:totNodes-1)  ! All node's boxes

#ifdef MPI
! Gather the boxes of all nodes
  call MPI_AllGather( box(1,1) , 6, MPI_Integer, &
                      boxes(1,1,0), 6, MPI_Integer, MPI_COMM_WORLD, MPIerror )
#else
! Copy the only node's box
  boxes(:,:,0) = box(:,:)
#endif

end subroutine gatherBoxes

!-------------------------------------------------------------------------------

integer function indexDistr( distrID )

! Returns the internal memory position in which a distribution is stored or
!   0 if distrID==0
!  -1 if there is no allocated distribution with that ID

  implicit none
  integer,intent(in):: distrID ! Distribution ID

  integer:: iDistr

  indexDistr = -1

  if (distrID==0) then
    indexDistr = 0
  else if (distrID>0) then
    do iDistr = 1,maxDistr
      if ( storedMeshDistr(iDistr)%defined .and. &
           any(storedMeshDistr(iDistr)%ID==distrID) ) then
        indexDistr = iDistr
        return
      end if
    end do ! iDistr
  end if

end function indexDistr

!-------------------------------------------------------------------------------

integer function indexTask( taskID )

! Returns the internal memory position in which a task is stored or
!  -1 if there is no allocated task with that ID

  implicit none
  integer,intent(in):: taskID ! Task ID

  integer:: iTask

  indexTask = -1

  if (taskID>0) then
    do iTask = 1,maxTasks
      if ( any(storedMeshTask(iTask)%ID==taskID) ) then
        indexTask = iTask
        return
      end if
    end do ! iTask
  end if

end function indexTask

!-------------------------------------------------------------------------------

subroutine initDistr( distrID, nMesh, firstNode, nNodes )

! Initializes a parallel distribution of mesh points among processor nodes.
! In serial execution (totNodes==1) it simply returns distrID=0.

  implicit none
  integer,intent(out):: distrID   ! ID assigned to the mesh distrib.
  integer,intent(in) :: nMesh(3)  ! Mesh divisions in each axis
                                  ! (including "subpoints")
  integer,intent(in) :: firstNode ! First node in the mesh distr.
  integer,intent(in) :: nNodes    ! Total nodes in the mesh distr
                                  ! If not present all nodes are used

  character(len=*),parameter:: myName = moduleName//'initDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  type(distrType),pointer:: distr
  integer:: iDistr, iID
  logical:: found

! Trap the serial case (totNodes available from module parallel)
  if (totNodes==1) then
    distrID = 0
    return
  end if

! Find an available internal distribution slot
  found = .false.
  do iDistr = 1,maxDistr
    if (all(storedMeshDistr(iDistr)%ID<0)) then  ! Available (empty) slot
      found = .true.
      exit ! iDistr loop
    end if
  end do ! iDistr
  if (.not.found) call die(errHead//'parameter maxDistr too small')

! Increase private distribution ID counter and assign it to new distribution
  nDistrID = nDistrID + 1
  distrID = nDistrID

! Set distribution parameters, except box
  distr => storedMeshDistr(iDistr)  ! Just a shorter name
  distr%defined = .true.
  distr%firstNode = firstNode
  distr%nNodes = nNodes
  distr%nMesh = nMesh
  distr%task = -1
  distr%ID(:) = -1
  distr%ID(1) = distrID

! Allocate box array
  allocate( distr%box(2,3,0:totNodes-1) )

! Initialize node boxes so that lower bound > upper bound (null box)
  distr%box(1,:,:) = 0
  distr%box(2,:,:) = -1

#ifdef DEBUG_XC
! write(udebug,'(a,2i6)') myName//'distrID,iDistr=', distrID, iDistr
#endif /* DEBUG_XC */

end subroutine initDistr

!-------------------------------------------------------------------------------

subroutine initTask( taskID )

! Initializes the ID of a communication task of mesh data

  implicit none
  integer,intent(out):: taskID

  character(len=*),parameter:: myName = moduleName//'initTask '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  logical:: found
  integer:: iTask
  type(taskType),pointer:: task

! Trap the serial case (totNodes available from module parallel)
  if (totNodes==1) then
    taskID = 0
    return
  end if

! Find an available internal distribution slot
  found = .false.
  do iTask = 1,maxTasks
    if (all(storedMeshTask(iTask)%ID<0)) then ! Available (empty) slot
      found = .true.
      exit ! iTask loop
    end if
  end do ! iTask
  if (.not.found) call die(errHead//'parameter maxTasks too small')

! Increase private distribution ID counter and assign it to new distribution
  nTaskID = nTaskID + 1
  taskID = nTaskID

! Set task parameters, except boxes
  task => storedMeshTask(iTask)  ! Just a shorter name
  task%defined = .false.
  task%optimized = .false.
  task%associated = .false.
  task%distr = -1
  task%nTrsf = 0
  task%ID = -1
  task%ID(1) = taskID

! Allocate transfer sequence (a send and a receive to/from each other node)
  allocate( task%trsfNode(0:2*(totNodes-1)) )
  allocate( task%trsfDir(0:2*(totNodes-1)) )

! Allocate box arrays
  allocate( task%srcBox(2,3,0:totNodes-1) )
  allocate( task%dstBox(2,3,0:totNodes-1) )

! Initialize boxes so that lower bound > upper bound (null box)
  task%srcBox(1,:,:) = 0
  task%srcBox(2,:,:) = -1
  task%dstBox(1,:,:) = 0
  task%dstBox(2,:,:) = -1

#ifdef DEBUG_XC
! write(udebug,'(a,2i6)') myName//'taskID, iTask=', taskID, iTask
#endif /* DEBUG_XC */

end subroutine initTask

!-------------------------------------------------------------------------------

subroutine myMeshBox( nMesh, distrID, box )

! Finds the mesh box belonging to a node in a parallel mesh distribution,
! If iDistr=0, it returns box(1,:)=0 and  box(2,:)=nMesh(:)-1 
! If node stores no mesh points, its (empty) box returns with box(1,:)>box(2,:)
! Note: although nMesh is available from the distr. data, it is also included
! in the routine arguments, to handle directly the case iDistr=0.

  implicit none
  integer, intent(in) :: nMesh(3)  ! Mesh divisions in each axis
  integer, intent(in) :: distrID   ! Mesh-distribution ID
  integer, intent(out):: box(2,3)  ! Mesh box: box(1,:) = lower bounds
                                   !           box(2,:) = upper bounds

  call nodeMeshBox( nMesh, distrID, myNode, box )

end subroutine myMeshBox

!-------------------------------------------------------------------------------

subroutine nodeMeshBox( nMesh, distrID, node, box )

! Finds the mesh box belonging to a node in a parallel mesh distribution,
! If iDistr=0, it returns box(1,:)=0 and  box(2,:)=nMesh(:)-1 
! If node stores no mesh points, its (empty) box returns with box(1,:)>box(2,:)
! Note: although nMesh is available from the distr. data, it is also included
! in the routine arguments, to handle directly the case iDistr=0.

  implicit none
  integer,intent(in) :: nMesh(3)  ! Mesh divisions in each axis
  integer,intent(in) :: distrID   ! Mesh-distribution ID
  integer,intent(in) :: node      ! Processor node whose box is needed
                                  ! If not present, node=myNode
  integer,intent(out):: box(2,3)  ! Mesh box: box(1,:) = lower bounds
                                  !           box(2,:) = upper bounds

  character(len=*),parameter:: myName = 'nodeMeshBox '
  character(len=*),parameter:: errHead = myName//'ERROR: '

  integer:: iDistr

  if (distrID==0) then ! Data not distributed
    box(1,:) = 0
    box(2,:) = nMesh(:)-1
  else
    iDistr = indexDistr( distrID )
    if (iDistr>0 .and. iDistr<=maxDistr) then
      if (any(nMesh/=storedMeshDistr(iDistr)%nMesh)) &
        call die(errHead//'nMesh/=distr%nMesh')
      box(:,:) = storedMeshDistr(iDistr)%box(:,:,node)
    else
      call die(errHead//'undefined mesh distribution')
    end if
  end if

end subroutine nodeMeshBox

!-------------------------------------------------------------------------------

subroutine optimizeNodeDistr( nMesh, nNodes, axisNodes )

! Finds a reasonably optimal distribution of nodes over the different axes.
! Optimal here means that it leads to most cubic mesh boxes for each node.

  implicit none
  integer,intent(in) :: nMesh(:)     ! Mesh points in each axis
  integer,intent(in) :: nNodes       ! Total number of processor nodes
  integer,intent(out):: axisNodes(:) ! Number of nodes in each axis

  character(len=*),parameter:: myName = moduleName//'optimizeNodeDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer,parameter:: maxFac = 1000
  integer :: fac(maxFac), i2, i3, i5, iFac1, iFac2, iRem, &
             n, n2, n3, n5, nAxes, nFac, nod1, nod2, nod3, nRem, remFac
  real(dp):: box1, box2, box3, boxAvg, boxErr, minErr

! Find number of axes and handle the trivial single-axis case
  nAxes = size(nMesh)
  if (nAxes<1) then
    return
  else if (nAxes==1) then
    axisNodes(1) = nNodes
    return
  else if (nAxes>3) then
    call die(errHead//'not prepared for nAxes>3')
  end if

! Find prime factors of nNodes (only 2, 3, and 5 coded presently) so that
!   nNodes = 2**n2 * 3**n3 * 5**n5 * remFac**nRem
  n = nNodes
  n2 = 0
  n3 = 0
  n5 = 0
  do
    if (mod(n,2)==0) then
      n2 = n2 + 1
      n = n / 2
    else if (mod(n,3)==0) then
      n3 = n3 + 1
      n = n / 3
    else if (mod(n,5)==0) then
      n5 = n5 + 1
      n = n / 5
    else
      exit ! do loop
    end if
  end do
  remFac = n
  if (remFac==1) then
    nRem = 0
  else
    nRem = 1
  end if

! Find all possible factors of nNodes
  nFac = 0
  do i2 = 0,n2
    do i3 = 0,n3
      do i5 = 0,n5
        do iRem = 0,nRem
          nFac = nFac + 1
          if (nFac>maxFac) call die(errHead//'parameter maxFac too small')
          fac(nFac) = 2**i2 * 3**i3 * 5**i5 * remFac**iRem
        end do ! iRem
      end do ! i5
    end do ! i3
  end do ! i2

! Try all possible combinations of factors on each axis and select that with
! most cubic boxes
  minErr = huge(minErr)
  do iFac1 = 1,nFac
    nod1 = fac(iFac1)
    box1 = nMesh(1) / real(nod1,dp)
    if (nAxes==2) then
      nod2 = nNodes / nod1
      box2 = nMesh(2) / real(nod2,dp)
      boxAvg = (box1 + box2) / 2
      boxErr = log(box1/boxAvg)**2 + log(box2/boxAvg)**2
      if (boxErr<minErr) then
        axisNodes(1) = nod1
        axisNodes(2) = nod2
        minErr = boxErr
      end if
    else ! (nAxes==3)
      do iFac2 = 1,nFac
        nod2 = fac(iFac2)
        if (mod(nNodes,nod1*nod2) /= 0) cycle ! iFac2 loop
        nod3 = nNodes / (nod1*nod2)
        box2 = nMesh(2) / real(nod2,dp)
        box3 = nMesh(3) / real(nod3,dp)
        boxAvg = (box1 + box2 + box3) / 3
        boxErr = log(box1/boxAvg)**2 + log(box2/boxAvg)**2 + log(box3/boxAvg)**2
        if (boxErr<minErr) then
          axisNodes(1) = nod1
          axisNodes(2) = nod2
          axisNodes(3) = nod3
          minErr = boxErr
        end if
      end do ! iFac2
    end if ! (nAxes==2)
  end do ! iFac1

end subroutine optimizeNodeDistr

!-------------------------------------------------------------------------------

subroutine optimizeTransferOrder( nNodes, node, nTrsf, trsfNode, trsfDir )

! Searches an optimal order of transfer communications


  implicit none
  integer,intent(in) :: nNodes   ! Number of parallel processor nodes. The
                                 ! number of transfers per node is 2*nNodes-1
  integer,intent(in) :: node     ! Local processor node, in range (0:nNodes-1)
  integer,intent(in) :: nTrsf    ! Number of required transfers
  integer,intent(inout):: trsfNode(nTrsf) ! To/from transfer node.
                                 ! First 'transfer' is local: trsfNode(0)=node
  integer,intent(inout):: trsfDir(nTrsf)  ! Transfer direction:
                                 ! trsfDir=+1 => from node to trsfNode
                                 ! trsfDir=-1 => from trsfNode to node

  character(len=*),parameter:: myName = moduleName//'optimizeTransferOrder '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  real(dp),        parameter:: incrFact = 2.0_dp ! Memory increase factor >1
  integer:: dir12, i1, i2, iter, iTrsf, maxTime, maxTrsf, MPIerror, &
            myNode, node1, node2, nTrsfNode1, time, trsf12
  integer:: orderNode1(nNodes), orderNode2(2*nNodes)
  logical:: found
  real(dp):: nTrsfNode(2*nNodes)
  integer,pointer:: inTrsf(:,:)=>null(), &
                    myTrsf(:)=>null(), outTrsf(:,:)=>null()

! In serial execution, do nothing
#ifdef MPI

! Find the maximun number of transfers of any node
  call MPI_AllReduce( nTrsf, maxTrsf, 1, MPI_Integer, &
                      MPI_Max, MPI_Comm_World, MPIerror )

! Trap a trivial case, in which there is nothing to optimize
  if (maxTrsf < 2) return

#ifdef DEBUG_XC
!  write(udebug,'(a,2i6)') myName//'nTrsf,maxTrsf=', nTrsf, maxTrsf
!  if (maxTrsf > 2*nNodes) call die(errHead//'too many transfers per node')
#endif /* DEBUG_XC */

! Allocate arrays for the input transfer sequences of all nodes
  call re_alloc( myTrsf, 1,maxTrsf, name=myName//'myTrsf' )
  call re_alloc( inTrsf, 1,maxTrsf, 1,nNodes, name=myName//'inTrsf' )

! Copy transfer sequence of my node to a different format
! Absolute value stores transfer node. Sign stores transfer direction.
! Internally, we use range (1:nNodes) rather than (0:nNodes-1)
  myNode = node + 1
  myTrsf(1:nTrsf) = trsfDir(:)*(trsfNode(:)+1)

! Gather initial transfer sequence of all nodes
  call MPI_AllGather( myTrsf, maxTrsf, MPI_Integer, &
                      inTrsf(1,1), maxTrsf, MPI_Integer, MPI_COMM_WORLD, MPIerror )

! Find the number of transfers of each node
  nTrsfNode(1:nNodes) = count( inTrsf(:,1:nNodes)/=0, dim=1 )

! Find the order of nodes by increasing number of transfers
  call ordix( nTrsfNode, 1, nNodes, orderNode1 )

! Allocate array for the optimized transfer times of all nodes
  maxTime = maxTrsf * incrFact
  call re_alloc( outTrsf, 1,maxTime, 1,nNodes, name=myName//'outTrsf' )

! Schedule the transfers of each node
  do i1 = nNodes,1,-1  ! Give priority to busiest nodes
    node1 = orderNode1(i1)

    ! Find the number of transfers of node 1
    nTrsfNode1 = count( inTrsf(:,node1)/=0 )

    ! Find the number of transfers of each node2, that transfers to/from node1
    do i2 = 1,nTrsfNode1
      node2 = abs( inTrsf(i2,node1) )
      nTrsfNode(i2) = count( inTrsf(:,node2)/=0 )
    end do

    ! Order node2 by increasing number of transfers
    call ordix( nTrsfNode, 1, nTrsfNode1, orderNode2 )

    ! Assign a time to each unassigned transfer of node1
    do iTrsf = nTrsfNode1,1,-1  ! Give priority to busiest node2
      i2 = orderNode2(iTrsf)

      trsf12 = inTrsf(i2,node1)
      node2 = abs( trsf12 )
      dir12 = sign( 1, trsf12 )

      ! Check that this transfer is not already assigned
      if (any(outTrsf(:,node1)==trsf12)) cycle ! iTrsf loop

      ! Look for the first available transfer time and assign it
      found = .false.
      realloc_loop: do iter = 1,2  ! Loop for case that maxTime is too small
        do time = 1,maxTime
          if (outTrsf(time,node1)==0 .and. outTrsf(time,node2)==0) then
            outTrsf(time,node1) =  dir12*node2
            outTrsf(time,node2) = -dir12*node1
            found = .true.
            exit realloc_loop
          end if
        end do ! iTime
        ! If no available time found, increase maxTime
        maxTime = maxTime*incrFact + 1
        call re_alloc( outTrsf, 1,maxTime, 1,nNodes, &
                       name=myName//'outTrsf', copy=.true. )
      end do realloc_loop

#ifdef DEBUG_XC
      if (.not.found) call die(errHead//'parameter incrFact too small')
#endif /* DEBUG_XC */
    end do ! iTrsf

    ! Copy transfer sequence to output arrays
    if (node1==myNode) then
      ! Pack transfer sequence, removing idle times for myNode
      myTrsf(1:nTrsf) = pack( outTrsf(:,myNode), outTrsf(:,myNode)/=0 )
      ! The -1 is to change back from node range (1:nNodes) to (0:nNodes-1)
      trsfNode = abs( myTrsf(1:nTrsf) ) - 1
      trsfDir = sign( 1, myTrsf(1:nTrsf) )
#ifdef DEBUG_XC
!      write(udebug,'(a,20i4)') myName//' inTrsf=', inTrsf(1:nTrsf,myNode)
!      write(udebug,'(a,20i4)') myName//'outTrsf=', myTrsf(1:nTrsf)
#endif /* DEBUG_XC */
      exit ! i1 loop
    end if

  end do ! i1

! Deallocate arrays
  call de_alloc( outTrsf, name=myName//'outTrsf' )
  call de_alloc( inTrsf, name=myName//'inTrsf' )
  call de_alloc( myTrsf, name=myName//'myTrsf' )

#endif

end subroutine optimizeTransferOrder

!-------------------------------------------------------------------------------

subroutine primeFactors( n, factor, power, nFactors )

! Finds prime factors of n

  implicit none
  integer,intent(in) :: n          ! Number to be factorized
  integer,intent(out):: factor(:)  ! Different prime factors in ascending order
  integer,intent(out):: power(:)   ! Power of each prime factor
  integer,intent(out):: nFactors   ! Number of different factors found

  integer:: i, m

  nFactors = 0
  m = n
  do i = 2,n  ! This is crude, but simple!
    if (mod(m,i)==0) then  ! i is a factor o m
      nFactors = nFactors + 1
      if (nFactors>size(factor)) &
        call die(moduleName//'primeFactors: ERROR: size(factor) too small')
      factor(nFactors) = i
      power(nFactors) = 0
      do ! loop to find power
        if (mod(m,i)==0) then
          power(nFactors) = power(nFactors) + 1
          m = m / i
        else
          exit ! power loop
        end if ! (mod(m,i)==0)
      end do ! power loop
    end if ! (mod(m,i)==0)
    if (m==1) exit ! i loop
  end do ! i

end subroutine primeFactors

!-------------------------------------------------------------------------------

subroutine projection( nMesh, srcBox, srcData, prjBox, prjData )

! Projects a box of values, defined in a 3D mesh of a periodic unit cell,
! onto each of the three mesh axes. The projected box must be entirely 
! inside the source box.

  implicit none
  integer, intent(in) :: nMesh(3)          ! Global mesh size in each axis
  integer, intent(in) :: srcBox(2,3)       ! Source data mesh box:
                                           ! srcBox(1,iAxis)=lower box limits
                                           ! srcBox(2,iAxis)=upper box limits
  real(gp),intent(in) :: srcData(0:,0:,0:) ! Source data values to be copied
                                           ! Its upper array bounds must be
                                           ! srcBox(2,iAxis)-srcBox(1,iAxis)
  integer, intent(in) :: prjBox(2,3)       ! Projection mesh box. Must be
                                           ! prjBox(1,:)>=srcBox(1,:) and
                                           ! prjBox(2,:)<=srcBox(2,:)
  real(gp),intent(out):: prjData(:,0:)     ! Projected data values:
                          ! prjData(iAxis,0:prjBox(2,iAxis)-prjBox(1,iAxis))

  character(len=*),parameter:: myName = moduleName//'projection '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: i1, i1max, i1min, i2, i2max, i2min, i3, i3max, i3min

! Check box limits
  if (any(prjBox(1,:)<srcBox(1,:)) .or. &
      any(prjBox(2,:)>srcBox(2,:)) )    &
    call die(errHead//'prjBox outside srcBox')

! Check array sizes
  if (any( shape(srcData) /= srcBox(2,:)-srcBox(1,:)+1 )) &
    call die(errHead//'shape of array srcData inconsistent with srcBox')
  if (size(prjData,2) < maxval(prjBox(2,:)-prjBox(1,:)+1)) &
    call die(errHead//'size of array prjData too small')

! Find prjBox limits relative to srcBox origin
  i1min = prjBox(1,1) - srcBox(1,1)
  i1max = prjBox(2,1) - srcBox(1,1)
  i2min = prjBox(1,2) - srcBox(1,2)
  i2max = prjBox(2,2) - srcBox(1,2)
  i3min = prjBox(1,3) - srcBox(1,3)
  i3max = prjBox(2,3) - srcBox(1,3)

! Find projection along each axis
  prjData = 0
  do i1 = i1min,i1max
    prjData(1,i1-i1min) = sum( srcData(i1,i2min:i2max,i3min:i3max) )
  end do
  do i2 = i2min,i2max
    prjData(2,i2-i2min) = sum( srcData(i1min:i1max,i2,i3min:i3max) )
  end do
  do i3 = i3min,i3max
    prjData(3,i3-i3min) = sum( srcData(i1min:i1max,i2min:i2max,i3) )
  end do

end subroutine projection

!-------------------------------------------------------------------------------

subroutine projectMeshData( nMesh, srcDistr, srcData, prjBox, prjData, task )

! Projects a box of values, defined in a 3D mesh of a periodic unit cell,
! onto each of the three mesh axes.
! The box limits may be partly or entirely outside the first unit cell.
! In parallel execution, the routine handles the interprocessor transfers
! needed to bring the requested box data from the nodes that store them.

  implicit none
  integer, intent(in) :: nMesh(3)          ! Global mesh size in each axis
  integer, intent(in) :: srcDistr          ! Distribution ID of source data
                                           ! or zero for serial execution
  real(gp),intent(in) :: srcData(0:,0:,0:) ! Source data values to be copied
                                           ! Its upper array bounds must be
                                           ! those assigned to the node by
                                           ! the srcDistr parallel distrib.
                                           ! or nMesh(iAxis)-1 in serial exe
  integer, intent(in) :: prjBox(2,3)       ! Projection mesh box:
                                           ! box(1,iAxis)=lower box limits
                                           ! box(2,iAxis)=upper box limits
                                           ! Box needs not be within the
                                           ! (0:nMesh-1) unit-cell range
  real(gp),intent(out):: prjData(:,0:)     ! Projected data values:
                          ! prjData(iAxis,0:prjBox(2,iAxis)-prjBox(1,iAxis))
  integer,intent(inout),optional:: task    ! ID of communication task

  character(len=*),parameter:: myName = moduleName//'projectMeshData '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer :: srcBox(2,3), srcMesh(3)
  real(gp),pointer:: mySrcData(:,:,:,:)=>null()

! Associate task to distributions
  if (present(task)) call associateMeshTask( task, srcDistr )

! Find box limits of source data in this node
  call myMeshBox( nMesh, srcDistr, srcBox )
  srcMesh = srcBox(2,:) - srcBox(1,:) + 1

! Check array shapes
!  if (any( shape(srcData) /= srcMesh )) &
!    call die( errHead//'incorrect shape of srcData array' )
  if (any( shape(srcData) /= srcMesh )) then
#ifdef DEBUG_XC
    write(udebug,'(a,i6,3(2x,2i4))') &
      errHead//'srcDistr,srcBox =', srcDistr, srcBox
    write(udebug,'(a,3i6,3x,3i6)') &
      errHead//'shape(srcData),srcMesh =', shape(srcData), srcMesh
#endif /* DEBUG_XC */
    call die( errHead//'incorrect shape of srcData array' )
  end if
  if (size(prjData,2) < maxval(prjBox(2,:)-prjBox(1,:)+1)) &
    call die( errHead//'size of prjData array too small' )

! Copy srcData into a rank-4 array to call reduceData
  call re_alloc( mySrcData, 0,srcMesh(1)-1, 0,srcMesh(2)-1, 0,srcMesh(3)-1, &
                 1,1, myName//'mySrcData' )
  mySrcData(0:srcMesh(1)-1,0:srcMesh(2)-1,0:srcMesh(3)-1,1) = &
    srcData(0:srcMesh(1)-1,0:srcMesh(2)-1,0:srcMesh(3)-1)

! Find projection
  call reduceData( nMesh, srcBox, mySrcData, prjBox, prjData=prjData, &
                       taskID=task )

! Deallocate temporary array
  call de_alloc( mySrcData, myName//'mySrcData' )

end subroutine projectMeshData

!-------------------------------------------------------------------------------

subroutine redistributeMeshData( srcDistr, srcData, dstDistr, dstData, task )

! Copies data values on the mesh with a new mesh distribution among nodes

  implicit none
  integer,        intent(in) :: srcDistr            ! Distrib. ID of srcData
  real(gp),target,intent(in) :: srcData(0:,0:,0:,:) ! Input source data
  integer,        intent(in) :: dstDistr            ! Distrib. ID of dstData
  real(gp),pointer           :: dstData(:,:,:,:)    ! Output destination array 
                                          ! (it must be different from srcData)
  integer,intent(inout),optional:: task   ! ID of communication task

  character(len=*),parameter:: myName = 'redistributeMeshData '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: dstBox(2,3), dstMesh(3), iDistr, nData, nMesh(3)
  integer:: i1, i2, i3, idat   ! provisional addition. JMS, June 2013

! For serial execution, avoid allocating a new array
  if (srcDistr==0 .and. dstDistr==0) then
    dstData => srcData
    return
  else if (srcDistr<=0 .or. dstDistr<=0) then
    call die(errHead//'invalid srcDistr or dstDistr')
  end if

! Find mesh box of node in new distribution
  iDistr = indexDistr( dstDistr )
  if (iDistr<=0 .or. iDistr>maxDistr) call die(errHead//'invalid dstDistr')
  nMesh = storedMeshDistr(iDistr)%nMesh
  call myMeshBox( nMesh, dstDistr, dstBox )
  dstMesh(:) = dstBox(2,:) - dstBox(1,:) + 1

! Allocate destination data array
  nData = size(srcData,4)
  call re_alloc( dstData, 0,dstMesh(1)-1, 0,dstMesh(2)-1, 0,dstMesh(3)-1, &
                 1,nData, name=myName//'dstData', copy=.false., shrink=.true. )

! Copy srcData to dstData
  if ( sameMeshDistr(srcDistr,dstDistr) ) then
    ! Provisional loops to circunvent an apparent compiler bug. JMS, June 2013.
    do idat = 1,nData
      do i3 = 0,dstMesh(3)-1
      do i2 = 0,dstMesh(2)-1
      do i1 = 0,dstMesh(1)-1
        dstData(i1,i2,i3,idat) = srcData(i1,i2,i3,idat)
      end do
      end do
      end do
    end do
!    dstData = srcData   ! original code
  else
    call copyMeshData( nMesh, srcDistr, srcData, dstBox, dstData, task )
  end if

end subroutine redistributeMeshData

!-------------------------------------------------------------------------------

subroutine reduceData( nMesh, srcBox, srcData, dstBox, dstData, prjData, &
                       taskID )

! Reduces (adds) a box of values defined in a 3D mesh of a periodic unit cell.
! The box limits may be entirely inside, or partially or totally outside 
! the unit cell, which is periodically repeated to fill the box.
! In parallel execution, the routine handles the interprocessor transfers
! needed to bring the requested box data from the nodes that store them.

  implicit none
  integer, intent(in) :: nMesh(3)          ! Global mesh size in each axis
  integer, intent(in) :: srcBox(2,3)       ! Source mesh box:
                                           ! box(1,iAxis)=lower box limits
                                           ! box(2,iAxis)=upper box limits
                                           ! Boxes may be larger than the
                                           ! (0:nMesh-1) unit-cell range
  real(gp),intent(in) :: srcData(0:,0:,0:,:) ! Source data values to be copied
                                           ! Its upper array bounds must be
                                           ! enough to contain the node box in
                                           ! the srcDistr parallel distrib.
                                           ! or nMesh(iAxis)-1 in serial exe
  integer, intent(in) :: dstBox(2,3)       ! Destination mesh box
  real(gp),intent(inout),optional:: dstData(0:,0:,0:,:)
                                           ! Destination data values.
                                           ! The upper array bounds must be 
                                           ! dstBox(2,iAxis)-dst(1,iAxis)
                                           ! dstData is NOT initialized
  real(gp),intent(out),optional:: prjData(:,0:) ! Projected data values:
                          ! prjData(iAxis,0:prjBox(2,iAxis)-prjBox(1,iAxis))
                          ! Only srcData(:,:,:,1) is projected
  integer,intent(inout),optional:: taskID  ! ID of communication task

  character(len=*),parameter:: myName = moduleName//'reduceData '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer,allocatable:: dstBoxes(:,:,:), srcBoxes(:,:,:), &
                        trsfDir(:), trsfNode(:)
  integer :: comBox(2,3), arrayS(4,2), &
             dstComBox(2,3,maxParts), dstNode, dstSize, &
             i, ic, ix, iNode, iPart, iTask, iTrsf, &
             maxMesh, nData, nParts, nTrsf, partSize, &
             sizeSum(0:maxParts), srcComBox(2,3,maxParts), &
             srcNode, srcSize, trsfSize, trueTrsf
  logical :: taskDefined, taskOptimized
  real(gp),allocatable:: prjBuff(:,:)
  real(gp),pointer:: trsfBuff(:)=>null()
  type(taskType),pointer:: task

#ifdef DEBUG_XC
!  write(udebug,'(a,i4,3(2x,2i4),2x,3(2x,2i4))') &
!    myName//'myNode+1,srcBox,dstBox=', myNode+1, srcBox, dstBox
#endif /* DEBUG_XC */

! Check that one of dstData or prjData is present
  if (.not.present(dstData) .and. .not.present(prjData)) return

! Check array shapes
  nData = size(srcData,4)
  do ic = 1,3
    if (size(srcData,ic) < srcBox(2,ic)-srcBox(1,ic)+1) &
      call die( errHead//'incorrect array shape of srcData' )
  end do
  if (present(dstData)) then
    do ic = 1,3
      if (size(dstData,ic) < dstBox(2,ic)-dstBox(1,ic)+1) &
        call die( errHead//'incorrect array shape of dstData' )
    end do
    if (size(dstData,4) /= nData) &
      call die( errHead//'inconsistent size of srcData and dstData' )
  end if
  if (present(prjData)) then
    if (size(prjData,2) < maxval(dstBox(2,:)-dstBox(1,:)+1)) &
      call die( errHead//'size of prjData array too small' )
  end if

! Allocate buffer for projection data transfer
  if (present(prjData)) then
    maxMesh = maxval(nMesh)
    allocate( prjBuff(3,0:maxMesh-1) )
  end if

! Find common intersection between srcBox and dstBox
  call commonBox( nMesh, srcBox, dstBox, srcComBox, dstComBox, sizeSum, nParts )

! Copy local data from each part of source box
  if (present(dstData)) then
    do iPart = 1,nParts
      dstData(dstComBox(1,1,iPart):dstComBox(2,1,iPart),     &
              dstComBox(1,2,iPart):dstComBox(2,2,iPart),     &
              dstComBox(1,3,iPart):dstComBox(2,3,iPart),:) = &
      dstData(dstComBox(1,1,iPart):dstComBox(2,1,iPart),     &
              dstComBox(1,2,iPart):dstComBox(2,2,iPart),     &
              dstComBox(1,3,iPart):dstComBox(2,3,iPart),:) + &
      srcData(srcComBox(1,1,iPart):srcComBox(2,1,iPart),     &
              srcComBox(1,2,iPart):srcComBox(2,2,iPart),     &
              srcComBox(1,3,iPart):srcComBox(2,3,iPart),:)
    end do ! iPart
  end if ! (present(dstData))

! Project local data onto axes of projection box
  if (present(prjData)) then
    prjData(:,:) = 0
    do iPart = 1,nParts
      ! Find common box relative to origin of unit cell
      comBox(1,:) = srcComBox(1,:,iPart) + srcBox(1,:)
      comBox(2,:) = srcComBox(2,:,iPart) + srcBox(1,:)
      ! Find projection of srcData in comBox
      call projection( nMesh, srcBox, srcData(:,:,:,1), comBox, prjBuff )
      ! Accumulate projection onto output array
      do ix = 1,3
        prjData(ix,dstComBox(1,ix,iPart):dstComBox(2,ix,iPart)) = &
        prjData(ix,dstComBox(1,ix,iPart):dstComBox(2,ix,iPart)) + &
          prjBuff(ix,0:comBox(2,ix)-comBox(1,ix))
      end do
    end do ! iPart
  end if ! (present(prjData))

! If totNodes=1 => serial case (all data are local) => all is done
  if (totNodes==1) return

! Else we need MPI code
#ifdef MPI

! Find the stored data about this communication task
  taskDefined = .false.
  taskOptimized = .false.
  if (present(taskID)) then
    iTask = indexTask( taskID )
    if (iTask<0) then  ! Create a new task
      call initTask( taskID )
      iTask = indexTask( taskID )
    end if ! (iTask<0)
    task => storedMeshTask(iTask)
    taskDefined = task%defined
    taskOptimized = task%optimized
  end if ! (present(task))

! Check task boxes
  if (taskDefined) then
    if (any( task%srcBox(:,:,myNode) /= srcBox )) &
      call die( errHead//'srcBox inconsistent with taskID' )
    if (any( task%dstBox(:,:,myNode) /= dstBox )) &
      call die( errHead//'dstBox inconsistent with taskID' )
  end if

! Find the box limits of all nodes
  allocate( srcBoxes(2,3,0:totNodes-1), dstBoxes(2,3,0:totNodes-1) )
  if (taskDefined) then
    srcBoxes = task%srcBox
    dstBoxes = task%dstBox
  else  ! Gather the boxes from all nodes
    call gatherBoxes( srcBox, srcBoxes )
    call gatherBoxes( dstBox, dstBoxes )
    if (present(taskID)) then  ! Store boxes for future calls
      task%srcBox = srcBoxes
      task%dstBox = dstBoxes
      task%defined = .true.
#ifdef DEBUG_XC
!      write(udebug,'(a,2i4,3(2x,2i4))') &
!        myName//'taskID,iTask,srcBox=', taskID, iTask, srcBox
!      write(udebug,'(31x,a,2i4,3(2x,2i4))') &
!                             'dstBox=', taskID, iTask, dstBox
#endif /* DEBUG_XC */
    end if
  end if

#ifdef DEBUG_XC
!  do iNode = 0,totNodes-1
!    call nodeMeshBox( nMesh, srcDistr, iNode, srcBox )
!    write(udebug,'(a,i4,2x,3(2i4,2x),3(2i4,2x))') &
!      myName//'node,srcBox,dstBox=', &
!      iNode, srcBoxes(:,:,iNode), dstBoxes(:,:,iNode)
!  end do
#endif /* DEBUG_XC */

! Find the order of inter-node transfers
  nTrsf = 2*(totNodes-1) ! One send and one receive to/from each other node
  allocate( trsfDir(0:nTrsf), trsfNode(0:nTrsf) )
  if (taskOptimized) then
    nTrsf = task%nTrsf
    trsfDir = task%trsfDir
    trsfNode = task%trsfNode
  else ! (.not.taskOptimized)
    call all2allTransferOrder( totNodes, myNode, nTrsf, trsfNode, trsfDir )
  end if ! (taskOptimized)

#ifdef DEBUG_XC
!  write(udebug,'(a,i4,2x,20i4)') myName//'myNode+1,Trsfs=', &
!    myNode+1, trsfDir(1:nTrsf)*(trsfNode(1:nTrsf)+1)
#endif /* DEBUG_XC */

! Find required size of buffer to transfer data
  trsfSize = 0
  if (present(dstData)) then
    do iTrsf = 1,nTrsf
      iNode = trsfNode(iTrsf)
      if (trsfDir(iTrsf)==+1) then ! Size required to send
        srcSize = product( srcBox(2,:)-srcBox(1,:)+1 )
        dstSize = product( dstBoxes(2,:,iNode)-dstBoxes(1,:,iNode)+1 )
        trsfSize = max( trsfSize, min(srcSize,dstSize) )
      else ! Size required to receive
        srcSize = product( srcBoxes(2,:,iNode)-srcBoxes(1,:,iNode)+1 )
        dstSize = product( dstBox(2,:)-dstBox(1,:)+1 )
        trsfSize = max( trsfSize, min(srcSize,dstSize) )
      end if ! (trsfDir(iTrsf)==+1)
    end do ! iTrsf
    trsfSize = trsfSize*nData
  end if
  if (present(prjData)) trsfSize = trsfSize + 3*maxMesh*maxParts

! Allocate buffer to transfer data
  call re_alloc( trsfBuff, 1,trsfSize, name=myName//'trsfBuff' )

! Loop on transfer sequence
  trueTrsf = 0
  do iTrsf = 1,nTrsf

    if (trsfDir(iTrsf)==+1) then ! My node will send data

      ! Find destination node
      dstNode = trsfNode(iTrsf)

      ! Find common intersection between srcBox and dstBox
      call commonBox( nMesh, srcBox, dstBoxes(:,:,dstNode), &
                      srcComBox, dstComBox, sizeSum, nParts )
      sizeSum = sizeSum * nData

      ! Select this transfer only if it is needed
      if (nParts>0) then
        trueTrsf = trueTrsf + 1
        trsfDir(trueTrsf) = trsfDir(iTrsf)
        trsfNode(trueTrsf) = trsfNode(iTrsf)
      else ! (nParts==0)
        cycle ! iTrsf loop
      end if ! (nParts>0)

#ifdef DEBUG_XC
      ! Check buffer size
      if (present(dstData) .and. sizeSum(nParts) > size(trsfBuff)) &
        call die(errHead//'size(trsfBuff) too small')
#endif /* DEBUG_XC */

#ifdef DEBUG_XC
!      if (nParts>0) then
!        write(udebug,'(a,2i6)') myName//'srcNode,dstNode=', myNode, dstNode
!        write(udebug,'(a,2i6)') myName//'        sizeSum=',sizeSum(0:nParts)
!        write(udebug,'(a,3(2i4,2x))') myName//'   srcBox=', srcBox
!        write(udebug,'(a,3(2i4,2x))') myName//'   dstBox=', &
!                                         dstBoxes(:,:,dstNode)
!      end if

!      do iPart = 1,nParts
!        write(udebug,'(a,3(2i4,2x))') myName//'srcComBox=', &
!          ((srcComBox(i,ix,iPart),i=1,2),ix=1,3)
!        write(udebug,'(a,3(2i4,2x))') myName//'dstComBox=', &
!          ((dstComBox(i,ix,iPart),i=1,2),ix=1,3)
!      end do ! iPart
#endif /* DEBUG_XC */

      ! Copy data from each part of common box to a transfer buffer
      if (present(dstData)) then
        do iPart = 1,nParts
          arrayS(1:3,1) = srcComBox(1,:,iPart) + 1
          arrayS(4,1) = 1
          arrayS(1:3,2) = srcComBox(2,:,iPart) + 1
          arrayS(4,2) = nData
          call array_copy(arrayS(:,1), arrayS(:,2), srcData(:,:,:,:), &
               sizeSum(iPart-1)+1, sizeSum(iPart), trsfBuff(:) )
          !partSize = sizeSum(iPart) - sizeSum(iPart-1)
          !trsfBuff(sizeSum(iPart-1)+1:sizeSum(iPart)) =                    &
          !  reshape( srcData(srcComBox(1,1,iPart):srcComBox(2,1,iPart),    &
          !                   srcComBox(1,2,iPart):srcComBox(2,2,iPart),    &
          !                   srcComBox(1,3,iPart):srcComBox(2,3,iPart),:), &
          !           (/partSize/) )
        end do ! iPart
        trsfSize = sizeSum(nParts)
      else ! (.not.present(dstData))
        trsfSize = 0
      end if ! (present(dstData))

      ! Project data from each part of common box onto the transfer buffer
      if (present(prjData)) then
        do iPart = 1,nParts
          ! Find common box relative to origin of first unit cell
          comBox(1,:) = srcComBox(1,:,iPart) + srcBox(1,:)
          comBox(2,:) = srcComBox(2,:,iPart) + srcBox(1,:)
          ! Find projection of srcData in comBox
          call projection( nMesh, srcBox, srcData(:,:,:,1), comBox, prjBuff )
          ! Copy projection data onto a transfer buffer
          call array_copy([1, 1], [3, maxMesh], prjBuff(:,:), &
               trsfSize+1, trsfSize+3*maxMesh, trsfBuff(:))
          ! AG: Work around gfortran bug: specify 0 lbound in reshape
          !trsfBuff(trsfSize+1:trsfSize+3*maxMesh) = &
          !   reshape( prjBuff(:,0:), (/3*maxMesh/) )
          trsfSize = trsfSize + 3*maxMesh
        end do ! iPart
      end if ! (present(prjData))

#ifdef DEBUG_XC
!      write(udebug,'(a,2i4,3(2x,2i4),2x,3(2x,2i4),i8)') &
!        myName//' myNode,dstNode, myBox,dstBox,trsfSize=', &
!        myNode+1, dstNode+1, srcBox, dstBoxes(:,:,dstNode), trsfSize
#endif /* DEBUG_XC */

      ! Send data buffer
      if (trsfSize>0) &
        call MPI_Send( trsfBuff(1:trsfSize), trsfSize, MPI_grid_real, &
                       dstNode, 0, MPI_COMM_WORLD, MPIerror )

    else ! My node will receive data

      ! Find source node and its mesh box
      srcNode = trsfNode(iTrsf)

      ! Find common intersection between srcBox and dstBox
      call commonBox( nMesh, srcBoxes(:,:,srcNode), dstBox, &
                      srcComBox, dstComBox, sizeSum, nParts )
      sizeSum = sizeSum * nData

      ! Select this transfer only if it is needed
      if (nParts>0) then
        trueTrsf = trueTrsf + 1
        trsfDir(trueTrsf) = trsfDir(iTrsf)
        trsfNode(trueTrsf) = trsfNode(iTrsf)
      else ! (nParts==0)
        cycle ! iTrsf loop
      end if ! (nParts>0)

#ifdef DEBUG_XC
      ! Check buffer size
      if (present(dstData) .and. sizeSum(nParts)>size(trsfBuff)) &
        call die(errHead//'size(trsfBuff) too small')
#endif /* DEBUG_XC */

#ifdef DEBUG_XC
!      if (nParts>0) then
!        write(udebug,'(a,2i6)') myName//'srcNode,dstNode=', srcNode, myNode
!        write(udebug,'(a,2i6)') myName//'        sizeSum=',sizeSum(0:nParts)
!        write(udebug,'(a,3(2i4,2x))') myName//'   srcBox=', srcBox
!        write(udebug,'(a,3(2i4,2x))') myName//'   dstBox=', dstBox
!      end if

!      do iPart = 1,nParts
!        write(udebug,'(a,3(2i4,2x))') myName//'srcComBox=', &
!          ((srcComBox(i,ix,iPart),i=1,2),ix=1,3)
!        write(udebug,'(a,3(2i4,2x))') myName//'dstComBox=', &
!          ((dstComBox(i,ix,iPart),i=1,2),ix=1,3)
!      end do ! iPart
#endif /* DEBUG_XC */

      ! Receive data buffer
      trsfSize = 0
      if (present(dstData)) trsfSize = trsfSize + sizeSum(nParts)
      if (present(prjData)) trsfSize = trsfSize + 3*maxMesh*nParts
      if (trsfSize>0) &
        call MPI_Recv( trsfBuff(1:trsfSize), trsfSize, MPI_grid_real, &
                       srcNode, 0, MPI_COMM_WORLD, MPIstatus, MPIerror )

#ifdef DEBUG_XC
!      write(udebug,'(a,2i4,3(2x,2i4),2x,3(2x,2i4),i8)') &
!        myName//'srcNode, myNode,srcBox, myBox,trsfSize=', &
!        srcNode+1, myNode+1, srcBoxes(:,:,srcNode), dstBox, trsfSize
#endif /* DEBUG_XC */

      ! Copy buffer data for each part of common box to destination array
      if (present(dstData)) then
        do iPart = 1,nParts
          arrayS(1:3,1) = dstComBox(1,:,iPart) + 1
          arrayS(4,1) = 1
          arrayS(1:3,2) = dstComBox(2,:,iPart) + 1
          arrayS(4,2) = nData
          call array_add(sizeSum(iPart-1)+1, sizeSum(iPart), trsfBuff(:), &
               arrayS(:,1), arrayS(:,2), dstData(:,:,:,:))
          !comShape(1:3) = dstComBox(2,:,iPart) - dstComBox(1,:,iPart) + 1
          !comShape(4) = nData
          !dstData(dstComBox(1,1,iPart):dstComBox(2,1,iPart),     &
          !        dstComBox(1,2,iPart):dstComBox(2,2,iPart),     &
          !        dstComBox(1,3,iPart):dstComBox(2,3,iPart),:) = &
          !dstData(dstComBox(1,1,iPart):dstComBox(2,1,iPart),     &
          !        dstComBox(1,2,iPart):dstComBox(2,2,iPart),     &
          !        dstComBox(1,3,iPart):dstComBox(2,3,iPart),:) + &
          !  reshape( trsfBuff(sizeSum(iPart-1)+1:sizeSum(iPart)), comShape )
        end do ! iPart
        trsfSize = sizeSum(nParts)
      else ! (.not.present(dstData))
        trsfSize = 0
      end if ! (present(dstData))

      ! Copy buffer data for each part of common box to destination array
      if (present(prjData)) then
        do iPart = 1,nParts
          ! Find common box relative to origin of first unit cell
          comBox(1,:) = srcComBox(1,:,iPart) + srcBox(1,:)
          comBox(2,:) = srcComBox(2,:,iPart) + srcBox(1,:)
          ! Copy transfer buffer onto a projection buffer
          call array_copy(trsfSize+1, trsfSize+3*maxMesh, trsfBuff(:), &
               [1, 1], [3, maxMesh], prjBuff(:,:))
          !prjBuff(:,:) = &
          !  reshape( trsfBuff(trsfSize+1:trsfSize+3*maxMesh), (/3,maxMesh/) )
          ! Add projection of common box to output array
          do ix = 1,3
            prjData(ix,dstComBox(1,ix,iPart):dstComBox(2,ix,iPart)) = &
            prjData(ix,dstComBox(1,ix,iPart):dstComBox(2,ix,iPart)) + &
              prjBuff(ix,0:comBox(2,ix)-comBox(1,ix))
          end do
          trsfSize = trsfSize + 3*maxMesh
        end do ! iPart
      end if ! (present(prjData))

    end if ! (iDir(iTrsf)==+1)

  end do ! iTrsf

  ! Optimize transfer sequence and store it
  if (present(taskID) .and. .not.taskOptimized) then
    nTrsf = trueTrsf
    call optimizeTransferOrder( totNodes, myNode, nTrsf, &
                                trsfNode(1:nTrsf), trsfDir(1:nTrsf) )
    task%nTrsf = nTrsf
    task%trsfNode = trsfNode
    task%trsfDir = trsfDir
    task%optimized = .true.
#ifdef DEBUG_XC
!    write(udebug,'(a,4i4,/,(20i4))') &
!      myName//'taskID,iTask,distr1,distr2,Trsfs=', &
!      taskID, iTask, task%distr, trsfDir(1:nTrsf)*trsfNode(1:nTrsf)
#endif /* DEBUG_XC */
  end if

  ! Deallocations
  call de_alloc( trsfBuff, name=myName//'trsfBuff' )
  deallocate( trsfDir, trsfNode, dstBoxes, srcBoxes )

#endif
! End of MPI code

  if (present(prjData)) deallocate( prjBuff )

end subroutine reduceData

!-------------------------------------------------------------------------------

subroutine reduceDistr( distrID )

! Find if this distribution was already defined, and resets it if so

  implicit none
  integer,intent(inout):: distrID  ! Distribution ID

  character(len=*),parameter:: myName = moduleName//'reduceDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  logical:: found
  integer:: iDistr, iID, jDistr
  type(distrType),pointer:: newDistr, oldDistr

! Find the new distribution
  iDistr = indexDistr( distrID )
  if (iDistr==0) then
    return
  else if (iDistr<0 .or. iDistr>maxDistr) then
    call die(errHead//'invalid distrID')
  end if
  newDistr => storedMeshDistr(iDistr)

! Compare it with all previous distributions
  do jDistr = 1,maxDistr

    ! Check that everything agrees, otherwise cycle
    if (jDistr==iDistr) cycle ! jDistr loop
    oldDistr => storedMeshDistr(jDistr)
    if (.not.oldDistr%defined) cycle
    if (oldDistr%nNodes/=newDistr%nNodes) cycle
    if (oldDistr%firstNode/=newDistr%firstNode) cycle
    if (any(oldDistr%nMesh(:)/=newDistr%nMesh(:))) cycle
    if (any(oldDistr%box(:,:,:)/=newDistr%box(:,:,:))) cycle
    ! If this point is reached, both distributions are equal.

#ifdef DEBUG_XC
!    write(udebug,'(a,3i4)') &
!      myName//'distrID,oldDistr,newDistr=', distrID, jDistr, iDistr
#endif /* DEBUG_XC */

    ! Free new distribution
    call freeMeshDistr( distrID )

    ! Assign new ID to oldDistr
    found = .false.
    do iID = 1,maxDistrID  ! Loop on ID slots of oldDistr
      if (oldDistr%ID(iID)<0) then ! Empty ID slot
        oldDistr%ID(iID) = distrID
        found = .true.
        exit ! iID loop
      end if
    end do ! iID
    if (.not.found) call die(errHead//'parameter maxDistrID too small')

    exit ! jDistr loop
  end do ! jDistr

end subroutine reduceDistr

!-------------------------------------------------------------------------------

logical function sameMeshDistr( ID1, ID2 )

  implicit none
  integer,intent(in):: ID1, ID2  ! Mesh distribution IDs

  integer:: i1, i2
  type(distrType),pointer:: distr1, distr2

  if (ID1==ID2) then           ! Same IDs => same distr.
    sameMeshDistr = .true.
  else                         ! Different IDs
    i1 = indexDistr(ID1)
    i2 = indexDistr(ID2)
    if (i1<0 .or. i2<0) then   ! One or both distr. are not defined
      sameMeshDistr = .false.
    else                       ! Both distr. are defined
      if (i1==i2) then         ! Different IDs but same distr.
        sameMeshDistr = .true.
      else if (i1==0 .or. i2==0) then ! Since 0 is valid only in serial mode
        sameMeshDistr = .false.
      else                     ! Different but possibly equivalent distr.
        distr1 => storedMeshDistr(i1)
        distr2 => storedMeshDistr(i2)
        if (all(distr1%box==distr2%box)) then
          sameMeshDistr = .true.
        else
          sameMeshDistr = .false.
        end if ! (all(distr1%box==distr2%box))
      end if ! (i1==i2)
    end if ! (i1>0 .and. i2>0)
  end if ! (ID1==ID2)

end function sameMeshDistr

!-------------------------------------------------------------------------------

subroutine setMeshDistr( distrID, nMesh, box, firstNode, nNodes, &
                         nNodesX, nNodesY, nNodesZ, nBlock, &
                         wlDistr, workload )

! Defines a parallel distribution of mesh points among processor nodes.
! In serial execution (totNodes==1) it simply returns distrID=0.
! If the input distribution ID is still valid, the same value is returned.
! Otherwise, it is freed. New IDs are never repeated.

  implicit none
  integer,       intent(inout):: distrID   ! ID assigned to the mesh distrib.
  integer,         intent(in) :: nMesh(3)  ! Mesh divisions in each axis
                                           ! (including "subpoints")
  integer,optional,intent(in) :: box(2,3)  ! Mesh box of my processor node
                                           ! If present(box), all other  
                                           ! optional arguments are ignored
  integer,optional,intent(in) :: firstNode ! First node in the mesh distr.
  integer,optional,intent(in) :: nNodes    ! Total nodes in the mesh distr
                                           ! If not present all nodes are used
  integer,optional,intent(in) :: nNodesX   ! Nodes in the X (first) axis. Must
                                           ! be present if present(nNodesYZ)
  integer,optional,intent(in) :: nNodesY   ! Nodes in the Y (second) axis
                                           ! Must be present if present(nNodesZ)
  integer,optional,intent(in) :: nNodesZ   ! Nodes in the Z (third) axis
                                           ! If nNodesXYZ are not present,
                                           ! the distribution of nodes on all
                                           ! axes is optimized by setMeshDistr.
                                           ! If only nNodesX is present, it is
                                           ! optimized in Y and Z axes
  integer,optional,intent(in) :: nBlock    ! Size of blocks of mesh points, in
                                           ! each axis, which are not splitted
                                           ! in different nodes
  integer,optional,intent(in) :: wlDistr   ! Distr. index of workload array
  real(gp),optional,intent(in):: workload(0:,0:,0:) ! Approx. relative workload
                                           ! of mesh points. Must be nonnegative
                                           ! at all points and have nonzero sum

  character(len=*),parameter:: myName = 'setMeshDistr '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer,         parameter:: maxFactors = 100  ! Max prime factors in nNodes
  type(distrType),pointer:: distr, newDistr, oldDistr
  integer,allocatable:: axisBox(:,:,:), partBox(:,:,:)
  integer:: axis, axisNodes(3), blockSize, boxSize, &
            factor(maxFactors), groupSize, &
            i1, i2, i3, iAxis, iBox, iDistr, iFac, iID, iNode, iPow, &
            jDistr, largerBoxSize, lastNode, &
            maxFactor, meshNodes, myGroup, myBox(2,3), myPart, &
            nFactors, node0, nParts, nRem, nBlocks, newDistrID, &
            oldDistrID, partSize, power(maxFactors)
  logical:: found
  character(len=5):: dat

! Make sure that myNode and totNodes are defined
  call parallel_init()

! Trap the serial case (totNodes available from module parallel)
  if (totNodes==1) then
    distrID = 0
    return
  end if

! Find number of mesh nodes
  if (present(nNodes)) then
    if (nNodes<1 .or. nNodes>totNodes) then
      call die(errHead//'nNodes out of range')
    else
      meshNodes = nNodes
    end if
  else ! Assume that all nodes participate in the mesh distribution
    meshNodes = totNodes
  end if

! Find first node in mesh distribution
  if (present(firstNode)) then
    lastNode = firstNode + meshNodes - 1
    if (firstNode<0 .or. lastNode>totNodes-1) then
      call die(errHead//'firstNode/nNodes out of range')
    else 
      node0 = firstNode
    end if
  else ! (.not.present(firstNode)) => use firstNode=0
    node0 = 0
  end if

! Find size of blocks of points and check that nMesh is a multiple of it
  if (present(nBlock)) then
    blockSize = nBlock
  else
    blockSize = 1
  end if
  if (any(mod(nMesh,blockSize)/=0)) &
    call die(errHead//'nMesh inconsistent with nBlock')

! Store input distribution ID for later comparison
  oldDistrID = distrID

! Initialize new distribution
  call initDistr( newDistrID, nMesh, node0, meshNodes )
  iDistr = indexDistr( newDistrID )
  distr => storedMeshDistr(iDistr)  ! Just a shorter name

! Handle box argument with priority
  if (present(box)) then
    ! Collect all node boxes and store them
    call gatherBoxes( box, distr%box )
    goto 999  ! Exit, since no other arguments must be considered in this case
  end if ! (present(box))

! Find box sizes
  if (present(workload)) then

    ! Check that the distribution of workload array is available
    if (.not.present(wlDistr)) call die(errHead//'workload without wlDistr')

    ! Check that myNode (available from module parallel) contains some points
    if (myNode>=node0 .and. myNode<node0+meshNodes) then

      ! Find prime factors of meshNodes, in ascending order
      call primeFactors( meshNodes, factor, power, nFactors )

      ! Allocate a small temporary array
      maxFactor = maxval( factor(1:nFactors) )
      allocate( partBox(2,3,0:maxFactor-1) )

      ! Iteratively divide full box and select my node's part
      myBox(1,:) = 0           ! Initialize all boxes with full box
      myBox(2,:) = nMesh(:)-1
      groupSize = meshNodes    ! Group = set of nodes assigned to present box
      do iFac = nFactors,1,-1  ! Use larger factors first (for larger boxes)
        do iPow = 1,power(iFac)
          nParts = factor(iFac) ! Number of parts to divide present box
          partSize = groupSize / nParts  ! Nodes that will be assigned to part
          myGroup = (myNode - node0) / groupSize  ! My group for present box
          myPart = (myNode-node0-myGroup*groupSize) / partSize ! My part within 
                                                               ! my group
          axis = 0  ! So that it will be chosen by divideBox3D
#ifdef DEBUG_XC
!          write(udebug,'(a,3(2x,2i4))') myName//'dividing box=', myBox
#endif /* DEBUG_XC */
          call divideBox3D( nMesh, wlDistr, workload, myBox, &
                              nParts, blockSize, axis, partBox )
          myBox(:,:) = partBox(:,:,myPart)  ! Select my part of the box
          groupSize = partSize              ! Reduce group size
#ifdef DEBUG_XC
!          write(udebug,'(a,3(2x,2i4))') myName//' my part box=', myBox
#endif /* DEBUG_XC */
        end do ! iPow
      end do ! iFac

      deallocate( partBox )

    else ! my node contains no mesh points
      myBox(1,:) = 0
      myBox(2,:) = -1
    end if ! (myNode>=node0 .and. myNode<node0+meshNodes)

    ! Gather the box limits of all nodes
    call gatherBoxes( myBox, distr%box )

  else ! (.not.present(workload)) => distribute nodes homogeneously in space

    ! Find distribution of nodes along each axis
    if (present(nNodesX)) then
      axisNodes(1) = nNodesX
      if (present(nNodesY)) then
        axisNodes(2) = nNodesY
        if (present(nNodesZ)) then
          axisNodes(3) = nNodesZ
        else
          axisNodes(3) = meshNodes / (nNodesX*nNodesY)
        end if
      else if (present(nNodesZ)) then
        call die(errHead//'nNodesZ present without nNodesY')
      else ! (.not.present(nNodesYZ)) => optimize distribution in Y and Z axes
        call optimizeNodeDistr( nMesh(2:3), meshNodes/nNodesX, axisNodes(2:3) )
      end if ! (present(nNodesY))
    else if (present(nNodesY).or.present(nNodesZ)) then
      call die(errHead//'nNodesY or nNodesZ present without nNodesX')
    else ! (.not.present(nNodesXYZ)) => optimize distribution in all three axes
      call optimizeNodeDistr( nMesh(1:3), meshNodes, axisNodes(1:3) )
    end if ! (present(nNodesX))

    ! Adjust number of nodes to that along the axes
    if (product(axisNodes) > meshNodes) then
      call die(errHead//'inconsistent nNodesXYZ')
    else
      meshNodes = product(axisNodes)
    end if

    ! Allocate a small temporary array
    allocate( axisBox(2,3,-1:maxval(axisNodes)-1) )

    ! Find mesh partition along each axis
    do iAxis = 1,3
      nBlocks = nMesh(iAxis) / blockSize         ! Blocks along axis
      boxSize = nBlocks / axisNodes(iAxis)       ! Blocks per node
      nRem = nBlocks - boxSize*axisNodes(iAxis)  ! Remaining blocks
      largerBoxSize = (boxSize+1)*blockSize      ! Mesh points per node
      boxSize = boxSize*blockSize                ! Mesh points per node
      axisBox(2,iAxis,-1) = -1             ! So that axisBox(1,ix,0)=0
      do iBox = 0,axisNodes(iAxis)-1
        axisBox(1,iAxis,iBox) = axisBox(2,iAxis,iBox-1) + 1
        if (iBox < nRem) then ! First nRem boxes of larger size
          axisBox(2,iAxis,iBox) = axisBox(2,iAxis,iBox-1) + largerBoxSize
        else                    ! Rest of boxes of normal size
          axisBox(2,iAxis,iBox) = axisBox(2,iAxis,iBox-1) + boxSize
        end if
      end do ! iBox
    end do ! iAxis

    ! Find box bounds of each node
    ! Innermost loop is i3 for compatibility with JDG distribution
    iNode = node0
    do i1 = 0,axisNodes(1)-1
    do i2 = 0,axisNodes(2)-1
    do i3 = 0,axisNodes(3)-1
      distr%box(:,1,iNode) = axisBox(:,1,i1)
      distr%box(:,2,iNode) = axisBox(:,2,i2)
      distr%box(:,3,iNode) = axisBox(:,3,i3)
      iNode = iNode + 1
    end do ! i3
    end do ! i2
    end do ! i1

    ! Deallocate temporary array
    deallocate( axisBox )
    
  end if ! (present(workload))

! Exit point
999 continue

! Find if the old (input) distrID is still valid
  if (sameMeshDistr(oldDistrID,newDistrID)) then  ! Old distr. is still valid
    call freeMeshDistr( newDistrID )
    distrID = oldDistrID
  else   ! Old distribution is no longer valid
    call freeMeshDistr( oldDistrID )
    ! Reset newDistrID if this distribution was already defined
    call reduceDistr( newDistrID )
    distrID = newDistrID
#ifdef DEBUG_XC
    iDistr = indexDistr( distrID )
    distr => storedMeshDistr(iDistr)
    if (present(box)) then
      dat = 'box'
    elseif (present(workload)) then
      dat = 'wkld'
    elseif (present(nNodesX)) then
      dat = 'xNod'
    else
      dat = 'none'
    endif
    write(udebug,'(a,a6,3i4,3(2x,2i4))') &
      myName//'dat,old/newDistrID,iDistr,myBox=', &
      dat, oldDistrID, newDistrID, iDistr, distr%Box(:,:,myNode)
#endif /* DEBUG_XC */
  end if

end subroutine setMeshDistr

!-------------------------------------------------------------------------------

subroutine unitCellParts( nMesh, fullBox, partBox, partBox0, nParts )

! Divides a mesh box into parts, each within a single periodic cell
! and returns the equivalent partial boxes within the first unit cell

  implicit none
  integer,intent(in) :: nMesh(3)        ! Number of mesh points in each axis
                                        ! of the periodic unit cell
  integer,intent(in) :: fullBox(2,3)    ! Box to be partitioned, possibly
                                        ! extending over several unit cells
  integer,intent(out):: partBox(:,:,:)  ! Partitioned boxes, each within a
                                        ! single unit cell:
                                        ! box(1,iAxis,iPart)=lower bounds
                                        ! box(2,iAxis,iPart)=upper bounds
  integer,intent(out):: partBox0(:,:,:) ! Equivalent part boxes in the first
                                        ! unit cell
  integer,intent(out):: nParts          ! Number of parts found

  character(len=*),parameter:: myName = moduleName//'unitCellParts '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: iDiv, iPart, ix, maxParts
  real(dp):: rMesh(3)

  maxParts = size(partBox,3)
  rMesh = nMesh   ! Real version of nMesh, so that partBox/rMesh is real
  partBox(:,:,1) = fullBox
  nParts = 1
  divideLoop: do ! Until no part can be further divided
    do iPart = 1,nParts  ! Loop on already existing parts
      do ix = 1,3        ! Axis loop
        ! iDiv is the highest surface, separating two unit cells, below the
        ! upper bound of partBox
        iDiv = floor(partBox(2,ix,iPart)/rMesh(ix)) * nMesh(ix)
        if (partBox(1,ix,iPart) < iDiv) then ! iDiv is between lower and
                                             ! upper bounds of partBox
          nParts = nParts + 1
          if (nParts > maxParts) call die(errHead//'maxParts too small')
          partBox(:,:,nParts) = partBox(:,:,iPart)
          partBox(1,ix,nParts) = iDiv       ! New part is upper one
          partBox(2,ix,iPart) = iDiv - 1    ! Old part reduced to lower one
          cycle divideLoop
        endif
      end do ! ix
    end do ! iPart
    exit divideLoop ! since no part can be further divided
  end do divideLoop

  ! Find equivalent boxes in first unit cell
  forall(ix=1:3) &
    partBox0(:,ix,1:nParts) = modulo( partBox(:,ix,1:nParts), nMesh(ix) )

end subroutine unitCellParts

!-------------------------------------------------------------------------------

END MODULE mesh3D

