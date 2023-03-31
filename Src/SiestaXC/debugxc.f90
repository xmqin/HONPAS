!!@LICENSE
!
MODULE debugXC

! Initializes output files for debug info. J.M.Soler. Jan'2008

USE parallel,         only: node, Nodes   ! This node number, Number of nodes
USE parallel,         only: parallel_init ! Initializes parallel variables
USE moreParallelSubs, only: nodeString    ! Returns a string with a node number
USE moreParallelSubs, only: copyFile      ! Copies a file to another node
USE m_io,             only: io_assign     ! Get and reserve an available IO unit

  implicit none

PUBLIC &
  setDebugOutputUnit,   &! Sets debug output unit and opens it for debug output
  closeDebugOutputFile, &! Close debug output file and copies it to node 0
  udebug                 ! Output file unit for debug info

PRIVATE ! Nothing is declared public beyond this point

  character(len=*),parameter:: filePrefix = 'debugXC' ! Prefix of file name
  character(len=32),    save:: fileName  ! Output file name
  integer,              save:: udebug=0  ! Output file unit for debug info

CONTAINS

!-----------------------------------------------------------------------------

subroutine setDebugOutputUnit()

  ! Sets debug output unit and opens file for debug output

  implicit none

  ! If already initialized, do nothing
  if (udebug>0) return

  ! Make sure that node and Nodes are set
  call parallel_init()

  ! Set output file name
  if (Nodes==1) then
    fileName = filePrefix//'.out'
  else
    fileName = filePrefix//'.node'//nodeString()
  end if

  ! Find an available output unit
  call io_assign( udebug )

  ! Open file and write header on it
  open( udebug, file=fileName )
  if (Nodes>1) &
    write(udebug,'(/,a)') 'Debug output for processor '//nodeString()

end subroutine setDebugOutputUnit

!-----------------------------------------------------------------------------

subroutine closeDebugOutputFile()

  ! Closes debug output file and copies it to node 0

  implicit none
  integer:: iNode
  character(len=32):: nodeFileName

  close( unit=udebug )

  ! Make sure that node and Nodes are set
  call parallel_init()

  ! Copy all debug.out* files to root node. Notice that in each call to
  ! copyFile, only node iNode will actually send its file to node=0
  do iNode = 1,Nodes-1
    nodeFileName = filePrefix//'.node'//nodeString(iNode)
    call copyFile( srcNode=iNode, srcFile=nodeFileName, &
                   dstNode=0,     dstFile=nodeFileName, &
                   writeOption='overwrite' )
  end do

end subroutine closeDebugOutputFile

END MODULE debugXC
