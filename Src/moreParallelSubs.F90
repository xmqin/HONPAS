! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
!******************************************************************************
! MODULE moreParallelSubs
! Provides some utility routines for parallel execution
! Written by J.M.Soler. Feb.2008-Jul.2009
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! copyFile      : Copies a formatted file from one node to another
! miscAllReduce : Reduces a miscellaneous set of variables and arrays
! nodeString    : Returns a string with a node number
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
! use m_io,       only: io_assign, io_close  ! Get and reserve an available IO unit
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
!
!   USED MPI procedures, interfaces and types:
! use mpi_siesta, only: MPI_AllReduce
! use mpi_siesta, only: MPI_Recv
! use mpi_siesta, only: MPI_Send
! use mpi_siesta, only: MPI_Integer
! use mpi_siesta, only: MPI_double_precision
! use mpi_siesta, only: MPI_character
! use mpi_siesta, only: MPI_Max
! use mpi_siesta, only: MPI_Min
! use mpi_siesta, only: MPI_Prod
! use mpi_siesta, only: MPI_Sum
! use mpi_siesta, only: MPI_COMM_WORLD
! use mpi_siesta, only: MPI_STATUS_SIZE
!
!   EXTERNAL procedures used:
! none
!
!******************************************************************************
!
! SUBROUTINE copyFile( srcNode, srcFile, dstNode, dstFile, writeOption )
!
! Copies a formatted file from one node to another
!------------------------------ INPUT -----------------------------------------
! integer,        :: srcNode   ! Source processor index in range (0:nNodes-1)
! character(len=*):: srcFile   ! Source file name
! integer,        :: dstNode   ! Destination processor index
! character(len=*):: dstFile   ! Destination file name
! character(len=*):: writeOption  ! ('append'|'overwrite')
!----------------------------- USAGE ------------------------------------------
! Sample usage to send each node's file into the file system of the root node
!   use mpi, only: MPI_Comm_Size, MPI_Comm_World
!   use moreParallelSubs, only: copyFile, nodeString
!   integer:: node, nNodes, MPIerror
!   character(len=32):: myFile
!   call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
!   myFile = 'node'//trim(nodeString())//'.out'  ! e.g. myFile='node021.out'
!   do node = 0,nNodes-1   ! The whole loop must be executed in each node
!     call copyFile( node, myFile, 0, myFile, 'overwrite' )
!   end do
!
! Sample usage to copy each node's file into the same destination file
!   use mpi, only: MPI_Comm_Size, MPI_Comm_World
!   use moreParallelSubs, only: copyFile
!   integer:: node, nNodes, MPIerror
!   call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
!   do node = 0,nNodes-1   ! The whole loop must be executed in each node
!     call copyFile( node, 'myNode.out', 0, 'allNodes.out', 'append' )
!   end do
!
! IMPORTANT: copyFile should be called from ALL nodes, even if they are neither
! the sending nor the receiving node, because there is an MPI_Barrier inside
!---------------------------- BEHAVIOUR ---------------------------------------
! If srcNode==dstNode and srcFile==dstFile, it does nothing
! If srcNode==dstNode but srcFile/=dstFile, it just copies the files, without
!   using MPI communications
! It allows copying several node's files sequentially into the same
!   destination file, in a given order, or to copy just one. But see
!   the note at the end of the USAGE section
!--------------------------- ALGORITHMS ---------------------------------------
! In srcNode, it reads srcFile into a character buffer, then sends it through
! MPI to dstNode, where it is writen into dstFile.
! An internal MPI_Barrier is used to impose that calls to copyFile from
! different nodes are executed in the required order.
!
!******************************************************************************
!
! SUBROUTINE miscAllReduce( op, a0, b0, c0, d0, e0, f0, &
!                               a1, b1, c1, &
!                               a2, b2, &
!                               a3 )
!
! Reduces a miscellaneous set of variables and arrays
!------------------------------ INPUT -----------------------------------------
! character(len=*) op  ! Operation to be applied ('sum'|'prod'|'max'|'min')
!-------------------- OPTIONAL INPUT and OUTPUT -------------------------------
! type                 :: a0, b0, c0, d0, e0, f0  ! Scalar arguments
! type,dimension(:)    :: a1, b1, c1              ! Vector arguments
! type,dimension(:,:)  :: a2, b2                  ! Rank-2 array arguments
! type,dimension(:,:,:):: a3                      ! Rank-3 array arguments
!   Where type=(integer|double precision)
!   All arguments must be of the same type
!----------------------------- UNITS ------------------------------------------
! Units of all arguments are arbitrary
!----------------------------- USAGE ------------------------------------------
! use precision,        only: dp
! use moreParallelSubs, only: miscAllReduce
! real(dp):: Ek, Eh, polarization(3), stress(3,3)
! real(dp),allocatable:: force(:,:)
! ... Find Ek, Eh, polarization, and stress. Allocate and find force.
! call miscAllReduce( 'sum', Ek, Eh, a1=polarization, a2=force, b2=stress )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution, it does nothing.
! If op/=('sum'|'prod'|'max'|'min'), it stops with an error message.
!--------------------------- ALGORITHMS ---------------------------------------
! All present variables and arrays are packed consecutively in a buffer vector
! The buffer is reduced with MPI_AllReduce and unpacked back into the present
!   variables and arrays.
!
!******************************************************************************
!
! character(len=6) function nodeString( node )
!
! Returns a string with a node number, or blank if nNodes<2, i.e. 00,01,...,63
! If node is not present, it assumes node=myNode
!------------------------------ INPUT -----------------------------------------
! integer,optional:: node  ! a valid node index: 0 <= node <= nNodes-1
!----------------------------- USAGE ------------------------------------------
! use moreParallelSubs, only: nodeString
! character(len=30):: fileName
! fileName = 'atoms.node'//trim(nodeString())//'.out'
!---------------------------- BEHAVIOUR ---------------------------------------
! If node is not present, it assumes node=myNode
! If node<0 or node>=nNodes, it stops with an error message
! In serial execution, it returns a blank string
! The string length is the same for all nodes (i.e. '0','1',...,'7', or 
!    '000','001',...,'131'), but as short as possible for the highest node
! The strings are left-justified, with right-hand blanks, i.e. '00    '
!
!******************************************************************************

MODULE moreParallelSubs

! Used module procedures
  use alloc,     only: de_alloc  ! De-allocation routine
  use alloc,     only: re_alloc  ! Re-allocation routine
  use sys,       only: die       ! Termination routine
  use m_io,      only: io_assign, io_close ! Get and reserve an available IO unit

! Used module parameters
  use precision, only: dp              ! Real double precision type
  use parallel,  only: myNode=>node    ! This process node
  use parallel,  only: totNodes=>nodes ! Total number of processor nodes

! Used MPI procedures and types
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

! All public procedures (there are no public types, parameters, or variables):
PUBLIC:: &
  copyFile,         &! Copies a formatted file from one node to another
  miscAllReduce,    &! Reduces a miscellaneous set of variables and arrays
  miscAllReduceInt, &! Integer version
  nodeString         ! Returns a character string with a node index

PRIVATE ! Nothing is declared public beyond this point

! Overload procedure name
! Note: a compiler warning may occur because calling miscAllReduce without any
!       optional arguments would not allow to resolve the right module proced.
!       But such a call would make no sense and it will never happen.
  interface miscAllReduce
    module procedure      &
!      miscAllReduceInt,   &! Integer version
      miscAllReduceDouble  ! Double precision real version
  end interface miscAllReduce

#ifdef MPI
  integer:: MPIerror, MPIstatus(MPI_STATUS_SIZE), MPItag
#endif

CONTAINS

!******************************************************************************

SUBROUTINE copyFile( srcNode, srcFile, dstNode, dstFile, writeOption )

  implicit none
  integer,         intent(in):: srcNode      ! Source processor index
                                             ! in range (0:nNodes-1)
  character(len=*),intent(in):: srcFile      ! Source file name
  integer,         intent(in):: dstNode      ! Destination processor index
  character(len=*),intent(in):: dstFile      ! Destination file name
  character(len=*),intent(in):: writeOption  ! ('append'|'overwrite')

! Internal parameters
  character(len=*),parameter:: myName  = 'copyFile '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer, parameter:: minLines =      10000  ! Initial number of file lines
  integer, parameter:: maxLines =     100000  ! Maximum number of file lines
  integer, parameter:: minBuffers =       10  ! Initial number of buffers
  integer, parameter:: maxBuffers =     1000  ! Maximum number of buffers
  integer, parameter:: bufferSize =    10000  ! Size of each buffer, in chars.
  real(dp),parameter:: incrFactor =   1.5_dp  ! Increment factor for re_alloc

! Internal variables and arrays
  integer,pointer:: lineBuf(:)=>null()! Which buffer contains each input line
  integer,pointer:: lineEnd(:)=>null()! Last character of each line in buffer
  character(len=bufferSize),pointer:: buffer(:)=>null()! Buffers to store input
  character(len=bufferSize):: line               ! One input line
  integer:: ib, iBuffer(2), il, iu, lineBegin, lineLength, &
            mBuffers, mLines, nBuffers, nLines, two

! Check whether source and destination are the same
  if (srcNode==dstNode .and. srcFile==dstFile) then
    if (writeOption=='APPEND' .or. writeOption=='append') then
      call die(errHead//'trying to append file to itself')
    else   ! (writeOption=='overwrite') => do nothing
      return
    end if ! (writeOption=='append')
  end if ! (srcNode==dstNode .and. srcFile==dstFile)

! Read source file
  if (myNode == srcNode) then

    ! Allocate buffers and line pointers
    mBuffers = minBuffers
    mLines   = minLines
    call re_alloc( buffer,  1,mBuffers, name=myName//'buffer'  )
    call re_alloc( lineBuf, 0,mLines,   name=myName//'lineBuf' )
    call re_alloc( lineEnd, 0,mLines,   name=myName//'lineEnd' )

    ! Find an available IO unit and open source file
    call io_assign( iu )
    open( unit=iu, file=srcFile, form='formatted', status='unknown' )

    ! Read source lines sequentially and store them in buffers
    lineBuf(0) = 1
    lineEnd(0) = 0
    do il = 1,maxLines

      ! Place source and destination file names in first two lines
      if (il==1) then
        line = srcFile
      else if (il==2) then
        line = dstFile
      else
        ! Read one input line. Exit il loop at end of file.
        read(iu,'(a)',end=1) line
      end if
      nLines = il

      ! Increase array sizes if necessary
      if (nLines > mLines) then
        if (mLines >= maxLines) &
            call die(errHead//'too many input lines')
        mLines = min(int(mLines*incrFactor),maxLines)
        call re_alloc( lineBuf, 0,mLines, name=myName//'lineBuff',copy=.true.)
        call re_alloc( lineEnd, 0,mLines, name=myName//'lineEnd' ,copy=.true.)
      end if ! (mLines >= maxLines)

      ! Find begining and end of line within present buffer
      lineLength  = len(trim(line))
      lineBegin   = lineEnd(il-1) + 1
      lineEnd(il) = lineEnd(il-1) + lineLength
      lineBuf(il) = lineBuf(il-1)

      ! Switch to new buffer if line does not fit within present one
      if (lineEnd(il) > bufferSize) then
        lineBegin   = 1
        lineEnd(il) = lineLength
        lineBuf(il) = lineBuf(il-1) + 1

        ! Increase number of buffers if necessary
        if (lineBuf(il) > mBuffers) then
          if (mBuffers >= maxBuffers) &
            call die(errHead//'parameter maxBuffers too small')
          mBuffers = min(int(mBuffers*incrFactor),maxBuffers)
          call re_alloc( buffer, 1,mBuffers, name=myName//'buffer',copy=.true.)
        end if ! (mBuffers >= maxBuffers)

      end if ! (lineEnd(il) > bufferSize)

      ! Copy input line into buffer
      buffer(lineBuf(il))(lineBegin:lineEnd(il)) = trim(line)
      nBuffers = lineBuf(il)

    end do ! il
    call die(errHead//'too many lines in source file')
1   continue  ! Jump here when reaching end of input file
    call io_close( iu )

  end if ! (myNode == srcNode)

#ifdef MPI
! Transmit data from source to destination nodes, if they are different
  if (myNode==srcNode .and. myNode/=dstNode) then  ! Send data

     iBuffer(1) = nLines
     iBuffer(2) = nBuffers
     call MPI_Send( iBuffer, 2, MPI_Integer, &
                    dstNode, 0, MPI_COMM_WORLD, MPIerror )
     call MPI_Send( lineBuf(0:nLines), nLines+1, MPI_Integer, &
                    dstNode, 0, MPI_COMM_WORLD, MPIerror )
     call MPI_Send( lineEnd(0:nLines), nLines+1, MPI_Integer, &
                    dstNode, 0, MPI_COMM_WORLD, MPIerror )
     call MPI_Send( buffer(1:nBuffers), nBuffers*bufferSize, MPI_Character, &
                    dstNode, 0, MPI_COMM_WORLD, MPIerror )

  else if (myNode/=srcNode .and. myNode==dstNode) then  ! Receive data

! DEBUG
!    print*, 'copyFile: receiving iBuffer from node ', srcNode
! END DEBUG

    ! Receive the sizes of the required arrays
    call MPI_Recv( iBuffer, 2, MPI_Integer, &
                   srcNode, MPItag, MPI_COMM_WORLD, MPIstatus, MPIerror )
    nLines   = iBuffer(1)
    nBuffers = iBuffer(2)

! DEBUG
!    print*, 'copyFile: nLines, nBuffers =', nLines, nBuffers
! END DEBUG

    ! Allocate arrays to receive data
    call re_alloc( buffer,  1,nBuffers, name=myName//'buffer'  )
    call re_alloc( lineBuf, 0,nLines,   name=myName//'lineBuf' )
    call re_alloc( lineEnd, 0,nLines,   name=myName//'lineEnd' )

    ! Receive data
    call MPI_Recv( lineBuf(0:nLines), nLines+1, MPI_Integer, &
                   srcNode, MPItag, MPI_COMM_WORLD, MPIstatus, MPIerror )
    call MPI_Recv( lineEnd(0:nLines), nLines+1, MPI_Integer, &
                   srcNode, MPItag, MPI_COMM_WORLD, MPIstatus, MPIerror )
! DEBUG
!    print'(a,/,3a8,/,(3i8))', &
!      'copyFile:', '    line', ' lineBuf', ' lineEnd', &
!      (il,lineBuf(il),lineEnd(il),il=0,nLines)
! END DEBUG
    call MPI_Recv( buffer(1:nBuffers), nBuffers*bufferSize, MPI_Character, &
                   srcNode, MPItag, MPI_COMM_WORLD, MPIstatus, MPIerror )
! DEBUG
!    do il = 1,nLines
!      print'(i4,2x,a)', il, buffer(1)(lineEnd(il-1)+1:lineEnd(il))
!    end do
! END DEBUG

  end if ! (myNode==srcNode .and. myNode/=dstNode)

  ! All nodes must wait here
  call MPI_Barrier( MPI_Comm_World, MPIerror )
#endif

! Write buffer into destination file
  if (myNode == dstNode) then

! DEBUG
!    if (srcNode==dstNode) then
!      print*, 'copyFile: nLines, nBuffers =', nLines, nBuffers
!      print'(a,/,3a8,/,(3i8))', &
!        'copyFile:', '    line', ' lineBuf', ' lineEnd', &
!        (il,lineBuf(il),lineEnd(il),il=0,nLines)
!      do il = 1,nLines
!        print'(i4,2x,a)', il, buffer(1)(lineEnd(il-1)+1:lineEnd(il))
!      end do
!    end if
! END DEBUG

    ! Find an available IO unit and open destination file
    call io_assign( iu )
    if (writeOption=='append' .or. writeOption=='APPEND') then
      open( unit=iu, file=dstFile, form='formatted', status='unknown', &
            position='append' )
    else if (writeOption=='overwrite' .or. writeOption=='OVERWRITE') then
      open( unit=iu, file=dstFile, form='formatted', status='unknown' )
    else
      call die(errHead//'unknown writeOption = '//writeOption)
    end if ! (writeOption == 'append')

    ! Write lines sequentially
    lineBuf(0) = 0
    do il = 1,nLines

      ! Check if switching to new buffer is required
      if (lineBuf(il)==lineBuf(il-1)) then
        lineBegin = lineEnd(il-1) + 1
      else
        lineBegin = 1
      end if

      ! Copy line from buffer
      lineLength = lineEnd(il) - lineBegin + 1
      line(1:lineLength) = buffer(lineBuf(il))(lineBegin:lineEnd(il))

      ! Check that source and destination files are equal in both nodes
      if (il==1) then
        if (line(1:lineLength) /= trim(srcFile)) &
          call die(errHead//'srcFile mismatch: ' &
                   //line(1:lineLength)//' '//trim(srcFile))
      else if (il==2) then
        if (line(1:lineLength) /= trim(dstFile)) &
          call die(errHead//'dstFile mismatch: ' &
                   //line(1:lineLength)//' '//trim(dstFile))
      else
        ! Copy line from buffer to destination file
        write(iu,'(a)') line(1:lineLength)
      end if

    end do ! il

    call io_close( iu )
  end if ! (myNode == dstNode)

! Deallocations
  if (myNode==srcNode .or. myNode==dstNode) then
    call de_alloc( buffer,  name=myName//'buffer'  )
    call de_alloc( lineBuf, name=myName//'lineBuf' )
    call de_alloc( lineEnd, name=myName//'lineEnd' )
  end if

END SUBROUTINE copyFile

!******************************************************************************

SUBROUTINE miscAllReduceInt( op, a0, b0, c0, d0, e0, f0, &
                                 a1, b1, c1, &
                                 a2, b2, &
                                 a3 )

  implicit none
  character(len=*),intent(in)    :: op  ! ('sum'|'prod'|'max'|'min')
  integer,optional,intent(inout)                 :: a0, b0, c0, d0, e0, f0
  integer,optional,intent(inout),dimension(:)    :: a1, b1, c1
  integer,optional,intent(inout),dimension(:,:)  :: a2, b2
  integer,optional,intent(inout),dimension(:,:,:):: a3

  character(len=*),parameter:: myName  = 'miscAllReduceInt '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: MPIerror, m, n
  integer,pointer:: recvBuff(:)=>null(), sendBuff(:)=>null()

! Do nothing unless execution is parallel
#ifdef MPI

! Find total size of variables and arrays to be gathered
  n = 0
  if (present(a0)) n = n + 1
  if (present(b0)) n = n + 1
  if (present(c0)) n = n + 1
  if (present(d0)) n = n + 1
  if (present(e0)) n = n + 1
  if (present(f0)) n = n + 1
  if (present(a1)) n = n + size(a1)
  if (present(b1)) n = n + size(b1)
  if (present(c1)) n = n + size(c1)
  if (present(a2)) n = n + size(a2)
  if (present(b2)) n = n + size(b2)
  if (present(a3)) n = n + size(a3)

! Allocate buffers to send and receive
  call re_alloc( sendBuff, 1,n, name=myName//'sendBuff' )
  call re_alloc( recvBuff, 1,n, name=myName//'recvBuff' )

! Pack all variables and arrays into sendBuff
  n = 0
  if (present(a0)) then
    sendBuff(n+1) = a0
    n = n + 1
  end if
  if (present(b0)) then
    sendBuff(n+1) = b0
    n = n + 1
  end if
  if (present(c0)) then
    sendBuff(n+1) = c0
    n = n + 1
  end if
  if (present(d0)) then
    sendBuff(n+1) = d0
    n = n + 1
  end if
  if (present(e0)) then
    sendBuff(n+1) = e0
    n = n + 1
  end if
  if (present(f0)) then
    sendBuff(n+1) = f0
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    sendBuff(n+1:n+m) = a1
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    sendBuff(n+1:n+m) = b1
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    sendBuff(n+1:n+m) = c1
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    sendBuff(n+1:n+m) = reshape( a2, (/m/) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    sendBuff(n+1:n+m) = reshape( b2, (/m/) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    sendBuff(n+1:n+m) = reshape( a3, (/m/) )
    n = n + m
  end if

! Reduce sendBuff into recvBuff
  if (op=='sum' .or. op=='Sum' .or. op=='SUM') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Sum, MPI_Comm_World, MPIerror )
  else if (op=='prod' .or. op=='Prod' .or. op=='PROD') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Prod, MPI_Comm_World, MPIerror )
  else if (op=='max' .or. op=='Max' .or. op=='MAX') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Max, MPI_Comm_World, MPIerror )
  else if (op=='min' .or. op=='Min' .or. op=='MIN') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_Integer, &
                        MPI_Min, MPI_Comm_World, MPIerror )
  else
    call die(errHead//'unknown operator')
  end if

! Unpack recvBuff
  n = 0
  if (present(a0)) then
    a0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(b0)) then
    b0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(c0)) then
    c0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(d0)) then
    d0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(e0)) then
    e0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(f0)) then
    f0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    a1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    b1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    c1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    a2 = reshape( recvBuff(n+1:n+m), shape(a2) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    b2 = reshape( recvBuff(n+1:n+m), shape(b2) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    a3 = reshape( recvBuff(n+1:n+m), shape(a3) )
    n = n + m
  end if

  call de_alloc( recvBuff, name=myName//'recvBuff' )
  call de_alloc( sendBuff, name=myName//'sendBuff' )

#endif

END SUBROUTINE miscAllReduceInt

!******************************************************************************

SUBROUTINE miscAllReduceDouble( op, a0, b0, c0, d0, e0, f0, &
                                    a1, b1, c1, &
                                    a2, b2, &
                                    a3 )

  implicit none
  character(len=*),intent(in)    :: op  ! ('sum'|'prod'|'max'|'min')
  real(dp),optional,intent(inout)                 :: a0, b0, c0, d0, e0, f0
  real(dp),optional,intent(inout),dimension(:)    :: a1, b1, c1
  real(dp),optional,intent(inout),dimension(:,:)  :: a2, b2
  real(dp),optional,intent(inout),dimension(:,:,:):: a3

  character(len=*),parameter:: myName  = 'miscAllReduceInt '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer:: MPIerror, m, n
  real(dp),pointer:: recvBuff(:)=>null(), sendBuff(:)=>null()

! Do nothing unless execution is parallel
#ifdef MPI

! Find total size of variables and arrays to be reduced
  n = 0
  if (present(a0)) n = n + 1
  if (present(b0)) n = n + 1
  if (present(c0)) n = n + 1
  if (present(d0)) n = n + 1
  if (present(e0)) n = n + 1
  if (present(f0)) n = n + 1
  if (present(a1)) n = n + size(a1)
  if (present(b1)) n = n + size(b1)
  if (present(c1)) n = n + size(c1)
  if (present(a2)) n = n + size(a2)
  if (present(b2)) n = n + size(b2)
  if (present(a3)) n = n + size(a3)

! Allocate buffers to send and receive
  call re_alloc( sendBuff, 1,n, name=myName//'sendBuff' )
  call re_alloc( recvBuff, 1,n, name=myName//'recvBuff' )

! Pack all variables and arrays into sendBuff
  n = 0
  if (present(a0)) then
    sendBuff(n+1) = a0
    n = n + 1
  end if
  if (present(b0)) then
    sendBuff(n+1) = b0
    n = n + 1
  end if
  if (present(c0)) then
    sendBuff(n+1) = c0
    n = n + 1
  end if
  if (present(d0)) then
    sendBuff(n+1) = d0
    n = n + 1
  end if
  if (present(e0)) then
    sendBuff(n+1) = e0
    n = n + 1
  end if
  if (present(f0)) then
    sendBuff(n+1) = f0
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    sendBuff(n+1:n+m) = a1
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    sendBuff(n+1:n+m) = b1
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    sendBuff(n+1:n+m) = c1
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    sendBuff(n+1:n+m) = reshape( a2, (/m/) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    sendBuff(n+1:n+m) = reshape( b2, (/m/) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    sendBuff(n+1:n+m) = reshape( a3, (/m/) )
    n = n + m
  end if

! Reduce sendBuff into recvBuff
  if (op=='sum' .or. op=='Sum' .or. op=='SUM') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Sum, MPI_Comm_World, MPIerror )
  else if (op=='prod' .or. op=='Prod' .or. op=='PROD') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Prod, MPI_Comm_World, MPIerror )
  else if (op=='max' .or. op=='Max' .or. op=='MAX') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Max, MPI_Comm_World, MPIerror )
  else if (op=='min' .or. op=='Min' .or. op=='MIN') then
    call MPI_AllReduce( sendBuff, recvBuff, n, MPI_double_precision, &
                        MPI_Min, MPI_Comm_World, MPIerror )
  else
    call die(errHead//'unknown operator')
  end if

! Unpack recvBuff
  n = 0
  if (present(a0)) then
    a0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(b0)) then
    b0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(c0)) then
    c0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(d0)) then
    d0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(e0)) then
    e0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(f0)) then
    f0 = recvBuff(n+1)
    n = n + 1
  end if
  if (present(a1)) then
    m = size(a1)
    a1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(b1)) then
    m = size(b1)
    b1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(c1)) then
    m = size(c1)
    c1 = recvBuff(n+1:n+m)
    n = n + m
  end if
  if (present(a2)) then
    m = size(a2)
    a2 = reshape( recvBuff(n+1:n+m), shape(a2) )
    n = n + m
  end if
  if (present(b2)) then
    m = size(b2)
    b2 = reshape( recvBuff(n+1:n+m), shape(b2) )
    n = n + m
  end if
  if (present(a3)) then
    m = size(a3)
    a3 = reshape( recvBuff(n+1:n+m), shape(a3) )
    n = n + m
  end if

  call de_alloc( recvBuff, name=myName//'recvBuff' )
  call de_alloc( sendBuff, name=myName//'sendBuff' )

#endif

END SUBROUTINE miscAllReduceDouble

!******************************************************************************

character(len=6) function nodeString( node )

! Returns a string with node index, or blank if nNodes<2

  implicit none
  integer,optional:: node  ! a valid node index 0 <= node <= nNodes-1

  character(len=20):: numName
  character(len=80):: errMsg
  integer:: fileLen, iNode, maxLen, numLen

  if (totNodes<2) then   ! Serial execution
    nodeString = ' '
  else                   ! Parallel execution
    if (present(node)) then
      if (0<=node .and. node<totNodes) then ! Valid node index
        iNode = node
      else                                  ! Invalid node index
        write(errMsg,*) 'nodeString: ERROR: invalid argument: node =', node
        call die( trim(errMsg) )
      end if ! (0<=node<totNodes)
    else ! (.not.present(node))
      iNode = myNode
    end if ! (present(node))

    ! Find maximum name length of node numbers
    write(numName,*) totNodes-1   ! 0 <= node <= Nodes-1
    numName = adjustl(numName)
    maxLen = len_trim(numName)

    ! Find name of this node's number, and its name length
    write(numName,*) iNode
    numName = adjustl(numName)
    numLen = len_trim(numName)

    ! Set node number string, e.g. 00, 01, ... 63
    nodeString = '000000000'
    nodeString(maxLen-numLen+1:maxLen) = trim(numName)
    nodeString(maxLen+1:) = ' '

  end if ! (totNodes<2)
  
end function nodeString

END MODULE moreParallelSubs
