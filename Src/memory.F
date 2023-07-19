! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine memory( Task, Type, NElements, CallingRoutine)
C 
C This subroutine keeps track of information relating to the use 
C of dynamic memory
C
C Input :
C
C character*1 Task  : job type = 'A' -> allocate
C                   :            'D' -> deallocate
C character*1 Type  : type of variable = 'I' = integer
C                   :                    'S' = single precision real
C                   :                    'D' = double precision real
C                   :                    'X' = grid precision real
C                   :                    'L' = logical
C                   :                    'C' = single precision complex
C                   :                    'Z' = double precision complex
C                   :                    'S' = character data (we assume takes one word)
C integer NElements : number of array elements being 
C                   : allocated/deallocated
C character         :
C   CallingRoutine  : string containing the name of the calling routine
C
C Created by J.D. Gale, October 1999
C
      use memoryinfo
      use parallel, only: Node

      implicit none

      integer, intent(in)                 :: NElements
      character(len=1), intent(in)        :: Task, Type
      character(len=*), intent(in)        :: CallingRoutine

C Local variables
      integer    Sign
CAG
!      if (Node.eq.0) then
!        write(66,'(a,a,a,a,i10)') "---memory: ", CallingRoutine // " ",
!     $                          task // " ",
!     $                          type // " ", nelements
!        call flush(66)
!      endif

C Assign sign of operation based on task
      if (Task .eq. 'A' .or. Task .eq. 'a') then
        Sign = +1
      else
        Sign = -1
      endif

C Add number of elements to appropriate pointer
      if (Type .eq. 'I' .or. Type .eq. 'i') then
        WordsInteger = WordsInteger + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(1)
      elseif (Type .eq. 'S' .or. Type .eq. 's') then
        WordsSP = WordsSP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(2)
      elseif (Type .eq. 'D' .or. Type .eq. 'd') then
        WordsDP = WordsDP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(3)
!
!     Support for flexible precision in grid-related arrays
!
      elseif (Type .eq. 'X' .or. Type .eq. 'x') then
#ifdef GRID_DP
        WordsDP = WordsDP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(3)
#else
        WordsSP = WordsSP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(2)
#endif

      elseif (Type .eq. 'L' .or. Type .eq. 'l') then
        WordsLogical = WordsLogical + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(4)
      elseif (Type .eq. 'C' .or. Type .eq. 'c') then
        WordsSC = WordsSC + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(5)
      elseif (Type .eq. 'Z' .or. Type .eq. 'z') then
        WordsDC = WordsDC + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(6)
      elseif (Type .eq. 'S' .or. Type .eq. 's') then
        WordsString = WordsString + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(7)
      endif

C Check whether memory is greater than peak memory so far
      if (CurrentMemory .gt. PeakMemory) then
        PeakMemory = CurrentMemory
        PeakRoutine = CallingRoutine
      endif

      end subroutine memory

      subroutine printmemory( Unit, Level )
C
C Outputs the information about the dynamic memory useage
C
C Created by J.D. Gale, October 1999
C
C Input :
C
C integer Unit  : channel number for output
C integer Level : number controlling the level of output
C               : 0 = peak memory details only
C               : 1 = peak memory on all Nodes
C
      use memoryinfo
      use precision
      use parallel,  only : Node, Nodes
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C Passed arguments
      integer
     .  Unit, Level

C Local arguments
      integer
     .  MaxMemory, n
      integer, allocatable :: maxmemlist(:)

#ifdef MPI
      integer
     .  MPIerror
#endif

C Find peak memory use over all Nodes
#ifdef MPI
      call MPI_AllReduce(PeakMemory,MaxMemory,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)
#else
      MaxMemory = PeakMemory
#endif

C Output memory use information
      if (Node.eq.0) then
        write(Unit,'(a)') ''
      endif
      if (Level.eq.0) then
        if (Node.eq.0) then
          write(Unit,'(''* Maximum dynamic memory allocated = '',i5,
     .      '' MB'')') (MaxMemory/1000000)+1
        endif
      else
        if (node==0) then
          allocate(maxmemlist(nodes))
        else
          allocate(maxmemlist(0))
        endif
#ifdef MPI
        call mpi_gather(peakmemory,1,mpi_integer,
     .                  maxmemlist,1,mpi_integer,
     .                  0,mpi_comm_world,mpierror)
#else
        maxmemlist=peakmemory
#endif
        if (node==0) then
          do n=0,Nodes-1
            write(Unit,'(''* Maximum dynamic memory allocated : Node ''
     .          ,i4,'' = '',i5,'' MB'')') n,(maxmemlist(n+1)/1000000)+1
          enddo
          write(Unit,'(a)') ''
          write(Unit,'(''* Maximum memory occured during '',a30)')
     .               PeakRoutine
        endif
        deallocate(maxmemlist)
      endif

      end subroutine printmemory
