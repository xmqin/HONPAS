! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

C This file contains module scheComm. It tries to create an
C optimal scheduling of communications.
C
C Written by Rogeli Grima (BSC) Dec.2007
C
C Includes the following subroutine:
C
C scheduleComm     = Schedule communications
C
C Includes the following data type:
C
C COMM_T = Stores all the communications
C 
      module scheComm

      TYPE COMM_T
      integer            :: np        ! Number of processors
      integer            :: ncol      ! Number of colours/steps
      integer, pointer   :: ind(:,:)  ! Communications of every process
      END TYPE COMM_T

      contains

C ==================================================================
C Schedule the communications trying to minimize the communication
C time. The idea is that in every time step we should make as many
C communications as possible (In a time step one process only can do
C one send or one receive). In every time step we should give
C preference to those processes with more remaining communications.
C ==================================================================
C SUBROUTINE scheduleComm( Ncom, src, dst, comm )
C
C INPUT:
C integer      Ncom             : Total number of communications
C integer      src              : Process ID of the source.
C integer      dst              : Process ID of the destiny.
C
C OUTPUT:
C type(COMM_T) comm             : Contains the communications of every process
C
C BEHAVIOR:
C We want to create an ordered list of communications for every process.
C
C We have a graph of communications, where every node is a process and
C every segment is a communication between two processes. First, we build
C its adjacency graph, where every node is a communication between two
C processes and every segment means that two communications involve the
C same process (they can't happen at the same time).
C
C In order to minimize the number of communication steps we use the
C graph coloring technique. We reorder the nodes of the adjacency graph
C according to its number of neightbours. We try to assign the smaller
C color to every node. Two neightbours can't share the same color.
C Every time that we assign a color to a communication, we fill the
C comm struct.
C
C ==================================================================
      subroutine scheduleComm( Ncom, src, dst, comm )
      use alloc,     only : re_alloc, de_alloc
      implicit none
C     Input variables
      integer              :: Ncom, src(Ncom), dst(Ncom)
      type(COMM_T)         :: comm
C     Local variables
      integer              :: i, j, com1, com2, P1, P2, nadj, Ncol, col
      integer,     pointer :: neights(:), xadj(:), adj(:),
     &                        perm(:), color(:)
      logical              :: found
!------------------------------------------------------------------------- BEGIN
      nullify( neights, xadj, adj, perm, color )
      call re_alloc( neights, 1, NCom, 'neights', 'scheComm' )
      neights = 0

C     Build the adjacency graph. In this new graph every node
C     means a communication. Two nodes are neighbours if they involve
C     the same node

C     neights(com) = number of neightbours of communication 'com'
      do com1= 1, Ncom
        P1 = src(com1)
        P2 = dst(com1)
        do com2= com1+1, Ncom
          if (P1.eq.src(com2) .or. P1.eq.dst(com2) .or.
     &        P2.eq.src(com2) .or. P2.eq.dst(com2) ) then
            neights(com1) = neights(com1) + 1
            neights(com2) = neights(com2) + 1
          endif
        enddo
      enddo

C     nadj= size of the adjacency graph
      call re_alloc( xadj, 1, NCom+1, 'xadj', 'scheComm' )
      nadj = 1
      do com1= 1, Ncom
        xadj(com1) = nadj
        nadj       = nadj + neights(com1)
      enddo
      xadj(Ncom+1) = nadj
      nadj         = nadj - 1

C     Fill the adjacency graph
      call re_alloc( adj, 1, nadj, 'adj', 'scheComm' )
      do com1= 1, Ncom
        P1 = src(com1)
        P2 = dst(com1)
        do com2= com1+1, Ncom
          if (P1.eq.src(com2) .or. P1.eq.dst(com2) .or.
     &        P2.eq.src(com2) .or. P2.eq.dst(com2) ) then
            adj(xadj(com1)) = com2
            adj(xadj(com2)) = com1
            xadj(com1) = xadj(com1) + 1
            xadj(com2) = xadj(com2) + 1
          endif
        enddo
      enddo

      do com1= 1, Ncom
        xadj(com1) = xadj(com1) - neights(com1)
      enddo

C     GRAPH COLORING: in order to optimize the communications
C     we assign a color to every node. We want to minimize the
C     number of colors. Two neightbours can't share the same
C     color.
      call re_alloc( perm, 1, NCom, 'perm', 'scheComm' )
      do com1= 1, Ncom
        perm(com1) = com1
      enddo

C     We order the nodes by its number of neightbours.
      call myQsort( Ncom, neights, perm )

      call re_alloc( color, 1, NCom, 'color', 'scheComm' )
      do com1= 1, Ncom
        color(com1) = 0
      enddo

C     The first node to get a color is the one with more neghtbours
      com1        = perm(1)
      color(com1) = 1
      Ncol        = 1
      do i= 2, Ncom
        com1 = perm(i)
C       For every node, try to find if we can assing one of the previous
C       colors.
        do col= 1, Ncol
          found = .false.
          do j= xadj(com1), xadj(com1+1)-1
            com2 = adj(j)
            if (color(com2).eq.col) then
C             If one of my neightbours have this color. I can't
C             use it.
              found = .true.
              exit
            endif
          enddo
          if (.not.found) then
C           None of my neightbours have this color. I can
C           use it.
            color(com1) = col
            exit
          endif
        enddo
        if (found) then
C         All the existings colors are used by any of my neightbours.
C         Create a new color
          Ncol        = Ncol + 1
          color(com1) = Ncol
        endif
      enddo

C     Fill the communication structure with the new scheduling
      comm%ncol = Ncol
      nullify( comm%ind )
      call re_alloc( comm%ind, 1, NCol, 1, comm%np,
     &              'comm%ind', 'scheComm' )

      comm%ind = 0
      do com1= 1, Ncom
        col = color(com1)
        P1  = src(com1)
        P2  = dst(com1)
        comm%ind(col,P1) = com1
        comm%ind(col,P2) = com1
      enddo

      call de_alloc(   color,   'color', 'scheComm' )
      call de_alloc(    perm,    'perm', 'scheComm' )
      call de_alloc(     adj,     'adj', 'scheComm' )
      call de_alloc(    xadj,    'xadj', 'scheComm' )
      call de_alloc( neights, 'neights', 'scheComm' )
!--------------------------------------------------------------------------- END
      end subroutine scheduleComm
      end module scheComm
