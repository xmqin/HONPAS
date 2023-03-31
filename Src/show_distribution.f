! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      subroutine show_distribution()
!     This is buggy for dd
      use parallel, only: Nodes
      use parallelsubs
      use siesta_geom, only: na_u
      use atomlist, only: lasto
      use sys, only: die

      implicit none

      integer  :: Node, nuotot, iu, localorbs
      logical  :: abort_flag

      nuotot = lasto(na_u)

      call io_assign(iu)
      open(unit=iu,file="PARALLEL_DIST",form="formatted",
     $     status="replace",
     $     action="write",position="rewind")

      abort_flag = .false.
      do Node = 0, Nodes-1
         call getNodeOrbs(nuotot,Node,Nodes,localorbs)
         write(iu,*) "Node ", Node, " handles ", localorbs, " orbitals."
         if (localorbs == 0) abort_flag = .true.
!         do io = 1, localorbs
!            call LocalToGlobalOrb(io,Node,Nodes,gio)
!            write(iu,*) gio
!         enddo
      enddo
      call io_close(iu)
      if (abort_flag) then
          write(6,"(a)") "Some processors are idle. Check PARALLEL_DIST"
          write(6,"(a)")
     $         "You have too many processors for the system size !!!"
          call die()
       endif

      end subroutine show_distribution

