! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE printMat
      CONTAINS
      subroutine printMatrix( nf, nc, ia, ja, filename, title )
      use alloc, only : re_alloc, de_alloc
      implicit none
      integer          :: nf, nc, ia(*), ja(*)
      character*(*)    :: filename, title
      integer, pointer :: iwa(:)
      integer          :: fd

      call re_alloc( iwa, 1, nc+1, name='iwa', routine='printmatrix' )

      call io_assign( fd )
      open( unit=fd, file=filename, status='unknown' )

      call pspltm1( nf, nc, 0, ia, ja, title, 0, 20.0, 'cm', 0, 0,
     &              iwa, 1, fd )

      call io_close( fd )
      call de_alloc( iwa, name='iwa', routine='printmatrix' )
      end subroutine printMatrix

      subroutine printPermMatrix( nuo, ia, ja, perm, invp,
     &                            filename, title )
      use alloc, only : re_alloc, de_alloc
      implicit none
      integer          :: nuo, ia(*), ja(0:nuo), perm(nuo), invp(nuo)
      character*(*)    :: filename, title

      integer, pointer :: adj(:), xadj(:), auxx(:)
      integer          :: nnz, nnzL, io, jo, ii, jj

      nnz = ja(nuo)-1
      call re_alloc( adj, 1, nnz, name='adj',
     &               routine='printPermMatrix' )
      call re_alloc( xadj, 0, nuo, name='xadj',
     &               routine='printPermMatrix' )
      call re_alloc( auxx, 1, nuo, name='auxx',
     &               routine='printPermMatrix' )
      auxx = 0
      ii   = 1
      do io= 1, nuo
        jo         = invp(io)
        xadj(io-1) = ii
        do jj= ja(jo-1), ja(jo)-1
          adj(ii) = perm(ia(jj))
          ii      = ii + 1
        enddo
        nnzL = ja(jo) - ja(jo-1)
        call myQsort( nnzL, adj(xadj(io-1)), auxx )
      enddo
      xadj(nuo) = ii
      call printMatrix( nuo, nuo, adj, xadj, filename, title )

      call de_alloc( auxx, name='auxx', routine='printPermMatrix' )
      call de_alloc( xadj, name='xadj', routine='printPermMatrix' )
      call de_alloc(  adj, name= 'adj', routine='printPermMatrix' )
      end subroutine printPermMatrix
      END MODULE printMat
