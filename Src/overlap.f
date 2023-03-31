! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_overlap

      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, orb_gindex
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb,
     &                          reset_neighbour_arrays
      use alloc,         only : re_alloc, de_alloc
      use m_new_matel,   only : new_matel
      use m_iodm_old,    only : write_dm
      use m_matio,       only : write_mat
      use atomlist, only: no_l
      use fdf

      implicit none

      public :: overlap
      private

      CONTAINS

      subroutine overlap(nua, na, no, scell, xa, indxua, rmaxo,
     &                   maxnh, lasto, iphorb, isa, 
     &                   numh, listhptr, listh, S)
C *********************************************************************
C Computes the overlap matrix
C Energies in Ry. Lengths in Bohr.
C Based on overfsm, but removing all references to Escf/Emat  (AG)
C Overfsm Writen by J.Soler and P.Ordejon, August-October'96, June'98
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer no               : Number of orbitals in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxnh            : First dimension of S and listh
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numh(nuotot)     : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuotot) : Pointer to start of rows (-1) of overlap
C                            matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix

      integer, intent(in)   ::  maxnh, na, no, nua
      integer, intent(in)   :: indxua(na), iphorb(no), isa(na),
     &                         lasto(0:na), listh(maxnh), numh(*),
     &                         listhptr(*)
      real(dp) , intent(in) :: scell(3,3), rmaxo, xa(3,na)
      real(dp), intent(out) :: S(maxnh)
C Internal variables ......................................................
      integer               :: ia, ind, io, ioa, is,  iio, j, ja, jn,
     &                         jo, joa, js, jua, nnia, ig, jg
      real(dp)              :: grSij(3) , rij, Sij
      real(dp),     pointer :: Si(:)
      external  timer

C     Start timer
      call timer( 'overlap', 1 )

C     Initialize neighb subroutine 
      call mneighb( scell, 2.d0*rmaxo, na, xa, 0, 0, nnia )

C     Allocate local memory
      nullify(Si)
      call re_alloc( Si, 1, no, 'Si', 'overlap' )

      do ia = 1,nua
        is = isa(ia)
        call mneighb( scell, 2.d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C         Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C           Valid orbital
            ioa = iphorb(io)
            ig = orb_gindex(is,ioa)
            do jn = 1,nnia
              ja = jna(jn)
              jua = indxua(ja)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)
                !
                ! Use global indexes for new version of matel
                !
                jg = orb_gindex(js,joa)
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                  call new_MATEL( 'S', ig, jg, xij(1:3,jn),
     &                        Sij, grSij )
                  Si(jo) = Si(jo) + Sij
                endif
              enddo
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio)+j
              jo = listh(ind)
              S(ind) = Si(jo)
              Si(jo) = 0.0d0
            enddo
          endif
        enddo
      enddo

C     Deallocate local memory
!      call new_MATEL( 'S', 0, 0, xij, Sij, grSij )
      call reset_neighbour_arrays( )
      call de_alloc( Si, 'Si', 'overlap' )

      if (fdf_get("Sonly",.false.)) then
         call write_dm(maxnh, no_l, 1,
     &               numh, listhptr, listh, S,
     $               userfile="SOLD")

         call write_mat(maxnh, no_l, 1,
     &               numh, listhptr, listh, S,
     $               userfile="SMAT")

         call timer("fastWriteMat",1)
         call write_mat(maxnh, no_l, 1,
     &               numh, listhptr, listh, S,
     $               userfile="SMATBS",compatible=.false.)
         call timer("fastWriteMat",2)
      endif


C     Finish timer
      call timer( 'overlap', 2 )
      end subroutine overlap
      end module m_overlap
