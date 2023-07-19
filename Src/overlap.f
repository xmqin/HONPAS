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
      module m_overlap

      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb
      use alloc,         only : re_alloc, de_alloc

      implicit none

      public :: overlap
      private

      CONTAINS

      subroutine overlap(nua, na, no, scell, xa, indxua, rmaxo,
     .                   maxnh, lasto, iphorb, isa, 
     .                   numh, listhptr, listh, S)
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

      integer, intent(in) ::  maxnh, na, no, nua

      integer, intent(in) ::
     . indxua(na), iphorb(no), isa(na), lasto(0:na), 
     . listh(maxnh), numh(*), listhptr(*)
      real(dp) , intent(in)     :: scell(3,3), rmaxo, xa(3,na)
      real(dp), intent(out)     :: S(maxnh)

C Internal variables ......................................................
  
      integer
     .  ia, ind, io, ioa, is,  iio, 
     .  j, ja, jn, jo, joa, js, jx, nnia

      real(dp) grSij(3) , rij, Sij, volcel, volume

      real(dp), dimension(:), pointer  :: Si

      external  timer
C ......................

C Start timer
      call timer( 'overlap', 1 )

      volume = nua * volcel(scell) / na

C Initialize neighb subroutine 
      call mneighb( scell, 2.d0*rmaxo, na, xa, 0, 0, nnia )

C Allocate local memory

      nullify(Si)
      call re_alloc(Si,1,no,name="Si",routine="overlap")

      do ia = 1,nua
        call mneighb( scell, 2.d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C Valid orbital 

            do jn = 1,nnia
              ja = jna(jn)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                ioa = iphorb(io)
                joa = iphorb(jo)
                is = isa(ia)
                js = isa(ja)
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                  call matel( 'S', is, js, ioa, joa, xij(1,jn),
     .                      Sij, grSij )
                  Si(jo) = Si(jo) + Sij
                endif
              enddo
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio)+j
              jo = listh(ind)
              S(ind) = Si(jo)
              Si(jo) = 0.0d0
!           if(node.eq.0) write(98,*) io, jo, ind
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
      call de_alloc(Si,name="Si",routine="overlap")

C Finish timer
      call timer( 'overlap', 2 )

      end subroutine overlap
      end module m_overlap




