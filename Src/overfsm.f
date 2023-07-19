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
      module m_overfsm

      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb
      use alloc,         only : re_alloc, de_alloc

      implicit none

      public :: overfsm
      private

      CONTAINS

      subroutine overfsm(nua, na, no, scell, xa, indxua, rmaxo, maxo,
     .                   maxnh, maxnd, lasto, iphorb, isa, 
     .                   numd, listdptr, listd, numh, listhptr, listh, 
     .                   nspin, Emat, fa, stress, S)
C *********************************************************************
C Overlap matrix and orthonormalization contribution to forces and stress.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, August-October'96, June'98
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer no               : Number of orbitals in supercell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
c integer indxua(na)       : Index of equivalent atom in unit cell
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxo             : Second dimension of Emat
C integer maxnh            : First dimension of S and listh
C integer maxnd            : First dimension of Escf and listd
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numd(nuotot)     : Number of nonzero elements of each row
C                            of Emat
C integer listdptr(nuotot) : Pointer to start of rows (-1) of Emat
C integer listd(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of Emat
C integer numh(nuotot)     : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuotot) : Pointer to start of rows (-1) of overlap
C                            matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C integer nspin            : Number of spin components of Emat
C integer Emat(maxnd,nspin): Energy-Density matrix
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix
C ********************** INPUT and OUTPUT *****************************
C real*8  fa(3,nua)        : Atomic forces (Orthog. part added to input)
C real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
C *********************************************************************

      integer, intent(in) ::
     . maxnd, maxnh, maxo, na, no, nspin, nua

      integer, intent(in) ::
     . indxua(na), iphorb(no), isa(na), lasto(0:na), 
     . listd(*), numd(*), listh(maxnh), numh(*), listdptr(*),
     . listhptr(*)

      real(dp) , intent(inout)  :: fa(3,nua), stress(3,3)
      real(dp) , intent(in)     :: scell(3,3), Emat(maxnd,nspin), 
     .                 rmaxo, xa(3,na)
      real(dp), intent(out)     :: S(maxnh)

C Internal variables ......................................................
  
      integer
     .  ia, ind, io, ioa, is, ispin, ix, iio, 
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      real(dp)
     .  fij(3), grSij(3) , rij, Sij, volcel, volume

      real(dp), dimension(:), pointer  ::  Di, Si

      external  timer
C ......................

C Start timer
      call timer( 'overfsm', 1 )

      volume = nua * volcel(scell) / na

C Initialize neighb subroutine 
      call mneighb( scell, 2.d0*rmaxo, na, xa, 0, 0, nnia)

C Allocate local memory
      nullify( Di )
      call re_alloc( Di, 1, no, name='Di', routine='overfsm' )
      nullify( Si )
      call re_alloc( Si, 1, no, name='Si', routine='overfsm' )

      Si(1:no) = 0.0d0  ! AG: Superfluous
      Di(1:no) = 0.0d0

      do ia = 1,nua
        call mneighb( scell, 2.d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then

C Valid orbital 
            do j = 1,numd(iio)
              ind = listdptr(iio)+j
              jo = listd(ind)
              do ispin = 1,nspin  
                Di(jo) = Di(jo) + Emat(ind,ispin)
              enddo
            enddo
            do jn = 1,nnia
              ja = jna(jn)
              jua = indxua(ja)
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
                  do ix = 1,3
                    fij(ix) = (- Di(jo)) * grSij(ix)
                    fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                    fa(ix,jua) = fa(ix,jua) - fij(ix)
                    do jx = 1,3
                      stress(jx,ix) = stress(jx,ix) +
     .                              xij(jx,jn) * fij(ix) / volume
                    enddo
                  enddo
                endif
              enddo
            enddo
            do j = 1,numd(iio)
              jo = listd(listdptr(iio)+j)
              Di(jo) = 0.0d0
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

C Deallocate local memory

      call de_alloc( Si, name='Si' )
      call de_alloc( Di, name='Di' )

C Finish timer
      call timer( 'overfsm', 2 )

      end subroutine overfsm
      end module m_overfsm




