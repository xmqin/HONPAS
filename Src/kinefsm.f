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
      subroutine kinefsm(nua, na, no, scell, xa, indxua, rmaxo, maxo,
     .                  maxnh, maxnd, lasto, iphorb, isa, 
     .                  numd, listdptr, listd, numh, listhptr, listh, 
     .                  nspin, Dscf, Ekin, fa, stress, H )
C *********************************************************************
C Kinetic contribution to energy, forces, stress and matrix elements.
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
C integer maxo             : Second dimension of Dscf and H
C integer maxnh            : First dimension of H and listh
C integer maxnd            : First dimension of Dscf and listd
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numd(nuo)        : Number of nonzero elements of each row
C                            of the density matrix
C integer listdptr(nuo)    : Pointer to the start of rows (-1) of
C                            the density matrix
C integer listd(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the density matrix
C integer numh(nuo)        : Number of nonzero elements of each row
C                            of the hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of rows (-1) of
C                            the hamiltonian matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the hamiltonian matrix
C integer nspin            : Number of spin components of Dscf and H
C integer Dscf(maxnd,nspin): Density matrix
C **************************** OUTPUT *********************************
C real*8 Ekin              : Kinetic energy in unit cell
C ********************** INPUT and OUTPUT *****************************
C real*8 fa(3,nua)         : Atomic forces (Kinetic part added to input)
C real*8 stress(3,3)       : Stress tensor (Kinetic part added to input)
C real*8 H(maxnh,nspin)    : Hamiltonian (Kinetic part added to input)
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb
      use alloc,         only : re_alloc, de_alloc

      implicit none

      integer ::  maxnd, maxnh, maxo, na, no, nspin, nua

      integer
     .  indxua(na), iphorb(no), isa(na), lasto(0:na), 
     .  listd(maxnd), listh(maxnh), numd(*), numh(*), listdptr(*),
     .  listhptr(*)

      real(dp)
     .  scell(3,3), Dscf(maxnd,nspin), Ekin, 
     .  fa(3,nua), H(maxnh,nspin), rmaxo, 
     .  stress(3,3), xa(3,na)

C Internal variables ..................................................
  
      integer
     .  ia, ind, io, iio, ioa, is, ispin, ix, 
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      real(dp)
     .  fij(3), grTij(3) , rij, Tij, volcel, volume

      real(dp), dimension(:), pointer :: Di, Ti

      external ::  volcel, timer
C ......................

C Start timer
      call timer( 'kinefsm', 1 )

C Allocate local memory
      nullify( Di )
      call re_alloc( Di, 1, no, name='Di', routine='kinefsm' )
      nullify( Ti )
      call re_alloc( Ti, 1, no, name='Ti', routine='kinefsm' )

      volume = nua * volcel(scell) / na

      call mneighb( scell, 2.0d0*rmaxo, na, xa, 0, 0, nnia )

      Ekin = 0.0_dp
      Di(1:no) = 0.0_dp
      Ti(1:no) = 0.0_dp

      do ia = 1,nua
        call mneighb( scell, 2.0d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then
            ioa = iphorb(io)
            is = isa(ia)

C Valid orbital 
            do j = 1,numd(iio)
              ind = listdptr(iio) + j
              jo = listd(ind)
              do ispin = 1,nspin
                Di(jo) = Di(jo) + Dscf(ind,ispin)
              enddo
            enddo
            do jn = 1,nnia
              ja = jna(jn)
              jua = indxua(ja)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                  call matel( 'T', is, js, ioa, joa, xij(1,jn),
     .                      Tij, grTij )
                  Ti(jo) = Ti(jo) + Tij
                  Ekin = Ekin + Di(jo) * Tij
                  do ix = 1,3
                    fij(ix) = Di(jo) * grTij(ix)
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
              Di(jo) = 0.0_dp
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio)+j
              jo = listh(ind)
              do ispin = 1,nspin
                H(ind,ispin) = H(ind,ispin) + Ti(jo)
              enddo
              Ti(jo) = 0.0_dp
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
      call de_alloc( Ti, name='Ti' )
      call de_alloc( Di, name='Di' )

C Finish timer
      call timer( 'kinefsm', 2 )

      end
