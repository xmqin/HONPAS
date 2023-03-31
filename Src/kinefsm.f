! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_kinefsm
      public :: kinefsm
      CONTAINS
      subroutine kinefsm(nua, na, no, scell, xa, indxua, rmaxo,
     .                  maxnh, maxnd, lasto, iphorb, isa, 
     .                  numd, listdptr, listd, numh, listhptr, listh, 
     .                  nspin, Dscf, Ekin, fa, stress, H,
     .                  matrix_elements_only )
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
C                            If computing only matrix elements, it
C                            can be set to 1.
C logical matrix_elements_only:
C integer Dscf(maxnd,nspin): Density matrix. Not touched if computing
C                            only matrix elements.
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
      use atmfuncs,      only : rcut, orb_gindex
      use neighbour,     only : jna=>jan, r2ij, xij, mneighb,
     &                          reset_neighbour_arrays 
      use m_new_matel,   only : new_matel
      use alloc,         only : re_alloc, de_alloc

      implicit none

      integer ::  maxnd, maxnh, na, no, nspin, nua

      integer
     .  indxua(na), iphorb(no), isa(na), lasto(0:na), 
     .  listd(maxnd), listh(maxnh), numd(*), numh(*), listdptr(*),
     .  listhptr(*)

      real(dp)
     .  scell(3,3), Dscf(maxnd,nspin), Ekin, 
     .  fa(3,nua), H(maxnh,nspin), rmaxo, 
     .  stress(3,3), xa(3,na)
      logical, intent(in)  :: matrix_elements_only

C Internal variables ..................................................
  
      integer
     .  ia, ind, io, iio, ioa, is, ispin, ix, ig, jg,
     .  j, ja, jn, jo, joa, js, jua, jx, nnia

      real(dp)
     .  fij(3), grTij(3) , rij, Tij, volcel, volume

      real(dp), dimension(:), pointer :: Di, Ti

      external ::  volcel, timer
C ......................

C Start timer
      call timer( 'kinefsm', 1 )

C Allocate local memory
      if (.not. matrix_elements_only) then
         nullify( Di )
         call re_alloc( Di, 1, no, 'Di', 'kinefsm' )
         Di(1:no) = 0.0_dp
      endif
      nullify( Ti )
      call re_alloc( Ti, 1, no, 'Ti', 'kinefsm' )
      Ti(1:no) = 0.0_dp

      volume = nua * volcel(scell) / na

      call mneighb( scell, 2.0d0*rmaxo, na, xa, 0, 0, nnia )

      Ekin = 0.0_dp

      do ia = 1,nua
        is = isa(ia)
        call mneighb( scell, 2.0d0*rmaxo, na, xa, ia, 0, nnia )
        do io = lasto(ia-1)+1,lasto(ia)

C Is this orbital on this Node?
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio.gt.0) then  ! Local orbital
            ioa = iphorb(io)
            ig = orb_gindex(is,ioa)

            if (.not. matrix_elements_only) then
               do j = 1,numd(iio)
                  ind = listdptr(iio) + j
                  jo = listd(ind)
                  do ispin = 1,nspin
                     Di(jo) = Di(jo) + Dscf(ind,ispin)
                  enddo
               enddo
            endif
            do jn = 1,nnia
              ja = jna(jn)
              jua = indxua(ja)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)
                jg = orb_gindex(js,joa)
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then
                  call new_MATEL( 'T', ig, jg, xij(1:3,jn),
     .                      Tij, grTij )
                  Ti(jo) = Ti(jo) + Tij

                  if (.not. matrix_elements_only) then
                     Ekin = Ekin + Di(jo) * Tij
                    do ix = 1,3
                      fij(ix) = Di(jo) * grTij(ix)
                      fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                      fa(ix,jua) = fa(ix,jua) - fij(ix)
                      do jx = 1,3
                        stress(jx,ix) = stress(jx,ix) +
     .                                  xij(jx,jn) * fij(ix) / volume
                      enddo
                    enddo
                  endif 
                endif

              enddo
            enddo
            if (.not. matrix_elements_only) then
               do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  Di(jo) = 0.0_dp
               enddo
            endif
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
!      call new_MATEL( 'T', 0, 0, 0, 0, xij, Tij, grTij )
      call reset_neighbour_arrays( )
      call de_alloc( Ti, 'Ti', 'kinefsm' )
      if (.not. matrix_elements_only) then
         call de_alloc( Di, 'Di', 'kinefsm' )
      endif

C Finish timer
      call timer( 'kinefsm', 2 )

      end subroutine kinefsm
      end module m_kinefsm
