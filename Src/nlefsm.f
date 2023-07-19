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
      module m_nlefsm

      implicit none

      public :: nlefsm

      private

      CONTAINS

      subroutine nlefsm( scell, nua, na, isa, xa, indxua,
     .                   maxnh, maxnd, lasto, lastkb, iphorb, 
     .                   iphKB, numd, listdptr, listd, numh, 
     .                   listhptr, listh, nspin, Dscf, Enl, 
     .                   fa, stress, H )
C *********************************************************************
C Calculates non-local (NL) pseudopotential contribution to total 
C energy, atomic forces, stress and hamiltonian matrix elements.
C Energies in Ry. Lengths in Bohr.
C Writen by J.Soler and P.Ordejon, June 1997.
C **************************** INPUT **********************************
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer isa(na)          : Species index of each atom
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer indxua(na)       : Index of equivalent atom in unit cell
C integer maxnh            : First dimension of H and listh
C integer maxnd            : Maximum number of elements of the
C                            density matrix
C integer lasto(0:na)      : Position of last orbital of each atom
C integer lastkb(0:na)     : Position of last KB projector of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom,
C                            where no=lasto(na)
C integer iphKB(nokb)      : Index of each KB projector in its atom,
C                            where nokb=lastkb(na)
C integer numd(nuo)        : Number of nonzero elements of each row of the
C                            density matrix
C integer listdptr(nuo)    : Pointer to the start of each row (-1) of the
C                            density matrix
C integer listd(maxnd)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer numh(nuo)        : Number of nonzero elements of each row of the
C                            hamiltonian matrix
C integer listhptr(nuo)    : Pointer to the start of each row (-1) of the
C                            hamiltonian matrix
C integer listh(maxnh)     : Nonzero hamiltonian-matrix element column 
C                            indexes for each matrix row
C integer nspin            : Number of spin components
C real*8  Dscf(maxnd,nspin): Density matrix
C ******************* INPUT and OUTPUT *********************************
C real*8 fa(3,na)          : NL forces (added to input fa)
C real*8 stress(3,3)       : NL stress (added to input stress)
C real*8 H(maxnh,nspin)    : NL Hamiltonian (added to input H)
C **************************** OUTPUT *********************************
C real*8 Enl               : NL energy
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs,  only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,  only : GlobalToLocalOrb
      use atmfuncs,      only : rcut, epskb
      use neighbour    , only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour    , only : mneighb
      use alloc,         only : re_alloc, de_alloc

      integer, intent(in) ::
     .   maxnh, na, maxnd, nspin, nua

      integer, intent(in)  ::
     .  indxua(na), iphKB(*), iphorb(*), isa(na),  
     .  lasto(0:na), lastkb(0:na), listd(maxnd), listh(maxnh),
     .  numd(*), numh(*), listdptr(*), listhptr(*)

      real(dp), intent(in) :: scell(3,3), Dscf(maxnd,nspin),
     .                        xa(3,na)
      real(dp), intent(inout) :: fa(3,nua), stress(3,3)
      real(dp), intent(inout) :: H(maxnh,nspin)
      real(dp), intent(out)   :: Enl

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector

      integer, save ::  maxno = 2000
  
      integer
     .  ia, ikb, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, 
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxkba

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, volume

      real(dp), dimension(:), pointer :: Di, Vi
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:,:), pointer :: grSki

      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall
C ......................

C Start time counter
      call timer( 'nlefsm', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range
      rmaxo = 0.0d0
      rmaxkb = 0.0d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once
      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory
      nullify( Di )
      call re_alloc( Di, 1, no, name='Di', routine='nlefsm' )
      nullify( Vi )
      call re_alloc( Vi, 1, no, name='Vi', routine='nlefsm' )
      nullify( listed )
      call re_alloc( listed, 1, no, name='listed', routine='nlefsm' )
      nullify( listedall )
      call re_alloc( listedall, 1, no, name='listedall',
     $               routine='nlefsm' )

      Enl = 0.0d0
      do jo = 1,no
        Di(jo) = 0.0d0         !AG: superfluous after initial re_alloc
        Vi(jo) = 0.0d0
        listed(jo) = .false.
        listedall(jo) = .false.
      enddo

C Make list of all orbitals needed for this node
      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Find maximum number of KB projectors of one atom = maxkba
      maxkba = 0
      do ka = 1,na
        nkb = lastkb(ka) - lastkb(ka-1)
        maxkba = max(maxkba,nkb)
      enddo

C Allocate local arrays that depend on saved parameters
      nullify( iano )
      call re_alloc( iano, 1, maxno, name='iano', routine='nlefsm' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, name='iono', routine='nlefsm' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &               routine='nlefsm' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxkba, 1, maxno, name='Ski', 
     &               routine='nlefsm' )
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno, name='grSki',
     &               routine='nlefsm' )

C Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

C Loop on atoms with KB projectors      
      do ka = 1,na
        kua = indxua(ka)
        ks = isa(ka)
        nkb = lastkb(ka) - lastkb(ka-1)

C Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )

C Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = iana(ina)
            is = isa(ia)
            rki = sqrt(r2ki(ina))
            do io = lasto(ia-1)+1,lasto(ia)

C Only calculate if needed locally
              if (listedall(io)) then

                ioa = iphorb(io)
          
C Find if orbital is within range
                within = .false.
                do ko = lastkb(ka-1)+1,lastkb(ka)
                  koa = iphKB(ko)
                  if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) 
     .              within = .true.
                enddo

C Find overlap between neighbour orbitals and KB projectors
                if (within) then
C Check maxno - if too small then increase array sizes
                  if (nno.eq.maxno) then
                    maxno = maxno + 100
                    call re_alloc( iano, 1, maxno, name='iano',
     $                   copy=.true., routine='nlefsm' )
                    call re_alloc( iono, 1, maxno, name='iono',
     $                   copy=.true., routine='nlefsm' )
                    call re_alloc( xno, 1, 3, 1, maxno, name='xno', 
     &                   copy=.true., routine='nlefsm' )
                    call re_alloc( Ski, 1, maxkba, 1, maxno, name='Ski', 
     &                   copy=.true., routine='nlefsm' )
                    call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno,
     $                   name='grSki', routine='nlefsm', copy=.true. )
                  endif
                  nno = nno + 1
                  iono(nno) = io
                  iano(nno) = ia
                  do ix = 1,3
                    xno(ix,nno) = xki(ix,ina)
                  enddo
                  ikb = 0
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ikb = ikb + 1
                    ioa = iphorb(io)
                    koa = iphKB(ko)
                    call matel( 'S', ks, is, koa, ioa, xki(1,ina),
     .                    Ski(ikb,nno), grSki(1,ikb,nno) )
                  enddo

                endif

              endif

            enddo

          enddo

C Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            ia = iano(ino)

            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio.gt.0) then
C Valid orbital

              if (ia .le. nua) then

C Scatter density matrix row of orbital io
                do j = 1,numd(iio)
                  ind = listdptr(iio)+j
                  jo = listd(ind)
                  do ispin = 1,nspin
                    Di(jo) = Di(jo) + Dscf(ind,ispin)
                  enddo
                enddo
          
C Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  jo = listh(listhptr(iio)+j)
                  listed(jo) = .true.
                enddo

C Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  if (listed(jo)) then

C Loop on KB projectors
                    ikb = 0
                    do ko = lastkb(ka-1)+1,lastkb(ka)
                      ikb = ikb + 1
                      koa = iphKB(ko)
                      epsk = epskb(ks,koa)
                      Sik = Ski(ikb,ino)
                      Sjk = Ski(ikb,jno)
                      Vi(jo) = Vi(jo) + epsk * Sik * Sjk
                      Cijk = Di(jo) * epsk
                      Enl = Enl + Cijk * Sik * Sjk
                      do ix = 1,3
                        fik = 2.d0 * Cijk * Sjk * grSki(ix,ikb,ino)
                        fa(ix,ia)  = fa(ix,ia)  - fik
                        fa(ix,kua) = fa(ix,kua) + fik
                        do jx = 1,3
                          stress(jx,ix) = stress(jx,ix) +
     .                                xno(jx,ino) * fik / volume
                        enddo
                      enddo
                    enddo

                  endif
                enddo

C Pick up contributions to H and restore Di and Vi
                do j = 1,numh(iio)
                  ind = listhptr(iio)+j
                  jo = listh(ind)
                  do ispin = 1,nspin
                    H(ind,ispin) = H(ind,ispin) + Vi(jo)
                  enddo
                  Vi(jo) = 0.0d0
                  listed(jo) = .false.
                enddo
                do j = 1,numd(iio)
                  jo = listd(listdptr(iio)+j)
                  Di(jo) = 0.0d0
                enddo

              endif

            endif

          enddo

      enddo

C Deallocate local memory

      call de_alloc( grSki, name='grSki' )
      call de_alloc( Ski, name='Ski' )
      call de_alloc( xno, name='xno' )
      call de_alloc( iono, name='iono' )
      call de_alloc( iano, name='iano' )

      call de_alloc( listedall, name='listedall' )
      call de_alloc( listed, name='listed' )
      call de_alloc( Vi, name='Vi' )
      call de_alloc( Di, name='Di' )

      call timer( 'nlefsm', 2 )

      end subroutine nlefsm

      end module m_nlefsm
