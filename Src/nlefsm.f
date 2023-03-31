! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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
     .                   fa, stress, H , matrix_elements_only)
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
C integer nspin            : Number of spin components of Dscf and H
C                            If computing only matrix elements, it
C                            can be set to 1.
C logical matrix_elements_only:
C integer Dscf(maxnd,nspin): Density matrix. Not touched if computing
C                            only matrix elements.
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
      use atm_types,     only : nspecies
      use atomlist,      only : in_kb_orb_u_range
      use atmfuncs,      only : rcut, epskb, orb_gindex, kbproj_gindex
      use atmfuncs,      only : nofis, nkbfis
      use chemical,      only : is_floating
      use neighbour,     only : iana=>jan, r2ki=>r2ij, xki=>xij
      use neighbour,     only : mneighb, reset_neighbour_arrays
      use alloc,         only : re_alloc, de_alloc
      use m_new_matel,   only : new_matel

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
      logical, intent(in)     :: matrix_elements_only

      real(dp) ::   volcel
      external ::   timer, volcel

C Internal variables ................................................
C maxno  = maximum number of basis orbitals overlapping a KB projector
! This should be estimated beforehand, to avoid checks,
! or a "guard" region implemented for a single check at the end
      
      integer, save ::  maxno = 500
  
      integer
     .  ia, ikb, ina, ind, ino,
     .  io, iio, ioa, is, ispin, ix, ig, kg,
     .  j, jno, jo, jx, ka, ko, koa, ks, kua,
     .  nkb, nna, nno, no, nuo, nuotot, maxkba
      integer :: natoms_k_over, max_nno_used

      integer, dimension(:), pointer :: iano, iono

      real(dp)
     .  Cijk, epsk, fik, rki, rmax, rmaxkb, rmaxo, 
     .  Sik, Sjk, volume

      real(dp), dimension(:), pointer :: Di, Vi
      real(dp), dimension(:,:), pointer :: Ski, xno
      real(dp), dimension(:,:,:), pointer :: grSki

      logical ::   within
      logical, dimension(:), pointer ::  listed, listedall
      !
      real(dp), allocatable :: rorbmax(:), rkbmax(:)
C ......................

C Start time counter
      call timer( 'nlefsm', 1 )

C Find unit cell volume
      volume = volcel( scell ) * nua / na

C Find maximum range and maximum number of KB projectors
      maxkba = 0
      
      allocate(rorbmax(nspecies),rkbmax(nspecies))
      do is = 1, nspecies

         ! Species orbital range
         rorbmax(is) = 0.0_dp
         do io = 1, nofis(is)
            rorbmax(is) = max(rorbmax(is), rcut(is,io))
         enddo
         
         ! Species KB range
         io = nkbfis(is)
         rkbmax(is) = 0.0_dp
         do ikb = 1, io
            rkbmax(is) = max(rkbmax(is), rcut(is,-ikb))
         enddo
         maxkba = max(maxkba,io)

      enddo
      rmaxo = maxval(rorbmax(1:nspecies))
      rmaxkb = maxval(rkbmax(1:nspecies))
      ! Calculate max extend
      rmax = rmaxo + rmaxkb

C Initialize arrays Di and Vi only once

      no = lasto(na)
      nuotot = lasto(nua)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)

C Allocate local memory

      nullify( Vi )
      call re_alloc( Vi, 1, no, 'Vi', 'nlefsm' )
      Vi(1:no) = 0.0_dp         ! OK. Later on, any non-zero elements
                                ! will be zero-ed out explicitly
      nullify( listed )
      call re_alloc( listed, 1, no, 'listed', 'nlefsm' )
      listed(1:no) = .false.
      nullify( listedall )
      call re_alloc( listedall, 1, no, 'listedall', 'nlefsm' )
      listedall(1:no) = .false.

      if (.not. matrix_elements_only) then
         nullify( Di )
         call re_alloc( Di, 1, no, 'Di', 'nlefsm' )
         Di(1:no) = 0.0_dp
      endif

      Enl = 0.0d0

!     Make list of all orbitals needed for this node
      
      do io = 1,nuo
        ! we need this process's orbitals...
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        listedall(iio) = .true.
        ! ... and those with which they interact
        do j = 1,numh(io)
          jo = listh(listhptr(io)+j)
          listedall(jo) = .true.
        enddo
      enddo

C Allocate local arrays that depend on saved parameters
      nullify( iano )
      call re_alloc( iano, 1, maxno, 'iano', 'nlefsm' )
      nullify( iono )
      call re_alloc( iono, 1, maxno, 'iono', 'nlefsm' )
      nullify( xno )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'nlefsm' )
      nullify( Ski )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski', 'nlefsm' )
      nullify( grSki )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno, 'grSki',
     &               'nlefsm' )

C     Initialize neighb subroutine
      call mneighb( scell, rmax, na, xa, 0, 0, nna )

! Loop on atoms with KB projectors
! All processes will be doing this loop over atoms.
! This is one reason for non-scalability
!
!        And what happens if there is no supercell?
!        How do we count out-of-unit-cell interactions?
!        ... they are automatically accounted for, in
!        the same way as the Hmu_nu terms themselves.
!
      natoms_k_over = 0
      max_nno_used = 0
      do ka = 1,na
!        Only the atoms within the proper
!        distance of a unit cell orbital (in our process) should
!        be considered, not the whole supercell.
!        This array was initialized in hsparse
         
        if (.not. in_kb_orb_u_range(ka)) CYCLE
        
        ks = isa(ka)
        ! Cycle also if ghost-orbital species...
        if (is_floating(ks)) CYCLE
        
        kua = indxua(ka)  ! Used only if forces and energies are comp.

C       Find neighbour atoms
        call mneighb( scell, rmax, na, xa, ka, 0, nna )

        nno = 0
        do ina = 1,nna
          rki = sqrt(r2ki(ina))
          ia = iana(ina)
          is = isa(ia)
          !     Early exit if too far
          !     This duplicates the test in hsparse...
          if (rki - rkbmax(ks) - rorbmax(is) > 0.d0) CYCLE

          ! Loop over orbitals close enough to overlap
          do io = lasto(ia-1)+1,lasto(ia)

C           Only calculate if needed locally in our MPI process
             if (.not. listedall(io)) CYCLE
             
             ioa = iphorb(io)
             ! rki_minus_rc_orb= rki - rcut(is,ioa)

             ! Find if orbital is within range
             ! This can be done with rkbmax(ks):
             within = (rki-rkbmax(ks)) < rcut(is,ioa)
             if (.not. within) CYCLE
              
!             Find overlap between neighbour orbitals and KB projectors

             if (nno.eq.maxno) call increase_maxno()
              
              nno = nno + 1  ! Update number of overlaps to keep
              iono(nno) = io 
              iano(nno) = ia
              do ix = 1,3
                 xno(ix,nno) = xki(ix,ina)
              enddo

! For each overlap family we keep the individual
! KB-orb matrix elements
! This will store some zeros sometimes, as some
! of the KBs might not actually overlap
! We could re-check the distances...
! Not worth it, as then we would have different
! numbers of matrix elements for different orbitals, and
! the bookeeping would get messy
              
              
              ikb = 0
              ig = orb_gindex(is,ioa)
              do ko = lastkb(ka-1)+1,lastkb(ka)
                 koa = iphKB(ko)
                 ! if ( rki_minus_rc_orb > rcut(ks,koa) CYCLE
                 ikb = ikb + 1
                 ! epsk_sqrt = sqrt(epskb(ks,koa))
                 kg = kbproj_gindex(ks,koa)
                 call new_MATEL( 'S', kg, ig, xki(1:3,ina),
     &                Ski(ikb,nno), grSki(1:3,ikb,nno) )
                 !  Maybe: Ski = epskb_sqrt * Ski
                 !         grSki = epskb_sqrt * grSki
              enddo

           enddo ! loop over orbitals

        enddo ! loop over neighbor atoms

!     Now we check which of the overlaps of our atom's KB's involve
!     two orbitals: one in the unit cell, and handled by our process,
!     and the other unrestricted
        
        max_nno_used = max(max_nno_used, nno)
        do ino = 1,nno    ! loop over overlaps
          ia = iano(ino)
          if (ia > nua) CYCLE  ! We want the 1st orb to be in the unit cell

          io = iono(ino)
          ! Note that if ia is in the unit cell, io is <= nuo,
          ! so that this call makes sense
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio == 0) CYCLE

          !  Scatter filter of desired matrix elements
          do j = 1,numh(iio)
             ind = listhptr(iio)+j
             jo = listh(ind)
             listed(jo) = .true.
             if (.not. matrix_elements_only) then
                do ispin = 1,nspin ! Both spins add up...
                   Di(jo) = Di(jo) + Dscf(ind,ispin)
                enddo
             endif
          enddo

! Find matrix elements with other neighbour orbitals
! Note that several overlaps might contribute to the
! same matrix element, hence the additions above (Dscf) and below (H)
          
          do jno = 1,nno
             jo = iono(jno)
             ! Check whether there is H_io_jo...
             if (.not. listed(jo)) CYCLE 

! Loop on KB projectors again. Note that ikb and ko run
! in step. ko is only needed for the Epskb factor.
! maybe we can store it with the value of the projector.

             ikb = 0
             do ko = lastkb(ka-1)+1,lastkb(ka)
                ikb = ikb + 1
                koa = iphKB(ko)
                epsk = epskb(ks,koa)
                Sik = Ski(ikb,ino)
                Sjk = Ski(ikb,jno)
                Vi(jo) = Vi(jo) + epsk * Sik * Sjk
                ! We should distinguish "energy-only" and
                ! "forces-and-stress"
                if (.not. matrix_elements_only) then
                   Cijk = Di(jo) * epsk
                   Enl = Enl + Cijk * Sik * Sjk
                   do ix = 1,3
                      fik = 2.d0 * Cijk * Sjk * grSki(ix,ikb,ino)
                      fa(ix,ia)  = fa(ix,ia)  - fik
                      fa(ix,kua) = fa(ix,kua) + fik
                      do jx = 1,3
                         stress(jx,ix) = stress(jx,ix) +
     &                        xno(jx,ino) * fik / volume
                      enddo
                   enddo
                endif

             enddo
          
          enddo ! loop over second orbitals

C         Pick up contributions to H and restore Di and Vi
          do j = 1,numh(iio)
             ind = listhptr(iio)+j
             jo = listh(ind)
             do ispin = 1,nspin
                H(ind,ispin) = H(ind,ispin) + Vi(jo)
             enddo
             Vi(jo) = 0.0d0     ! See initial zero-out at top
             listed(jo) = .false.
             if (.not. matrix_elements_only) Di(jo) = 0.0d0
          enddo

       enddo  ! loop over 1st orbitals
       natoms_k_over = natoms_k_over + 1
      enddo   ! loop over atoms holding KB projectors

      if (Node == 0) then
         ! For future diagnostics
         ! Currently only the root process outputs info
         write(6,"(a,2i8)")
     $     "No. of atoms with KB's overlaping orbs in proc 0." //
     $     " Max # of overlaps:", natoms_k_over, max_nno_used
      endif
      
C     Deallocate local memory
!      call new_MATEL( 'S', 0, 0, 0, 0, xki, Ski, grSki )
      call reset_neighbour_arrays( )
      call de_alloc( grSki, 'grSki', 'nlefsm' )
      call de_alloc( Ski, 'Ski', 'nlefsm' )
      call de_alloc( xno, 'xno', 'nlefsm' )
      call de_alloc( iono, 'iono', 'nlefsm' )
      call de_alloc( iano, 'iano', 'nlefsm' )
      call de_alloc( listedall, 'listedall', 'nlefsm' )
      call de_alloc( listed, 'listed', 'nlefsm' )
      call de_alloc( Vi, 'Vi', 'nlefsm' )
      if (.not. matrix_elements_only) then
         call de_alloc( Di, 'Di', 'nlefsm' )
      endif
      
      deallocate(rkbmax,rorbmax)

      call timer( 'nlefsm', 2 )
      
      CONTAINS
      
      subroutine increase_maxno()
      
! if too small then increase array sizes
      maxno = maxno + 10
      call re_alloc( iano, 1, maxno, 'iano', 'nlefsm',
     &        .true. )
      call re_alloc( iono, 1, maxno, 'iono', 'nlefsm',
     &        .true. )
      call re_alloc( xno, 1, 3, 1, maxno, 'xno',  'nlefsm',
     &        .true. )
      call re_alloc( Ski, 1, maxkba, 1, maxno, 'Ski',
     &        'nlefsm', .true. )
      call re_alloc( grSki, 1, 3, 1, maxkba, 1, maxno,
     &        'grSki', 'nlefsm', .true. )

      end subroutine increase_maxno
      
      end subroutine nlefsm

      end module m_nlefsm
