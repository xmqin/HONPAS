! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module atomlist

      use precision
      use alloc
      use parallel, only: IOnode
      use  atmfuncs, only: nofis, nkbfis, izofis, massfis,
     $                     rcut, atmpopfio, zvalfis, floating
      use atm_types, only: species
      use siesta_geom, only: na_u, na_s, xa, isa, xa_last
      implicit none

      private
      public :: initatomlists, superc
      public :: superx          ! for backwards compatibility
      public :: reset_atomlists

!
!     Instead of "generic" na, no, and nokb, we use:
!
! For "supercell" (intended for k-point calcs)
      integer, save, public          :: no_s    ! Number of orbitals
      integer, save, public          :: nokb_s  ! Number of KB projs

! Same for "unit", or "real" cell:
      integer, save, public          :: no_u, nokb_u
      integer, save, public          :: no_l=1  !      Local to node

! Here 'na' is a generic number. It could be na_u or na_s, depending
! on whether we need a supercell or not.

C character*2 elem(na)      : Element name of each atom.
C integer lasto(0:na)       : Position of last orbital of each atom
C integer lastkb(0:na)      : Position of last KB proj. of each atom
C integer iza(na)           : Atomic number of each atom
C real*8 amass(na)          : Atomic mass of each atom
C real*8 qa(na)             : Neutral atom charge of each atom

      character(len=2), pointer, save, public :: elem(:) => null()
! elem will contain element names, so is 2 chars in length

      integer, pointer, save, public  :: iza(:) => null()
      integer, pointer, save, public  :: lasto(:) => null()
      integer, pointer, save, public  :: lastkb(:) => null()
      real(dp), pointer, save, public :: amass(:) => null()
      real(dp), pointer, save, public :: qa(:) => null()
      integer, pointer, save, public  :: indxua(:) => null()

      ! This array depends on the actual geometry, so
      ! it does not properly belong here. It is used only
      ! for communication between nlefsm and hsparse.
      ! For safety, it is initialized in each invokation
      ! of hsparse.
      logical, pointer, save, public  :: in_kb_orb_u_range(:) => null()

!     Index of equivalent atom in "u" cell
      real(dp), save, public          :: rmaxv    ! Max cutoff for Vna
      real(dp), save, public          :: rmaxo    ! Max cuoff for at. orb.
      real(dp), save, public          :: rmaxkb   ! Max cuoff for KB proj.
      real(dp), save, public          :: rmaxdftu ! Max cuoff for DFTU proj.

      real(dp), save, public          :: qtot ! Total number of elect.
      real(dp), save, public          :: qtots(2) ! Total number of electrons per spin

      real(dp), save, public          :: zvaltot
                                         ! Total number of pseudoprotons
                                         ! (excluding those of ghost atoms)


      integer, pointer, save, public  :: iaorb(:) => null()
                                ! Atomic index of each orbital
      integer, pointer, save, public  :: iphorb(:)  => null()
                         ! Orbital index of each  orbital in its atom
      real(dp), pointer, save, public :: Datm(:) => null()
                         !  Neutral atom charge 
                         !  of each orbital
      real(dp), pointer, save, public :: rco(:) => null()
                         ! Cutoff radius of each orbital

      integer, pointer, save, public :: indxuo(:) => null()
                   !        Index of equivalent orbital in "u" cell

      integer, pointer, save, public :: iakb(:) => null()
!         Atomic index of each KB projector
      integer, pointer, save, public :: iphKB(:) => null()
!         Index of each KB projector in its atom (negative)
      real(dp), pointer, save, public :: rckb(:) => null()
!         Cutoff radius of each KB projector
!

      CONTAINS

!=======================================================
      subroutine initatomlists()

      use atm_types, only: species_info
      use radial, only: rad_func
      use dftu_specs, only: switch_dftu
      
C Routine to initialize the atomic lists.
C
      integer  ia, io, is, nkba, noa, nol, nokbl, ioa, ikb
      type(species_info), pointer :: spp
      type(rad_func), pointer :: pp

      call re_alloc( indxua, 1, na_u, 'indxua', 'atomlist' )
      call re_alloc( lastkb, 0, na_u, 'lastkb', 'atomlist' )
      call re_alloc( lasto, 0, na_u, 'lasto', 'atomlist' )
      call re_alloc( qa, 1, na_u, 'qa', 'atomlist' )
      call re_alloc(xa_last,1,3,1,na_u,'xa_last','atomlist')
      call re_alloc( amass, 1, na_u, 'amass', 'atomlist' )
      call re_alloc( in_kb_orb_u_range, 1, na_u, 'in_kb', 'atomlist' )
!
!     Find number of orbitals and KB projectors in cell
!
      no_u = 0
      nokb_u = 0
      do ia = 1,na_u
        is = isa(ia)
        noa  = species(is)%norbs
        nkba = species(is)%nprojs
        no_u = no_u + noa
        nokb_u = nokb_u + nkba
      enddo

      na_s = na_u
      no_s = no_u
      nokb_s = nokb_u

      call re_alloc( iaorb, 1, no_u, 'iaorb', 'atomlist' )
      call re_alloc( indxuo, 1, no_u, 'indxuo', 'atomlist' )
      call re_alloc( iphorb, 1, no_u, 'iphorb', 'atomlist' )
      call re_alloc( Datm, 1, no_u, 'Datm', 'atomlist' )
      call re_alloc( rco, 1, no_u, 'rco', 'atomlist' )
      call re_alloc( iaKB, 1, nokb_u, 'iaKB', 'atomlist' )
      call re_alloc( iphKB, 1, nokb_u, 'iphKB', 'atomlist' )
      call re_alloc( rckb, 1, nokb_u, 'rckb', 'atomlist' )

c Initialize atomic lists
      nol = 0
      nokbl = 0
      qtot = 0._dp
      rmaxv  = 0._dp
      rmaxo  = 0._dp
      rmaxkb = 0._dp
      rmaxdftu = 0._dp
      lasto(0) = 0
      lastkb(0) = 0
      zvaltot = 0.0_dp
      do ia = 1,na_u
        is = isa(ia)
        if (.not. floating(is)) then
           zvaltot = zvaltot + zvalfis(is)
        endif
        noa  = nofis(is)
        nkba = nkbfis(is)
        lasto(ia)  = lasto(ia-1)  + noa
        lastkb(ia) = lastkb(ia-1) + nkba
        rmaxv = max( rmaxv, rcut(is,0) )
        iza(ia) = izofis(is)
        amass(ia) = massfis(is)
        qa(ia) = 0.0_dp
        do io = 1,noa
          nol = nol + 1
          rmaxo = max( rmaxo, rcut(is,io) )
          iaorb(nol) = ia
          iphorb(nol) = io
          Datm(nol) = atmpopfio(is,io)
          qa(ia) = qa(ia) + Datm(nol)
          qtot = qtot + Datm(nol)
        enddo
        do io = 1,nkba
          nokbl = nokbl + 1
          rmaxkb = max( rmaxkb, rcut(is,-io) ) 
          iaKB(nokbl) = ia
          iphKB(nokbl) = -io
        enddo
        if( switch_dftu ) then
           spp => species(is)
           do io = 1, spp%n_pjdftunl
              pp => spp%pjdftu(io)
              rmaxdftu = max( rmaxdftu, pp%cutoff )
           enddo
        endif
      enddo

! Find rco and rckb .............................

      do ia = 1,na_u
        is = isa(ia)
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rco(io) = rcut(is,ioa)
        enddo
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rckb(ikb) = rcut(is,ioa)
        enddo
      enddo
      
      if (IOnode)
     $   write(6,'(/a,3(1x,i5))')
     $   'initatomlists: Number of atoms, orbitals, and projectors: ',
     $     na_u, no_u, nokb_u

      end subroutine initatomlists

      subroutine reset_atomlists()

      use alloc, only: de_alloc

      call de_alloc( indxua, 'indxua', 'atomlist' )
      call de_alloc( isa, 'isa', 'atomlist')
      call de_alloc( iza, 'iza', 'atomlist')
      call de_alloc( lastkb, 'lastkb', 'atomlist' )
      call de_alloc( lasto, 'lasto', 'atomlist' )
      call de_alloc( qa, 'qa', 'atomlist' )
      call de_alloc(xa_last,'xa_last','atomlist')
      call de_alloc( amass, 'amass', 'atomlist' )
      call de_alloc( in_kb_orb_u_range, 'in_kb', 'atomlist' )

      call de_alloc( iaorb, 'iaorb', 'superc' )
      call de_alloc( indxuo, 'indxuo', 'superc' )
      call de_alloc( iphorb, 'iphorb', 'superc' )
      call de_alloc( Datm, 'Datm', 'superc' )
      call de_alloc( rco, 'rco', 'superc' )
      call de_alloc( iaKB, 'iaKB', 'superc' )
      call de_alloc( iphKB, 'iphKB', 'superc' )
      call de_alloc( rckb, 'rckb', 'superc' )

      end subroutine reset_atomlists

      subroutine superc( ucell, scell, nsc)

C Finds the supercell required to avoid multiple image overlaps,
C and expands arrays from unit cell to supercell
C Written by J.M.Soler. August 1998.
! Rewritten Alberto Garcia, May 2000.

      implicit none
      integer, intent(in)  :: nsc(3)      ! Diagonal elements of mscell
      real(dp), intent(in) :: ucell(3,3)  ! Unit cell vectors
      real(dp), intent(out) :: scell(3,3) ! Supercell vectors

C Internal variables
      integer           ia, io, iua, iuo, ja, ncells,
     $                  na, no, nokb

!
!      Find number of cells, atoms and orbitals in supercell
      ncells = nsc(1) * nsc(2) * nsc(3)
      na    = na_u   * ncells
      no    = no_u   * ncells
      nokb  = nokb_u * ncells

!
!     Reallocate arrays if needed
!
      if (na.gt.na_s) then
        call re_alloc( indxua, 1, na, 'indxua', 'atomlist', .true. )
        call re_alloc( isa, 1, na, 'isa', 'atomlist', .true. )
        call re_alloc( iza, 1, na, 'iza', 'atomlist', .true. )
        call re_alloc( lastkb, 0, na, 'lastkb', 'atomlist', .true. )
        call re_alloc( lasto, 0, na, 'lasto', 'atomlist', .true. )
        call re_alloc( in_kb_orb_u_range, 1, na,
     $                 'in_kb', 'atomlist', .true. )
        call re_alloc( qa, 1, na, 'qa', 'atomlist', .true. )
        call re_alloc( xa, 1, 3, 1, na, 'xa', 'atomlist', .true. )
        call re_alloc(xa_last, 1,3, 1,na, 'xa_last', 'superc',
     &                copy=.true. )
      endif

      na_s  = na

C Find supercell vectors and atomic coordinates in supercell 
      call superx( ucell, nsc, na_u, na_s, xa, scell )

C Find indxua and expand isa, iza, lasto and lastkb to supercell 
      do ia = 1,na_s
        ja = mod(ia-1,na_u) + 1
        indxua(ia) = ja
        isa(ia)    = isa(ja)
        iza(ia)    = iza(ja)
        lasto(ia)  = lasto(ia-1)  + lasto(ja)  - lasto(ja-1)
        lastkb(ia) = lastkb(ia-1) + lastkb(ja) - lastkb(ja-1)
      enddo

! Reallocate orbital arrays

      if (no.gt.no_s) then
        call re_alloc(iaorb, 1,no, routine='superc',copy=.true.)
        call re_alloc(indxuo, 1,no, routine='superc',copy=.true.)
        call re_alloc(iphorb, 1,no, routine='superc',copy=.true.)
        call re_alloc(Datm, 1,no, routine='superc',copy=.true.)
        call re_alloc(rco, 1,no, routine='superc',copy=.true.)
      endif

      no_s = no

C Find indxuo and expand iaorb, iphorb, and rco 
      do io = 1,no_s
        indxuo(io) = mod(io-1,no_u) + 1
      enddo
      do ia = 1,na_s
        do io = lasto(ia-1)+1,lasto(ia)
          iuo = indxuo(io)
          iaorb(io)  = ia
          iphorb(io) = iphorb(iuo)
          rco(io)    = rco(iuo)
        enddo
      enddo

! Reallocate projector arrays

      if (nokb .gt. nokb_s) then
        call re_alloc(iaKB, 1,nokb, routine='superc',copy=.true.)
        call re_alloc(iphKB, 1,nokb, routine='superc',copy=.true.)
        call re_alloc(rckb, 1,nokb, routine='superc',copy=.true.)
      endif

      nokb_s = nokb

C Expand iakb and iphKB and rckb

      do ia = 1,na_s
        iua = indxua(ia)
        iuo = lastkb(iua-1)
        do io = lastkb(ia-1)+1,lastkb(ia)
          iuo = iuo + 1
          iakb(io)  = ia
          iphKB(io) = iphKB(iuo)
          rckb(io)  = rckb(iuo)
        enddo
      enddo

      if (IOnode .and. ncells.gt.1) then
         write(6,'(/,a,i6,a,i6,a,i6,a,i8)')
     .    'superc: Internal auxiliary supercell:',
     .     nsc(1), ' x', nsc(2), ' x', nsc(3), '  =', ncells

         write(6,'(a,1x,i5,2(1x,i6))')
     $     'superc: Number of atoms, orbitals, and projectors: ',
     $     na_s, no_s, nokb_s
           
      endif

      end subroutine superc

      SUBROUTINE SUPERX( UCELL, NSC, NA, MAXA, XA, SCELL )

C **********************************************************************
C Generates supercell vectors and atomic positions.
C Written by J.M.Soler, August 1998
C *************** Input ************************************************
C Real*8  UCELL(3,3)  : Unit cell vectors UCELL(Ixyz,Ivector)
C Integer NSC(3)      : Number of cells in each supercell direction:
C                         SCELL(ix,i) = UCELL(ix,i) * NSC(i)
C Integer NA          : Number of atoms in unit cell
C Integer MAXA        : Second dimension of XA
C *************** Input and output *************************************
C Real*8  XA(3,MAXA)  : Atomic positions in unit cell (input) and
C                       in supercell (output), in cartesian coord.
C Real*8  SCELL(3,3)  : Supercell vectors
C *********** Units ****************************************************
C Units of CELL and XA are arbitrary but must be the same
C *********** Behavior *************************************************
C - If NA*NCELLS > MAXA (where NCELLS is the total number of cells),
C   the supercell atomic coordinates are not generated.
C - The first supercell atoms are those of the initial unit cell, i.e.
C   the positions XA(i,ia) for (ia.le.NA) are not modified.
C - The remaining atoms are ordered by unit cells, i.e. the atom ia
C   is equivalent to the unit-cell atom ja=MOD(ia-1,NA)+1
C **********************************************************************

      IMPLICIT          NONE
      INTEGER           MAXA, NA, NSC(3)
      DOUBLE PRECISION  SCELL(3,3), UCELL(3,3), XA(3,MAXA)

C Internal variables
      INTEGER           I, I1, I2, I3, IA, IX, JA, NCELLS
      DOUBLE PRECISION  XC(3)

C Find supercell vectors
      DO 10 I = 1,3
        DO 5 IX = 1,3
          SCELL(IX,I) = UCELL(IX,I) * NSC(I)
    5   CONTINUE
   10 CONTINUE

C Expand atomic positions to supercell
      NCELLS = NSC(1) * NSC(2) * NSC(3)
      IF (NA*NCELLS .LE. MAXA) THEN
        IA = 0
        DO 60 I3 = 0,NSC(3)-1
        DO 50 I2 = 0,NSC(2)-1
        DO 40 I1 = 0,NSC(1)-1
          DO 15 IX = 1,3
            XC(IX) = UCELL(IX,1)*I1 + UCELL(IX,2)*I2 + UCELL(IX,3)*I3
   15     CONTINUE
          DO 30 JA = 1,NA
            IA = IA + 1
            DO 20 IX = 1,3
              XA(IX,IA) = XA(IX,JA) + XC(IX)
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
   50   CONTINUE
   60   CONTINUE
      ENDIF

      END subroutine superx

      end module atomlist



