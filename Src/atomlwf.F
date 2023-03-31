! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine cspa(ioptlwf,iopt,natoms,no_u,no_l,lasto,isa,
     .                qa,rcoor,rh,cell,xa,nhmax,numh,listh,listhptr,
     .                maxnc,ncmax,nctmax,nfmax,nftmax,nhijmax,nbands,
     .                no_cl,nspin,Node)
C ******************************************************************************
C This subroutine builds the Localized Wave Functions, centered
C on ATOMS (searching for atoms within a cutoff radius rcoor), and 
C assigns a RANDOM initial guess to start the CG minimization.
C
C Criterion to build LWF's: 
C 1) Method of Kim et al: use more localized orbitals than 
C    occupied orbitals.
C    We assign the minimum number of orbitals so that there
C    is place for more electrons than those in the system;
C    for instance:
C      H:        1 LWF
C      C,Si:     3 LWF's
C      N:        3 LWF's
C      O:        4 LWF's
C      ...
C 2) Method of Ordejon et al: number of localized orbitals 
C    equal to number of occupied orbitals. 
C    For the initial assignment of LWF centers to atoms, atoms
C    with even number of electrons, n, get n/2 LWFs. Odd atoms
C    get (n+1)/2 and (n-1)/2 in an alternating sequence, ir order
C    of appearance (controlled by the input in the atomic coor block). 
C
C Written by P.Ordejon, 1993. 
C Re-written by P.Ordejon, November'96.
C Corrected by P.Ordejon, April'97,  May'97  
C lmax, lmaxs and nzls erased from the input by DSP, Aug 1998.
C Alternating sequence for odd species in Ordejon, by E.Artacho, Aug 2008.
C ******************************* INPUT ***************************************
C integer ioptlwf           : Build LWF's according to:
C                             0 = Read blindly from disk
C                             1 = Functional of Kim et al.
C                             2 = Functional of Ordejon-Mauri
C integer iopt              : 0 = Find structure of sparse C matrix and
C                                  build initial guess
C                             1 = Just find structure of C matrix
C integer natoms            : Number of atoms
C integer no_u            : Number of basis orbitals
C integer lasto(0:natoms)   : Index of last orbital of each atom
C integer isa(natoms)       : Species index of each atom
C real*8 qa(natoms)         : Neutral atom charge
C real*8 rcoor              : Cutoff radius of Localized Wave Functions
C real*8 rh                 : Maximum cutoff radius of Hamiltonian matrix
C real*8 cell(3,3)          : Supercell vectors
C real*8 xa(3,natoms)       : Atomic coordinates
C integer maxnc             : First dimension of C matrix, and maximum
C                             number of nonzero elements of each row of C
C integer nspin             : Number of spins
C ****************************** OUTPUT **************************************
C real*8 c(ncmax,no_u)    : Initial guess for sparse coefficients 
C                             of LWF's  (only if iopt = 0)
C integer numc(no_u)      : Control vector of C matrix
C integer listc(ncmax,no_u): Control vector of C matrix
C integer ncmax             : True value for ncmax, 
C                             If ncmax it is too small, then
C                             c, numc and listc are NOT initialized!!
C integer nctmax            : Maximum number of nonzero elements
C                             of eaxh column of C
C integer nfmax             : Maximum number of nonzero elements 
C                             of each row of F
C integer nftmax            : Maximum number of nonzero elements 
C                             of eaxh column of F
C integer nhijmax           : Maximum number of nonzero elements 
C                             of each row of Hij
C integer nbands            : Number of LWF's
C ****************************************************************************

      use precision, only: dp
      use alloc,     only: re_alloc, de_alloc
      use on_main,   only: numc, listc, c, cold, listcold
      use on_main,   only: ncg2l, ncl2g, nct2p, ncp2t
      use neighbour, only: jan, r2ij, xij, mneighb, maxnna,
     &                     reset_neighbour_arrays
      use sys,      only : die
      use spatial,  only : nL2G

      implicit none

      integer, intent(in) :: 
     .  iopt,ioptlwf,natoms,no_u,nhmax,Node,no_l,nspin

      integer, intent(inout) :: maxnc

      integer, intent(out) ::
     $     nbands,ncmax,nctmax, nfmax,nftmax,nhijmax,
     $     no_cl

      integer, intent(in) ::
     .  isa(natoms),lasto(0:natoms),
     .  numh(no_l),listh(nhmax),listhptr(no_l)

      real(dp), intent(in) ::
     .  cell(3,3),qa(natoms),rcoor,rh,xa(3,natoms)

C  Internal variables .......................................................

      integer, dimension(:), pointer ::  alist, ilist, numft
      integer, pointer               ::  indexloc(:) => null()

      integer
     . i,ia,imu,in,index,indexa,indexb,indexi,indexj,iorb,is,j,ja,
     . jj,k,mu,nct,nelectr,nf,nm,norb,nqtot,nu,numloc,nhij,indmu,
     . mul, mull, ind, nelectr1, nelectr2

      integer,                            save :: iseed = 17
      integer,                            save :: nna = 200

      logical, dimension(:), pointer ::  lneeded

      real(dp) ::
     .  cg,fact(2),qtot,r,rmax,rr(3),rrmod,
     .  snor,tiny,randomg,cgval

      external :: randomg

      logical, save :: firstcall = .true., secondodd = .false.

      if (firstcall) then

C Initialise random number generator
        cgval = randomg(-iseed)

        firstcall = .false.

      endif

      tiny = 1.d-10

C Check that iopt is correct ................................................
      if (iopt .ne. 0 .and. iopt .ne. 1) then
        call die('cspa: Wrong iopt in cspa')
      endif


C Allocate local arrays
      nullify( alist, ilist, numft, lneeded )
      call re_alloc( alist,   1, natoms, 'alist',   'cspa' )
      call re_alloc( ilist, 1, no_u, name='ilist', routine='cspa' )
      call re_alloc( numft, 1, no_u, name='numft', routine='cspa' )
      call re_alloc( lneeded, 1, no_u, name='lneeded',
     &               routine='cspa' )

C Work out which basis functions are required in local C copy
!     These will be those which interact with the local orbitals

      lneeded(1:no_u) = .false.
      do mul = 1,no_l
        ind = listhptr(mul)
        lneeded(nL2G(mul,Node+1)) = .true.
        do i = 1,numh(mul)
          mu = listh(ind+i)
          lneeded(mu) = .true.
        enddo
      enddo
!
!     Now build up the indexes ncG2L and ncL2G
!
      no_cl = 0
      ncG2L(1:no_u) = 0
      ncT2P(1:no_u) = 0
      ! First set of orbitals: those handled by 
      ! the local node
      do mu = 1,no_l
        no_cl = no_cl + 1  ! = mu
        ncL2G(no_cl) = nL2G(mu,Node+1)
        ncG2L(nL2G(mu,Node+1)) = no_cl  ! = mu
      enddo
      ! Second set of orbitals: those interacting
      ! with the orbitals handled by the local node
      ! (and not already counted)
      do mu = 1,no_u
        if (lneeded(mu).and.ncG2L(mu).eq.0) then
          no_cl = no_cl + 1
          ncL2G(no_cl) = mu
          ncG2L(mu) = no_cl
        endif
      enddo
      do mul = 1,no_l
        ! These are just identity mappings,
        ! since the first no_l entries are
        ! the same in the no_l and the no_cl sets

        ! nct2p is just a filter: if >0, its argument
        ! belongs to the (partial) set of indexes
        ! that go from 1 to no_l. It could really be
        ! a logical mask. It is always used as
        ! if (nct2p(i)>0) then...

        ! ncp2t is just the identity over 1:no_l,
        ! where it is used in the rest of the program

        ! Both arrays are dimensioned to no_u, even though
        ! only the first no_l entries are used for ncP2T.

        ncP2T(mul) = ncG2L(nL2G(mul,Node+1)) ! ncP2T(mul) = mul
        ! Note also that the range of ncG2L is 1:no_cl, 
        ncT2P(ncG2L(nL2G(mul,Node+1))) = mul ! ncT2P(mul) = mul
      enddo

C Resize arrays that depend on no_cl
      call re_alloc(listc,1,maxnc,1,no_cl,
     .              shrink=.false.,name='listc')
      call re_alloc(listcold,1,maxnc,1,no_l,
     .              shrink=.false.,name='listcold')
      call re_alloc(c,1,maxnc,1,no_cl,1,nspin,
     .              shrink=.false.,name='c')
      call re_alloc(cold,1,maxnc,1,no_l,1,nspin,
     .              shrink=.false.,name='cold')

C Initialize some stuff

      do ia = 1,natoms
        alist(ia) = 0
      enddo

      do mul = 1,no_cl
        numc(mul) = 0
      enddo

      do mu = 1,no_u
        numft(mu) = 0
      enddo

      if (iopt .eq. 0) then
        do is = 1,nspin
          do mu = 1,no_cl
            do i = 1,maxnc
              c(i,mu,is) = 0.0d0
            enddo
          enddo
        enddo
      endif

      ncmax   = 0
      nctmax  = 0
      nfmax   = 0
      nftmax  = 0
      nhijmax = 0
C ........................

C Calculate maximum length in unit cell ......................................
C determine maximum cell length
      !AG: Possible bug: cell vectors are given by the columns of cell!
      rmax = 0.0d0
      do i = -1,1
        do j = -1,1
          do k = -1,1
            !AG It should be
            ! rr(1) = i*cell(1,1) + j*cell(1,2) + k*cell(1,3)
            rr(1) = i*cell(1,1) + j*cell(2,1) + k*cell(3,1)
            rr(2) = i*cell(1,2) + j*cell(2,2) + k*cell(3,2)
            rr(3) = i*cell(1,3) + j*cell(2,3) + k*cell(3,3)
            rrmod = sqrt( rr(1)**2 + rr(2)**2 + rr(3)**2 )
            if (rrmod .gt. rmax) rmax = rrmod
          enddo
        enddo
      enddo
C ........................

C Check that there is an integer number of electrons 
      qtot = 0.0d0
      do ia = 1,natoms
        qtot = qtot + qa(ia)
      enddo
      qtot = qtot + tiny
      nqtot = nint(qtot)
      if (abs(nqtot-qtot+tiny) .gt. 1e-3) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Wrong total charge; non integer:',qtot
        endif
        call die()
      endif
C ..................

C Build control vectors of sparse LWF 
C loop over the localized wave funcions (centered at atoms)..................

C Initialize routine for neighbour search
      ! Note that this will not update nna!!
      call mneighb(cell,rcoor,natoms,xa,0,0,nna) 

C Allocate local memory
      call re_alloc(indexloc,1,max(maxnna,natoms),name='indexloc')

      index = 0   ! Counter for total number of LWF
      loop_ia: do ia = 1,natoms

        if (2.0*rcoor .lt. rmax) then
          ! Localization size smaller than max cell size
          ! Look for neighbours of atom ia
          call mneighb(cell,rcoor,natoms,xa,ia,0,nna)
          ! Reallocate in case nna has increased  (Tight form)
          call re_alloc( indexloc, 1, nna, 'indexloc', 'cspa' )
        else
          ! Localization size greater than max cell size
          ! We can treat all atoms in the cell as neighbors
          ! Make sure we have enough space
          ! (Cleaner alternative to initial allocation with  
          ! max(maxnna,natoms)
          !call sizeup_neighbour_arrays(natoms)
          !call re_alloc( indexloc, 1, natoms, name='indexloc')
          nna = natoms
          do  jj = 1,natoms
            jan(jj) = jj
          enddo
        endif

        ! All procs have information about the number of LWFs
        ! but they keep only the coefficients associated to
        ! the orbitals they manage and to those that interact with
        ! them

        ! loop over LWF's centered on atom ia
        call get_number_of_lwfs_on_atom()  ! gets indexi, nelectr
        loop_lwfs_on_ia: do indexb = 1,indexi
          index = index + 1  ! global counter for number of LWFs

c clear list of atoms considered within loc. range
          do indexa = 1,nna
            indexloc(indexa) = 0
          enddo
          numloc = 0

C initialize stuff...
          nct = 0     ! total number of coeffs for this LWF in this node
          snor = 0.0d0
          nf = 0
          do nu = 1,no_u
            ilist(nu) = 0
          enddo
 
C loop over the neighbors of ia within rcoor
          loop_neighbor_atoms: do  jj = 1,nna  
            ja = jan(jj)

            ! Check if ja has already been included in current lwf
            ! (an image atom, maybe?)
            do indexa = 1,numloc
              if (ja .eq. indexloc(indexa)) cycle loop_neighbor_atoms
            enddo
            numloc = numloc + 1
            indexloc(numloc) = ja

            !  Loop over orbitals of ja
            norb = lasto(ja) - lasto(ja-1)
            loop_orbs_on_neighbor: do iorb = 1,norb
              mu = iorb + lasto(ja-1)

              ! Get random number here for reproducibility over numbers of processors
              get_random: do
                 cgval = (randomg(iseed) - 0.5d0) * 2.0d0
                 if (abs(cgval) .ge. 1d-5) exit get_random
              enddo get_random

              ! If this orbital (global index mu) is part of our
              ! c-list, include it in the indexes

              !! alternative: 
              !! if (ncG2L(mu) == 0) cycle loop_orbs_on_neighbor

              if (ncG2L(mu).ne.0) then
                mul = ncG2L(mu)
                nm = numc(mul)
                numc(mul) = nm + 1
                if (numc(mul) .gt. maxnc) then
                  ! Reallocate all arrays that depend on maxnc
                  maxnc = maxnc + 50   ! this is not optimal
                                       ! but the arrays could be shrunk
                                       ! in routine ordern
                  call re_alloc(listc,1,maxnc,1,no_cl,name='listc')
                  call re_alloc(c,1,maxnc,1,no_cl,1,nspin,name='c')
                endif

                listc(nm+1,mul) = index  ! LWF index
                nct = nct + 1            ! number of coefficients so far

                ! Find out structure of F and Ft matrices
                ! F = S*C
                !   find orbitals (nu, global index) which interact with mu

                if (ncT2P(mul).gt.0) then
                  ! mul is one of 1:no_l
                  mull = ncT2P(mul)
                  indmu = listhptr(mull)  
                  do imu = 1,numh(mull)
                    nu = listh(indmu+imu)
                    if (ilist(nu) .eq. 0) then  ! Avoid overcounting
                      ilist(nu) = 1
                      numft(nu) = numft(nu) + 1
                      nf = nf + 1
                    endif
                  enddo
                endif
                ! At the end of this process, numft(nu) will contain 
                ! the number of locally managed orbitals (indexes 1:no_l)
                ! that interact with nu. nf will be the total number of such interactions.
                ! It is somehow the "sparsity" of the transpose of the S and H matrices
                ! But note that numft is deallocated here, and only nftmax = max(numft(:))
                ! is retained.

                ! Assign random guess for orbitals in atom ia if iopt = 0
                ! Note: only those on atom ia, to avoid duplication of work
                if (iopt .eq. 0) then
                  if (ja .eq. ia) then
                    call initguess(ia,iorb,isa(ia),nelectr,
     .                cg,cgval,Node)
                    do is = 1,nspin
                      c(nm+1,mul,is) = cg
                    enddo
                    snor = snor + cg**2
                  endif
                endif

              endif  ! if iorb is part of our c-list

            enddo  loop_orbs_on_neighbor ! iorb
          enddo     loop_neighbor_atoms ! ja

          ! Normalize LWF's if iopt = 0  .............................................
          !   (normalize to one if functions are expected to be occupied, 
          !   0.5 if half occupied 
          !   0.1 if empty)
          !
          if (iopt .eq. 0) then
            fact(1:2) = 1.0d0
            if (ioptlwf .eq. 1) then
              if (indexb.eq.indexi) then
                if (nspin.eq.1) then
                  if (2*(nelectr/2) .eq. nelectr) fact(1:2) = sqrt(0.1)
                  if (2*(nelectr/2) .ne. nelectr) fact(1:2) = sqrt(0.5)
                else
                  if (indexb.gt.nelectr1) fact(1) = sqrt(0.1)
                  if (indexb.gt.nelectr2) fact(2) = sqrt(0.1)
                endif
              endif
            endif

            do is = 1,nspin
              do mu = lasto(ia-1)+1, lasto(ia)
                mul = ncG2L(mu)
                if (mul.gt.0) then     ! mu is in our c-list
                  do in = 1,numc(mul)
                    if (listc(in,mul) .eq. index) then
                      !AG: can put spin loop here
                      c(in,mul,is) = c(in,mul,is) * fact(is)/sqrt(snor)
                    endif
                  enddo
                endif
              enddo
            enddo
          endif

          nctmax = max ( nctmax , nct )
          nfmax = max ( nfmax , nf )
          
        enddo loop_lwfs_on_ia
      enddo loop_ia

      do mul = 1,no_cl
        ncmax  = max ( ncmax  , numc(mul) )
      enddo

      do mu = 1,no_u
        nftmax = max ( nftmax , numft(mu) )
      enddo

      nbands = index

      if (index .gt. no_u) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Number of LWFs larger than  basis set size'
          write(6,*) '      Increase basis set, or use less LWFs'
        endif
        call die()
      endif

      if ((ioptlwf .eq. 2) .and. (nbands .ne. nqtot/2)) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Number of LWFs incorrectly calculated'
          write(6,*) '      Something went wrong in generating the'
          write(6,*) '      LWFs for the Ordejon-Mauri functional'
        endif
        call die()
      endif


C Find out sparse dimensions of Hij
C loop over the localized wave funcions (centered at atoms)..................

C Maximum interacion range between LWF's
      r = 2.0 * rcoor + rh

      if (2*r .ge. rmax) then
        nhijmax = nbands
      else
        call mneighb(cell,r,natoms,xa,0,0,nna)   !! ,nnmax)

        index = 0
        do ia = 1,natoms
  
C Look for neighbours of atom ia within maximum interaction range
          call mneighb(cell,r,natoms,xa,ia,0,nna) !! ,nnmax)

          nhij = 0
C Loop over the neighbors of ia within rcoor
          do  jj = 1,nna
            ja = jan(jj)
            alist(ja) = 0
          enddo
          do  jj = 1,nna
            ja = jan(jj)
            if (alist(ja) .eq. 1) goto 20
            alist(ja) = 1

C  determine how many LWF's centered in ja, depending on the atomic species 
C  (or the number of electrons)
            nelectr = qa (ja) + tiny
            if (ioptlwf .eq. 1) then
              indexj = ( ( nelectr + 2 ) / 2 )
              nelectr2 = nelectr/2
              nelectr1 = nelectr - nelectr2
            else if (ioptlwf .eq. 2) then
              !AG: Is this right?
              if ( (nelectr/2)*2 .ne. nelectr) then
                 call
     $        die("Check Ordejon-Mauri Func. indexj assignment in cspa")
              endif
              indexj = ( ( nelectr / 2 ) )
            else
              call die('cspa: Wrong functional option in cspa')
            endif

            nhij = nhij + indexj
20          continue
          enddo
          do  jj = 1,nna
            ja = jan(jj)
            alist(ja) = 0
          enddo
          nhijmax = max ( nhijmax , nhij )
        enddo  ! ia
      endif

C Deallocate local memory
      call reset_neighbour_arrays( )
      call de_alloc( lneeded, name='lneeded' )
      call de_alloc( numft, name='numft' )
      call de_alloc( ilist, name='ilist' )
      call de_alloc( alist, name='alist' )
      call de_alloc( indexloc, name='indexloc' )

      CONTAINS

      !---------------------------------------------
      subroutine get_number_of_lwfs_on_atom()

      ! Determine how many LWF's depending on the atomic species 
      ! (or the number of electrons)
        nelectr = qa(ia) + tiny 
        if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
          if (Node.eq.0) then
            write(6,*) 'cspa: Wrong atomic charge for atom ',ia
            write(6,*) '      qa = ',qa(ia),' must be an integer'
          endif
          call die()
        endif
        if (ioptlwf .eq. 1) then
          indexi = ( ( nelectr + 2 ) / 2 )
          nelectr2 = nelectr/2
          nelectr1 = nelectr - nelectr2
        else if (ioptlwf .eq. 2) then
          if ( (nelectr/2)*2 .ne. nelectr) then
c EA: Instead of dying if there is any atom of odd species, use
c an alternating scheme for assigning 1/2 more or 1/2 less, in
c strict order of appearance. The user controls where the odd LWFs
c go by defining the order of atoms in the AtomicSpeciesAndAtomicCoor... block
c OLD------ if (Node.eq.0) then
c             write(6,*) 'cspa: Wrong Order-N functional option in ',
c    .                   'cspa.'
c             write(6,*) '      You can only use the functional of'
c             write(6,*) '      Ordejon-Mauri for atoms with an even'
c             write(6,*) '      number of electrons.'
c           endif
c OLD------ call die()
c give one-extra/one-less LWF to odd species in turn with flag secondodd
            if (secondodd) then
               indexi = ( ( nelectr - 1 ) / 2 )
               secondodd = .false.
            else
               indexi = ( ( nelectr + 1 ) / 2 )
               secondodd = .true.
            endif
          else
            indexi = ( ( nelectr ) / 2 )
          endif
        else
          call die('cspa: Wrong functional option in cspa')
        endif
      end subroutine get_number_of_lwfs_on_atom

      end subroutine cspa

      subroutine initguess(ia,iorb,is,ne,cg,cgval,Node)
C *****************************************************************************
C Routine to assign an initial guess for an atomic orbital iorb in a localized
C wave function centered on atom ia. 
C Assigns a random guess if the orbital belongs to the first 'zeta' of the
C atom (lowest energy shell of its angular momentum), and if the angular
C momentum is populated in the free atom. Otherwise, sets coefficient to cero.
C
C Written by P.Ordejon, November'96 
C lmax, lmaxs and nzls erased from the input by DSP, Aug. 1998.
C ******************************* INPUT ***************************************
C integer ia                   : Atom to which orbital belongs
C integer mu                   : Orbital index within atom ia
C integer is                   : Species index of atom ia
C integer ne                   : Number of electrons of atom ia
C real*8  cgval                : Random value to set to cg if needed
C integer Node                 : Node number
C ***************************** OUTPUT ****************************************
C real*8 cg                    : Initial guess for WF coefficient
C ***************************************************************************** 
C The following functions must exist:
C
C INTEGER FUNCTION LOMAXFIS(IS)
C    Returns the maximum angular momentum of orbitals
C Input:
C     INTEGER IS : Species index
C
C INTEGER FUNCTION NZTFL(IS,L)
C    Returns the number of different basis functions with the
C    same angular momentum L.
C Input:
C     INTEGER IS : Species index
C     INTEGER L  : Angular momentum
C
C ***************************************************************************** 

      use precision

      use atmfuncs, only : lomaxfis, nztfl
      use sys,      only : die

      implicit none

      integer 
     .  ia,iorb,is,ne,Node

      real(dp) ::
     .  cg, cgval

C Internal variables .........................................................
      integer
     .  index,iz,l,lmaxp,m
      
C Initialize cg to cero
      cg = 0.0d0

C Find out angular momentum of orbital iorb

      index = 0
      do l = 0,lomaxfis(is)
        do iz = 1,nztfl(is,l)
          do m = -l,l
            index = index + 1
          enddo
        if (index .ge. iorb) goto 10
        enddo
      enddo
      call die('cspa: Error in orbital indexing in initguess')
10    continue

C Return if orbital is not the first zeta of its shell
      if (iz .ne. 1) return

C Assign initial guess.
C If 2 or less electrons, populate lowest s orbital
C If 8 or less electrons, populate lowest s and p  orbitals
C If 18 or less electrons, populate lowest s, p and d orbitals
C If 32 or less electrons, populate lowest s, p, d and f orbitals

      lmaxp = 0
      if (ne .le. 32) lmaxp = 3
      if (ne .le. 18) lmaxp = 2
      if (ne .le. 8) lmaxp = 1
      if (ne .le. 2) lmaxp = 0
      if (ne .gt. 32) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Cannot build initial guess in initguess.'
          write(6,*) '      Reason: Too many electrons for this routine'
        endif
        call die()
      endif 

      if (lmaxp .gt. lomaxfis(is)) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Cannot build initial guess in initguess.'
          write(6,*) '      Reason: Max. angular moment for atom ',ia,
     .               '      is not large enough'
        endif
        call die()
      endif 

      if (ne .gt. 32) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Cannot build initial guess in initguess.'
          write(6,*) '      Too many valence electrons in atom ',ia
        endif
        call die()
      endif 

      if (l .le. lmaxp) then
        cg = cgval
      endif

      end subroutine initguess
