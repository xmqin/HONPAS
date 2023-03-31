! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_ordern
!
!     We need an explicit interface since assumed-shape arrays are used
!
      public   :: ordern
      private
      ! These variables were moved here by Rogeli Grima
      integer, save :: iterm = 0
      integer, save :: maxnc = 1
      integer, save :: no_cl            ! set in cspa

      CONTAINS

      subroutine ordern(usesavelwf,ioptlwf,natoms,no_u,no_l,
     .                  lasto,isa,qa,rcoor,rh,cell,xa,iscf,istep,
     .                  itmax,ftol,eta,enum,nhmax,numh,listhptr,
     .                  listh,h,s,chebef,noeta,rcoorcp,beta,ipcheb,
     .                  dm,edm,Ecorrec,nspin,qs)
C ****************************************************************************
C Order-N solver of the Hamiltonian problem.
C It calls the appropriate routines, and returns the
C density matrix to the calling program
C Uses the funcional of Kim et al. (PRB 52, 1640 (95))
C
C Written by P.Ordejon, October'96
C ****************************** INPUT *************************************
C logical usesavelwf           : True = try to use saved lwf from disk
C integer ioptlwf              : Build LWF's according to:
C                                0 = Read blindly from disk
C                                1 = Functional of Kim et al.
C                                2 = Functional of Ordejon-Mauri
C integer natoms               : Number of atoms
C integer no_u                 : Number of basis orbitals globally
C integer no_l                 : Number of basis orbitals locally
C integer lasto(0:natoms)      : Index of last orbital of each atom
C integer isa(natoms)          : Species index of each atom
C real*8 qa(natoms)            : Neutral atom charge
C real*8 rcoor                 : Cutoff radius of Localized Wave Functions
C real*8 rh                    : Maximum cutoff radius of Hamiltonian matrix
C real*8 cell(3,3)             : Supercell vectors
C real*8 xa(3,natoms)          : Atomic coordinates
C integer iscf                 : SCF Iteration cycle (used to find 
C                                control vectors only if iscf=1)
C integer istep                : MD step
C integer itmax                : Maximum number of CG iterations
C real*8 ftol                  : Relative tolerance in CG minimization
C                                  (recomended: 1e-8)
C real*8 eta(2)                : Fermi level parameter of Kim et al.
C real*8 enum                  : Total number of electrons
C integer nhmax                : First dimension of listh and H, and maximum
C                                number of nonzero elements of each row of H
C integer numh(no_u)           : Control vector of H matrix
C                                (number of nonzero elements of each row of H)
C integer listh(nhmax)         : Control vector of H matrix
C                                (list of nonzero elements of each row of H)
C integer listhptr(no_u)       : Control vector of H matrix
C                                (pointer to start of row in listh / H)
C real*8 h(nhmax,nspin)        : Hamiltonian matrix (sparse)
C real*8 s(nhmax)              : Overlap matrix (sparse)
C logical chebef               : Compute the chemical potential
C logical noeta                : Use computed Chem.pot. instead of eta
C real*8 rcoorcp               : Cutoff radius to compute Fermi level by
C                                projection.
C integer ipcheb               : Order of Chebishev expansion to compute Ef
C real*8 beta                  : Inverse Temperature for Chebishev expansion
C integer nspin                : Number of spins (1 or 2)
C real*8  qs(2)                : Number of electrons per spin
C **************************** OUTPUT **************************************
C real*8 dm(nhmax,nspin)       : Density matrix (sparse)
C real*8 edm(nhmax,nspin)      : Energy Density matrix (sparse)
C real*8 Ecorrec               : Energy correction of Kim functional:
C                                eta * (etot-qs) , where qs is the charge
C                                of the Order-N solution
C **************************** BEHAVIOUR ***********************************
C If istep=1 and iscf=1, an initial guess is build. Otherwise, the LWF's of
C the former time steps are used,  extrapolated between MD steps.
C **************************************************************************** 

C  Modules
      use precision, only: dp
      use alloc
      use on_main,    only: c, cold, listc, listcold, numc, numcold
      use on_main,    only: ncg2l, ncl2g, ncp2t, nct2p
      use on_main,    only: xi, g, hg

      use parallel,  only : IONode, Node, NOdes
      use siesta_cml
      use sys,       only : die
      use units,     only : eV
      use m_chempot, only : chempot
      use m_cgwf,    only : cgwf
#ifdef MPI
      use globalise,  only: setglobalise, globaliseC
      use m_mpi_utils, only: globalize_sum, globalize_max  
#endif

      implicit none

      integer, intent(in) :: 
     .  ioptlwf, ipcheb, iscf, istep, itmax, natoms, 
     .  no_u, no_l, nhmax, nspin

      integer, intent(in) ::
     .  isa(natoms), lasto(0:natoms), listh(nhmax), 
     .  listhptr(no_l), numh(no_l)

      real(dp) , intent(inout) :: eta(2)
      real(dp) , intent(in) ::
     .  beta, cell(3,3),
     .  enum, ftol, h(nhmax,nspin), qa(natoms), rcoor, rcoorcp,
     .  rh, s(nhmax), xa(3,natoms), qs(nspin)

      real(dp) , intent(out) :: dm(nhmax,nspin), edm(nhmax,nspin)
      real(dp) , intent(out) :: Ecorrec

      logical, intent(in)    :: chebef, noeta, usesavelwf

C Internal variables ..........................................................
      integer
     .  in, io, iopt, iord, is, iter, ncmax, nctmax, 
     .  nfmax, nftmax, nhijmax, nbands, nh, gen_istep

#ifdef MPI
      integer :: ntmp
      real(dp) ::
     .  rtmp(2)
#endif

      real(dp), pointer, save ::
     .  aux(:,:)

      real(dp) ::
     .  fe, qtot(nspin)

      real(dp) ::
     .  chpot,emax,emin

      logical
     .  itest, found

      logical, save :: frstme = .true.

      save 
     .  nbands
     
C ...................
#ifdef DEBUG
      call write_debug( '    PRE ordern' )
#endif

*     call timer( 'ordern', 1 )

      if (frstme) then
C Nullify pointers
        nullify(aux)
        nullify(c)
        nullify(cold)
        nullify(listc)
        nullify(listcold)
        nullify(numc)
        nullify(numcold)
        nullify(ncG2L)
        nullify(ncL2G)
        nullify(ncP2T)
        nullify(ncT2P)
        nullify(xi)
        nullify(g)
        nullify(hg)

C Set default sizes for order N arrays
        call re_alloc(listc,1,maxnc,1,no_l,name='listc')
        call re_alloc(listcold,1,maxnc,1,no_l,name='listcold')
        call re_alloc(numc,1,no_u,name='numc')
        call re_alloc(numcold,1,no_u,name='numcold')
        call re_alloc(ncG2L,1,no_u,name='ncG2L')
        call re_alloc(ncL2G,1,no_u,name='ncL2G')
        call re_alloc(ncP2T,1,no_u,name='ncP2T')
        call re_alloc(ncT2P,1,no_u,name='ncT2P')
        call re_alloc(c,1,maxnc,1,no_l,1,nspin,name='c')
        call re_alloc(cold,1,maxnc,1,no_l,1,nspin,name='cold')
        call re_alloc(xi,1,maxnc,1,no_l,1,nspin,name='xi')
        call re_alloc(g,1,maxnc,1,no_l,1,nspin,name='g')
        call re_alloc(hg,1,maxnc,1,no_l,1,nspin,name='hg')
        call re_alloc(aux,1,2,1,no_u,name='aux')

      endif

C ............................

      if (IONode) then
        write(6,"(/a,f12.4)") 'ordern: enum =', enum
      endif
      if (cml_p) call cmlAddProperty(xf=mainXML, value=enum, 
     .     dictref='siesta:enum', title='Number of electrons',
     .     units='cmlUnits:countable')

C  Check if options are compatible
      if (ioptlwf.eq.0.and.(.not.usesavelwf)) then
        if (IONode) then
          write(6,"(/a)") 'ordern: ERROR: You must use LWF files.'
          write(6,"(a)") '        If you choose ON.functional = Files'
          write(6,"(a)") '        Please set ON.UseSaveLWF = True'
        endif
        call die('LWF file error')
      endif

C  If iscf = 1 (that is, if we are in a new MD step), find out initial
C  structure of localized wave functions, and initial guess .............
      if (iscf.eq.1) then
        if (istep.eq.1) then
          iopt = 0
        else
          iopt = 1
        endif
C  Call cspa to initialise Wannier function related data structures
        call cspa(ioptlwf,iopt,natoms,no_u,no_l,lasto,isa,
     .            qa,rcoor,rh,cell,xa,nhmax,numh,listh,
     .            listhptr,maxnc,ncmax,nctmax,nfmax,nftmax,
     .            nhijmax,nbands,no_cl,nspin,Node)
      endif

C  Resize arrays in case maxnc has changed
!!!! AG: Shrink if needed (maxnc increased in chunks of 50 in cspa)
!!
      call re_alloc(listc,1,maxnc,1,no_cl,name='listc')
      call re_alloc(listcold,1,maxnc,1,no_l,name='listcold')
      call re_alloc(c,1,maxnc,1,no_cl,1,nspin,name='c')
      call re_alloc(cold,1,maxnc,1,no_l,1,nspin,name='cold')
      call re_alloc(xi,1,maxnc,1,no_cl,1,nspin,name='xi')
      call re_alloc(g,1,maxnc,1,no_cl,1,nspin,name='g')
      call re_alloc(hg,1,maxnc,1,no_cl,1,nspin,name='hg')

      if (iscf.eq.1.and.istep.eq.1) then
        if (usesavelwf) then
C  Read Wannier functions if file is present
          call iolwf( 'read', no_u, no_cl, no_l, maxnc, 
     .                found, nspin)
          if (found) then
C  Find out number of bands
            nbands = 0
            do io = 1,no_l
              do in = 1,numcold(io)
                nbands = max( nbands, listcold(in,io) )
              enddo
            enddo
#ifdef MPI
            call Globalize_max(nbands,ntmp)
            nbands = ntmp
#endif
C  Copy coefficients read in from Cold to C ready for extrapolon
            do is = 1,nspin
              do io = 1,no_l
                do in = 1,numcold(io)
                  c(in,io,is) = cold(in,io,is)
                enddo
              enddo
            enddo

          endif
        endif
      endif

C Calculate Chemical Potential, Max and Min eigenvalues, energy gap
C and HOMO and LUMO levels ..........
      if (chebef) then

        call chempot(h,s,numh,listhptr,listh,rcoorcp,ipcheb,beta,
     .               lasto,cell,xa,enum,no_u,no_l,natoms,
     .               nhmax,chpot,emax,emin)

        if (ionode) then
          write(6,'(a,f8.4,a)') 'ordern:   Chemical Potential = ',
     .                         chpot/eV,' eV'
          write(6,'(a,f8.4,a)') 'ordern:   Maximum Eigenvalue = ',
     .                         emax/eV,' eV'
          write(6,'(a,f8.4,a)') 'ordern:   Minimum Eigenvalue = ',
     .                         emin/eV,' eV'
        endif
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, value=chpot/eV, 
     .          dictref='siesta:chemical_pot',
     .          units='siestaUnits:ev', title='Chemical potential')
          call cmlAddProperty(xf=mainXML, value=emax/eV,
     .         dictref='siesta:max_eigenv',
     .     units='siestaUnits:ev', title='Maximum eigenvalue')
          call cmlAddProperty(xf=mainXML, value=emin/eV, 
     .         dictref='siesta:min_eigenv',
     .     units='siestaUnits:ev', title='Minimum eigenvalue')
        endif

        if (noeta) eta = chpot
      endif

C  Extrapolate wave funcions from those of former time steps if 
C  iscf = 1 (that is, if we are in a new MD step) ...................
      if (iscf .eq. 1)  then
        if (ionode) then
          write(6,"(/a,i3)") 'ordern: ioptlwf =',ioptlwf
        endif
        if (ioptlwf.eq.0) then
          do io = 1,no_l
            numc(io) = numcold(io)
            do in = 1,numcold(io)
              listc(in,io) = listcold(in,io)
            enddo
          enddo
        endif
        iord = 1
        if (iterm .gt. 50) iord = 0
C If LWF's have just been read from disk, 
C call extrapol with istep = 2 and iord = 1
C to make it update the structure of c, if needed
        gen_istep = istep
        if (istep.eq.1 .and. usesavelwf .and. found) then
          gen_istep = 2
          iord = 0
        endif
        call extrapolon(gen_istep,iord,nspin,no_u,no_cl,maxnc,
     .                  numc,listc,aux,numcold,listcold,cold,c,
     .                  no_l)

#ifdef MPI
        if (gen_istep.ne.1) then

C Zero parts of C that current Node does not have primary responsibility for
          do is = 1,nspin
            do io = no_l+1,no_cl
              do in = 1,numc(io)
                c(in,io,is) = 0.0d0
              enddo
            enddo
          enddo

C Globalise C to get values for rows from other nodes
          call setglobalise(no_u,no_l,no_cl,nhmax,numh,
     .                    listh,listhptr,Node,Nodes)
          call globaliseC(no_cl,maxnc,numc,c,nspin,Node)

        endif
#endif

      endif

C Call the CG routines ...............................................
      if (iscf .eq. 1) iterm = 0
      call cgwf(iscf,itmax,ftol,eta,qs,no_u,no_l,nbands,
     .          nhmax,numh,listh,listhptr,maxnc,h,s,fe,iter,dm,
     .          edm,no_cl,nspin,Node,Nodes)
      iterm = max(iterm,iter)

C Calculate correction to the total energy from Kim's functional .....
C First calculate total charge of the solution
      qtot(1:nspin) = 0.0d0
      if (no_l.gt.0) then
        nh = listhptr(no_l) + numh(no_l)
        do is = 1,nspin
          do in = 1,nh
            qtot(is) = qtot(is) + dm(in,is) * s(in)
          enddo
        enddo
      endif
#ifdef MPI
C Global reduction of qtot
      call Globalize_sum(qtot(1:nspin),rtmp(1:nspin))
      qtot(1:nspin) = rtmp(1:nspin)
#endif

      Ecorrec = 0.0d0
      do is = 1,nspin
        Ecorrec = Ecorrec + eta(is) * (qs(is) - qtot(is))
      enddo
      if (ionode) then
        write(6,'(a,f12.4)') 
     .       'ordern: qtot (after  DM normalization) = ', sum(qtot)
      endif
      if (cml_p) then
        call cmlAddProperty(xf=mainXML, value=sum(qtot), 
     .       dictref='siesta:qtot', title='qtot',
     .       units='siestaUnits:e')
      endif

C Save LWF's to disk .................................................
      call iolwf( 'write', no_u, no_cl, no_l, maxnc, 
     .            found, nspin)
      frstme = .false.

#ifdef DEBUG
      call write_debug( '    POS ordern' )
#endif
      end subroutine ordern

      end module m_ordern
