! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine denmatlomem(c,eta,h,s,qs,nbasis,nbasisloc,nbands,
     .                       ncmax,nctmax,nfmax,nftmax,nhmax,nhijmax,
     .                       numc,listc,numct,listct,cttoc,numf,listf,
     .                       numft,listft,numh,listh,listhptr,numhij,
     .                       listhij,dm,edm,nbasisCloc,nspin,Node)
C *******************************************************************
C Subroutine to compute the Density and Energy Density matrices
C for the Order-N functional of Kim et al. (PRB 52, 1640 (95))
C (generalization of that proposed by Mauri et al, and Ordejon et al)
C
C Density Matrix:
C  D_mu,nu = 2 * C_i,mu * ( 2 * delta_i,j - S_i,j) * C_j,nu
C
C Energy Density Matrix:
C  E_mu,nu = 2 * C_i,mu * ( H_i,j + 2 * eta * (delta_i,j - S_i,j) ) * C_j,nu
C
C (The factor 2 is for spin)
C
C (See Ordejon et al, PRB 51, 1456 (95))
C
C The DM is normalized to the exact number of electrons!!!
C
C Written by P.Ordejon, Noviembre'96
C Modified by J.M.Soler, May'97
C Multiplication for cHc x c / cSc x c changed to suit parallel
C version and spatial decomposition added by J.D. Gale December '04
C Lower memory version written by J.D. Gale June '05
C ************************** INPUT **********************************
C real*8 c(ncmax,nbasis)      : Localized Wave Functions (sparse)
C real*8 eta                  : Fermi level parameter of Kim et al.
C real*8 h(nhmax)             : Hamiltonian matrix (sparse)
C real*8 s(nhmax)             : Overlap matrix (sparse)
C real*8 qs(nspin)            : Total number of electrons
C integer nbasis              : Global number of atomic orbitals
C integer nbasisloc           : Local number of atomic orbitals
C integer nbands              : Number of Localized Wave Functions
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer nctmax              : Max num of <>0 elements of each col of C
C integer nfmax               : Max num of <>0 elements of each row of 
C                               F = Ct x H
C integer nftmax              : Max num of <>0 elements of each col of F
C integer nhmax               : First dimension of listh and H, and maximum
C                               number of nonzero elements of each row of H
C integer nhijmax             : Maximum number of non-zero elements of each
C                               row of Hij
C integer numc(nbasis)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,nbasis) : Control vector of C matrix
C                              (list of nonzero elements of each row of C)
C integer numct(nbands)       : Control vector of C transpose matrix
C                              (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands): Control vector of C transpose matrix
C                              (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands) : Map from Ct to C indexing
C integer numf(nbands)        : Control vector of F matrix
C                               (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands) : Control vector of F matrix
C                               (list of <>0  elements of each row of F)
C integer numft(nbasis)       : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasis): Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer numh(nbasis)        : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listh(nhmax)        : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C integer listhptr(nbasis)    : Control vector of H matrix
C                               (pointer to start of row in listh/h/s)
C integer numhij(nbands)      : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C ************************* OUTPUT **********************************
C real*8 dm(nhmax,nbasis)     : Density Matrix
C real*8 edm(nhmax,nbasis)    : Energy density matrix
C *******************************************************************

      use precision
      use on_main,   only : ncG2L, ncP2T, ncT2P, nbL2G, nbandsloc
      use alloc,     only : re_alloc, de_alloc

#ifdef MPI
      use globalise, only: setglobalise, maxft2, numft2, listft2,
     .                     globalinitb, globalloadb1, globalcommb,
     .                     globalreloadb1, globalrezerob1,
     .                     globalisef
      use mpi_siesta, only: mpi_comm_world, mpi_sum
      use mpi_siesta, only: mpi_double_precision
#endif

      implicit none

      integer
     .  nbasis,nbasisloc,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,
     .  nbasisCloc,Node,nspin

      integer
     .  cttoc(nctmax,nbandsloc),listc(ncmax,nbasisCloc),
     .  listct(nctmax,nbandsloc),listf(nfmax,nbandsloc),
     .  listh(nhmax),listhptr(nbasisloc),listhij(nhijmax,nbandsloc),
     .  numc(nbasisCloc),numct(nbandsloc),numf(nbandsloc),
     .  numh(nbasisloc),numhij(nbandsloc)
     
      integer
     .  nftmax, listft(nftmax,nbasisCloc), numft(nbasisCloc)

      real(dp) ::
     .  c(ncmax,nbasisCloc,nspin),dm(nhmax,nspin),edm(nhmax,nspin),
     .  qs(nspin),eta(nspin),h(nhmax,nspin),s(nhmax)
     
      external
     .  timer

C Internal variales ..................................................
C   Notation hints:
C     m,n : basis orbital inexes (mu,nu)
C     i,j : band (and LWF) indexes
C     im  : index for LWF's of basis orbital m
C     mi  : index for basis orbitals of LWF i
C     nm  : index for basis orbitals connected to basis orbital m

      integer 
     .  i, il, in, indm, indn, im, is, j, jn, m, mi, mm, mn, 
     .  n, nh, ni, nm, nn

      real(dp), pointer :: cHrow(:), cSrow(:), chscrow(:), chccCol(:),
     &                     csccCol(:)

      real(dp) ::
     .  cim, cnj, chin, csin, cchccmn, ccsccmn,
     .  Hmn, Smn, qout(2), fact, spinfct

#ifdef MPI
      real(dp) ::
     .  qtmp(2)
      real(dp), pointer :: ctmp(:), chscrowsave(:,:), ftG(:,:),
     &                     fstG(:,:)
      integer ::
     .  MPIerror
#else
       real(dp), pointer ::   ft(:,:), fst(:,:)
#endif

C Start time counter
      call timer('denmat',1)

      if (nspin.eq.1) then
        spinfct = 2.0_dp
      else
        spinfct = 1.0_dp
      endif

C Allocate local arrays

      nullify( cHrow )
      call re_alloc( cHrow, 1, nbasis, name='cHrow',
     &               routine='denmatlomem' )
      nullify( cSrow )
      call re_alloc( cSrow, 1, nbasis, name='cSrow',
     &               routine='denmatlomem' )
      nullify( chscrow )
      call re_alloc( chscrow, 1, nbands, name='chscrow',
     &               routine='denmatlomem' )
      nullify( chccCol )
      call re_alloc( chccCol, 1, nbands, name='chccCol',
     &               routine='denmatlomem' )
      nullify( csccCol )
      call re_alloc( csccCol, 1, nbands, name='csccCol',
     &               routine='denmatlomem' )

#ifdef MPI
      nullify( chscrowsave )
      call re_alloc( chscrowsave, 1, nhijmax, 1, nbandsloc,
     &               name='chscrowsave',
     &               routine='denmatlomem' )
      nullify( ctmp )
      call re_alloc( ctmp, 1, nbands, name='ctmp',
     &               routine='denmatlomem' )
      nullify( ftG )
      call re_alloc( ftG, 1, MaxFt2, 1, nbasisCloc, name='ftG',
     &               routine='denmatlomem' )
      nullify( fstG )
      call re_alloc( fstG, 1, MaxFt2, 1, nbasisCloc, name='fstG',
     &               routine='denmatlomem' )
#else
      nullify( ft )
      call re_alloc( ft, 1, nftmax, 1, nbasisCloc, name='ft',
     &               routine='denmatlomem' )
      nullify( fst )
      call re_alloc( fst, 1, nftmax, 1, nbasisCloc, name='fst',
     &               routine='denmatlomem' )
#endif

C Initialize temporary arrays
      do m = 1,nbasis
        cHrow(m) = 0.0_dp
        cSrow(m) = 0.0_dp
      enddo
      do i = 1,nbands
        csccCol(i) = 0.0_dp
        chccCol(i) = 0.0_dp
        chscrow(i) = 0.0_dp
#ifdef MPI
        ctmp(i) = 0.0_dp
#endif
      enddo
#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          chscrowsave(j,i) = 0.0_dp
        enddo
      enddo
      do i = 1,nbasisCloc
        do j = 1,MaxFT2
          ftG(j,i) = 0.0_dp
          fstG(j,i) = 0.0_dp
        enddo
      enddo
#else
      ft = 0.0_dp
      fst = 0.0_dp
#endif

C Loop over spins
      do is = 1,nspin

#ifdef MPI
C Initialise communication arrays
        call globalinitB(1)
#endif

C Find cscc=(2-ct*S*c)*ct and chcc=(ct*H*c+2eta(1-ct*S*c))*ct.
        do il = 1,nbandsloc
          i = nbL2G(il)
      
C Find row i of cS=ct*S and cH=ct*H
          do mi = 1,numct(il)
            m = listct(mi,il)
            mm = ncT2P(m)
            if (mm.gt.0) then
              im = cttoc(mi,il)
              cim = c(im,m,is)
              indm = listhptr(mm)
              do nm = 1,numh(mm)
                n = listh(indm+nm)
                Hmn = H(indm+nm,is)
                Smn = S(indm+nm)
                cHrow(n) = cHrow(n) + cim*(Hmn-2.0_dp*eta(is)*Smn)
              enddo
            endif
          enddo

C Find row i of csc=2-ct*S*c and chc=ct*H*c+2eta(1-ct*S*c)

C Now use the list of nonzero elements of f=ct*H
          do ni = 1,numf(il)
            n = listf(ni,il)
            nn = ncG2L(n)
            chin = cHrow(n)
            do jn = 1,numc(nn)
              j = listc(jn,nn)
              cnj = c(jn,nn,is)
              chscrow(j) = chscrow(j) + chin*cnj
            enddo

C Restore cSrow and cHrow for next row
            cHrow(n) = 0.0_dp
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB1(il,nbands,chscrow)

C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)       
            chscrowsave(jn,il) = chscrow(j)
            chscrow(j) = 0.0_dp
          enddo

        enddo

C Globalise chc/csc
        call globalcommB(Node)

C Restore local chc/csc terms
        do il = 1,nbandsloc
          i = nbL2G(il)
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chscrow(j) = chscrowsave(jn,il)
          enddo
          
C Reload data after globalisation
          call globalreloadB1(il,nbands,chscrow,ctmp)

C Add diagonal terms - 2 and 2eta
          ctmp(i) = ctmp(i) + 2.0_dp * eta(is)
#else
          chscrow(i) = chscrow(i) + 2.0_dp * eta(is)
#endif

C Find row i of cscc=csc*ct and chcc=chc*ct. 
C Only the nonzero elements of f=cH will be required.
          do mi = 1,numct(il)
            m = listct(mi,il)
#ifdef MPI
            if (ncT2P(m).gt.0) then
              im = cttoc(mi,il)
              do in = 1,numft2(m)
                j = listft2(in,m)
                ftG(in,m)  = ftG(in,m)  + ctmp(j)*c(im,m,is)
              enddo
            endif
#else
            im = cttoc(mi,il)
            do in = 1,numft(m)
              j = listft(in,m)
              ft(in,m)  = ft(in,m)  + chscrow(j)*c(im,m,is)
            enddo
#endif
          enddo

C Reinitialise chcrow
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chscrow(j) = 0.0_dp
          enddo
#ifdef MPI
C Reinitialise ctmp
          call globalrezeroB1(il,nbands,ctmp)
#endif
        
        enddo

#ifdef MPI
C Initialise communication arrays
        call globalinitB(1)
#endif

C Find cscc=(2-ct*S*c)*ct and chcc=(ct*H*c+2eta(1-ct*S*c))*ct.
        do il = 1,nbandsloc
          i = nbL2G(il)
      
C Find row i of cS=ct*S and cH=ct*H
          do mi = 1,numct(il)
            m = listct(mi,il)
            mm = ncT2P(m)
            if (mm.gt.0) then
              im = cttoc(mi,il)
              cim = c(im,m,is)
              indm = listhptr(mm)
              do nm = 1,numh(mm)
                n = listh(indm+nm)
                Smn = S(indm+nm)
                cSrow(n) = cSrow(n) + cim * Smn
              enddo
            endif
          enddo

C Find row i of csc=2-ct*S*c and chc=ct*H*c+2eta(1-ct*S*c)

C Now use the list of nonzero elements of f=ct*H
          do ni = 1,numf(il)
            n = listf(ni,il)
            nn = ncG2L(n)
            csin = cSrow(n)
            do jn = 1,numc(nn)
              j = listc(jn,nn)
              cnj = c(jn,nn,is)
              chscrow(j) = chscrow(j) - csin * cnj
            enddo

C Restore cSrow and cHrow for next row
            cSrow(n) = 0.0_dp
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB1(il,nbands,chscrow)

C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)       
            chscrowsave(jn,il) = chscrow(j)
            chscrow(j) = 0.0_dp
          enddo

        enddo

C Globalise chc/csc
        call globalcommB(Node)

C Restore local chc/csc terms
        do il = 1,nbandsloc
          i = nbL2G(il)
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chscrow(j) = chscrowsave(jn,il)
          enddo
          
C Reload data after globalisation
          call globalreloadB1(il,nbands,chscrow,ctmp)

C Add diagonal terms - 2 and 2eta
          ctmp(i) = ctmp(i) + 2.0_dp
#else
          chscrow(i) = chscrow(i) + 2.0_dp
#endif

C Find row i of cscc=csc*ct and chcc=chc*ct. 
C Only the nonzero elements of f=cH will be required.
          do mi = 1,numct(il)
            m = listct(mi,il)
#ifdef MPI
            if (ncT2P(m).gt.0) then
              im = cttoc(mi,il)
              do in = 1,numft2(m)
                j = listft2(in,m)
                fstG(in,m) = fstG(in,m) + ctmp(j)*c(im,m,is)
              enddo
            endif
#else
            im = cttoc(mi,il)
            do in = 1,numft(m)
              j = listft(in,m)
              fst(in,m) = fst(in,m) + chscrow(j)*c(im,m,is)
            enddo
#endif
          enddo

C Reinitialise cscrow
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chscrow(j) = 0.0_dp
          enddo
#ifdef MPI
C Reinitialise ctmp
          call globalrezeroB1(il,nbands,ctmp)
#endif
        
        enddo

C Find dm=c*cscc and edm=c*chcc. Only the nonzero elements of H.
        do n = 1,nbasisloc
          nn = ncP2T(n)
      
C Use listft to expand a column of cscc
#ifdef MPI
          do in = 1,numft2(nn)
            i = listft2(in,nn)
            chccCol(i) = ftG(in,nn)
            csccCol(i) = fstG(in,nn)
          enddo
#else
          do in = 1,numft(nn)
            i = listft(in,nn)
            chccCol(i) = ft(in,nn)
            csccCol(i) = fst(in,nn)
          enddo
#endif

C Find column n of c*cscc and c*chcc
C Use that H is symmetric to determine required elements
          indn = listhptr(n)
          do mn = 1,numh(n)
            m = listh(indn+mn)
            mm = ncG2L(m)
            if (mm.gt.0) then
C Find element (m,n) of c*cscc and c*chcc
              ccsccmn = 0.0_dp
              cchccmn = 0.0_dp
              do im = 1,numc(mm)
                i = listc(im,mm)
                ccsccmn = ccsccmn + c(im,mm,is)*csccCol(i)
                cchccmn = cchccmn + c(im,mm,is)*chccCol(i)
              enddo
C Use that dm and edm are symmetric
              dm(indn+mn,is)  = spinfct*ccsccmn
              edm(indn+mn,is) = spinfct*cchccmn
            endif
          enddo
        
C Restore csccCol and chccCol for next column
          do in = 1,numft(nn)
            i = listft(in,nn)
            csccCol(i) = 0.0_dp
            chccCol(i) = 0.0_dp
          enddo
        enddo

C End of loop over spins
      enddo

C Normalize DM to exact charge .........................
C Calculate total output charge ...
      qout(1:2) = 0.0_dp
      if (nbasisloc.gt.0) then
        nh = listhptr(nbasisloc) + numh(nbasisloc)
      else
        nh = 0
      endif
      do is = 1,nspin
        do in = 1,nh
          qout(is) = qout(is) + dm(in,is) * s(in)
        enddo
      enddo

#ifdef MPI
C Globalise total charge
      qtmp(1:2) = qout(1:2)
      call MPI_AllReduce(qtmp,qout,nspin,MPI_double_precision,
     .  MPI_sum,MPI_Comm_World,MPIerror)
#endif

      if (Node.eq.0) then
        write(6,"(/a,2f12.4)") 
     .    'denmat: qtot (before DM normalization) = ',qout(1:nspin)
      endif

C Normalize ...
      do is = 1,nspin
C        if (dabs(qs(is)-qout(is)) .gt. 0.05_dp) then
          fact = qs(is) / qout(is)
          dm(1:nh,is) = fact*dm(1:nh,is)
          edm(1:nh,is) = fact*edm(1:nh,is)
C        endif
      enddo

C Deallocate local arrays
#ifdef MPI
      call de_alloc( fstG, name='fstG' )
      call de_alloc( ftG, name='ftG' )
      call de_alloc( ctmp, name='ctmp' )
      call de_alloc( chscrowsave, name='chscrowsave' )
#else
      call de_alloc( fst, name='fst' )
      call de_alloc( ft, name='ft' )
#endif
      call de_alloc( csccCol, name='csccCol' )
      call de_alloc( chccCol, name='chccCol' )
      call de_alloc( chscrow, name='chscrow' )
      call de_alloc( cSrow, name='cSrow' )
      call de_alloc( cHrow, name='cHrow' )

C Stop time counter and return ..................
      call timer('denmat',2)
      end
