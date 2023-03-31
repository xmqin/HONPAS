! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_denmat
      public :: denmat
      CONTAINS
      subroutine denmat(c,eta,h,s,qs,no_u,no_l,nbands,ncmax,
     .                  nctmax,nfmax,nftmax,nhmax,nhijmax,numc,listc,
     .                  numct,listct,cttoc,numf,listf,numft,listft,
     .                  numh,listh,listhptr,numhij,listhij,dm,edm,
     .                  no_cl,nspin,Node)
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
C ************************** INPUT **********************************
C real*8 c(ncmax,no_u)      : Localized Wave Functions (sparse)
C real*8 eta                  : Fermi level parameter of Kim et al.
C real*8 h(nhmax)             : Hamiltonian matrix (sparse)
C real*8 s(nhmax)             : Overlap matrix (sparse)
C real*8 qs(nspin)            : Total number of electrons
C integer no_u              : Global number of atomic orbitals
C integer no_l           : Local number of atomic orbitals
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
C integer numc(no_u)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,no_u) : Control vector of C matrix
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
C integer numft(no_u)       : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,no_u): Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer numh(no_u)        : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listh(nhmax)        : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C integer listhptr(no_u)    : Control vector of H matrix
C                               (pointer to start of row in listh/h/s)
C integer numhij(nbands)      : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C ************************* OUTPUT **********************************
C real*8 dm(nhmax,no_u)     : Density Matrix
C real*8 edm(nhmax,no_u)    : Energy density matrix
C *******************************************************************

      use precision, only : dp
      use on_main,   only : ncG2L, ncP2T, ncT2P, nbL2G, nbandsloc
      use m_mpi_utils, only : globalize_sum
      use alloc,     only : re_alloc, de_alloc

#ifdef MPI
      use globalise, only: setglobalise, maxft2, numft2, listft2,
     .                     globalinitb, globalloadb2, globalcommb,
     .                     globalreloadb2, globalrezerob2,
     .                     globalisef
#endif

      implicit none

      integer, intent(in) ::
     .  no_u,no_l,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,
     .  no_cl,Node,nspin

      integer, intent(in) ::
     .  cttoc(nctmax,nbandsloc),listc(ncmax,no_cl),
     .  listct(nctmax,nbandsloc),listf(nfmax,nbandsloc),
     .  listh(nhmax),listhptr(no_l),listhij(nhijmax,nbandsloc),
     .  numc(no_cl),numct(nbandsloc),numf(nbandsloc),
     .  numh(no_l),numhij(nbandsloc)
     
      integer, intent(in) ::
     .  nftmax,
     .  listft(nftmax,no_cl), numft(no_cl)

      real(dp), intent(in) ::
     .  c(ncmax,no_cl,nspin),
     .  qs(nspin),eta(nspin),h(nhmax,nspin),s(nhmax)

      real(dp), intent(out) :: dm(nhmax,nspin),edm(nhmax,nspin)

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

      real(dp), pointer ::
     .  cHrow(:), cSrow(:), chcrow(:), cscrow(:), chccCol(:),
     .  csccCol(:)

      real(dp) ::
     .  cim, cnj, chin, csin, cchccmn, ccsccmn,
     .  Hmn, Smn, qout(2), fact, spinfct

#ifdef MPI
      real(dp) ::  qtmp(nspin)
      real(dp), pointer ::
     .  ctmp(:,:), chcrowsave(:,:), cscrowsave(:,:), 
     .  ftG(:,:), fstG(:,:)
#else
       real(dp), pointer ::  ft(:,:), fst(:,:)
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
      call re_alloc( cHrow, 1, no_u, name='cHrow', routine='denmat' )
      nullify( cSrow )
      call re_alloc( cSrow, 1, no_u, name='cSrow', routine='denmat' )
      nullify( chcrow )
      call re_alloc( chcrow, 1, nbands, name='chcrow',
     &               routine='denmat' )
      nullify( cscrow )
      call re_alloc( cscrow, 1, nbands, name='cscrow',
     &               routine='denmat' )
      nullify( chccCol )
      call re_alloc( chccCol, 1, nbands, name='chccCol',
     &               routine='denmat' )
      nullify( csccCol )
      call re_alloc( csccCol, 1, nbands, name='csccCol',
     &               routine='denmat' )

#ifdef MPI
      nullify( chcrowsave )
      call re_alloc( chcrowsave, 1, nhijmax, 1, nbandsloc,
     &               name='chcrowsave', routine='denmat' )
      nullify( cscrowsave )
      call re_alloc( cscrowsave, 1, nhijmax, 1, nbandsloc,
     &               name='cscrowsave', routine='denmat' )
      nullify( ctmp )
      call re_alloc( ctmp, 1, nbands, 1, 2,
     &               name='ctmp', routine='denmat' )
      nullify( ftG )
      call re_alloc( ftG, 1, MaxFt2, 1, no_cl,
     &               name='ftG', routine='denmat' )
      nullify( fstG )
      call re_alloc( fstG, 1, MaxFt2, 1, no_cl,
     &               name='fstG', routine='denmat' )
#else
      nullify( ft )
      call re_alloc( ft, 1, nftmax, 1, no_cl,
     &               name='ft', routine='denmat' )
      nullify( fst )
      call re_alloc( fst, 1, nftmax, 1, no_cl,
     &               name='fst', routine='denmat' )
#endif

C Initialize temporary arrays
      do m = 1,no_u
        cHrow(m) = 0.0_dp
        cSrow(m) = 0.0_dp
      enddo
      do i = 1,nbands
        csccCol(i) = 0.0_dp
        chccCol(i) = 0.0_dp
        cscrow(i) = 0.0_dp
        chcrow(i) = 0.0_dp
#ifdef MPI
        ctmp(i,1) = 0.0_dp
        ctmp(i,2) = 0.0_dp
#endif
      enddo
#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          cscrowsave(j,i) = 0.0_dp
          chcrowsave(j,i) = 0.0_dp
        enddo
      enddo
      do i = 1,no_cl
        do j = 1,MaxFT2
          ftG(j,i) = 0.0_dp
          fstG(j,i) = 0.0_dp
        enddo
      enddo
#else
      ft=0.0_dp
      fst=0.0_dp
#endif

C Loop over spins
      do is = 1,nspin

#ifdef MPI
C Initialise communication arrays
        call globalinitB(2)
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
                cHrow(n) = cHrow(n) + cim * Hmn
                cSrow(n) = cSrow(n) + cim * Smn
              enddo
            endif
          enddo

C Find row i of csc=2-ct*S*c and chc=ct*H*c+2eta(1-ct*S*c)

C Now use the list of nonzero elements of f=ct*H
          do ni = 1,numf(il)
            n = listf(ni,il)
            nn = ncG2L(n)
            csin = - cSrow(n)
            chin = cHrow(n) - 2.0_dp*eta(is)*cSrow(n)
            do jn = 1,numc(nn)
              j = listc(jn,nn)
              cnj = c(jn,nn,is)
              chcrow(j) = chcrow(j) + chin * cnj
              cscrow(j) = cscrow(j) + csin * cnj
            enddo

C Restore cSrow and cHrow for next row
            cSrow(n) = 0.0_dp
            cHrow(n) = 0.0_dp
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB2(il,nbands,chcrow,cscrow)

C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)       
            chcrowsave(jn,il) = chcrow(j)
            cscrowsave(jn,il) = cscrow(j)
            chcrow(j) = 0.0_dp
            cscrow(j) = 0.0_dp
          enddo

        enddo

C Globalise chc/csc
        call globalcommB(Node)

C Restore local chc/csc terms
        do il = 1,nbandsloc
          i = nbL2G(il)
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chcrow(j) = chcrowsave(jn,il)
            cscrow(j) = cscrowsave(jn,il)
          enddo
          
C Reload data after globalisation
          call globalreloadB2(il,nbands,chcrow,cscrow,ctmp)

C Add diagonal terms - 2 and 2eta
          ctmp(i,1) = ctmp(i,1) + 2.0_dp * eta(is)
          ctmp(i,2) = ctmp(i,2) + 2.0_dp
#else
          chcrow(i) = chcrow(i) + 2.0_dp * eta(is)
          cscrow(i) = cscrow(i) + 2.0_dp
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
                ftG(in,m)  = ftG(in,m)  + ctmp(j,1)*c(im,m,is)
                fstG(in,m) = fstG(in,m) + ctmp(j,2)*c(im,m,is)
              enddo
            endif
#else
            im = cttoc(mi,il)
            do in = 1,numft(m)
              j = listft(in,m)
              ft(in,m)  = ft(in,m)  + chcrow(j)*c(im,m,is)
              fst(in,m) = fst(in,m) + cscrow(j)*c(im,m,is)
            enddo
#endif
          enddo

C Reinitialise chcrow/cscrow
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            chcrow(j) = 0.0_dp
            cscrow(j) = 0.0_dp
          enddo
#ifdef MPI
C Reinitialise ctmp
          call globalrezeroB2(il,nbands,ctmp)
#endif
        
        enddo

C Find dm=c*cscc and edm=c*chcc. Only the nonzero elements of H.
        do n = 1,no_l
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
      if (no_l.gt.0) then
        nh = listhptr(no_l) + numh(no_l)
      else
        nh = 0
      endif
      do is = 1,nspin
        do in = 1,nh
          qout(is) = qout(is) + dm(in,is) * s(in)
        enddo
      enddo

C Globalize total charge
#ifdef MPI
      qtmp(1:nspin) = qout(1:nspin)
      call Globalize_sum(qtmp(1:nspin),qout(1:nspin))
#endif

!! AG  #ifdef MPI
!! AG      qtmp(1:2) = qout(1:2)
!! AG      call MPI_AllReduce(qtmp,qout,nspin,MPI_double_precision,
!! AG     .  MPI_sum,MPI_Comm_World,MPIerror)
!! AG #endif

      if (Node.eq.0) then
        write(6,"(/a,2f12.4)") 
     .    'denmat: qtot (before DM normalization) = ',qout(1:nspin)
      endif

C Normalize ...
      do is = 1,nspin
C        if (dabs(qs(is)-qout(is)) .gt. 0.05_dp) then
          fact = qs(is) / qout(is)
          call dscal(nh,fact,dm(1,is),1)
          call dscal(nh,fact,edm(1,is),1)
C        endif
      enddo

C Deallocate local arrays
#ifdef MPI
      call de_alloc( fstG, name='fstG' )
      call de_alloc( ftG, name='ftG' )
      call de_alloc( ctmp, name='ctmp' )
      call de_alloc( cscrowsave, name='cscrowsave' )
      call de_alloc( chcrowsave, name='chcrowsave' )
#else
      call de_alloc( fst, name='fst' )
      call de_alloc( ft, name='ft' )
#endif
      call de_alloc( csccCol, name='csccCol' )
      call de_alloc( chccCol, name='chccCol' )
      call de_alloc( cscrow, name='cscrow' )
      call de_alloc( chcrow, name='chcrow' )
      call de_alloc( cSrow, name='cSrow' )
      call de_alloc( cHrow, name='cHrow' )

C Stop time counter and return ..................
      call timer('denmat',2)
      end subroutine denmat
      end module m_denmat
