! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine gradientlomem(c,eta,qs,h,s,nbasis,nbasisloc,nbands,
     .                         ncmax,nctmax,nfmax,nftmax,nhmax,nhijmax,
     .                         numc,listc,numct,listct,cttoc,numf,listf,
     .                         numft,listft,fttof,numh,listh,listhptr,
     .                         numhij,listhij,f,fs,grad,ener,nbasisCloc,
     .                         nspin,Node)

C ************************************************************************
C Finds the energy and gradient at point C.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C Modified: J.M.Soler. 30/04/97
C Adapted from gradient by J.D. Gale. June 2005
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasis)       : Current point (wave function coeff.
C                                  in sparse form)
C real*8 eta(nspin)            : Fermi level parameter of Kim et al.
C real*8 qs(nspin)             : Total number of electrons
C real*8 h(nhmax)              : Hamiltonian matrix (sparse)
C real*8 s(nhmax)              : Overlap matrix (sparse)
C integer nbasis               : Global number of basis orbitals
C integer nbasisloc            : Local number of basis orbitals
C integer nbands               : Number of LWF's
C integer ncmax                : Max num of <>0 elements of each row of C
C integer nctmax               : Max num of <>0 elements of each col of C
C integer nfmax                : Max num of <>0 elements of each row of 
C                                   F = Ct x H
C integer nftmax               : Max num of <>0 elements of each col of F
C integer nhmax                : Max num of <>0 elements of each row of H
C integer nhijmax              : Max num of <>0 elements of each row of 
C                                   Hij=Ct x H x C
C integer numc(nbasis)         : Control vector of C matrix
C                                (number of <>0  elements of each row of C)
C integer listc(ncmax,nbasis)  : Control vector of C matrix 
C                               (list of <>0  elements of each row of C)
C integer numct(nbands)        : Control vector of C transpose matrix
C                               (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands) : Control vector of C transpose matrix
C                               (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands)  : Map from Ct to C indexing
C integer numf(nbands)         : Control vector of F matrix
C                                (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands)  : Control vector of F matrix 
C                                (list of <>0  elements of each row of F)
C integer numft(nbasis)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasis) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasis)  : Map from Ft to F indexing
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer listhptr(nbasis)     : Control vector of H matrix 
C                               (pointer to start of rows in listh/h/s)
C integer numhij(nbands)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C real*8 f(nfmax,nbands)       : Auxiliary space
C real*8 fs(nfmax,nbands)      : Auxiliary space
C ***************************** OUTPUT ***********************************
C real*8 ener                  : Energy at point C
C real*8 grad(ncmax,nbasis)    : Gradient of functional (sparse)
C ************************************************************************

      use precision
      use on_main,   only : ncG2L, ncT2P, nbL2G, nbG2L, nbandsloc
      use alloc,     only : re_alloc, de_alloc
      use globalise
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world, mpi_sum
      use mpi_siesta, only: mpi_double_precision
#endif

      implicit none

      integer
     .  nbasis,nbasisloc,nbands,ncmax,nctmax,nfmax,nftmax,nhmax,
     .  nhijmax,nbasisCloc,nspin,Node

      integer
     .  cttoc(nctmax,nbandsloc),fttof(nftmax,nbasisCloc),
     .  listc(ncmax,nbasisCloc),listct(nctmax,nbandsloc),
     .  listf(nfmax,nbandsloc),listft(nftmax,nbasisCloc),
     .  listh(nhmax),listhptr(nbasisloc),listhij(nhijmax,nbandsloc),
     .  numc(nbasisCloc),numct(nbandsloc),numf(nbandsloc),
     .  numft(nbasisCloc),numh(nbasisloc),numhij(nbandsloc)

      real(dp) ::
     .  c(ncmax,nbasisCloc,nspin),ener,eta(nspin),qs(nspin),
     .  grad(ncmax,nbasisCloc,nspin),h(nhmax,nspin),s(nhmax),
     .  f(nfmax,nbandsloc),fs(nfmax,nbandsloc)

C Internal variables ......................................................

      integer
     .  i,ik,il,in,indk,is,j,jk,jl,jn,k,kk,kn,mu,muk

#ifdef MPI
      integer ::
     .  MPIerror
      real(dp) ::
     .  ftmp(2), ftmp2(2)
      real(dp), pointer, save ::
     .  buxg(:), bux1save(:,:), bux2save(:,:), ftG(:,:), fstG(:,:)
#endif

      real(dp), pointer, save ::
     .  aux1(:), aux2(:), bux(:), bux1(:), ft(:,:), fst(:,:)

      real(dp) ::
     .  a0,b0,p0,func1,func2,spinfct

C..................

      call timer('gradient',1)
      
C Allocate local arrays
      call re_alloc( aux1, 1, nbasis, 'aux1', 'gradient' )
      call re_alloc( aux2, 1, nbasis, 'aux2', 'gradient' )
      call re_alloc( bux, 1, nbands, 'bux', 'gradient' )
      call re_alloc( bux1, 1, nbands, 'bux1', 'gradient' )
      call re_alloc( ft, 1, nftmax, 1, nbasisCloc, 'ft', 'gradient' )
      call re_alloc( fst, 1, nftmax, 1, nbasisCloc, 'fst', 'gradient' )
#ifdef MPI
      call re_alloc( ftG, 1, MaxFt2, 1, nbasisCloc, 'ftG', 'gradient' )
      call re_alloc( fstG, 1, MaxFt2, 1, nbasisCloc, 'fstG', 'gradient')
      call re_alloc( bux1save, 1, nhijmax, 1, nbandsloc, 'bux1save',
     &               'gradient' )
      call re_alloc( bux2save, 1, nhijmax, 1, nbandsloc, 'bux2save',
     &               'gradient' )
      call re_alloc( buxg, 1, nbands, 'buxg', 'gradient')
#endif

C Initialize output and auxiliary varialbles ...............................

      if (nspin.eq.1) then
        spinfct = 2.0d0
      else
        spinfct = 1.0d0
      endif

      ener = 0.0d0

      do i = 1,nbasis
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux(i) = 0.0d0
#ifdef MPI
        buxg(i) = 0.0d0
#endif
      enddo

#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          bux1save(j,i) = 0.0_dp
          bux2save(j,i) = 0.0_dp
        enddo
      enddo
#endif

C Initialise gradient vector
      grad = 0.0_dp

C Calculate Functional .....................................................

C Loop over spin
      do is = 1,nspin

C Initialise variables
        func1 = 0.0_dp
        func2 = 0.0_dp

C F=CtH  ---> JMS: F=Ct*(H-eta*S)
C Fs=CtS
        do il = 1,nbandsloc
          do in = 1,numct(il)
            k = listct(in,il)
            kk = ncT2P(k)
            if (kk.gt.0) then
              ik = cttoc(in,il)
              p0 = c(ik,k,is)
              indk = listhptr(kk)
              do kn = 1,numh(kk)
                aux1(listh(indk+kn)) = aux1(listh(indk+kn)) + 
     .            p0*( h(indk+kn,is) - eta(is)*s(indk+kn) )
                aux2(listh(indk+kn)) = aux2(listh(indk+kn)) + 
     .            p0*s(indk+kn)
              enddo
            endif
          enddo
          do in = 1,numf(il)
            k = listf(in,il)
            f(in,il) = aux1(k)
            fs(in,il) = aux2(k)
            aux1(k) = 0.0_dp
            aux2(k) = 0.0_dp
          enddo
        enddo

C-JMS Find transpose of F and Fs
        do mu = 1,nbasisCloc
          do muk = 1,numft(mu)
            j = listft(muk,mu)
            jl = nbG2L(j)
            if (jl.gt.0) then
              jk = fttof(muk,mu)
              ft(muk,mu) = f(jk,jl)
              fst(muk,mu) = fs(jk,jl)
            endif
          enddo
        enddo

#ifdef MPI
C Globalise ft/fst
        call globaliseF(nbasisCloc,nbands,nftmax,numft,listft,
     .                  ft,ftG,Node)
        call globaliseF(nbasisCloc,nbands,nftmax,numft,listft,
     .                  fst,fstG,Node)
#endif

C------------
C H - eta*S -
C------------
#ifdef MPI
C Initialise communication arrays
        call globalinitB(1)
#endif

C Sij=CtSC
C multiply FxC and FsxC row by row
        do il = 1,nbandsloc
          do in = 1,numf(il)
            k = listf(in,il)
            kk = ncG2L(k)
            a0 = f(in,il)
            do kn = 1,numc(kk)
              bux(listc(kn,kk)) = bux(listc(kn,kk)) + a0*c(kn,kk,is)
            enddo
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB1(il,nbands,bux)
#endif

#ifdef MPI
C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1save(jn,il) = bux(j)
            bux(j) = 0.0_dp
          enddo

        enddo

C Globalise Hij-eta*Sij
        call globalcommB(Node)

C Restore local bux1/2 terms
        do il = 1,nbandsloc
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux(j) = bux1save(jn,il)
          enddo

C Reload data after globalisation
          call globalreloadB1(il,nbands,bux,buxg)
#endif

C Multiply Hij x Fs row by row
C (only products of neccesary elements)
          do ik = 1,numct(il)
            mu = listct(ik,il)
            kk = ncT2P(mu)
            if (kk.gt.0) then
              a0 = 0.0d0
#ifdef MPI
              do muk = 1,numft2(mu)
                j = listft2(muk,mu)
                a0 = a0 + buxg(j)*fstG(muk,mu) 
#else
              do muk = 1,numft(mu)
                j = listft(muk,mu)
                a0 = a0 + bux(j)*fst(muk,mu) 
#endif
              enddo
              grad(cttoc(ik,il),mu,is) = - spinfct*a0
            endif
          enddo

C Calculate first energy term 
          i = nbL2G(il)
          func1 = func1 + bux(i)

#ifdef MPI
C Reinitialise buxg
          call globalrezeroB1(il,nbands,buxg)
#endif
C Reinitialise bux1 / bux2
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux(j) = 0.0_dp
          enddo

        enddo
C----
C S -
C----
#ifdef MPI
C Initialise communication arrays
        call globalinitB(1)
#endif

        do il = 1,nbandsloc
          do in = 1,numf(il)
            k = listf(in,il)
            kk = ncG2L(k)
            b0 = fs(in,il)
            do kn = 1,numc(kk)
              bux(listc(kn,kk)) = bux(listc(kn,kk)) + b0*c(kn,kk,is)
            enddo
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB1(il,nbands,bux)
#endif

#ifdef MPI
C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux2save(jn,il) = bux(j)
            bux(j) = 0.0d0
          enddo

        enddo

C Globalise Sij
        call globalcommB(Node)

C Restore local bux terms
        do il = 1,nbandsloc
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1(j) = bux1save(jn,il)
            bux(j)  = bux2save(jn,il)
          enddo

C Reload data after globalisation
          call globalreloadB1(il,nbands,bux,buxg)
#endif

C Multiply Sij x F row by row
C (only products of neccesary elements)
          do ik = 1,numct(il)
            mu = listct(ik,il)
            kk = ncT2P(mu)
            if (kk.gt.0) then
              a0 = 0.0d0
#ifdef MPI
              do muk = 1,numft2(mu)
                j = listft2(muk,mu)
                a0 = a0 + buxg(j)*ftG(muk,mu)
#else
              do muk = 1,numft(mu)
                j = listft(muk,mu)
                a0 = a0 + bux(j)*ft(muk,mu)
#endif
              enddo
              grad(cttoc(ik,il),mu,is) = grad(cttoc(ik,il),mu,is) - 
     .          spinfct*a0
            endif
          enddo

C Calculate second energy term 
          do jn = 1,numhij(il)
            j = listhij(jn,il)
#ifdef MPI
            func2 = func2 + bux1(j)*buxg(j)
#else
            func2 = func2 + bux1(j)*bux(j)
#endif
          enddo

#ifdef MPI
C Reinitialise buxg
          call globalrezeroB1(il,nbands,buxg)
#endif
C Reinitialise bux1 / bux2
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1(j) = 0.0d0
            bux(j) = 0.0d0
          enddo

        enddo

#ifdef MPI
C Global reduction of func1/2
        ftmp(1) = func1
        ftmp(2) = func2
        call MPI_AllReduce(ftmp,ftmp2,2,MPI_double_precision,MPI_sum,
     .    MPI_Comm_World,MPIerror)
        func1 = ftmp2(1)
        func2 = ftmp2(2)
#endif

        bux1(1:nbands) = 0.0_dp
        do k = 1,nbasisCloc
          do ik = 1,numft(k)
            j = listft(ik,k)
            bux1(j) = 2.0_dp*spinfct*ft(ik,k) + bux1(j)
          enddo
          do ik = 1,numc(k)
            j = listc(ik,k)
            grad(ik,k,is) = bux1(j) + grad(ik,k,is)
          enddo
          do ik = 1,numft(k)
            j = listft(ik,k)
            bux1(j) = 0.0_dp
          enddo
        enddo

        ener = ener + spinfct*(func1 - 0.5_dp*func2) 
     .              + 0.5_dp*eta(is)*qs(is)

C End loop over spins
      enddo

#ifdef MPI
      call globaliseC(nbasisCloc,ncmax,numc,grad,nspin,Node)
#endif

C Deallocate local arrays
#ifdef MPI
      call de_alloc( buxg, 'buxg', 'gradient' )
      call de_alloc( bux2save, 'bux2save', 'gradient' )
      call de_alloc( bux1save, 'bux1save', 'gradient' )
      call de_alloc( fstG, 'fstG', 'gradient' )
      call de_alloc( ftG, 'ftG', 'gradient' )
#endif

      call de_alloc( fst, 'fst', 'gradient' )
      call de_alloc( ft, 'ft', 'gradient' )
      call de_alloc( bux1, 'bux1', 'gradient' )
      call de_alloc( bux, 'bux', 'gradient' )
      call de_alloc( aux2, 'aux2', 'gradient' )
      call de_alloc( aux1, 'aux1', 'gradient' )

      call timer('gradient',2)

      return
      end
