! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_gradient
      public :: gradient
      CONTAINS
      subroutine gradient(c,eta,qs,h,s,no_u,no_l,nbands,ncmax,
     .                    nctmax,nfmax,nftmax,nhmax,nhijmax,
     .                    numc,listc,numct,listct,cttoc,numf,listf,
     .                    numft,listft,fttof,numh,listh,
     .                    listhptr,numhij,listhij,f,fs,grad,ener,
     .                    no_cl,nspin,Node)

C ************************************************************************
C Finds the energy and gradient at point C.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C Last modified: J.M.Soler. 30/04/97
C ****************************** INPUT ***********************************
C real*8 c(ncmax,no_u)       : Current point (wave function coeff.
C                                  in sparse form)
C real*8 eta(nspin)            : Fermi level parameter of Kim et al.
C real*8 qs(nspin)             : Total number of electrons
C real*8 h(nhmax)              : Hamiltonian matrix (sparse)
C real*8 s(nhmax)              : Overlap matrix (sparse)
C integer no_u               : Global number of basis orbitals
C integer no_l            : Local number of basis orbitals
C integer nbands               : Number of LWF's
C integer ncmax                : Max num of <>0 elements of each row of C
C integer nctmax               : Max num of <>0 elements of each col of C
C integer nfmax                : Max num of <>0 elements of each row of 
C                                   F = Ct x H
C integer nftmax               : Max num of <>0 elements of each col of F
C integer nhmax                : Max num of <>0 elements of each row of H
C integer nhijmax              : Max num of <>0 elements of each row of 
C                                   Hij=Ct x H x C
C integer numc(no_u)         : Control vector of C matrix
C                                (number of <>0  elements of each row of C)
C integer listc(ncmax,no_u)  : Control vector of C matrix 
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
C integer numft(no_u)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,no_u) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,no_u)  : Map from Ft to F indexing
C integer numh(no_u)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer listhptr(no_u)     : Control vector of H matrix 
C                               (pointer to start of rows in listh/h/s)
C integer numhij(nbands)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C real*8 f(nfmax,nbands)       : Auxiliary space
C real*8 fs(nfmax,nbands)      : Auxiliary space
C ***************************** OUTPUT ***********************************
C real*8 ener                  : Energy at point C
C real*8 grad(ncmax,no_u)    : Gradient of functional (sparse)
C ************************************************************************

      use precision
      use on_main,     only : ncG2L, ncT2P, nbL2G, nbG2L, nbandsloc
      use alloc,       only : re_alloc, de_alloc
#ifdef MPI
      use globalise,   only : globalinitB, globalisef, globalisec
      use globalise,   only : globalloadb2, globalcommb
      use globalise,   only : globalreloadb2, globalrezerob2
      use globalise,   only : maxft2, numft2, listft2
#endif
      use m_mpi_utils, only : globalize_sum

      implicit none

      integer, intent(in) ::
     .  no_u,no_l,nbands,ncmax,nctmax,nfmax,nftmax,nhmax,
     .  nhijmax,no_cl,nspin,Node

      integer, intent(in) ::
     .  cttoc(nctmax,nbandsloc),fttof(nftmax,no_cl),
     .  listc(ncmax,no_cl),listct(nctmax,nbandsloc),
     .  listf(nfmax,nbandsloc),listft(nftmax,no_cl),
     .  listh(nhmax),listhptr(no_l),listhij(nhijmax,nbandsloc),
     .  numc(no_cl),numct(nbandsloc),numf(nbandsloc),
     .  numft(no_cl),numh(no_l),numhij(nbandsloc)

      real(dp), intent(in) ::
     .  c(ncmax,no_cl,nspin),eta(nspin),qs(nspin),
     .  h(nhmax,nspin),s(nhmax)

      real(dp), intent(inout) ::
     .     f(nfmax,nbandsloc),fs(nfmax,nbandsloc)

      real(dp), intent(out) :: grad(ncmax,no_cl,nspin), ener

C Internal variables ......................................................

      integer
     .  i,ik,il,in,indk,is,j,jk,jl,jn,k,kk,kn,mu,muk

#ifdef MPI
      real(dp) ::
     .  ftmp(2), ftmp2(2)
      real(dp), pointer, save ::
     .  buxg(:,:),bux1save(:,:),bux2save(:,:), ftG(:,:), fstG(:,:)
#endif
      real(dp), pointer, save ::
     .  aux1(:), aux2(:), bux1(:), bux2(:), ft(:,:), fst(:,:)

      real(dp) ::
     .  a0,b0,p0,func1,func2,spinfct

C..................

      call timer('gradient',1)
      
C Allocate local arrays
      call re_alloc( aux1, 1, no_u, 'aux1', 'gradient' )
      call re_alloc( aux2, 1, no_u, 'aux2', 'gradient' )
      call re_alloc( bux1, 1, nbands, 'bux1', 'gradient' )
      call re_alloc( bux2, 1, nbands, 'bux2', 'gradient' )
      call re_alloc( ft, 1, nftmax, 1, no_cl, 'ft', 'gradient' )
      call re_alloc( fst, 1, nftmax, 1, no_cl, 'fst', 'gradient' )
#ifdef MPI
      call re_alloc( ftG, 1, MaxFt2, 1, no_cl, 'ftG', 'gradient' )
      call re_alloc( fstG, 1, MaxFt2, 1, no_cl, 'fstG', 'gradient')
      call re_alloc( bux1save, 1, nhijmax, 1, nbandsloc, 'bux1save',
     &               'gradient' )
      call re_alloc( bux2save, 1, nhijmax, 1, nbandsloc, 'bux2save',
     &               'gradient' )
      call re_alloc( buxg, 1, nbands, 1, 2, 'buxg', 'gradient')
#endif

C Initialize output and auxiliary varialbles ...............................

      if (nspin.eq.1) then
        spinfct = 2.0d0
      else
        spinfct = 1.0d0
      endif

      ener = 0.0d0

      do i = 1,no_u
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux1(i) = 0.0d0
        bux2(i) = 0.0d0
#ifdef MPI
        buxg(i,1) = 0.0d0
        buxg(i,2) = 0.0d0
#endif
      enddo

#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          bux1save(j,i) = 0.0d0
          bux2save(j,i) = 0.0d0
        enddo
      enddo
#endif

C Initialise gradient vector
      grad = 0.0d0

C Calculate Functional .....................................................

C Loop over spin
      do is = 1,nspin

#ifdef MPI
C Initialise communication arrays
        call globalinitB(2)
#endif

C Initialise variables
        func1 = 0.0d0
        func2 = 0.0d0

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
            aux1(k) = 0.0d0
            aux2(k) = 0.0d0
          enddo
        enddo

C-JMS Find transpose of F and Fs
        do mu = 1,no_cl
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
        call globaliseF(no_cl,nbands,nftmax,numft,listft,
     .                  ft,ftG,Node)
        call globaliseF(no_cl,nbands,nftmax,numft,listft,
     .                  fst,fstG,Node)
#endif

C Sij=CtSC
C multiply FxC and FsxC row by row
        do il = 1,nbandsloc
          do in = 1,numf(il)
            k = listf(in,il)
            kk = ncG2L(k)
            a0 = f(in,il)
            b0 = fs(in,il)
            do kn = 1,numc(kk)
              bux1(listc(kn,kk)) = bux1(listc(kn,kk)) + a0 * c(kn,kk,is)
              bux2(listc(kn,kk)) = bux2(listc(kn,kk)) + b0 * c(kn,kk,is)
            enddo
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB2(il,nbands,bux1,bux2)
#endif

#ifdef MPI
C Reinitialise bux1/2 and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1save(jn,il) = bux1(j)
            bux2save(jn,il) = bux2(j)
            bux1(j) = 0.0d0
            bux2(j) = 0.0d0
          enddo

        enddo

C Globalise Hij/Sij
        call globalcommB(Node)

C Restore local bux1/2 terms
        do il = 1,nbandsloc
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1(j) = bux1save(jn,il)
            bux2(j) = bux2save(jn,il)
          enddo

C Reload data after globalisation
          call globalreloadB2(il,nbands,bux1,bux2,buxg)
#endif

C Calculate energy terms & reinitialise bux1/2
          i = nbL2G(il)
          func1 = func1 + bux1(i)
  
          do jn = 1,numhij(il)
            j = listhij(jn,il)
#ifdef MPI
            func2 = func2 + bux1(j)*buxg(j,2)
#else
            func2 = func2 + bux1(j)*bux2(j)
#endif
          enddo

C Multiply Hij x Fs and Sij x F row by row
C (only products of neccesary elements)
          do ik = 1,numct(il)
            mu = listct(ik,il)
            kk = ncT2P(mu)
            if (kk.gt.0) then
              a0 = 0.0d0
#ifdef MPI
              do muk = 1,numft2(mu)
                j = listft2(muk,mu)
                a0 = a0 + buxg(j,1)*fstG(muk,mu) 
     .                  + buxg(j,2)*ftG(muk,mu)
#else
              do muk = 1,numft(mu)
                j = listft(muk,mu)
                a0 = a0 + bux1(j)*fst(muk,mu) 
     .                  + bux2(j)*ft(muk,mu)
#endif
              enddo
              grad(cttoc(ik,il),mu,is) = - spinfct*a0
            endif
          enddo

#ifdef MPI
C Reinitialise buxg
          call globalrezeroB2(il,nbands,buxg)
#endif
C Reinitialise bux1 / bux2
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux1(j) = 0.0d0
            bux2(j) = 0.0d0
          enddo

        enddo

#ifdef MPI
C Global reduction of func1/2
        ftmp(1) = func1
        ftmp(2) = func2
        call Globalize_sum(ftmp(1:2),ftmp2(1:2))
        func1 = ftmp2(1)
        func2 = ftmp2(2)
#endif

        bux1(1:nbands) = 0.0d0
        do k = 1,no_cl
          do ik = 1,numft(k)
            j = listft(ik,k)
            bux1(j) = 2.0d0*spinfct*ft(ik,k) + bux1(j)
          enddo
          do ik = 1,numc(k)
            j = listc(ik,k)
            grad(ik,k,is) = bux1(j) + grad(ik,k,is)
          enddo
          do ik = 1,numft(k)
            j = listft(ik,k)
            bux1(j) = 0.0d0
          enddo
        enddo

        ener = ener + spinfct*(func1 - 0.5d0*func2) 
     .              + 0.5d0*eta(is)*qs(is)

C End loop over spins
      enddo

#ifdef MPI
      call globaliseC(no_cl,ncmax,numc,grad,nspin,Node)
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
      call de_alloc( bux2, 'bux2', 'gradient' )
      call de_alloc( bux1, 'bux1', 'gradient' )
      call de_alloc( aux2, 'aux2', 'gradient' )
      call de_alloc( aux1, 'aux1', 'gradient' )

      call timer('gradient',2)

      end subroutine gradient
      end module m_gradient
