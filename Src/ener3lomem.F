! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine ener3lomem(c,grad,lam,eta,qs,h,s,nbasis,nbasisloc,
     .                      nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,
     .                      numc,listc,numct,listct,cttoc,numf,listf,
     .                      numh,listh,listhptr,numhij,listhij,ener,
     .                      nbasisCloc,nspin,Node)

C ************************************************************************
C Finds the energy at three points of the line passing thru C in the
C direction of GRAD. LAM is the distance (in units of GRAD) between 
C points.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasis)       : Current point (wave function coeff.
C                                  in sparse form)
C real*8 grad(ncmax,nbasis)    : Direction of search (sparse)
C real*8 lam                   : Length of step
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
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer listhptr(nbasis)     : Control vector of H matrix 
C                               (pointer to start of row in listh/h/s)
C integer numhij(nbands)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C ***************************** OUTPUT ***********************************
C real*8 ener(3)               : Energy at the three points:
C                                     C +     lam * GRAD
C                                     C + 2 * lam * GRAD
C                                     C + 3 * lam * GRAD
C ************************************************************************

      use precision
      use globalise
      use on_main,     only : ncG2L,ncT2P,nbG2L,nbL2G,nbandsloc
      use m_mpi_utils, only : globalize_sum
      use alloc,       only : re_alloc, de_alloc

      implicit none

      integer
     .  nbasis,nbasisloc,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,
     .  nbasisCloc,nspin,Node

      integer
     .  cttoc(nctmax,nbandsloc),listc(ncmax,nbasisCloc),
     .  listct(nctmax,nbandsloc),listf(nfmax,nbandsloc),
     .  listh(nhmax),listhptr(nbasisloc),listhij(nhijmax,nbandsloc),
     .  numc(nbasisCloc),numct(nbandsloc),numf(nbandsloc),
     .  numh(nbasisloc),numhij(nbandsloc)

      real(dp) ::
     .  c(ncmax,nbasisCloc,nspin),ener(3),eta(nspin),qs(nspin),
     .  grad(ncmax,nbasisCloc,nspin),h(nhmax,nspin),lam,s(nhmax)

C Internal variables ......................................................

      integer
     .  i,il,ilam,in,is,j,jn,k,kk,kn,lc,indk,mu

      real(dp), pointer, save :: aux1(:), aux2(:), bux1(:), bux2(:)

      real(dp) ::
     .  a,b,c1,func1(3),func2(3),
     .  lam123(3),pp,hs,ss, spinfct

#ifdef MPI
      real(dp) ::
     .  ftmp(3,2),ftmp2(3,2)
      real(dp), pointer, save :: bux2g(:), bux1s(:,:), bux2s(:,:)
#endif
C
      call timer('ener3',1)

C Allocate workspace arrays
      call re_alloc( aux1, 1, nbasis, 'aux1', 'ener3' )
      call re_alloc( aux2, 1, nbasis, 'aux2', 'ener3' )
      call re_alloc( bux1, 1, nbasis, 'bux1', 'ener3' )
      call re_alloc( bux2, 1, nbasis, 'bux2', 'ener3' )
#ifdef MPI
      call re_alloc( bux1s, 1, nhijmax, 1, nbandsloc, 'bux1s', 'ener3' )
      call re_alloc( bux2s, 1, nhijmax, 1, nbandsloc, 'bux2s', 'ener3' )
      call re_alloc( bux2g, 1, nbasis, 'bux2g', 'ener3' )
#endif

C..................

C Initialize output and auxiliary varialbles ...............................

      if (nspin.eq.1) then
        spinfct = 2.0d0
      else
        spinfct = 1.0d0
      endif

      do i = 1,3
        ener(i) = 0.0d0
      enddo

      do i = 1,nbasis
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux1(i) = 0.0d0
        bux2(i) = 0.0d0
#ifdef MPI
        bux2g(i) = 0.0d0
#endif
      enddo

#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          bux1s(j,i) = 0.0d0
          bux2s(j,i) = 0.0d0
        enddo
      enddo
#endif

C Define points to compute energy ..........................................
      lam123(1) = lam
      lam123(2) = lam*2.0d0
      lam123(3) = lam*3.0d0

C Loop over spin states
      do is = 1,nspin

C Loop over lambda values
        do ilam = 1,3

#ifdef MPI
C Initialise communication arrays
          call globalinitB(1)
#endif

C Initialise variables
          func1(ilam) = 0.0d0
          func2(ilam) = 0.0d0

C Calculate Functional
C F=CtH
C Fs=CtS

          do il = 1,nbandsloc
            i = nbL2G(il)

            do in = 1,numct(il)
              k = listct(in,il)
              kk = ncT2P(k)
              if (kk.gt.0) then
                pp = c(cttoc(in,il),k,is) + 
     .               lam123(ilam)*grad(cttoc(in,il),k,is)
  
                indk = listhptr(kk)

                do kn = 1,numh(kk)
                  mu = listh(indk+kn)
                  hs = h(indk+kn,is) - eta(is)*s(indk+kn)
                  ss = s(indk+kn)
                  aux1(mu) = aux1(mu) + pp*hs
                  aux2(mu) = aux2(mu) + pp*ss
                enddo
              endif
            enddo

            do in = 1,numf(il)
              k = listf(in,il)
              a = aux1(k)
              b = aux2(k)
              aux1(k) = 0.0d0
              aux2(k) = 0.0d0
  
C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
              kk = ncG2L(k)
              do kn = 1,numc(kk)
                lc = listc(kn,kk)
                c1 = c(kn,kk,is) + lam123(ilam)*grad(kn,kk,is)
                bux1(lc) = bux1(lc) + a * c1
                bux2(lc) = bux2(lc) + b * c1
              enddo
            enddo

#ifdef MPI
C Load data into globalisation arrays
            call globalloadB1(il,nbands,bux2)
#endif

C First energy contribution
            func1(ilam) = func1(ilam) + bux1(i)

#ifdef MPI
C Reinitialise buxs and save for later in sparse form
            do jn = 1,numhij(il)
              j = listhij(jn,il) 
              bux1s(jn,il) = bux1(j)
              bux2s(jn,il) = bux2(j)
              bux1(j) = 0.0d0
              bux2(j) = 0.0d0
            enddo

          enddo

C Global sum of relevant bux values
          call globalcommB(Node)

C Restore local buxs terms
          do il = 1,nbandsloc
            do jn = 1,numhij(il)
              j = listhij(jn,il)
              bux2(j) = bux2s(jn,il)
            enddo

C Reload data after globalisation
            call globalreloadB1(il,nbands,bux2,bux2g)
#endif

C Second energy contribution & reinitialise bux2/4/6
            do jn = 1,numhij(il)
              j = listhij(jn,il)
#ifdef MPI
              func2(ilam) = func2(ilam) + bux1s(jn,il)*bux2g(j)
              bux2(j) = 0.0d0
#else
              func2(ilam) = func2(ilam) + bux1(j)*bux2(j)
              bux1(j) = 0.0d0
              bux2(j) = 0.0d0
#endif
            enddo

#ifdef MPI 
C Reinitialise bux246
            call globalrezeroB1(il,nbands,bux2g)
#endif

          enddo

        enddo

#ifdef MPI
C Global reduction of func1/2
        ftmp(1:3,1) = func1(1:3)
        call Globalize_sum(ftmp(1:3,1), ftmp2(1:3,1))
        func1(1:3) = ftmp2(1:3,1)
        ftmp(1:3,2) = func2(1:3)
        call Globalize_sum(ftmp(1:3,2), ftmp2(1:3,2))
        func2(1:3) = ftmp2(1:3,2)
#endif

C This is valid for an spin-unpolarized sytem
        do i = 1,3
          ener(i) = ener(i) + spinfct*(func1(i) - 0.5d0*func2(i))
     .                      + 0.5d0*eta(is)*qs(is)
        enddo

C End loop over spin states
      enddo

C Dellocate workspace arrays
#ifdef MPI
      call de_alloc( bux2g, 'bux2g', 'ener3' )
      call de_alloc( bux2s, 'bux2s', 'ener3' )
      call de_alloc( bux1s, 'bux1s', 'ener3' )
#endif
      call de_alloc( bux2, 'bux2', 'ener3' )
      call de_alloc( bux1, 'bux1', 'ener3' )
      call de_alloc( aux2, 'aux2', 'ener3' )
      call de_alloc( aux1, 'aux1', 'ener3' )

      call timer('ener3',2)

      return
      end
