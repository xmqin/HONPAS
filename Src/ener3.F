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
      module m_ener3
      public :: ener3
      CONTAINS
      subroutine ener3(c,grad,lam,eta,qs,h,s,no_u,no_l,nbands,
     .                 ncmax,nctmax,nfmax,nhmax,nhijmax,numc,listc,
     .                 numct,listct,cttoc,numf,listf,numh,listh,
     .                 listhptr,numhij,listhij,ener,no_cl,nspin,
     .                 Node)

C ************************************************************************
C Finds the energy at three points of the line passing thru C in the
C direction of GRAD. LAM is the distance (in units of GRAD) between 
C points.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C ****************************** INPUT ***********************************
C real*8 c(ncmax,no_u)       : Current point (wave function coeff.
C                                  in sparse form)
C real*8 grad(ncmax,no_u)    : Direction of search (sparse)
C real*8 lam                   : Length of step
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
C integer numh(no_u)         : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listh(nhmax)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer listhptr(no_u)     : Control vector of H matrix 
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

      use precision, only : dp
      use on_main,   only : ncG2L,ncT2P,nbG2L,nbL2G,nbandsloc
      use globalise
      use m_mpi_utils, only : globalize_sum

      implicit none

      integer, intent(in) ::
     .  no_u,no_l,nbands,ncmax,nctmax,nfmax,nhmax,nhijmax,
     .  no_cl,nspin,Node

      integer, intent(in) ::
     .  cttoc(nctmax,nbandsloc),listc(ncmax,no_cl),
     .  listct(nctmax,nbandsloc),listf(nfmax,nbandsloc),
     .  listh(nhmax),listhptr(no_l),listhij(nhijmax,nbandsloc),
     .  numc(no_cl),numct(nbandsloc),numf(nbandsloc),
     .  numh(no_l),numhij(nbandsloc)

      real(dp), intent(in) ::
     .  c(ncmax,no_cl,nspin),eta(nspin),qs(nspin),
     .  grad(ncmax,no_cl,nspin),h(nhmax,nspin),lam,s(nhmax)

      real(dp), intent(out) :: ener(3)

C Internal variables ......................................................

      integer
     .  i,il,in,is,j,jn,k,kk,kn,lc,indk,mu

      real(dp), dimension(:), allocatable, save ::
     .  aux1, aux2, aux3, aux4, aux5, aux6

      real(dp), dimension(:), allocatable, save ::
     .  bux1, bux2, bux3, bux4, bux5, bux6

      real(dp) ::
     .  a1,a2,a3,b1,b2,b3,c1,c2,c3,func1(3),func2(3),
     .  lam1,lam2,lam3,pp1,pp2,pp3,hs,ss, spinfct

#ifdef MPI
      real(dp) ::
     .  ftmp(3,2),ftmp2(3,2)
      real(dp), dimension(:,:), allocatable, save ::
     .  bux246, bux1s, bux2s, bux3s, bux4s, bux5s, bux6s
#endif
C
      call timer('ener3',1)

C Allocate workspace arrays
      allocate(aux1(no_u))
      call memory('A','D',size(aux1),'ener3')
      allocate(aux2(no_u))
      call memory('A','D',size(aux2),'ener3')
      allocate(aux3(no_u))
      call memory('A','D',size(aux3),'ener3')
      allocate(aux4(no_u))
      call memory('A','D',size(aux4),'ener3')
      allocate(aux5(no_u))
      call memory('A','D',size(aux5),'ener3')
      allocate(aux6(no_u))
      call memory('A','D',size(aux6),'ener3')
      allocate(bux1(nbands))
      call memory('A','D',size(bux1),'ener3')
      allocate(bux2(nbands))
      call memory('A','D',size(bux2),'ener3')
      allocate(bux3(nbands))
      call memory('A','D',size(bux3),'ener3')
      allocate(bux4(nbands))
      call memory('A','D',size(bux4),'ener3')
      allocate(bux5(nbands))
      call memory('A','D',size(bux5),'ener3')
      allocate(bux6(nbands))
      call memory('A','D',size(bux6),'ener3')
#ifdef MPI
      allocate(bux1s(nhijmax,nbandsloc))
      call memory('A','D',size(bux1s),'ener3')
      allocate(bux2s(nhijmax,nbandsloc))
      call memory('A','D',size(bux2s),'ener3')
      allocate(bux3s(nhijmax,nbandsloc))
      call memory('A','D',size(bux3s),'ener3')
      allocate(bux4s(nhijmax,nbandsloc))
      call memory('A','D',size(bux4s),'ener3')
      allocate(bux5s(nhijmax,nbandsloc))
      call memory('A','D',size(bux5s),'ener3')
      allocate(bux6s(nhijmax,nbandsloc))
      call memory('A','D',size(bux6s),'ener3')
      allocate(bux246(nbands,3))
      call memory('A','D',size(bux246),'ener3')
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

      do i = 1,no_u
        aux1(i) = 0.0d0
        aux2(i) = 0.0d0
        aux3(i) = 0.0d0
        aux4(i) = 0.0d0
        aux5(i) = 0.0d0
        aux6(i) = 0.0d0
      enddo

      do i = 1,nbands
        bux1(i) = 0.0d0
        bux2(i) = 0.0d0
        bux3(i) = 0.0d0
        bux4(i) = 0.0d0
        bux5(i) = 0.0d0
        bux6(i) = 0.0d0
#ifdef MPI
        bux246(i,1) = 0.0d0
        bux246(i,2) = 0.0d0
        bux246(i,3) = 0.0d0
#endif
      enddo

#ifdef MPI
      do i = 1,nbandsloc
        do j = 1,nhijmax
          bux1s(j,i) = 0.0d0
          bux2s(j,i) = 0.0d0
          bux3s(j,i) = 0.0d0
          bux4s(j,i) = 0.0d0
          bux5s(j,i) = 0.0d0
          bux6s(j,i) = 0.0d0
        enddo
      enddo
#endif

C Define points to compute energy ..........................................
      lam1 = lam
      lam2 = lam*2.0d0
      lam3 = lam*3.0d0

C Loop over spin states
      do is = 1,nspin

#ifdef MPI
C Initialise communication arrays
        call globalinitB(3)
#endif

C Initialise variables
        do i = 1,3
          func1(i) = 0.0d0
          func2(i) = 0.0d0
        enddo

C Calculate Functional
C F=CtH
C Fs=CtS

        do il = 1,nbandsloc
          i = nbL2G(il)

          do in = 1,numct(il)
            k = listct(in,il)
            kk = ncT2P(k)
            if (kk.gt.0) then
              pp1 = c(cttoc(in,il),k,is) + lam1*grad(cttoc(in,il),k,is)
              pp2 = c(cttoc(in,il),k,is) + lam2*grad(cttoc(in,il),k,is)
              pp3 = c(cttoc(in,il),k,is) + lam3*grad(cttoc(in,il),k,is)
  
              indk = listhptr(kk)

              do kn = 1,numh(kk)
                mu = listh(indk+kn)
                hs = h(indk+kn,is) - eta(is)*s(indk+kn)
                ss = s(indk+kn)
                aux1(mu) = aux1(mu) + pp1*hs
                aux2(mu) = aux2(mu) + pp2*hs
                aux3(mu) = aux3(mu) + pp3*hs
                aux4(mu) = aux4(mu) + pp1*ss
                aux5(mu) = aux5(mu) + pp2*ss
                aux6(mu) = aux6(mu) + pp3*ss
              enddo
            endif
          enddo

          do in = 1,numf(il)
            k = listf(in,il)
            a1 = aux1(k)
            a2 = aux2(k)
            a3 = aux3(k)
            b1 = aux4(k)
            b2 = aux5(k)
            b3 = aux6(k)
            aux1(k) = 0.0d0
            aux2(k) = 0.0d0
            aux3(k) = 0.0d0
            aux4(k) = 0.0d0
            aux5(k) = 0.0d0
            aux6(k) = 0.0d0

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
            kk = ncG2L(k)
            do kn = 1,numc(kk)
              lc = listc(kn,kk)
              c1 = c(kn,kk,is) + lam1 * grad(kn,kk,is)
              c2 = c(kn,kk,is) + lam2 * grad(kn,kk,is)
              c3 = c(kn,kk,is) + lam3 * grad(kn,kk,is)
              bux1(lc) = bux1(lc) + a1 * c1
              bux2(lc) = bux2(lc) + b1 * c1
              bux3(lc) = bux3(lc) + a2 * c2
              bux4(lc) = bux4(lc) + b2 * c2
              bux5(lc) = bux5(lc) + a3 * c3
              bux6(lc) = bux6(lc) + b3 * c3
            enddo
          enddo

#ifdef MPI
C Load data into globalisation arrays
          call globalloadB3(il,nbands,bux2,bux4,bux6)
#endif

C First energy contribution
          func1(1) = func1(1) + bux1(i)
          func1(2) = func1(2) + bux3(i)
          func1(3) = func1(3) + bux5(i)

#ifdef MPI
C Reinitialise buxs and save for later in sparse form
          do jn = 1,numhij(il)
            j = listhij(jn,il) 
            bux1s(jn,il) = bux1(j)
            bux2s(jn,il) = bux2(j)
            bux3s(jn,il) = bux3(j)
            bux4s(jn,il) = bux4(j)
            bux5s(jn,il) = bux5(j)
            bux6s(jn,il) = bux6(j)
            bux1(j) = 0.0d0
            bux2(j) = 0.0d0
            bux3(j) = 0.0d0
            bux4(j) = 0.0d0
            bux5(j) = 0.0d0
            bux6(j) = 0.0d0
          enddo

        enddo

C Global sum of relevant bux values
        call globalcommB(Node)

C Restore local buxs terms
        do il = 1,nbandsloc
          do jn = 1,numhij(il)
            j = listhij(jn,il)
            bux2(j) = bux2s(jn,il)
            bux4(j) = bux4s(jn,il)
            bux6(j) = bux6s(jn,il)
          enddo

C Reload data after globalisation
          call globalreloadB3(il,nbands,bux2,bux4,bux6,bux246)
#endif

C Second energy contribution & reinitialise bux2/4/6
          do jn = 1,numhij(il)
            j = listhij(jn,il)
#ifdef MPI
            func2(1) = func2(1) + bux1s(jn,il)*bux246(j,1)
            func2(2) = func2(2) + bux3s(jn,il)*bux246(j,2)
            func2(3) = func2(3) + bux5s(jn,il)*bux246(j,3)
            bux2(j) = 0.0d0
            bux4(j) = 0.0d0
            bux6(j) = 0.0d0
#else
            func2(1) = func2(1) + bux1(j)*bux2(j)
            func2(2) = func2(2) + bux3(j)*bux4(j)
            func2(3) = func2(3) + bux5(j)*bux6(j)
            bux1(j) = 0.0d0
            bux2(j) = 0.0d0
            bux3(j) = 0.0d0
            bux4(j) = 0.0d0
            bux5(j) = 0.0d0
            bux6(j) = 0.0d0
#endif
          enddo

#ifdef MPI 
C Reinitialise bux246
          call globalrezeroB3(il,nbands,bux246)
#endif

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
      call memory('D','D',size(bux246),'ener3')
      deallocate(bux246)
      call memory('D','D',size(bux6s),'ener3')
      deallocate(bux6s)
      call memory('D','D',size(bux5s),'ener3')
      deallocate(bux5s)
      call memory('D','D',size(bux4s),'ener3')
      deallocate(bux4s)
      call memory('D','D',size(bux3s),'ener3')
      deallocate(bux3s)
      call memory('D','D',size(bux2s),'ener3')
      deallocate(bux2s)
      call memory('D','D',size(bux1s),'ener3')
      deallocate(bux1s)
#endif
      call memory('D','D',size(bux6),'ener3')
      deallocate(bux6)
      call memory('D','D',size(bux5),'ener3')
      deallocate(bux5)
      call memory('D','D',size(bux4),'ener3')
      deallocate(bux4)
      call memory('D','D',size(bux3),'ener3')
      deallocate(bux3)
      call memory('D','D',size(bux2),'ener3')
      deallocate(bux2)
      call memory('D','D',size(bux1),'ener3')
      deallocate(bux1)
      call memory('D','D',size(aux6),'ener3')
      deallocate(aux6)
      call memory('D','D',size(aux5),'ener3')
      deallocate(aux5)
      call memory('D','D',size(aux4),'ener3')
      deallocate(aux4)
      call memory('D','D',size(aux3),'ener3')
      deallocate(aux3)
      call memory('D','D',size(aux2),'ener3')
      deallocate(aux2)
      call memory('D','D',size(aux1),'ener3')
      deallocate(aux1)

      call timer('ener3',2)

      end subroutine ener3
      end module m_ener3

