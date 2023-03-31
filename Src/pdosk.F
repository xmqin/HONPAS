! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine pdosk( nspin, nuo, no, maxspn, maxnh, 
     .                  maxo, numh, listhptr, listh, H, S,
     .                  E1, E2, nhist, sigma, 
     .                  xij, indxuo, nk, kpoint, wk, eo, 
     .                  Haux, Saux, psi, dtot, dpr, nuotot )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C ****  INPUT  *********************************************************
C INTEGER nspin             : Number of spin components (1 or 2)
C INTEGER nuo               : Number of atomic orbitals in the unit cell
C INTEGER NO                : Number of atomic orbitals in the supercell
C INTEGER maxspn            : Second dimension of eo and qo 
C                             (maximum number of differents spin polarizations)
C INTEGER maxnh             : Maximum number of orbitals interacting
C                             with any orbital
C INTEGER maxo              : First dimension of eo
C INTEGER numh(nuo)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C INTEGER listhptr(nuo)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C INTEGER listh(maxnh)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C REAL*8  H(maxnh,nspin)    : Hamiltonian in sparse format
C REAL*8  S(maxnh)          : Overlap in sparse format
C REAL*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C INTEGER nhist             : Number of the subdivisions of the histogram
C REAL*8  sigma             : Width of the gaussian to expand the eigenvectors
C REAL*8  xij(3,maxnh)      : Vectors between orbital centers (sparse)
C                             (not used if only gamma point)
C INTEGER indxuo(no)        : Index of equivalent orbital in unit cell
C INTEGER NK                : Number of k points
C REAL*8  kpoint(3,nk)      : k point vectors
C REAL*8  WK(nk)            : Weights for k points
C REAL*8  EO(maxo,maxspn,nk): Eigenvalues
C INTEGER nuotot            : Total number of orbitals per unit cell
C ****  AUXILIARY  *****************************************************
C REAL*8  Haux(2,nuo,nuo)   : Auxiliary space for the hamiltonian matrix
C REAL*8  Saux(2,nuo,nuo)   : Auxiliary space for the overlap matrix
C REAL*8  psi(2,nuo,nuo)    : Auxiliary space for the eigenvectors
C ****  OUTPUT  ********************************************************
C REAL*8  dtot(nhist,2)   : Total density of states
C REAL*8  dpr(nhist,nuo,2): Proyected density of states
C **********************************************************************

      use precision
      use parallel,     only : Node, Nodes, BlockSize
      use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
      use units,        only : pi
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif
      use sys,          only : die

      implicit none

      integer
     .  nspin, nuo, no, maxspn, maxnh, NK, 
     .  maxo, nhist, nuotot

      integer
     .  numh(nuo), listhptr(nuo), listh(maxnh),
     .  indxuo(no)

      real(dp)
     .  H(maxnh,nspin), S(maxnh), E1, E2, sigma, 
     .  xij(3,maxnh), kpoint(3,nk), eo(maxo,maxspn,nk),
     .  Haux(2,nuotot,nuo), Saux(2,nuotot,nuo), psi(2,nuotot,nuo),
     .  dtot(nhist,2), dpr(nhist,nuotot,2), wk(nk)

C Internal variables ---------------------------------------------------
      integer
     .  ik, ispin, iuo, juo, J, JO, ihist, iband, ind, ierror

      real(dp)
     .  kxij, Ckxij, Skxij, delta, ener, diff, pipj1, pipj2, 
     .  pipjS1, gauss, norm, wksum

#ifdef MPI
      integer ::
     .  BNode, Bnuo, ibandg, maxnuo, MPIerror
      real(dp), dimension(:,:,:), pointer :: Sloc
#endif

      external  cdiag

#ifdef DEBUG
      call write_debug( '    PRE pdosk' )
#endif

C Initialize some variables
      delta = (E2 - E1)/nhist

C Solve eigenvalue problem for each k-point
      do ispin = 1, nspin

        do IK = 1, NK

C Initialize auxiliary variables 
          do iuo = 1,nuo
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
              Haux(1,juo,iuo) = 0.0d0
              Haux(2,juo,iuo) = 0.0d0
            enddo
          enddo

          do iuo = 1, nuo
            do j = 1, numh(iuo)
              ind = listhptr(iuo) + j
              jo  = listh(ind)
              juo = indxuo(jo)
C Calculate the phases k*r_ij
              kxij = kpoint(1,IK) * xij(1,ind) +
     .               kpoint(2,IK) * xij(2,ind) +
     .               kpoint(3,IK) * xij(3,ind) 
              Ckxij = cos(kxij)
              Skxij = sin(kxij)
C Calculate the Hamiltonian and the overlap in k space
C H(k) = Sum(R) exp(i*k*R) * H(R)
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind) * Ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind) * Skxij
              Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin) * Ckxij
              Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin) * Skxij
            enddo
          enddo

C Diagonalize for each k point
          call cdiag( Haux, Saux, nuotot, nuo, nuotot,
     .                eo(1,ispin,IK), psi, nuotot, 1, ierror,
     $                BlockSize)

C Check error flag and take appropriate action
          if (ierror.gt.0) then
            call die('Terminating due to failed diagonalisation')
          elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
            do iuo = 1,nuo
              do juo = 1,nuotot
                Saux(1,juo,iuo) = 0.0d0
                Saux(2,juo,iuo) = 0.0d0
                Haux(1,juo,iuo) = 0.0d0
                Haux(2,juo,iuo) = 0.0d0
              enddo
            enddo
            do iuo = 1, nuo
              do j = 1, numh(iuo)
                ind = listhptr(iuo) + J
                jo  = listh(ind)
                juo = indxuo(jo)
                kxij = kpoint(1,IK) * xij(1,ind) +
     .                 kpoint(2,IK) * xij(2,ind) +
     .                 kpoint(3,IK) * xij(3,ind) 
                Ckxij = cos(kxij)
                Skxij = sin(kxij)
                Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind) * Ckxij
                Saux(2,juo,iuo) = Saux(2,juo,iuo) - S(ind) * Skxij
                Haux(1,juo,iuo) = Haux(1,juo,iuo) + H(ind,ispin) * Ckxij
                Haux(2,juo,iuo) = Haux(2,juo,iuo) - H(ind,ispin) * Skxij
              enddo
            enddo
            call cdiag( Haux, Saux, nuotot, nuo, nuotot,
     .                  eo(1,ispin,IK), psi, nuotot, 1, ierror,
     $                  BlockSize )
          endif

C Recalculate again the overlap matrix in k-space
          do iuo = 1, nuo
            do juo = 1, nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
            enddo
          enddo

          do iuo = 1,nuo
            do  j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo  = listh(ind)
              juo = indxuo(jo)
C Calculates the phases k*r_ij
              kxij = kpoint(1,IK) * xij(1,ind) +
     .               kpoint(2,IK) * xij(2,ind) +
     .               kpoint(3,IK) * xij(3,ind) 
              ckxij = cos(kxij)
              skxij = sin(kxij)
! Since we are doing element wise multiplications (and not dot-products)
! we might as well setup the transpose S(k)^T == S(-k) because this will
! mean that we can do a simpler multiplication further down
              Saux(1,juo,iuo) = Saux(1,juo,iuo) + S(ind) * ckxij
              Saux(2,juo,iuo) = Saux(2,juo,iuo) + S(ind) * skxij
            enddo
          enddo

#ifdef MPI
C Find maximum number of orbitals per node
          call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max,
     .      MPI_Comm_World,MPIerror)

C Allocate workspace array for broadcast overlap matrix
          nullify( Sloc )
          call re_alloc( Sloc, 1, 2, 1, nuotot, 
     &                   1, maxnuo, name='Sloc', routine='pdosk' )

C Loop over nodes broadcasting overlap matrix
          do BNode = 0,Nodes-1

C Find out how many orbitals there are on the broadcast node
            call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

C Transfer data
            if (Node.eq.BNode) then
              Sloc(1:2,1:nuotot,1:Bnuo) = Saux(1:2,1:nuotot,1:Bnuo)
            endif
            call MPI_Bcast(Sloc(1,1,1),2*nuotot*Bnuo,
     .        MPI_double_precision,BNode,MPI_Comm_World,MPIerror)

C Loop over all the energy range
            do ihist = 1, nhist
              ener = E1 + (ihist - 1) * delta
              do 170 iband = 1, nuo
                call LocalToGlobalOrb(iband,Node,Nodes,ibandg)
                diff = (ener - EO(ibandg,ispin,IK))**2 / (sigma ** 2)
                if (diff .gt. 15.0D0) then
                  cycle
                else
                  gauss = exp(-diff) * wk(ik)
                  if (Node.eq.BNode) then
C Only add once to dtot - not everytime loop over processors is executed
                    dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
                  endif
                  do jo = 1, Bnuo
                    call LocalToGlobalOrb(jo,BNode,Nodes,juo)
                    do iuo = 1, nuotot
!                     This is:  psi(iuo) * psi(juo)^*
                      pipj1 = psi(1,iuo,iband) * psi(1,juo,iband) +
     .                        psi(2,iuo,iband) * psi(2,juo,iband)
                      pipj2 = - psi(1,iuo,iband) * psi(2,juo,iband) +
     .                        psi(2,iuo,iband) * psi(1,juo,iband)
                      pipjS1= pipj1*Sloc(1,iuo,JO)-pipj2*Sloc(2,iuo,JO)
                      dpr(ihist,juo,ispin)= dpr(ihist,juo,ispin) +
     .                                       pipjS1*gauss
                    enddo
                  enddo
                endif
 170          enddo

            enddo

C End loop over broadcast nodes
          enddo

C Free workspace array for overlap
          call de_alloc( Sloc, name='Sloc' )

#else
C Loop over all the energy range
          do ihist = 1, nhist
            ener = E1 + (ihist - 1) * delta
            do 170 iband = 1,nuo
              diff = (ener - EO(iband,ispin,IK))**2 / (sigma ** 2)
              if (diff .gt. 15.0d0) then
                cycle
              else
                gauss = exp(-diff) * wk(ik)
                dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
                do iuo = 1, nuotot
                  do juo = 1, nuotot
                    pipj1 = psi(1,iuo,iband) * psi(1,juo,iband) +
     .                      psi(2,iuo,iband) * psi(2,juo,iband)
                    pipj2 = - psi(1,iuo,iband) * psi(2,juo,iband) +
     .                      psi(2,iuo,iband) * psi(1,juo,iband)
                    pipjS1= pipj1*Saux(1,iuo,juo)-pipj2*Saux(2,iuo,juo)
                    dpr(ihist,juo,ispin)= dpr(ihist,juo,ispin) + 
     .                                     pipjS1*gauss
                  enddo
                enddo
              endif
 170        enddo

          enddo
#endif

        enddo

      enddo

#ifdef MPI
C Allocate workspace array for global reduction
      nullify( Sloc )
      call re_alloc( Sloc, 1, nhist, 1, max(nuotot,nspin), 
     &               1, nspin, name='Sloc', routine='pdosk' )

C Global reduction of dpr matrix
      Sloc(1:nhist,1:nuotot,1:nspin) = 0.0d0
      call MPI_AllReduce(dpr(1,1,1),Sloc(1,1,1),nhist*nuotot*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dpr(1:nhist,1:nuotot,1:nspin) = Sloc(1:nhist,1:nuotot,1:nspin)

C Global reduction of dtot matrix
      Sloc(1:nhist,1:nspin,1) = 0.0d0
      call MPI_AllReduce(dtot(1,1),Sloc(1,1,1),nhist*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dtot(1:nhist,1:nspin) = Sloc(1:nhist,1:nspin,1)

C Free workspace array for global reduction
      call de_alloc( Sloc, name='Sloc' )
#endif

      wksum = 0.0d0
      do IK = 1,NK
        wksum = wksum + WK(IK)
      enddo

      norm = sigma * sqrt(PI) * wksum

      do ihist = 1, nhist
        do ispin = 1, nspin
          dtot(ihist,ispin) = dtot(ihist,ispin) / norm
          do iuo = 1, nuotot
            dpr(ihist,iuo,ispin) = dpr(ihist,iuo,ispin) /norm
          enddo
        enddo
      enddo

#ifdef DEBUG
      call write_debug( '    POS pdosk' )
#endif

      return
      end
