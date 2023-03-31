! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine pdosg( nspin, nuo, no, maxspn, maxnh, 
     .                  maxo, numh, listhptr, listh, H, S,
     .                  E1, E2, nhist, sigma, indxuo, eo, 
     .                  haux, saux, psi, dtot, dpr, nuotot )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C Gamma point version adapted from PDOSK by Julian Gale. Feb' 03
C ****  INPUT  *********************************************************
C integer nspin             : Number of spin components (1 or 2)
C integer nuo               : Number of atomic orbitals in the unit cell
C integer no                : Number of atomic orbitals in the supercell
C integer maxspn            : Second dimension of eo and qo 
C                             (maximum number of differents spin polarizations)
C integer maxnh             : Maximum number of orbitals interacting
C                             with any orbital
C integer maxo              : First dimension of eo
C integer numh(nuo)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C integer listhptr(nuo)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C integer listh(maxnh)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C real*8  H(maxnh,nspin)    : Hamiltonian in sparse format
C real*8  S(maxnh)          : Overlap in sparse format
C real*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C integer nhist             : Number of the subdivisions of the histogram
C real*8  sigma             : Width of the gaussian to expand the eigenvectors
C integer indxuo(no)        : Index of equivalent orbital in unit cell
C real*8  eo(maxo,maxspn)   : Eigenvalues
C integer nuotot            : Total number of orbitals per unit cell
C ****  AUXILIARY  *****************************************************
C real*8  haux(nuo,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8  saux(nuo,nuo)     : Auxiliary space for the overlap matrix
C real*8  psi(nuo,nuo)      : Auxiliary space for the eigenvectors
C ****  OUTPUT  ********************************************************
C real*8  dtot(nhist,2)   : Total density of states
C real*8  dpr(nhist,nuo,2): Proyected density of states
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
     .  nspin, nuo, no, maxspn, maxnh, 
     .  maxo, nhist, nuotot

      integer
     .  numh(nuo), listhptr(nuo), listh(maxnh),
     .  indxuo(no)

      real(dp)
     .  H(maxnh,nspin), S(maxnh), E1, E2, sigma, eo(maxo,maxspn),
     .  haux(nuotot,nuo), saux(nuotot,nuo), psi(nuotot,nuo),
     .  dtot(nhist,2), dpr(nhist,nuotot,2) 

C Internal variables ---------------------------------------------------
      integer
     .  ispin, iuo, juo, j, jo, ihist, iband, ind, ierror

      real(dp)
     .  delta, ener, diff, pipj1, gauss, norm

#ifdef MPI
      integer ::
     .  BNode, Bnuo, ibandg, maxnuo, MPIerror
      real(dp), dimension(:,:), pointer ::  Sloc
      real(dp), dimension(:,:,:), pointer :: tmp
#endif

      external rdiag

C Initialize some variables
      delta = (E2 - E1)/nhist

C Solve eigenvalue problem for each k-point
      do ispin = 1, nspin

C Initialize auxiliary variables 
        do iuo = 1, nuo
          do juo = 1, nuotot
            saux(juo,iuo) = 0.0d0
            haux(juo,iuo) = 0.0d0
          enddo
        enddo

        do iuo = 1, nuo
          do j = 1, numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo= indxuo(jo)

C Build the full Hamiltonian and overlap matrix
            saux(juo,iuo) = saux(juo,iuo) + S(ind)
            haux(juo,iuo) = haux(juo,iuo) + H(ind,ispin)
          enddo
        enddo

C Diagonalize at the Gamma point
        call rdiag( haux, saux, nuotot, nuo, nuotot, eo(1,ispin), 
     .              psi, nuotot, 1, ierror, BlockSize)

C Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
          do iuo = 1, nuo
            do juo = 1, nuotot
              saux(juo,iuo) = 0.0d0
              haux(juo,iuo) = 0.0d0
            enddo
          enddo
          do iuo = 1, nuo
            do j = 1, numh(iuo)
              ind = listhptr(iuo) + j
              jo  = listh(ind)
              juo = indxuo(jo)
              saux(juo,iuo) = saux(juo,iuo) + S(ind)
              haux(juo,iuo) = haux(juo,iuo) + H(ind,ispin)
            enddo
          enddo
          call rdiag( haux, saux, nuotot, nuo, nuotot, eo(1,ispin),
     .                psi, nuotot, 1, ierror, BlockSize)
        endif

C Rebuild the full overlap matrix
        do iuo = 1, nuo
          do juo = 1, nuotot
            saux(juo,iuo) = 0.0d0
          enddo
          do j = 1, numh(iuo)
            ind = listhptr(iuo) + j
            jo  = listh(ind)
            juo = indxuo(jo)
            saux(juo,iuo) = saux(juo,iuo) + S(ind)
          enddo
        enddo

#ifdef MPI
C Find maximum number of orbitals per node
        call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)

C Allocate workspace array for broadcast overlap matrix
        nullify( Sloc )
        call re_alloc( Sloc, 1, nuotot, 1, maxnuo, 
     &                name='Sloc', routine='pdosg' )

C Loop over nodes broadcasting overlap matrix
        do BNode = 0,Nodes-1

C Find out how many orbitals there are on the broadcast node
          call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

C Transfer data
          if (Node.eq.BNode) then
            Sloc(1:nuotot,1:Bnuo) = Saux(1:nuotot,1:Bnuo)
          endif
          call MPI_Bcast(Sloc(1,1),nuotot*Bnuo,
     .      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)

C Loop over all the energy range
          do ihist = 1, nhist
            ener = E1 + (ihist - 1) * delta
            do 170 iband = 1, nuo
              call LocalToGlobalOrb(iband,Node,Nodes,ibandg)
              diff = (ener - eo(ibandG,ispin))**2 / (sigma ** 2)
              if (diff .gt. 15.0d0) then
                cycle
              else
                gauss = exp(-diff)
                if (Node.eq.BNode) then
C Only add once to dtot - not everytime loop over processors is executed
                  dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
                endif
                do jo = 1, Bnuo
                  call LocalToGlobalOrb(jo,BNode,Nodes,juo)
                  do iuo = 1, nuotot
C Solo para los Juo que satisfagan el criterio del record...
                    pipj1 = psi(iuo,iband) * psi(juo,iband)
                    dpr(ihist,juo,ispin) = dpr(ihist,juo,ispin) + 
     .                                     pipj1*gauss*Sloc(iuo,jo)
                  enddo
                enddo
              endif
 170        enddo

          enddo

C End loop over broadcast nodes
        enddo

C Free workspace array for overlap
        call de_alloc( Sloc, name='Sloc' )

#else
C Loop over all the energy range
        do ihist = 1, nhist
          ener = E1 + (ihist - 1) * delta
          do 170 iband = 1, nuo
            diff = (ener - eo(iband,ispin))**2 / (sigma ** 2)
            if (diff .gt. 15.0d0) then
              cycle
            else
              gauss = exp(-diff)
              dtot(ihist,ispin) = dtot(ihist,ispin) + gauss
              do iuo = 1, nuotot
                do juo = 1, nuotot
                  pipj1 = psi(iuo,iband) * psi(juo,iband)
                  dpr(ihist,juo,ispin) = dpr(ihist,juo,ispin) + 
     .                                   pipj1*gauss*saux(iuo,juo)
                enddo
              enddo
            endif
 170      enddo

        enddo
#endif

      enddo

#ifdef MPI
C Allocate workspace array for global reduction
      nullify( tmp )
      call re_alloc( tmp, 1, nhist, 1, max(nuotot,nspin), 
     &               1, nspin, name='tmp', routine='pdosg' )

C Global reduction of dpr matrix
      tmp(1:nhist,1:nuotot,1:nspin) = 0.0d0
      call MPI_AllReduce(dpr(1,1,1),tmp(1,1,1),nhist*nuotot*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dpr(1:nhist,1:nuotot,1:nspin) = tmp(1:nhist,1:nuotot,1:nspin)

C Global reduction of dtot matrix
      tmp(1:nhist,1:nspin,1) = 0.0d0
      call MPI_AllReduce(dtot(1,1),tmp(1,1,1),nhist*nspin,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      dtot(1:nhist,1:nspin) = tmp(1:nhist,1:nspin,1)

C Free workspace array for global reduction
      call de_alloc( tmp, name='tmp' )
#endif

      norm = sigma * sqrt(pi)

      do ihist = 1, nhist
        do ispin = 1, nspin
          dtot(ihist,ispin) = dtot(ihist,ispin) / norm
          do iuo = 1, nuotot
            dpr(ihist,iuo,ispin) = dpr(ihist,iuo,ispin) /norm
          enddo
        enddo
      enddo

      end
