! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine pdos2g( nuo, no, maxuo, maxnh,
     .                  maxo, numh, listhptr, listh, H, S,
     .                  E1, E2, nhist, sigma, indxuo, eo,
     .                  haux, saux, psi, dtot, dpr, nuotot )

C **********************************************************************
C Find the density of states projected onto the atomic orbitals
C     D_mu(E) = Sum(n,k,nu) C(mu,n,k) C(nu,n,k) S(mu,nu,k) Delta(E-E(n,k))
C where n run over all the bands between two given energies
C Written by J. Junquera and E. Artacho. Nov' 99
C Gamma point version adapted from PDOSK by Julian Gale. Feb' 03
C Noncollinear spin version by J. Ferrer, October 2007
C ****  INPUT  *********************************************************
C integer nuo               : Number of atomic orbitals in the unit cell
C integer no                : Number of atomic orbitals in the supercell
C                             (maximum number of differents spin polarizations)
C integer maxuo             : Maximum number of atomic orbitals in the unit cell
C integer maxnh             : Maximum number of orbitals interacting
C                             with any orbital
C integer maxo              : First dimension of eo
C integer numh(nuo)         : Number of nonzero elements of each row
C                             of hamiltonian matrix
C integer listhptr(nuo)     : Pointer to each row (-1) of the
C                             hamiltonian matrix
C integer listh(maxnh)      : Nonzero hamiltonian-matrix element
C                             column indexes for each matrix row
C real*8  H(maxnh,4)    : Hamiltonian in sparse format
C real*8  S(maxnh)          : Overlap in sparse format
C real*8  E1, E2            : Energy range for density-matrix states
C                             (to find local density of states)
C                             Not used if e1 > e2
C integer nhist             : Number of the subdivisions of the histogram
C real*8  sigma             : Width of the gaussian to expand the eigenvectors
C integer indxuo(no)        : Index of equivalent orbital in unit cell
C real*8  eo(maxo,2)   : Eigenvalues
C integer nuotot            : Total number of orbitals per unit cell
C ****  AUXILIARY  *****************************************************
C real*8  haux(nuo,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8  saux(nuo,nuo)     : Auxiliary space for the overlap matrix
C complex*16 psi(2,nuo,nuo) : Auxiliary space for the eigenvectors
C ****  OUTPUT  ********************************************************
C real*8  dtot(nhist,4)   : Total density of states
C real*8  dpr(nhist,nuo,4): Proyected density of states
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
     .  nuo, no, maxuo, maxnh,
     .  maxo, nhist, nuotot

      integer
     .  numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)

      real(dp)
     .  H(maxnh,4), S(maxnh), E1, E2, sigma, eo(maxo*2),
     .  dtot(nhist,4), dpr(nhist,nuotot,4)
      complex(dp), target :: psi(2,nuotot,2*nuo)
      
      complex(dp) Haux(2,nuotot,2,nuo), Saux(2,nuotot,2,nuo)
      complex(dp), pointer :: caux(:,:)
      external  cdiag

C Internal variables ---------------------------------------------------
      integer
     .  ispin, iuo, juo, io, j, jo, ihist, iband, ind, ierror

      real(dp)
     .  delta, ener, diff, rpipj, ipipj, gauss, norm

      real(dp), pointer :: Spr(:,:) => null()

#ifdef MPI
      integer :: BNode, Bnuo, ibandg, maxnuo, MPIerror
      real(dp), pointer :: Sloc(:,:)
      real(dp)          :: tmp(nhist,nuotot,4)
#endif

C Initialize some variables
      delta = (E2 - E1)/nhist

      call re_alloc(Spr, 1, nuotot, 1, nuo, name='Spr',
     &     routine='pdos2g')

C Initialize auxiliary variables
      do io = 1,nuo
        Saux(:,:,:,io) = 0._dp
        Haux(:,:,:,io) = 0._dp
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          jo = indxuo(jo)
          Saux(1,jo,1,io) = Saux(1,jo,1,io) + dcmplx( S(ind), 0.0_dp)
          Saux(2,jo,2,io) = Saux(2,jo,2,io) + dcmplx( S(ind), 0.0_dp)
          Haux(1,jo,1,io) = Haux(1,jo,1,io) + dcmplx(H(ind,1),   0.0_dp)
          Haux(2,jo,2,io) = Haux(2,jo,2,io) + dcmplx(H(ind,2),   0.0_dp)
          Haux(2,jo,1,io) = Haux(2,jo,1,io) + dcmplx(H(ind,3), H(ind,4))
          Haux(1,jo,2,io) = Haux(1,jo,2,io) + dcmplx(H(ind,3),-H(ind,4))
        enddo
      enddo

C Diagonalize at the Gamma point
      call cdiag( Haux, Saux, 2*nuotot, 2*nuo, 2*nuotot,
     .            eo, psi, 2*nuotot, 1, ierror, 2*BlockSize )
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror .lt. 0) then
C Repeat diagonalisation with increased memory to handle clustering
        do io = 1,nuo
          Saux(:,:,:,io) = 0._dp
          Haux(:,:,:,io) = 0._dp
          do j = 1,numh(io)
           ind = listhptr(io) + j
           jo = listh(ind)
           jo = indxuo(jo)
           Saux(1,jo,1,io) = Saux(1,jo,1,io) +dcmplx( S(ind), 0.0_dp)
           Saux(2,jo,2,io) = Saux(2,jo,2,io) +dcmplx( S(ind), 0.0_dp)
           Haux(1,jo,1,io) = Haux(1,jo,1,io) +dcmplx(H(ind,1),   0.0_dp)
           Haux(2,jo,2,io) = Haux(2,jo,2,io) +dcmplx(H(ind,2),   0.0_dp)
           Haux(2,jo,1,io) = Haux(2,jo,1,io) +dcmplx(H(ind,3), H(ind,4))
           Haux(1,jo,2,io) = Haux(1,jo,2,io) +dcmplx(H(ind,3),-H(ind,4))
         enddo
        enddo
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo,psi,
     .                         2*nuotot,1,ierror,2*BlockSize)
      endif

C Rebuild the full overlap matrix
      Spr = 0.0_dp
      do io = 1, nuo
        do j = 1, numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          jo = indxuo(jo)
          Spr(jo,io) = Spr(jo,io) + S(ind)
        enddo
      enddo

#ifdef MPI
C Find maximum number of orbitals per node
        call MPI_AllReduce(nuo,maxnuo,1,MPI_integer,MPI_max,
     .                   MPI_Comm_World,MPIerror)

C Allocate workspace array for broadcast overlap matrix
        nullify( Sloc )
        call re_alloc( Sloc, 1, nuotot, 1, maxnuo,
     .                 name='Sloc', routine='pdos2g' )

C Loop over nodes broadcasting overlap matrix
        do BNode = 0,Nodes-1

C Find out how many orbitals there are on the broadcast node
          call GetNodeOrbs(nuotot,BNode,Nodes,Bnuo)

C Transfer data
          if (Node.eq.BNode) then
            Sloc(1:nuotot,1:Bnuo) = Spr(1:nuotot,1:Bnuo)
          endif
          call MPI_Bcast(Sloc(1,1),nuotot*Bnuo,
     .      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)

          do ihist = 1, nhist
            ener = E1 + (ihist - 1) * delta
            do iband = 1, nuo*2
              call LocalToGlobalOrb((iband+1)/2,Node,Nodes,ibandg)
              ibandg = ibandg * 2 - mod(iband, 2) 
              diff = (ener - eo(ibandg))**2 / (sigma ** 2)
              if (diff .gt. 15.0d0) cycle
              gauss = exp(-diff)
              caux => psi(:,:,iband) ! c_{up,j}, c_{down,j}
              do jo = 1, Bnuo
                call LocalToGlobalOrb(jo,BNode,Nodes,juo)
                do io = 1, nuotot
                  rpipj = dreal(caux(1,io)*dconjg(caux(1,juo)))*gauss
                  dpr(ihist,juo,1)=dpr(ihist,juo,1)+rpipj*Sloc(io,jo)

                  rpipj = dreal(caux(2,io)*dconjg(caux(2,juo)))*gauss
                  dpr(ihist,juo,2)=dpr(ihist,juo,2)+rpipj*Sloc(io,jo)

                  rpipj = dreal(caux(1,io)*dconjg(caux(2,juo)))*gauss
                  dpr(ihist,juo,3)=dpr(ihist,juo,3)+rpipj*Sloc(io,jo)

                  ipipj = dimag(caux(1,io)*dconjg(caux(2,juo)))*gauss
                  dpr(ihist,juo,4)=dpr(ihist,juo,4)-ipipj*Sloc(io,jo)
                enddo
              enddo
            enddo ! inband loop
          enddo ! ihist loop
        enddo ! BNode loop

C Free workspace array for overlap
        call de_alloc(Sloc, 'Sloc', 'pdos2g')

#else
C Loop over all the energy range

        do ihist = 1, nhist
          ener = E1 + (ihist - 1) * delta
          do iband = 1, nuo*2
            diff = (ener - eo(iband))**2 / (sigma ** 2)
            if (diff .gt. 15.0d0) cycle
            gauss = exp(-diff)
            caux => psi(:,:,iband) ! c_{up,j}, c_{down,j}
            do jo = 1, nuotot
              do io = 1, nuotot
                rpipj = dreal(caux(1,io)*dconjg(caux(1,jo)))*gauss
                dpr(ihist,jo,1)=dpr(ihist,jo,1)+rpipj*Spr(io,jo)

                rpipj = dreal(caux(2,io)*dconjg(caux(2,jo)))*gauss
                dpr(ihist,jo,2)=dpr(ihist,jo,2)+rpipj*Spr(io,jo)

                rpipj = dreal(caux(1,io)*dconjg(caux(2,jo)))*gauss
                dpr(ihist,jo,3)=dpr(ihist,jo,3)+rpipj*Spr(io,jo)

                ipipj = dimag(caux(1,io)*dconjg(caux(2,jo)))*gauss
                dpr(ihist,jo,4)=dpr(ihist,jo,4)-ipipj*Spr(io,jo)
              enddo
            enddo
          enddo
        enddo
#endif

      call de_alloc(Spr, 'Spr', 'pdos2g')
        
#ifdef MPI

C Global reduction of dpr matrix
      do io= 1, nuotot
        do ihist = 1, nhist
          tmp(ihist,io,1) = 0.0d0
          tmp(ihist,io,2) = 0.0d0
          tmp(ihist,io,3) = 0.0d0
          tmp(ihist,io,4) = 0.0d0
        enddo
      enddo
      call MPI_AllReduce(dpr(1,1,1),tmp(1,1,1),nhist*nuotot*4,
     .  MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
      do io= 1, nuotot
        do ihist = 1, nhist
          dpr(ihist,io,1) = tmp(ihist,io,1)
          dpr(ihist,io,2) = tmp(ihist,io,2)
          dpr(ihist,io,3) = tmp(ihist,io,3)
          dpr(ihist,io,4) = tmp(ihist,io,4)
        enddo
      enddo

#endif

      norm = sigma * sqrt(pi)
      dpr = dpr / norm
      
      do ihist = 1, nhist
        dtot(ihist,1) = sum(dpr(ihist,:,1))
        dtot(ihist,2) = sum(dpr(ihist,:,2))
        dtot(ihist,3) = sum(dpr(ihist,:,3))
        dtot(ihist,4) = sum(dpr(ihist,:,4))
      enddo

      return
      end
