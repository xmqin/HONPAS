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
      subroutine diag2g( nuo, no, maxuo, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr,
     .                   listd, H, S, getD, qtot, temp, e1, e2,
     .                   eo, qo, Dnew, Enew, ef, Entropy, psi, 
     .                   nuotot, occtol, iscf, neigwanted)

!
!     Contributed by Volodymyr Maslyuk
!      
      use precision
      use sys
      use parallel,        only : Node, Nodes, BlockSize
      use parallelsubs,    only : LocalToGlobalOrb,GlobalToLocalOrb
      use m_fermid,        only : fermid, stepf
      use m_spin,          only : MColl
#ifdef MPI
      use mpi_siesta
#endif
      
      implicit none
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin at gamma point.
C Writen by J.Soler, May and August 1998.
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals on local node
C integer no                  : Number of basis orbitals
C integer maxuo               : Maximum number of basis  orbitals
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of Dnew and Enew
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,4)          : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer iscf                : SCF cycle number
C integer neigwanted          : Number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo*2)           : Eigenvalues
C real*8 qo(maxo*2)           : Occupations of eigenstates
C real*8 Dnew(maxnd,4)        : Output Density Matrix
C real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,2,nuotot,2,nuo): Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,2,nuotot,2,nuo): Auxiliary space for the overlap matrix
C real*8 psi(2,2,nuotot,2*no) : Auxiliary space for the eigenvectors
C real*8 aux(5,2*nuotot)      : Extra auxiliary space
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
! *********************************************************************
!     Haux(js,juo,is,iuo) = <js,juo|H|is,iuo>
!     Indices is and js are for spin components 
!     Indices iuo and juo are for orbital components:
! *********************************************************************
      integer maxo, maxuo, maxnd, maxnh
      integer no, nuo, nuotot, iscf, neigwanted
      
      integer listh(maxnh), numh(nuo), listhptr(nuo)
      integer listd(maxnd), numd(nuo), listdptr(nuo)
      
      real(dp) Dnew(maxnd,4), e1, e2, ef, Enew(maxnd,4)
      real(dp) Entropy, eo(maxo*2), H(maxnh,4), qo(maxo*2)
      real(dp) qtot, S(maxnh), temp, occtol
      
      logical               getD
      
!     Internal variables .............................................
      real(dp)                       psi(2,2,nuotot,2*nuo)
      complex(dp),dimension(:,:,:,:), allocatable :: Haux, Saux
      complex(dp),dimension(:,:),     allocatable :: caux
      
      integer           BNode, BTest, ie, ierror, iie, iio
      integer           ind, io, j, jo, nd, iuo, juo
      real(dp)          ee, pipj, qe, t
      complex(dp)       cicj, D11, D22, D12, D21
#ifdef MPI
      integer MPIerror
#endif
      external              cdiag
!***********************************************************************
!     B E G I N
!***********************************************************************
!     define the arrays
      allocate(Haux(MColl,nuotot,MColl,nuo))
      call memory('A','Z',MColl*nuotot*MColl*nuo,'diag2g')
      allocate(Saux(MColl,nuotot,MColl,nuo))
      call memory('A','Z',MColl*nuotot*MColl*nuo,'diag2g')
      allocate(caux(MColl,nuotot))
      call memory('A','Z',MColl*nuotot,'diag2g')
      
      Saux = dcmplx(0.0_dp, 0.0_dp)
      Haux = dcmplx(0.0_dp, 0.0_dp)
      do io = 1,nuo
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
           Saux(1,jo,1,io) = dcmplx(  S(ind),   0.0_dp)
           Saux(2,jo,2,io) = dcmplx(  S(ind),   0.0_dp)
           Haux(1,jo,1,io) = dcmplx(H(ind,1),   0.0_dp)
           Haux(2,jo,2,io) = dcmplx(H(ind,2),   0.0_dp)
           Haux(2,jo,1,io) = dcmplx(H(ind,3), H(ind,4))
           Haux(1,jo,2,io) = dcmplx(H(ind,3),-H(ind,4))
        enddo
      enddo
      
!     Solve the eigenvalue problem
      call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo,psi,
     .                       MColl*neigwanted,iscf,ierror)
      
!     Check error flag and take appropriate action
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror.lt.0) then
!       Repeat diagonalisation with increased memory to handle clustering
        Saux = dcmplx(0.0_dp, 0.0_dp)
        Haux = dcmplx(0.0_dp, 0.0_dp)
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
             Saux(1,jo,1,io) = dcmplx(  S(ind),   0.0_dp)
             Saux(2,jo,2,io) = dcmplx(  S(ind),   0.0_dp)
             Haux(1,jo,1,io) = dcmplx(H(ind,1),   0.0_dp)
             Haux(2,jo,2,io) = dcmplx(H(ind,2),   0.0_dp)
             Haux(2,jo,1,io) = dcmplx(H(ind,3), H(ind,4))
             Haux(1,jo,2,io) = dcmplx(H(ind,3),-H(ind,4))
          enddo
        enddo
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo,psi,
     .                         MColl*neigwanted,iscf,ierror)
      endif
      
!     Check if we are done
      if (.not.getD) return
      
!     Find new Fermi energy and occupation weights 
      call fermid(1, 1, 1, (/0.5_dp/), MColl*maxo, MColl*neigwanted, 
     .                     eo, temp, qtot, qo, ef, Entropy)
      
!     Find weights for local density of states ............................
      if (e1 .lt. e2) then
        t = max( temp, 1.d-6 )
        do io = 1, nuotot*MColl
          qo(io) = ( stepf((eo(io)-e2)/t) - stepf((eo(io)-e1)/t))
        enddo
      endif
      
!     New density and energy-density matrices of unit-cell orbitals .......
      Dnew(:,:) = 0.0_dp
      Enew(:,:) = 0.0_dp
      
      BNode = 0
      iie   = 0
      
      do ie = 1,MColl*nuotot
        qe = qo(ie)
        if (Node.eq.BNode) then
          iie = iie + 1
        endif
      
        caux(:,:) = dcmplx(0.0_dp,0.0_dp)
        if (abs(qe).gt.occtol) then
          if (Node.eq.BNode) then
            do j = 1,nuotot
              caux(1,j)=dcmplx( psi(1,1,j,iie), psi(2,1,j,iie)) ! c_{i,up}
              caux(2,j)=dcmplx( psi(1,2,j,iie), psi(2,2,j,iie)) ! c_{i,dn}
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(caux(1,1),MColl*nuotot,MPI_double_complex,
     .         BNode,MPI_Comm_World,MPIerror)
#endif
          ee = qo(ie) * eo(ie)
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do j = 1,numd(io)
              ind = listdptr(io) + j
              jo = listd(ind)
      
!                 | ------- 1,1 -------     ------- 2,1 ------- |
!                 | c_{j,up} c_{i,up}^*     c_{j,dn} c_{i,up)^* |
!     D_{j,i} =   |                                             |
!                 | ------- 1,2 -------     ------- 2,2 ------- |
!                 | c_{j,up} c_{i,dn}^*     c_{j,dn} c_{i,dn)^* |
!
!------- 1,1 -----------------------------------------------------------
              D11  = dconjg(caux(1,iio)) * caux(1,jo)
!------- 2,2 -----------------------------------------------------------
              D22  = dconjg(caux(2,iio)) * caux(2,jo)
!------- 2,1 -----------------------------------------------------------
              D12  = dconjg(caux(1,iio)) * caux(2,jo)
!------- 1,2 -----------------------------------------------------------
              D21  = dconjg(caux(2,iio)) * caux(1,jo)
      
!------------ Density matrix has to be hermitian ----------------------
              D11 = 0.5_dp * (D11 + dconjg(D11))
              D22 = 0.5_dp * (D22 + dconjg(D22))
              D12 = 0.5_dp * (D12 + dconjg(D21))
              D21 = dconjg(D12)
      
!       Add contribution to density matrices of unit-cell orbitals
!       ----------------------------------------------------------------
!       | D11 = D_{j,i}(1)               D21 = D_{j,i}(3)-i D_{j,i}(4) |
!       | D12 = D_{i,j}(3)+i D_{i,j}(4)  D22 = D_{j,i}(2)              |
!       ----------------------------------------------------------------
!       ----------------------------------------------------------------
!       | D11 = Dnew(1)                  D21 = Dnew(3)-i Dnew(4)       |
!       | D12 = Dnew(3)+i Dnew(4)        D22 = Dnew(2)                 |
!       ----------------------------------------------------------------
              Dnew(ind,1) = Dnew(ind,1) + dreal(D11) * qe
              Dnew(ind,2) = Dnew(ind,2) + dreal(D22) * qe
              Dnew(ind,3) = Dnew(ind,3) + dreal(D12) * qe
              Dnew(ind,4) = Dnew(ind,4) + dimag(D12) * qe
      
              Enew(ind,1) = Enew(ind,1) + dreal(D11) * ee
              Enew(ind,2) = Enew(ind,2) + dreal(D22) * ee
              Enew(ind,3) = Enew(ind,3) + dreal(D12) * ee
              Enew(ind,4) = Enew(ind,4) + dimag(D12) * ee
      
            enddo
          enddo
        endif
        BTest = ie/(MColl*BlockSize)
        if (BTest*MColl*BlockSize.eq.ie) then
          BNode = BNode + 1
          if (BNode .gt. Nodes-1) BNode = 0
        endif
      enddo
      
!     deallocate the arrays
      call memory('D','Z',MColl*nuotot*MColl*nuo,'diag2g')
      deallocate(Haux)
      call memory('D','Z',MColl*nuotot*MColl*nuo,'diag2g')
      deallocate(Saux)
      call memory('D','Z',MColl*nuotot,'diag2g')
      deallocate(caux)
!***********************************************************************
      return
      end
!**********************************************************************
      
