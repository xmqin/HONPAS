!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diag3g( nuo, no, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr,
     .                   listd, H, S, getD, getPSI, qtot, temp, e1, e2,
     .                   eo, qo, Dnew, Enew, ef, Entropy,
     .                   Haux, Saux, psi, caux, 
     .                   nuotot, occtol, iscf, neigwanted)

      use precision
      use sys
      use parallel,        only : Node, Nodes, BlockSize
      use parallelsubs,    only : LocalToGlobalOrb,GlobalToLocalOrb
      use writewave,       only : writew
      use m_fermid,        only : fermid, stepf
      use intrinsic_missing, only: MODP
      use iso_c_binding,   only : c_loc, c_f_pointer
#ifdef MPI
      use mpi_siesta
#endif
      
      implicit none
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for spin-orbit at gamma point.
C Writen by J.Soler, May and August 1998.
C Reduced memory requirements, Nick, Aug. 2017
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals on local node
C integer no                  : Number of basis orbitals
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
C complex*16 Haux(2,nuotot,2,nuo): Auxiliary space for the hamiltonian matrix
C complex*16 Saux(2,nuotot,2,nuo): Auxiliary space for the overlap matrix
C complex*16 psi(2,nuotot,2*nuo) : Auxiliary space for the eigenvectors
C complex*16 caux(2,nuotot)      : Extra auxiliary space
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
      integer maxo, maxnd, maxnh
      integer no, nuo, nuotot, iscf, neigwanted
      
      integer listh(maxnh), numh(nuo), listhptr(nuo)
      integer listd(maxnd), numd(nuo), listdptr(nuo)
      
      real(dp) Dnew(maxnd,8), e1, e2, ef, Enew(maxnd,4)
      real(dp) Entropy, eo(maxo*2), H(maxnh,8), qo(maxo*2)
      real(dp) qtot, S(maxnh), temp, occtol
      
      complex(dp), dimension(2,nuotot,2*nuo), target :: psi
      real(dp), pointer                              :: psi_real_1d(:)
      complex(dp), dimension(2,nuotot,2,nuo) :: Haux, Saux
      complex(dp), dimension(2,nuotot) :: caux
      logical               getD, getPSI
      
!     Internal variables .............................................
      
      integer           BNode, BTest, ie, ierror, iie, iio
      integer           ind, io, j, jo, nd, iuo, juo
      real(dp)          ee, qe, t, k(3)
      complex(dp)       D11, D22, D12, D21
#ifdef MPI
      integer MPIerror
#endif
      external              cdiag
!***********************************************************************
!     B E G I N
!***********************************************************************

!***********************************************************************
! BUILD HAMILTONIAN
!***********************************************************************
! The different subroutines build H_{j,i}^{js,is} = <js,j|H|is,i>
! The spin notation is as follows:
!
!            | H_{j,i}^{u,u}  H_{j,i}^{u,d} |
!  H_(j,i} = |                              |
!            | H_{j,i}^{d,u}  H_{j,i}^{d,d} |
!
!            | H(ind,1) + i H(ind,5)   H(ind,3) - i H(ind,4) |
!          = |                                               |
!            | H(ind,7) + i H(ind,8)   H(ind,2) + i H(ind,6) |
!
! 1. Hermiticity imposes H_{i,j}^{is,js}=H_{j,i}^{js,is}^*
! 2. Since wave functions are real, if there are no single P or L operators:
!      (a) H_{i,j}^{is,js}=H_{j,i}^{is,js}
!      (b) These imply spin-box hermiticity: H_{i,j}^{is,js}=H_{i,j}^{js,is}^*
!
!      (c) Hence
!
!                               | H(ind,1)                H(ind,3) - i H(ind,4) |
!           H_{j,i} = H_{i,j} = |                                               |
!                               | H(ind,3) + i H(ind,4)   H(ind,2)              |
!
!
! The spin-orbit interaction and the orbital part of the Zeeman interaction
! break the  "spin-box hermiticity" of H. Hence the full format of H is needed
! in those cases. 
!
! Spin-box symmetrization of D:
! D_{i,j}(1,2) = 0.5 ( D_{i,j}(1,2) + D_{i,j}(2,1)^* )
! can not be enforced here.

      do io = 1,nuo
        Saux(:,:,:,io) = 0._dp
        Haux(:,:,:,io) = 0._dp
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          jo = MODP(jo,nuotot) ! To allow auxiliary supercells
          Saux(1,jo,1,io) = Saux(1,jo,1,io) + dcmplx( S(ind), 0.0_dp)
          Saux(2,jo,2,io) = Saux(2,jo,2,io) + dcmplx( S(ind), 0.0_dp)
          Haux(1,jo,1,io) = Haux(1,jo,1,io) + dcmplx(H(ind,1), H(ind,5))
          Haux(2,jo,2,io) = Haux(2,jo,2,io) + dcmplx(H(ind,2), H(ind,6))
          Haux(1,jo,2,io) = Haux(1,jo,2,io) + dcmplx(H(ind,3),-H(ind,4))
          Haux(2,jo,1,io) = Haux(2,jo,1,io) + dcmplx(H(ind,7), H(ind,8))
        enddo
      enddo
      
!     Solve the eigenvalue problem
      call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo,psi,
     .           2*neigwanted,iscf,ierror,2*BlockSize)
      
!     Check error flag and take appropriate action
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror.lt.0) then
!       Repeat diagonalisation with increased memory to handle clustering
      do io = 1,nuo
        Saux(:,:,:,io) = 0._dp
        Haux(:,:,:,io) = 0._dp
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          jo = MODP(jo,nuotot) ! To allow auxiliary supercells
          Saux(1,jo,1,io) = Saux(1,jo,1,io) + dcmplx( S(ind), 0.0_dp)
          Saux(2,jo,2,io) = Saux(2,jo,2,io) + dcmplx( S(ind), 0.0_dp)
          Haux(1,jo,1,io) = Haux(1,jo,1,io) + dcmplx(H(ind,1), H(ind,5))
          Haux(2,jo,2,io) = Haux(2,jo,2,io) + dcmplx(H(ind,2), H(ind,6))
          Haux(1,jo,2,io) = Haux(1,jo,2,io) + dcmplx(H(ind,3),-H(ind,4))
          Haux(2,jo,1,io) = Haux(2,jo,1,io) + dcmplx(H(ind,7), H(ind,8))
        enddo
      enddo
      
        call cdiag(Haux,Saux,2*nuotot,2*nuo,2*nuotot,eo,psi,
     .             2*neigwanted,iscf,ierror,2*BlockSize)
      endif

      if (getPSI) then
         k(1:3) = 0.0d0

          ! We do not have info about occupation
          ! In the current implementation, neigneeded
          ! is ALWAYS nuotot.
          ! This is because writewave expects the full spectrum
         call c_f_pointer( c_loc(psi), psi_real_1d, [size(psi)*2] )
         call writew(nuotot,nuo,1,k,1,
     &        eo,psi_real_1d,gamma=.true.,
     $        non_coll=.true.,blocksize=2*BlockSize)
      endif

!     Check if we are done
      if (.not.getD) return
      
!     Find new Fermi energy and occupation weights 

      call fermid(2, 1, 1, (/1._dp/), 2*maxo, 2*neigwanted,
     .                     eo, temp, qtot, qo, ef, Entropy)

!      write(6,'(a,f12.6)') ' Ef = ', ef
      
!     Find weights for local density of states ............................
      if (e1 .lt. e2) then
        t = max( temp, 1.d-6 )
        do io = 1, nuotot*2
          qo(io) = ( stepf((eo(io)-e2)/t) - stepf((eo(io)-e1)/t))
        enddo
      endif
      
!***********************************************************************
! BUILD NEW DENSITY MATRIX
!***********************************************************************
!
!                 | ------- 1,1 -------     ------- 1,2 --------- |
!                 | c_{i,up} c_{j,up}^*     c_{i,up} c_{j,down)^* |
!     D_{i,j} =   |                                               |
!                 | ------- 2,1 -------     ------- 2,2 --------- |
!                 | c_{i,down} c_{j,up}^*     c_{i,dn} c_{j,dn)^* |
!
!             =   | D_{i,j}(1)              D_{i,j}(3)-i D_{i,j}(4) |
!                 | D_{i,j}(7)+i D_{i,j}(8) D_{i,j}(2)              |

! The Energy is computed as E = Tr [ D_{i,j} H_{j,i} ]
!
! The Density Matrix is not "spin box hermitian" even if H does not contain P or
! L operators. 
!
! If H contains P or L operators, then spin-box "symetrization" of D leads to erroneous
! results for E

      Dnew(:,:) = 0.0_dp
      Enew(:,:) = 0.0_dp
      
      BNode = 0
      iie   = 0
      
      do ie = 1,2*nuotot
        qe = qo(ie)
        if (Node.eq.BNode) then
          iie = iie + 1
        endif
      
        caux(:,:) = dcmplx(0.0_dp,0.0_dp)
        if (abs(qe).gt.occtol) then
          if (Node.eq.BNode) then
            do j = 1,nuotot
              caux(1,j) = psi(1,j,iie) ! c_{i,up}
              caux(2,j) = psi(2,j,iie) ! c_{i,dn}
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(caux(1,1),2*nuotot,MPI_double_complex,
     .         BNode,MPI_Comm_World,MPIerror)
#endif
          ee = qo(ie) * eo(ie)
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            do j = 1,numd(io)
              ind = listdptr(io) + j
              jo = listd(ind)
              jo = MODP(jo,nuotot) ! To allow auxiliary supercells
!------- 1,1 -----------------------------------------------------------
              D11  = caux(1,iio) * dconjg(caux(1,jo))
!------- 2,2 -----------------------------------------------------------
              D22  = caux(2,iio) * dconjg(caux(2,jo))
!------- 1,2 -----------------------------------------------------------
              D12  = caux(1,iio) * dconjg(caux(2,jo))
!------- 2,1 -----------------------------------------------------------
              D21  = caux(2,iio) * dconjg(caux(1,jo))
      
              Dnew(ind,1) = Dnew(ind,1) + dreal(D11) * qe
              Dnew(ind,2) = Dnew(ind,2) + dreal(D22) * qe
              Dnew(ind,3) = Dnew(ind,3) + dreal(D12) * qe
              Dnew(ind,4) = Dnew(ind,4) - dimag(D12) * qe
              Dnew(ind,5) = Dnew(ind,5) + dimag(D11) * qe
              Dnew(ind,6) = Dnew(ind,6) + dimag(D22) * qe
              Dnew(ind,7) = Dnew(ind,7) + dreal(D21) * qe
              Dnew(ind,8) = Dnew(ind,8) + dimag(D21) * qe
      
              Enew(ind,1) = Enew(ind,1) + dreal(D11) * ee
              Enew(ind,2) = Enew(ind,2) + dreal(D22) * ee
              Enew(ind,3) = Enew(ind,3) + dreal(D12) * ee
              Enew(ind,4) = Enew(ind,4) - dimag(D12) * ee
      
            enddo
          enddo
        endif
        BTest = ie/(2*BlockSize)
        if (BTest*2*BlockSize.eq.ie) then
          BNode = BNode + 1
          if (BNode .gt. Nodes-1) BNode = 0
        endif
      enddo

!***********************************************************************
      return
      end
!**********************************************************************
      
