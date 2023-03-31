! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diagg( h_spin_dim, nuo, maxuo, maxnh, maxnd, 
     &                  maxo, numh, listhptr, listh, numd, 
     &                  listdptr, listd, H, S,
     &                  getD, getPSI, fixspin, qtot, qs, temp, e1, e2,
     &                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     &                  Haux, Saux, psi, nuotot, occtol,
     &                  iscf, neigwanted )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer h_spin_dim    : Number of spin components of H and D
C integer spinor_dim          : # of spin components of eo and qo, qs, efs
C integer e_spin_dim          : Number of spin components of E_dm  
C integer nuo                 : Number of basis orbitals local to node
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,h_spin_dim) : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical getPSI              : Find and print wavefunctions?
C logical fixspin             : Fix the spin of the system?
C real*8  qtot                : Number of electrons in unit cell
C real*8  qs(spinor_dim)      : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C integer iscf                : SCF cycle number
C real*8  occtol              : Occupancy threshold for DM build
C integer neigwanted          : Number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo,spinor_dim)  : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,spinor_dim)   : Occupations of eigenstates
C real*8 Dnew(maxnd,h_spin_dim): Output Density Matrix
C real*8 Enew(maxnd,e_spin_dim): Output Energy-Density Matrix
C real*8 ef                    : Fermi energy
C real*8 efs(spinor_dim)       : Fermi energy for each spin
C                                 (for fixed spin calculations)
C real*8 Entropy               : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
C real*8 psi(nuotot,maxuo,spinor_dim) : Auxiliary space for the eigenvectors
C real*8 aux(nuotot)          : Extra auxiliary space
C *************************** UNITS ***********************************
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave,     only : writew
      use m_fermid,      only : fermid, fermispin, stepf
      use m_spin,        only : spinor_dim, e_spin_dim
      use alloc
      use intrinsic_missing, only: MODP
      use iso_c_binding, only : c_f_pointer, c_loc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     &  MPIerror
#endif

      integer
     &  iscf, maxnd, maxnh, maxuo, maxo, nuo, nuotot,
     &  neigwanted, h_spin_dim

      integer 
     &  listh(maxnh), numh(nuo), listhptr(nuo),
     &  listd(maxnd), numd(nuo), listdptr(nuo)

      real(dp)
     &  Dnew(maxnd,h_spin_dim), e1, e2, ef, Enew(maxnd,e_spin_dim), 
     &  Entropy, eo(maxo,spinor_dim), H(maxnh,h_spin_dim), 
     &  qo(maxo,spinor_dim), qtot, qs(spinor_dim), S(maxnh), temp, 
     &  efs(spinor_dim), occtol
     
      real(dp) Haux(nuotot,nuo), Saux(nuotot,nuo)
      real(dp), target :: psi(nuotot,maxuo,spinor_dim) 
      real(dp), pointer :: psi_real_1d(:)
      real(dp), dimension(1), parameter :: wk = (/ 1.0_dp /)
      integer, parameter                :: nk = 1

      logical
     .  getD, getPSI, fixspin

      external rdiag

C  Internal variables .............................................
      integer           ie, io, iio, ispin, ix, j, jo, BNode, iie, ind,
     &                  ierror, nd
      real(dp)          ee, eei, qe, qei, rt, t, k(3)

      integer           :: maxnuo, mm
      integer, pointer  :: nuo_LOC(:)
      real(dp), pointer :: psi_tmp(:), paux(:)

#ifdef DEBUG
      call write_debug( '    PRE diagg' )
#endif
      if (h_spin_dim /= spinor_dim) then
         call die("Spin size mismatch in diagg")
      endif

C Solve eigenvalue problem
      
      do ispin = 1,spinor_dim

        call timer( 'r-eigvec', 1 )
        call timer( 'r-buildHS', 1 )
!$OMP parallel do default(shared), private(io,j,ind,jo)
        do io = 1,nuo
          Saux(:,io) = 0._dp
          Haux(:,io) = 0._dp
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            jo = MODP(jo,nuotot)   ! To allow auxiliary supercells
            Saux(jo,io) = Saux(jo,io) + S(ind)
            Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
          enddo
        enddo
!$OMP end parallel do
        call timer( 'r-buildHS', 2 )
        call rdiag( Haux, Saux, nuotot, nuo, nuotot,
     $       eo(1,ispin), 
     .       psi(1,1,ispin), neigwanted, iscf, ierror,
     .       BlockSize)
        call timer( 'r-eigvec', 2 )


C Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
!$OMP parallel do default(shared), private(io,j,ind,jo)
          do io = 1,nuo
            Saux(:,io) = 0._dp
            Haux(:,io) = 0._dp
            do j = 1,numh(io)
              ind = listhptr(io) + j
              jo = listh(ind)
              jo = MODP(jo,nuotot) ! To allow auxiliary supercells
              Saux(jo,io) = Saux(jo,io) + S(ind)
              Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
            enddo
          enddo
!$OMP end parallel do
          call rdiag( Haux, Saux, nuotot, nuo, nuotot, eo(1,ispin), 
     &         psi(1,1,ispin), neigwanted, iscf, ierror,
     .         BlockSize )
        endif

        if (getPSI) then
          ! We do not have info about occupation
          ! In the current implementation, neigneeded
          ! is ALWAYS nuotot.
          ! This is because writewave expects the full spectrum
          do ix = 1,3
            k(ix) = 0.0d0
          enddo
          call c_f_pointer(c_loc(psi(:,:,ispin)),psi_real_1d,
     $                     [size(psi,1)*size(psi,2)])
          call writew(nuotot,nuo,1,k,ispin,
     &                eo(1,ispin),psi_real_1d,gamma=.true.,
     $                non_coll=.false.,blocksize=BlockSize)
        endif

      enddo
C Check if we are done ................................................
      if (.not.getD) then
#ifdef DEBUG
         call write_debug( '    POS diagg' )
#endif
         return
      end if

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( spinor_dim, spinor_dim, nk, wk, maxo,
     &                  neigwanted, eo, temp, qs, qo, efs, Entropy )
      else
        call fermid(spinor_dim, spinor_dim, nk, wk, maxo, neigwanted,
     &               eo, temp, qtot, qo, ef, Entropy )
      endif


C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        rt = 1.0d0/t
!$OMP parallel do default(shared), private(ispin,io)
        do ispin = 1,spinor_dim
          do io = 1,nuotot
            qo(io,ispin) = ( stepf((eo(io,ispin)-e2)*rt) -
     &                       stepf((eo(io,ispin)-e1)*rt)) * 
     &                       2.0d0/dble(spinor_dim)
          enddo
        enddo
!$OMP end parallel do
      endif
      
      if (nuo.gt.0) then
C New density and energy-density matrices of unit-cell orbitals .......
        nd = listdptr(nuo) + numd(nuo)
        Dnew(1:nd,1:h_spin_dim) = 0.0d0
        Enew(1:nd,1:e_spin_dim) = 0.0d0
      endif

      call timer( 'r-buildD', 1 )
C Global operation to form new density matrix
      nullify( nuo_LOC )
      call re_alloc( nuo_LOC, 0, Nodes-1, 'nuo_LOC', 'diagg' )
#ifdef MPI
      call MPI_Allgather( nuo, 1, MPI_INTEGER, nuo_LOC, 1,
     &                    MPI_INTEGER, MPI_Comm_World, MPIerror )
#else
      nuo_LOC(0) = nuo
#endif
      maxnuo = 0
      do BNode=0, Nodes-1
        maxnuo = max(nuo_LOC(BNode),maxnuo)
      enddo
      nullify( PSI_TMP )
      call re_alloc( PSI_TMP, 1, nuotot*maxnuo*spinor_dim, 'PSI_TMP',
     &               'diagg' )

      do Bnode=0, Nodes-1
        if (BNode.eq.Node) then
          mm = 0
          do ispin= 1, spinor_dim
            do io= 1, nuo
              PSI_TMP(mm+1:mm+nuotot) = PSI(1:nuotot,io, ispin)
              mm = mm + nuotot
            enddo
          enddo
        endif
#ifdef MPI
        call MPI_Bcast( PSI_TMP, nuotot*nuo_LOC(Bnode)*spinor_dim,
     &                  MPI_double_precision, BNode, MPI_Comm_World,
     &                  MPIerror )
#endif
        do ispin = 1,h_spin_dim
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            mm   = (ispin-1)*nuo_LOC(Bnode)*nuotot
            do iie = 1,nuo_LOC(Bnode)
              call  LocalToGlobalOrb(iie,BNode,Nodes,ie)
              qe = qo(ie,ispin)
              if (abs(qe).gt.occtol) then
                paux => PSI_TMP(mm+1:mm+nuotot)
                qei = qe*paux(iio)
                eei = qei*eo(ie,ispin)
                do j = 1,numd(io)
                  ind = listdptr(io) + j
                  jo  = listd(ind)
                  jo = MODP(jo,nuotot) ! To allow auxiliary supercells
                  Dnew(ind,ispin) = Dnew(ind,ispin) + qei*paux(jo)
                  Enew(ind,ispin) = Enew(ind,ispin) + eei*paux(jo)
                enddo
              endif
              mm = mm + nuotot
            enddo
          enddo
        enddo
      enddo

      call timer( 'r-buildD', 2 )
      call de_alloc( nuo_LOC, 'nuo_LOC', 'diagg' )
      call de_alloc( PSI_TMP, 'PSI_TMP', 'diagg' )

#ifdef DEBUG
      call write_debug( '    POS diagg' )
#endif
      end 
