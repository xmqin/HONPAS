! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

#ifdef MPI
subroutine diag3kp( spin, no_l, no_u, no_s, nnz, &
    ncol, ptr, col, H, S, getD, &
    qtot, temp, e1, e2, xij, &
    nk, kpoint, wk, eo, qo, DM, EDM, ef, &
    Entropy, Hk, Sk, psi, aux, &
    occtol, iscf, neigwanted)
  
  ! *********************************************************************
  ! Subroutine to calculate the eigenvalues and eigenvectors, density
  ! and energy-density matrices, and occupation weights of each 
  ! eigenvector, for given Hamiltonian and Overlap matrices with spin-orbit
  ! K-sampling version.
  ! Created by Nick Papior, 2020
  ! Uses parallelisation over K points instead of parallelisation 
  ! within them.
  ! **************************** INPUT **********************************
  ! type(t_spin) spin           : Spin type
  ! integer no_l                : Number of basis orbitals in unit cell
  !                               local to this processor
  ! integer no_u                : Number of basis orbitals in unit cell
  ! integer no_s                : Number of basis orbitals in supercell
  ! integer nnz                 : Maximum number of orbitals interacting
  ! integer ncol(no_l)          : Number of nonzero elements of each row 
  !                               of hamiltonian/density matrix locally
  ! integer ptr(no_l)           : Pointer to each row (-1) of the
  !                               hamiltonian matrix locally
  ! integer col(nnz)            : Nonzero hamiltonian-matrix element  
  !                               column indexes for each matrix row
  ! real*8  H(nnz,spin%H)       : Hamiltonian in sparse form
  ! real*8  S(nnz)              : Overlap in sparse form
  ! logical getD                : Find occupations and density matrices?
  ! real*8  qtot                : Number of electrons in unit cell
  ! real*8  temp                : Electronic temperature 
  ! real*8  e1, e2              : Energy range for density-matrix states
  !                               (to find local density of states)
  !                               Not used if e1 > e2
  ! real*8  xij(3,nnz)          : Vectors between orbital centers (sparse)
  !                               (not used if only gamma point)
  ! integer nk                  : Number of k points
  ! real*8  kpoint(3,nk)        : k point vectors
  ! real*8  wk(nk)              : k point weights (must sum one)
  ! real*8  occtol              : Occupancy threshold for DM build
  ! integer neigwanted          : maximum number of eigenvalues wanted
  ! *************************** OUTPUT **********************************
  ! real*8 eo(maxo*spin%spinor,nk)   : Eigenvalues
  ! ******************** OUTPUT (only if getD=.true.) *******************
  ! real*8 qo(maxo*spin%spinor,nk)   : Occupations of eigenstates
  ! real*8 Dnew(maxnd,spin%DM)    : Output Density Matrix
  ! real*8 Enew(maxnd,spin%EDM)    : Output Energy-Density Matrix
  ! real*8 ef                   : Fermi energy
  ! real*8 Entropy              : Electronic entropy
  ! *************************** AUXILIARY *******************************
  ! complex*16 Hk(2,no_u,2,no_u)   : Auxiliary space for the hamiltonian matrix
  ! complex*16 Sk(2,no_u,2,no_u)   : Auxiliary space for the overlap matrix
  ! complex*16 psi(2,no_u,no_u*2)  : Auxiliary space for the eigenvectors
  ! real*8     aux(2*no_u)         : Extra auxiliary space
  ! *************************** UNITS ***********************************
  ! xij and kpoint must be in reciprocal coordinates of each other.
  ! temp and H must be in the same energy units.
  ! eo, Enew and ef returned in the units of H.
  ! *************************** PARALLEL ********************************
  ! The auxiliary arrays are now no longer symmetric and so the order
  ! of referencing has been changed in several places to reflect this.
  ! Note : It is assumed in a couple of places that the sparsity of
  ! H/S and Dscf are the same (which is the case in the code at present)
  ! and therefore the only one set of pointer arrays are globalised.
  ! *********************************************************************
  !
  !  Modules
  !
  use precision
  use sys
  use parallel,     only : Node, Nodes
  use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
  use mpi_siesta
  use m_fermid,     only : fermid, stepf
  use alloc,        only : re_alloc, de_alloc
  use t_spin, only: tSpin
  use intrinsic_missing, only: modp

  implicit none

  integer :: MPIerror

  type(tSpin), intent(in) :: spin

  ! Pass system size information
  ! Local, unit-cell, supercell
  integer, intent(in) :: no_l, no_u, no_s

  ! K-point information
  integer, intent(in) :: nk
  real(dp), intent(in) :: kpoint(3,nk), wk(nk)

  ! Now pass sparse patterns
  integer, intent(in) :: nnz
  integer, intent(in) :: ncol(no_l), ptr(no_l), col(nnz)

  ! Matrices
  real(dp), intent(in) :: H(nnz,spin%H), S(nnz), xij(3,nnz)
  real(dp), intent(inout) :: DM(nnz,spin%DM), EDM(nnz,spin%EDM)

  ! Sought charges per spin
  real(dp), intent(in) :: Qtot
  ! Fermi-levels
  real(dp), intent(inout) :: Ef
  ! Energy range of DM
  real(dp), intent(in) :: E1, E2
  
  ! Eigenvalues and occupations (charges)
  real(dp), intent(inout) :: qo(no_u*spin%spinor,nk), eo(no_u*spin%spinor,nk)
  real(dp), intent(in) :: Occtol, Temp
  ! Calculated quantities
  real(dp), intent(inout) :: Entropy

  ! Control whether DM/EDM is calculated
  logical, intent(in) :: getD

  ! Current SCF step
  integer, intent(in) :: iSCF

  ! Number of calculated eigenstates
  integer, intent(in) :: neigwanted

  ! Auxiliary arrays
  complex(dp), intent(inout), target :: Hk(2,no_u,2,no_u), Sk(2,no_u,2,no_u)
  complex(dp), intent(inout) :: psi(2,no_u,no_u*2)
  real(dp), intent(inout) :: aux(2*no_u)

  ! Internal variables
  integer :: BNode, ie, ierror, ik, is
  integer :: io, iio, jo, ind, neigneeded
  integer :: no_u2, neigwanted2
  real(dp) :: kxij, t
  complex(dp) :: kph, D11, D22, D12, D21, cp

  ! Globalized matrices
  integer :: g_nnz
  integer, pointer :: g_ncol(:), g_ptr(:), g_col(:)
  real(dp), pointer :: g_H(:,:), g_S(:), g_xij(:,:)
  real(dp), pointer :: g_DM(:,:), g_EDM(:,:)

  ! Not allocated, only pointers
  complex(dp), pointer :: Dk(:,:,:,:), Ek(:,:,:,:)

#ifdef DEBUG
  call write_debug( '    PRE diag3kp' )
#endif

  no_u2 = no_u * 2
  neigwanted2 = neigwanted * 2

  ! Globalise sparsity pattern
  call MPI_AllReduce(nnz, g_nnz, 1, MPI_Integer, MPI_Sum, &
      MPI_Comm_World, MPIerror)

  ! Nullify arrays
  nullify(g_ncol, g_ptr, g_col, g_H, g_S, g_xij, g_DM, g_EDM)
  
  ! Allocate local memory for global list arrays
  call re_alloc( g_ncol, 1, no_u, name='g_ncol', routine= 'diag3kp' )
  call re_alloc( g_ptr, 1, no_u, name='g_ptr', routine= 'diag3kp' )
  call re_alloc( g_col, 1, g_nnz, name='g_col', routine= 'diag3kp' )
  call re_alloc( g_H, 1, g_nnz, 1, spin%H, name='g_H', routine= 'diag3kp' )
  call re_alloc( g_S, 1, g_nnz, name='g_S', routine= 'diag3kp' )
  call re_alloc( g_xij, 1, 3, 1, g_nnz, name='g_xij', routine= 'diag3kp' )

  ! Create pointers for d* arrays
  ! Note that these arrays may then not be inter-used
  Dk => Hk
  Ek => Sk

  ! Globalize arrays
  g_ptr(1) = 0
  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)
    
    if ( Node == BNode ) then
      call GlobalToLocalOrb(io,Node,Nodes,iio)
      g_ncol(io) = ncol(iio)
      do jo = 1, ncol(iio)
        g_col(g_ptr(io)+jo) = col(ptr(iio)+jo)
        g_H(g_ptr(io)+jo,:) = H(ptr(iio)+jo,:)
        g_S(g_ptr(io)+jo) = S(ptr(iio)+jo)
        g_xij(:,g_ptr(io)+jo) = xij(:,ptr(iio)+jo)
      end do
    end if
    
    call MPI_Bcast(g_ncol(io),1,MPI_Integer, BNode, &
        MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_col(g_ptr(io)+1),g_ncol(io),MPI_Integer, &
        BNode,MPI_Comm_World,MPIerror)
    do is = 1, spin%H
      call MPI_Bcast(g_H(g_ptr(io)+1,is),g_ncol(io),MPI_Double_Precision, &
          BNode,MPI_Comm_World,MPIerror)
    end do
    call MPI_Bcast(g_S(g_ptr(io)+1),g_ncol(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
    call MPI_Bcast(g_xij(1,g_ptr(io)+1),3*g_ncol(io),MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)

    ! Update list-pointer
    if ( io < no_u ) g_ptr(io+1) = g_ptr(io) + g_ncol(io)
    
  end do

  ! Perform diagonalization loop
  do ik = 1 + Node, nk, Nodes
    
    call setup_k(kpoint(:,ik))

    ! Since we want to calculate the velocities as well we do need the eigenstates
    call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eo(1,ik),psi, &
        -neigwanted2,iscf,ierror, -1)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    else if ( ierror < 0 ) then
      call setup_k(kpoint(:,ik))
      call cdiag(Hk,Sk,no_u2,no_u2,no_u2,eo(1,ik),psi, &
          -neigwanted2,iscf,ierror, -1)
    end if

  end do

  ! Globalise eigenvalues
  do ik = 1, nk
    BNode = mod(ik-1, Nodes)
    call MPI_Bcast(eo(1,ik), neigwanted2, MPI_Double_Precision, &
        BNode,MPI_Comm_World,MPIerror)
  end do

  ! Check if we are done ................................................
  if ( .not. getD ) goto 999

  ! Find new Fermi energy and occupation weights ........................
  call fermid(2, 1, nk, wk, no_u2, &
      neigwanted2, eo, Temp, qtot, qo, ef, Entropy )

  ! Allocate globalized DM and EDM
  call re_alloc( g_DM, 1, g_nnz, 1, spin%DM, name='g_DM', routine= 'diag3kp' )
  call re_alloc( g_EDM, 1, g_nnz, 1, spin%EDM, name='g_EDM', routine= 'diag3kp' )


  ! Find weights for local density of states ............................
  if ( e1 < e2 ) then
    
    t = max( temp, 1.d-6 )
!$OMP parallel do default(shared), private(ik,io), firstprivate(t)
    do ik = 1,nk
      do io = 1, neigwanted2
        qo(io,ik) = wk(ik) * &
            ( stepf((eo(io,ik)-e2)/t) - stepf((eo(io,ik)-e1)/t) )
      end do
    end do
!$OMP end parallel do

  end if

  ! Initialize to 0
  g_DM(:,:) = 0._dp
  g_EDM(:,:) = 0._dp

  do ik = 1 + Node, nk, Nodes

    ! Find maximum eigenvector that is required for this k point and spin
    ! Note that since eo has averaged out the degeneracy eigenvalues this below
    ! block will also group *all* degenerate eigenvalues!
    neigneeded = 1
    do ie = neigwanted2, 1, -1
      if ( abs(qo(ie,ik)) > occtol ) then
        neigneeded = ie
        exit
      end if
    end do

    ! Find eigenvectors
    call setup_k(kpoint(:,ik))
    call cdiag(Hk,Sk,no_u2,no_u2,no_u2,aux,psi,neigneeded,iscf,ierror, -1)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    else if ( ierror < 0 ) then
      call setup_k(kpoint(:,ik))
      call cdiag(Hk,Sk,no_u2,no_u2,no_u2,aux,psi,neigneeded,iscf,ierror, -1)
    end if

    ! Expand the eigenvectors to the density matrix

    Dk = cmplx(0._dp, 0._dp, dp)
    Ek = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared), &
!$OMP&private(ie,io,jo,ind), &
!$OMP&private(kxij,kph,D11,D22,D12,D21,cp)

    ! Add contribution to density matrices of unit-cell orbitals
    ! Global operation to form new density matrix
    do ie = 1, neigneeded
        
!$OMP do
      do io = 1, no_u
        D11 = qo(ie,ik) * psi(1,io,ie)
        D22 = qo(ie,ik) * psi(2,io,ie)
        D12 = D11 * aux(ie)
        D21 = D22 * aux(ie)
        do jo = 1, no_u
          
          cp = conjg(psi(1,jo,ie))
          Dk(1,jo,1,io) = Dk(1,jo,1,io) + D11 * cp
          Ek(1,jo,1,io) = Ek(1,jo,1,io) + D12 * cp          
          Dk(2,jo,1,io) = Dk(2,jo,1,io) + D22 * cp
          Ek(2,jo,1,io) = Ek(2,jo,1,io) + D21 * cp

          cp = conjg(psi(2,jo,ie))
          Dk(1,jo,2,io) = Dk(1,jo,2,io) + D11 * cp
          Ek(1,jo,2,io) = Ek(1,jo,2,io) + D12 * cp
          Dk(2,jo,2,io) = Dk(2,jo,2,io) + D22 * cp
          Ek(2,jo,2,io) = Ek(2,jo,2,io) + D21 * cp

        end do
      end do
!$OMP end do

    end do
      
!$OMP do
    do io = 1, no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = kpoint(1,ik) * g_xij(1,ind) + kpoint(2,ik) * g_xij(2,ind) + kpoint(3,ik) * g_xij(3,ind)
        kph = exp(cmplx(0._dp, - kxij, dp))

        D11 = Dk(1,jo,1,io) * kph
        D22 = Dk(2,jo,2,io) * kph
        D12 = Dk(1,jo,2,io) * kph
        D21 = Dk(2,jo,1,io) * kph
          
        g_DM(ind,1) = g_DM(ind,1) + real(D11, dp)
        g_DM(ind,2) = g_DM(ind,2) + real(D22, dp)
        g_DM(ind,3) = g_DM(ind,3) + real(D12, dp)
        g_DM(ind,4) = g_DM(ind,4) - aimag(D12)
        g_DM(ind,5) = g_DM(ind,5) + aimag(D11)
        g_DM(ind,6) = g_DM(ind,6) + aimag(D22)
        g_DM(ind,7) = g_DM(ind,7) + real(D21, dp)
        g_DM(ind,8) = g_DM(ind,8) + aimag(D21)

        D11 = Ek(1,jo,1,io) * kph
        D22 = Ek(2,jo,2,io) * kph
        D12 = Ek(1,jo,2,io) * kph
        D21 = Ek(2,jo,1,io) * kph

        g_EDM(ind,1) = g_EDM(ind,1) + real(D11, dp)
        g_EDM(ind,2) = g_EDM(ind,2) + real(D22, dp)
        g_EDM(ind,3) = g_EDM(ind,3) + real(D12, dp)
        g_EDM(ind,4) = g_EDM(ind,4) - aimag(D12)
          
      end do
    end do
!$OMP end do nowait
    
!$OMP end parallel

  end do

  ! Now do bcast of the respective arrays

  do io = 1, no_u
    
    call WhichNodeOrb(io,Nodes,BNode)
    call GlobalToLocalOrb(io,BNode,Nodes,iio)

    if ( BNode == Node ) then
      
      do is = 1, spin%DM
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            DM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      do is = 1, spin%EDM
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            EDM(ptr(iio)+1,is),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do

    else

      ! This is because this is *NOT* the receiving node
      ! and hence the recv_buff is not important
      do is = 1, spin%DM
        call MPI_Reduce(g_DM(g_ptr(io)+1,is), &
            aux(1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      do is = 1, spin%EDM
        call MPI_Reduce(g_EDM(g_ptr(io)+1,is), &
            aux(1),g_ncol(io),MPI_double_precision, &
            MPI_sum,BNode,MPI_Comm_World,MPIerror)
      end do
      
    end if
  end do

  call de_alloc( g_DM, name='g_DM', routine= 'diag3kp' )
  call de_alloc( g_EDM, name='g_EDM', routine= 'diag3kp' )

  ! Exit point 
999 continue

  ! Clean up memory
  call de_alloc( g_ncol, name='g_ncol', routine= 'diag3kp' )
  call de_alloc( g_ptr, name='g_ptr', routine= 'diag3kp' )
  call de_alloc( g_col, name='g_col', routine= 'diag3kp' )
  call de_alloc( g_H, name='g_H', routine= 'diag3kp' )
  call de_alloc( g_S, name='g_S', routine= 'diag3kp' )
  call de_alloc( g_xij, name='g_xij', routine= 'diag3kp' )

#ifdef DEBUG
  call write_debug( '    POS diag3kp' )
#endif

contains
  
  subroutine setup_k(k)
    real(dp), intent(in) :: k(3)

    Sk = cmplx(0._dp, 0._dp, dp)
    Hk = cmplx(0._dp, 0._dp, dp)

!$OMP parallel do default(shared), private(io,jo,ind,kxij,kph)
    do io = 1,no_u
      do ind = g_ptr(io) + 1, g_ptr(io) + g_ncol(io)
        jo = modp(g_col(ind), no_u)
        kxij = k(1) * g_xij(1,ind) + k(2) * g_xij(2,ind) + k(3) * g_xij(3,ind)

        ! Calculate the complex phase
        kph = exp(cmplx(0._dp, - kxij, dp))
        
        Sk(1,jo,1,io) = Sk(1,jo,1,io) + g_S(ind) * kph
        Sk(2,jo,2,io) = Sk(2,jo,2,io) + g_S(ind) * kph
        Hk(1,jo,1,io) = Hk(1,jo,1,io) + cmplx(g_H(ind,1), g_H(ind,5), dp) * kph
        Hk(2,jo,2,io) = Hk(2,jo,2,io) + cmplx(g_H(ind,2), g_H(ind,6), dp) * kph
        Hk(1,jo,2,io) = Hk(1,jo,2,io) + cmplx(g_H(ind,3), -g_H(ind,4), dp) * kph
        Hk(2,jo,1,io) = Hk(2,jo,1,io) + cmplx(g_H(ind,7), g_H(ind,8), dp) * kph

      end do
    end do
!$OMP end parallel do

  end subroutine setup_k
  
end subroutine diag3kp


#else

subroutine diag3kp_dummy()
end subroutine diag3kp_dummy

#endif
