! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine diag2kspiral( nuo, no, maxnh, maxnd, maxo, &
    numh, listhptr, listh, numd, listdptr, &
    listd, H, S, getD, qtot, temp, e1, e2, &
    xij, indxuo, nk, kpoint, wk, &
    eo, qo, Dnew, Enew, ef, Entropy, qw, &
    Haux, Saux, psi, Dk, Ek, aux, nuotot, &
    occtol, iscf, neigwanted)
! *********************************************************************
! Calculates the eigenvalues and eigenvectors, density
! and energy-density matrices, and occupation weights of each 
! eigenvector, for given Hamiltonian and Overlap matrices.
! This version is for non-colinear spin with k-sampling and spiral
! arrangement of spins.
! Written by V. M. Garcia-Suarez. June 2002
! Modified to reduce eigenvector computation by J.D. Gale Nov. 2004
! Lots of bugs fixed and changed to complex by Nick Papior, June 2020
! **************************** INPUT **********************************
! integer nuo                 : Number of basis orbitals in unit cell
! integer no                  : Number of basis orbitals in supercell
! integer maxnh               : Maximum number of orbitals interacting  
! integer maxnd               : First dimension of listd / DM
! integer maxo                : First dimension of eo and qo
! integer numh(nuo)           : Number of nonzero elements of each row 
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to each row (-1) of the
!                               hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
!                               column indexes for each matrix row
! integer numd(nuo)           : Number of nonzero elements of each row 
!                               of density matrix
! integer listdptr(nuo)       : Pointer to each row (-1) of the
!                               density matrix
! integer listd(maxnd)        : Nonzero density-matrix element column 
!                               indexes for each matrix row
! real*8  H(maxnh,4)          : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! logical getD                : Find occupations and density matrices?
! real*8  qtot                : Number of electrons in unit cell
! real*8  temp                : Electronic temperature 
! real*8  e1, e2              : Energy range for density-matrix states
!                               (to find local density of states)
!                               Not used if e1 > e2
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not used if only gamma point)
! integer indxuo(no)          : Index of equivalent orbital in unit cell
!                               Unit cell orbitals must be the first in
!                               orbital lists, i.e. indxuo.le.nuo, with
!                               nuo the number of orbitals in unit cell
! real*8 qw(3)                : Wave vector for spiral configuration
! integer nk                  : Number of k points
! real*8  kpoint(3,nk)        : k point vectors
! real*8  wk(nk)              : k point weights (must sum one)
! integer nuotot              : total number of orbitals per unit cell
!                               over all processors
! real*8  occtol              : Occupancy threshold for DM build
! integer iscf                : SCF cycle number
! integer neigwanted          : Number of eigenvalues wanted
! *************************** OUTPUT **********************************
! real*8 eo(maxo*2,nk)        : Eigenvalues
! real*8 qo(maxo*2,nk)        : Occupations of eigenstates
! real*8 Dnew(maxnd,4)        : Output Density Matrix
! real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
! real*8 ef                   : Fermi energy
! real*8 Entropy              : Electronic entropy
! *************************** AUXILIARY *******************************
! complex*16 Haux(2,nuotot,2,nuo) : Aux. space for the hamiltonian matrix
! complex*16 Saux(2,nuotot,2,nuo) : Aux. space for the overlap matrix
! complex*16 psi(2,nuotot,2*nuo)  : Aux. space for the eigenvectors
! complex*16 aux(2,nuotot)        : Extra auxiliary space
! complex*16 Dk(2,nuotot,2,nuo)   : Aux. space that may be the same as Haux
! complex*16 Ek(2,nuotot,2,nuo)   : Aux. space that may be the same as Saux
! *************************** UNITS ***********************************
! xij and kpoint must be in reciprocal coordinates of each other.
! temp and H must be in the same energy units.
! eo, Enew and ef returned in the units of H.
! *************************** PARALLEL ********************************
! The auxiliary arrays are now no longer symmetry and so the order
! of referencing has been changed in several places to reflect this.
! *********************************************************************
!
!  Modules
!
  use precision
  use sys
  use parallel,     only : Node, Nodes, BlockSize
  use parallelsubs, only : LocalToGlobalOrb
  use m_fermid,     only : fermid, stepf
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

#ifdef MPI
  integer :: MPIerror
#endif

  integer :: maxnd, maxnh, maxo, nk, no, nuo, nuotot, iscf, neigwanted

  integer :: indxuo(no), listh(maxnh), numh(nuo), listd(maxnd), numd(nuo), &
      listhptr(nuo), listdptr(nuo)

  real(dp) :: Dnew(maxnd,4), &
      e1, e2, ef, Enew(maxnd,4), Entropy, eo(maxo*2,nk), &
      H(maxnh,4), kpoint(3,nk), qo(maxo*2,nk), qtot, &
      S(maxnh), temp, wk(nk), xij(3,maxnh), qw(3), occtol

  complex(dp), dimension(2,nuotot,2,nuo) :: Haux, Saux, Dk, Ek
  complex(dp), dimension(2,nuotot,2*nuo) :: psi
  complex(dp), dimension(2,nuotot) :: aux

  logical :: getD

  external :: cdiag

!  Internal variables .............................................
  integer :: BNode, BTest, ie, ierror, iie, ik, ind, io, iio, &
      iuo, j, jo, juo, neigneeded
  integer :: nuotot2, nuo2, BlockSize2, maxo2, neigwanted2
  complex(dp) :: kphs, qphs
  complex(dp) :: cicj, D11, D22, D12, D21
  real(dp) ::  kxij, qxij, ee, qe, t, qwh(3)

#ifdef DEBUG
  call write_debug( '    PRE diag2kspiral' )
#endif

  ! Easier
  nuotot2 = nuotot * 2
  neigwanted2 = neigwanted * 2
  nuo2 = nuo * 2
  maxo2 = maxo * 2
  BlockSize2 = BlockSize * 2
  qwh(:) = qw(:) * 0.5_dp

! Find eigenvalues at every k point ...............................
  do ik = 1,nk

    call setup_k()

    ! Find eigenvalues
    ! Possible memory optimization: equivalence Haux and psi
    call cdiag(Haux,Saux,nuotot2,nuotot2,nuo2,eo(1,ik),psi, &
        -neigwanted2,iscf,ierror,BlockSize2)
    if ( ierror /= 0 ) then
      call die('Terminating due to failed diagonalisation')
    end if
  end do

  ! Check if we are done ................................................
  if (.not.getD) then
#ifdef DEBUG
    call write_debug( '    POS diag2kspiral' )
#endif
    return
  end if

  ! Find new Fermi energy and occupation weights ........................
  call fermid(2, 1, nk, wk, maxo2, neigwanted2, eo, &
      temp, qtot, qo, ef, Entropy )

  ! Find weights for local density of states ............................
  if ( e1 < e2 ) then
    t = max( temp, 1.d-6 )
    do ik = 1,nk
      do io = 1, neigwanted2
        qo(io,ik) = wk(ik) * &
            ( stepf( (eo(io,ik)-e2)/t ) - stepf( (eo(io,ik)-e1)/t ) ) 
      end do
    end do
  end if

! New density and energy-density matrices of unit-cell orbitals .......
  Dnew(:,:) = 0._dp
  Enew(:,:) = 0._dp

  do ik = 1, nk

    ! Find maximum eigenvector that is required for this k point
    neigneeded = 0
    ie = neigwanted2
    do while ( ie > 0 .and. neigneeded == 0 )
      qe = qo(ie,ik)
      if ( abs(qe) > occtol ) neigneeded = ie
      ie = ie - 1
    end do

    call setup_k()

    call cdiag(Haux,Saux,nuotot2,nuotot2,nuo2,eo(1,ik),psi, &
        neigneeded,iscf,ierror,BlockSize2)

    ! Check error flag and take appropriate action
    if ( ierror > 0 ) then
      call die('Terminating due to failed diagonalisation')
    else if ( ierror < 0 ) then

      call setup_k()
      call cdiag(Haux,Saux,nuotot2,nuotot2,nuo2,eo(1,ik),psi, &
          neigneeded,iscf,ierror,BlockSize2)
 
    end if
    
! Store the products of eigenvectors in matrices Dk and Ek
! WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
    Dk = cmplx(0.0_dp, 0.0_dp, dp)
    Ek = cmplx(0.0_dp, 0.0_dp, dp)

    BNode = 0
    iie = 0
    do ie = 1, neigneeded
      if ( Node == BNode ) then
        iie = iie + 1
        do j = 1, nuotot
          aux(1,j) = psi(1,j,iie) ! c_i,up
          aux(2,j) = psi(2,j,iie) ! c_i,dn
        end do
      end if
#ifdef MPI
      call MPI_Bcast(aux(1,1),nuotot2,MPI_double_complex,BNode, &
          MPI_Comm_World,MPIerror)
#endif
      qe = qo(ie,ik)
      ee = qe * eo(ie,ik)
      do iuo = 1,nuo
        call LocalToGlobalOrb(iuo,Node,Nodes,iio)
        do juo = 1, nuotot
!------- 1,1 -----------------------------------------------------------
          cicj = conjg(aux(1,juo)) * aux(1,iio)
          Dk(1,juo,1,iuo) = Dk(1,juo,1,iuo) + qe * cicj
          Ek(1,juo,1,iuo) = Ek(1,juo,1,iuo) + ee * cicj
!------- 2,2 -----------------------------------------------------------
          cicj = conjg(aux(2,juo)) * aux(2,iio)
          Dk(2,juo,2,iuo) = Dk(2,juo,2,iuo) + qe * cicj
          Ek(2,juo,2,iuo) = Ek(2,juo,2,iuo) + ee * cicj
!------- 1,2 -----------------------------------------------------------
          cicj = conjg(aux(2,juo)) * aux(1,iio)
          Dk(1,juo,2,iuo) = Dk(1,juo,2,iuo) + qe * cicj
          Ek(1,juo,2,iuo) = Ek(1,juo,2,iuo) + ee * cicj
!------- 2,1 -----------------------------------------------------------
          cicj = conjg(aux(1,juo)) * aux(2,iio)
          Dk(2,juo,1,iuo) = Dk(2,juo,1,iuo) + qe * cicj
          Ek(2,juo,1,iuo) = Ek(2,juo,1,iuo) + ee * cicj
        end do
      end do
      BTest = ie/BlockSize2
      if ( BTest*BlockSize2 == ie ) then
        BNode = BNode + 1
        if ( BNode >= Nodes ) BNode = 0
      end if
    end do

! Add contribution to density matrices of unit-cell orbitals
    do iuo = 1,nuo
      do j = 1,numd(iuo)
        ind = listdptr(iuo) + j
        jo = listd(ind)
        juo = indxuo(jo)
        kxij = kpoint(1,ik) * xij(1,ind) + &
            kpoint(2,ik) * xij(2,ind) + &
            kpoint(3,ik) * xij(3,ind)
        kphs = exp(cmplx(0.0_dp, - kxij, dp))
        qxij = qwh(1) * xij(1,ind) + &
            qwh(2) * xij(2,ind) + &
            qwh(3) * xij(3,ind)
        qphs = exp(cmplx(0.0_dp, - qxij, dp))

        D11 = Dk(1,juo,1,iuo) * kphs * qphs
        D21 = Dk(2,juo,1,iuo) * kphs * conjg(qphs)
        D12 = Dk(1,juo,2,iuo) * kphs * qphs
        D22 = Dk(2,juo,2,iuo) * kphs * conjg(qphs)

        ! Make DM Spin-box hermitian 
        D12 = 0.5_dp * (D12 + conjg(D21))

        Dnew(ind,1) = Dnew(ind,1) + real(D11, dp)
        Dnew(ind,2) = Dnew(ind,2) + real(D22, dp)
        Dnew(ind,3) = Dnew(ind,3) + real(D12, dp)
        Dnew(ind,4) = Dnew(ind,4) - aimag(D12)

        D11 = Ek(1,juo,1,iuo) * kphs * qphs
        D21 = Ek(2,juo,1,iuo) * kphs * conjg(qphs)
        D12 = Ek(1,juo,2,iuo) * kphs * qphs
        D22 = Ek(2,juo,2,iuo) * kphs * conjg(qphs)

        ! Make EDM Spin-box hermitian 
        D12 = 0.5_dp * (D12 + conjg(D21))

        Enew(ind,1) = Enew(ind,1) + real(D11, dp)
        Enew(ind,2) = Enew(ind,2) + real(D22, dp)
        Enew(ind,3) = Enew(ind,3) + real(D12, dp)
        Enew(ind,4) = Enew(ind,4) - aimag(D12)

      end do
    end do
    
  end do

#ifdef DEBUG
  call write_debug( '    POS diag2ksipral' )
#endif

contains

  subroutine setup_k()

    ! Initialize Hamiltonian and overlap matrices in full format
    ! Index i is for real/imag parts
    ! Indices is and js are for spin components
    ! Indices iuo and juo are for orbital components:
    ! Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
    Saux = cmplx(0.0_dp, 0.0_dp, dp)
    Haux = cmplx(0.0_dp, 0.0_dp, dp)

    ! Transfer S,H matrices from sparse format in supercell to
    ! full format in unit cell
    ! Convention: ispin=1 => H11, ispin=2 => H22, 
    !             ispin=3 => Real(H12), ispin=4 => Imag(H12)
    do iuo = 1,nuo
      do j = 1,numh(iuo)
        ind = listhptr(iuo) + j
        jo = listh(ind)
        juo = indxuo(jo)
        kxij = kpoint(1,ik) * xij(1,ind) + &
            kpoint(2,ik) * xij(2,ind) + &
            kpoint(3,ik) * xij(3,ind)
        kphs = exp(cmplx(0.0_dp, - kxij, dp))
        qxij = qwh(1) * xij(1,ind) + &
            qwh(2) * xij(2,ind) + &
            qwh(3) * xij(3,ind)
        qphs = exp(cmplx(0.0_dp, - qxij, dp))
        D11 = kphs * qphs
        D22 = kphs * conjg(qphs)

        Saux(1,juo,1,iuo) = Saux(1,juo,1,iuo) + S(ind) * D11
        Saux(2,juo,2,iuo) = Saux(2,juo,2,iuo) + S(ind) * D22

        Haux(1,juo,1,iuo) = Haux(1,juo,1,iuo) + H(ind,1) * D11
        Haux(2,juo,1,iuo) = Haux(2,juo,1,iuo) + cmplx(H(ind,3), H(ind,4),dp) * D22
        Haux(1,juo,2,iuo) = Haux(1,juo,2,iuo) + cmplx(H(ind,3), -H(ind,4),dp) * D11
        Haux(2,juo,2,iuo) = Haux(2,juo,2,iuo) + H(ind,2) * D22

      end do
    end do

  end subroutine setup_k

end subroutine diag2kspiral
