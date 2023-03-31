! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!> See [this page](|page|/datastructures/2-sparse.html)
!> for background information on sparsity in SIESTA
module sparse_matrices

  use precision
  use class_dSpData1D
  use class_dSpData2D
  use class_Sparsity
  use class_OrbitalDistribution
  use class_Fstack_Pair_Geometry_dSpData2D

  implicit none
  
  private
  save

  ! The sparse matrix in Siesta is an intrinsic part of the
  ! entire code.
  !
  ! Here is a brief discussion of how it is handled.
  !
  ! The sparse matrix is in CSR/CSC matrix format with some slight
  ! differences from the elsewhere found CSR format.
  ! Because Siesta (primarily) deals with Hermitian matrices, CSR and
  ! CSC are equivalent.
  ! In Siesta the CSR format is composed of the following components.
  ! Note that these matrices are 1D block cyclic distributed and the
  ! routines found in class_OrbitalDistribution or module parallelsubs
  ! are required to translate local rows to global rows (columns
  ! are not distributed).
  ! The discussion below infers that only a single MPI node is used
  ! such that no_l == no_u, where the former is the number of distributed
  ! rows on the MPI node.
  !   integer maxnh
  !      number of total non-zero elements (sum of numh array)
  !   integer numh(no_u)
  !      number of non-zero elements each orbital connects too
  !      NOTES:
  !       This is a specific to Siesta only array, since many
  !       sparse implementations of CSR matrices in Fortran uses
  !       only the listhptr array with one more element.
  !   integer listhptr(no_u)
  !      index pointer to listh such that listh(listhptr(1) + 1)
  !      is the first non-zero element of orbital 1.
  !      listh(listhptr(io) + 1) is thus the first non-zero element
  !      of orbital 'io' while listh(listhptr(io) + numh(io)) is the
  !      last non-zero element of orbital 'io'.
  !      NOTES:
  !       listhptr is 0-based (listhptr(1) == 0, ALWAYS)
  !       while many implementations in Fortran would have used listhptr(1) == 1
  !       this is not so in Siesta.
  !   integer listh(maxnh)
  !      the column indices for the non-zero elements.
  !      listh(listhptr(io)+1)
  !      would correspond to the element M(io, listh(listhptr(io)+1))
  !      NOTES:
  !       See below for the exact handling of the listh elements when
  !       the auxiliary supercell is in effect.
  !
  ! The sparse pattern is contained in two forms:
  !   1. The old array form corresponding to the above mentioned
  !      arrays.
  !   2. type(Sparsity) hosted in sparse_pattern.
  !
  ! The old arrays are *actually* pointing to the arrays in the variable
  ! sparse_pattern.  
  !
  ! A specific handling of the sparse matrices in Siesta is that the column
  ! index is *also* an index for the supercell picture of the auxiliary supercell.
  ! The supercells in Siesta is governed by siesta_geom::nsc and siesta_geom::isc_off
  ! which describe, number of supercells along each lattice vector and the direct
  ! supercell offsets for each supercell index, i.e. isc_off(:, is) is the supercell
  ! offset for the 'is'th supercell.
  ! To retrieve the supercell index and the correct column index from the listh array,
  ! the following should be done:
  !    is = (listh(ind) - 1) / no_u + 1
  !    col = mod(listh(ind) - 1, no_u) + 1
  !
  ! The above destribes the sparse matrix format and thus to construct the
  ! Hamiltonian matrix one should do this simple loop:
  !  PLEASE CHECK IN DIAGK FOR SPECIFIC DETAILS!
  !  THIS IS A SAMPLE CODE NOT TO BE USED!
  !
  !  do io = 1, no_u
  !    do ind = listhptr(io) + 1, lishptr(io) + numh(io)
  !      is = (listh(ind) - 1) / no_u + 1
  !      col = mod(listh(ind) - 1, no_u) + 1
  !      H(io, col) = H(io, col) + H_sparse(ind) * exp(-1j * sum(R(:, is) * k(:)))
  !    end do
  !  end do
  !
  !
  ! NOTES:
  !  The suffix 'h' (numh, maxh, ...) is a legacy name and should be avoided
  !  in subsequent usages (i.e. when new routines are being passed the above arrays
  !  or the sparse_pattern.
  !  Currently Siesta does not distinguish between H/DM/S/EDM/... sparsity patterns
  !  and thus the sparse pattern is not specific to any array, but is more a global
  !  definition of all sparse matrices used in Siesta.

  !> Global sparse pattern used for H/DM/EDM/S matrices, see [sparse pattern](|page|/datastructures/2-sparse.html) for details
  type(Sparsity), public :: sparse_pattern
  
  !> Number of non-zero elements in the sparse matrix (local to MPI Node)
  integer, public :: maxnh = 0
  !> Column indices in the CSR matrix (local to MPI Node)
  integer, public, pointer :: listh(:) => null()
  !> Index pointer to `listh` for each row in the CSR matrix, 0-based (local to MPI Node)
  integer, public, pointer :: listhptr(:) => null()
  !> Number of non-zero elements per row in the CSR matrix (local to MPI Node)
  integer, public, pointer :: numh(:) => null()

  ! Density matrix for in and out
  real(dp), public, pointer :: Dold(:,:) => null(), Dscf(:,:) => null()
  !> Current SCF step density matrix
  type(dSpData2D), public :: DM_2D
  ! Energy density matrix for in and out
  real(dp), public, pointer :: Eold(:,:) => null(), Escf(:,:) => null()
  !> Current SCF step energy density matrix
  type(dSpData2D), public :: EDM_2D
  ! Hamiltonian matrix for in and out
  real(dp), public, pointer :: Hold(:,:) => null(), H(:,:) => null()
  !> Current SCF step Hamiltonian matrix
  type(dSpData2D), public :: H_2D

  ! Overlap matrix (constant for complete SCF loop, changes per MD)
  real(dp), public, pointer :: S(:) => null()
  !> Current MD step overlap matrix
  type(dSpData1D), public :: S_1D
  
  ! Orbital distance matrix (constant for complete SCF loop, changes per MD)
  real(dp), public, pointer :: xijo(:,:) => null()
  !> Inter-orbital [vector](|page|/implementation/1-auxiliary-supercell.html)
  type(dSpData2D), public :: xij_2D
      
  ! Pieces of H that do not depend on the SCF density matrix
  ! Formerly there was a single array H0 for this
  type(dSpData1D), public :: H_vkb_1D, H_kin_1D
  ! LDA+U and spin-orbit coupling Hamiltonian
  type(dSpData2D), public :: H_dftu_2D, H_so_2D

  !> Geometry density matrix history, see [here](|page|/implementation/1-auxiliary-supercell.html)
  type(Fstack_Pair_Geometry_dSpData2D), public :: DM_history

  !> MPI block distribution of the orbitals in the [sparse matrices](|page|/datastructures/2-sparse.html)
  type(OrbitalDistribution), public :: block_dist
  !> Dummy local distribution
  type(OrbitalDistribution), public :: single_dist

  public :: resetSparseMatrices

contains

  subroutine resetSparseMatrices( )
    use alloc, only : de_alloc

    implicit none

    call delete( block_dist )
    call delete( sparse_pattern )
    nullify(numh,listhptr,listh)
    maxnh = 0

    call delete( H_kin_1D )
    call delete( H_vkb_1D )
    call delete( H_dftu_2D )
    call delete( H_so_2D )

    call delete( DM_2D )
    nullify(Dscf)
    call delete( EDM_2D )
    nullify(Escf)
    call delete( S_1D )
    nullify(S)
    call delete( H_2D )
    nullify(H)
    call delete( xij_2D )
    nullify(xijo)

    call delete( DM_history )

    ! If never allocated de_alloc returns immediately.
    ! Currently these matrices are not stored in the new
    ! sparse matrix format.
    call de_alloc( Hold, 'Hold', 'sparseMat' )
    call de_alloc( Dold, 'Dold', 'sparseMat' )
    call de_alloc( Eold, 'Eold', 'sparseMat' )

  end subroutine resetSparseMatrices

end module sparse_matrices
