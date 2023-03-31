title: Sparse Matrices in SIESTA

Sparse matrices are widespread in Siesta due to the use of a basis set
of finite-range atomic orbitals. Matrix elements between orbitals
which are "too far" from each other are zero, leading to sparse
matrices for all but the smallest systems.

The following discussion is focused on the global data found in [[sparse_matrices]].

Take the Hamiltonian. If `no_u` is the number of orbitals, its matrix
elements \(H_{\mu\nu}\) can be simply arranged in a square matrix with
`no_u` rows and columns.

In the presence of a sizable number of zeros, it pays to represent
this matrix in a more compact form.

The sparse matrix is in CSR/CSC matrix format with some slight
differences from the elsewhere found CSR format. Because Siesta
(primarily) deals with Hermitian matrices, CSR and CSC are equivalent.
In Siesta the CSR format is composed of the following (integer) components.

*  `maxnh`: 
     number of total non-zero elements (sum of numh array)
*  `numh(no_u)`:
    number of non-zero elements each orbital connects too
    @note
      This is a specific to Siesta, since many
      sparse implementations of CSR matrices in Fortran use
      only the `listhptr` array with one more element.
     @endnote
*  `listhptr(no_u)`
     index pointer to listh such that `listh(listhptr(1) + 1)`
     is the first non-zero element of orbital 1.
     `listh(listhptr(io) + 1)` is thus the first non-zero element
     of orbital 'io' while `listh(listhptr(io) + numh(io))` is the
     last non-zero element of orbital 'io'.
     @note
      `listhptr` is 0-based (`listhptr(1) == 0`, ALWAYS)
      while many implementations in Fortran would have used listhptr(1) == 1
      this is not so in Siesta. Secondly, in other typical CSR implementations
	  this array would be of size `no_u + 1` with `listhptr(no_u+1) == maxnh + 1`.
     @endnote
*  `listh(maxnh)`:
     the column indices for the non-zero elements.
     `listh(listhptr(io)+1)`
     would correspond to the element `M(io, listh(listhptr(io)+1))`

Note that these matrices are 1D block cyclic distributed and the
routines found in [[class_OrbitalDistribution]] or module [[parallelsubs]]
are required to translate local rows to global rows (columns
are not distributed).

The discussion above assumes that only a single MPI node is used
such that `no_l == no_u`, where the former is the number of distributed
rows on the MPI node.

See the discussion of the [auxiliary supercell](|page|/implementation/1-auxiliary-supercell.html) for the modifications
introduced by the support for periodic systems with k-point sampling.

The sparse pattern information is handled by the program in two forms:

* Through a set of arrays, as mentioned above.
* Through a [reference-counted](|page|/datastructures/1-buds.html) derived type [[class_Sparsity:sparsity]] and its instantiation in ([[sparse_matrices:sparse_pattern]]).
  The old arrays are *actually* pointing to the arrays in the variable `sparse_pattern`.


@note
 The suffix 'h' (`numh`, `maxh`, ...) is a legacy name and should be avoided
 in subsequent usages (i.e. when new routines are being passed the above arrays
 or the sparse_pattern.
 Currently Siesta does not distinguish between H/DM/S/EDM/... sparsity patterns
 and thus the sparse pattern is not specific to any physical array, it is a global
 definition of all sparse matrices used in Siesta.
@endnote
