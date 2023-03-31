title: Auxiliary Supercell

For periodic systems SIESTA uses an auxiliary supercell as an indexing
tool to represent all the orbital interactions present in the
Hamiltonian (and other matrices). These interactions can extend to
orbitals beyond the unit cell, so extra bookeeping is needed
beyond the simple \(H_{\mu\nu}\) notation used previously.

Consider the simplified system depicted in the figure below. The unit
cell contains two atoms (open and black circles) with a single orbital
each (orbitals are not marked to avoid clutter). The orbitals of the
atoms in the unit cell at the center of the figure overlap with those
of other atoms in the crystal whenever the distance between atoms is
less than the sum of the orbital ranges. In the figure these
interactions extend to the first-neighboring cells only. In this
two-dimensional example, the atoms in the unit cell interact with
those in nine cells, and this geometrical construct is fit to
represent the complete set of non-zero Hamiltonian matrix elements
\(H_{\mu\nu}({\mathbf R}) \), where \(\mu\) and \(\nu\) are orbital
indexes (ranging over the unit cell only) and \({\mathbf R}\) is a
relative vector from the unit cell. This vector could just be a cell
lattice vector, but in practice SIESTA uses the relative vector
between the atoms involved.

![Orbital Interactions](|media|/Interactions.png "Orbital interactions")
{: style="text-align: center" ; width="50%" }

There are multiple image interactions for the same pair, each
associated with a different \({\mathbf R}\).  The set of nine cells is the
auxiliary supercell. In this case it is a 3x3 repetition of the unit
cell, and the "3" can be seen as "2x1+1": the central cell plus one
neighboring cell in either side.  We can also represent the
interactions with a rectangular matrix:

![RectangularMatrix](|media|/RectangularMatrix.png "Interaction
Matrix")


which has a row for each unit-cell orbital, and as many columns as
orbitals in the nine cells involved in the interactions. Any matrix
element $$H_{\mu\nu}({\mathbf R}) = \langle\phi_\mu({\mathbf R=0})|H|\phi_\nu({\mathbf
R})\rangle$$ can be stored by itself in the appropriate slot.

When it it time to build the Hamiltonian matrix for a given k-point \({\mathbf k}\):

$$ H_{\mu\nu}({\mathbf k}) =
\sum_{\mathbf R} { H_{\mu\nu}({\mathbf R}) e^{i{\mathbf k}\cdot{\mathbf R}}} $$

every slot's contribution, with the appropriate phase, is folded back
and reduced into the left-most square matrix, as the arrows in the
figure indicate.

For periodic systems with large unit cells, k-point sampling is not
really necessary, and the phases involved are all 1 (formally only the
\(\Gamma\) point \({\mathbf k=0}\) is used). In this case the auxiliary supercell
is not strictly necessary: the matrix elements can be reduced on the
fly to the unit-cell square interaction matrix, and other operations
can be similarly folded automatically throughout. It is possible,
however, to request that an auxiliary supercell be used, since the
extra level of bookeeping can be useful for other purposes (e.g. for
COHP analysis of the orbital-pairs contributions to the energy).

In the program, the auxiliary supercell is handled in several key routines:

* The need for a supercell is assessed in [[siesta_init]], by checking whether k-points are
  going to be used anywhere in the program.
* The size of the supercell needed is stored in the [[siesta_geom:nsc(variable)]] array
* The supercell offsets for each supercell index are stored in the [[siesta_geom:isc_off(variable)]]
  array.
* Indexing arrays live in [[sparse_matrices]] and are initialized in [[state_init]]

When an auxiliary supercell is used, the "column" stored in the
`listh` array defined in the [sparsity](|page|/datastructures/2-sparse.html) page ranges
over the long dimension of the rectangular interaction matrix. The
"unit-cell" column to which the folding is done is recorded in the
array `indxuo`. Its contents, however, can be really computed on the
fly, since the column indexes in the rectangular matrix are just
juxtapositions of blocks of size `no_u`.

The following idiom can be used to go through the arrays:
```fortran
do io = 1, no_u
  do ind = listhptr(io) + 1, lishptr(io) + numh(io)
    is = (listh(ind) - 1) / no_u + 1
    col = mod(listh(ind) - 1, no_u) + 1
	! or equivalently (found in some parts of the code):
    !   col = indxuo(listh(ind))
    H(io, col) = H(io, col) + H_sparse(ind) * exp(-1j * sum(xijo(:, ind) * k(:)))
  end do
end do
```

## MD or geometry optimizations

The auxiliary supercell is a fixed size for any given atomic configuration due to the
constant distances between atoms and neighbouring cell atoms.
In MD or geometry optimization simulations where atoms are being displaced the auxiliary cell
may change. In the following MD refers to _any_ kind of geometry movement, be it actual MD or
geometry relaxations.

<!---
Consider two atoms currently at a distance \(R>R_\mathrm{max}\)
with \(R_\mathrm{max}\) being the sum of orbital range belonging to the atoms. In a following geometry
iteration the atoms are brought closer to each other such that \(R<R_\mathrm{max}\). This will lead
to an increase in the number of non-zero elements in the sparse matrix.
Concintly for the reverse action where the atoms are moved farther apart one finds a reduction in the
number of non-zero elements.
 --->

In the following we will use `nsc` to refer to the variable [[siesta_geom:nsc(variable)]].

Generically one may find that the following list of actions are carried out:

1. [[siesta_init]]: figure out initial number of supercells and call [[atomlist:superc]]
2. Begin MD iterations

	1. [[state_init]]: figure out number of supercells for this geometry configuration
	2. If any of the supercells has changed we need to copy the old information to the
       new supercell information
	3. Perform SCF, move atoms and go to 2.1.
  
3. End Siesta and analyse

The step 2.1 and 2.2 requries some further explanation. Siesta keeps a history of SCF converged
density matrices from previous MD steps. When 2.1 is initiated the new initial density matrix is
an extrapolation of these previous density matrices. Since all density matrices are stored using
the [sparsity matrix](|page|/datastructures/2-sparse.html) definitions we need to make a conversion
between different auxiliary supercells. The conversion (in-place) between two different auxiliary
supercells may be found in [[m_handle_sparse:correct_supercell_SpD]].
Additionally, since the number of non-zero elements may also change during the MD cycle one also
needs to remove/add new non-zero elements in the sparse matrix. This is done
in [[m_restruct_SpData2D:restruct_dSpData2D]].

