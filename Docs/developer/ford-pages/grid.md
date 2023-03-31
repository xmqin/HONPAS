title: The real-space grid

SIESTA uses a real-space grid:

* To compute and represent the charge density in real-space.
* To solve Poisson's equation to compute Hartree potentials from charge densities.
* To compute matrix elements of the part of the Hamiltonian which is not amenable to
the standard ``matel`` treatment, basically the Hartree, XC, and neutral-atom potentials.

All these operations are carried out in module [[m_dhscf(module)]].

This grid is similar in spirit to the ``density`` grid used in
plane-wave codes (*not* the ``wave-function`` grid). Hence its
associated cutoff is typically from 100 Ryd to a few-hundred Ryd.  But
since SIESTA's basis set is made of atomic-like orbitals, a number of
rather more technical details are needed to provide the right indexing and
mapping.

In no particular order:

* Some arrays are needed to determine which orbitals have a given grid
  point in their support (the sphere in which they are confined). This
  is needed, for example, to compute the charge density. They live in
  module [[meshphi]], and are initialized in [[meshsubs:PhiOnMesh]].

* In normal operation, the actual values of the orbitals on every grid
  point are stored in array [[meshphi:phi]] for faster execution of
  some operations. This array is kept in single precision by default,
  and still uses a significant fraction of the total memory of the
  program. Orbital values can instead be computed on-the-fly if the
  ``Direct-Phi`` option is used.

* SIESTA uses a two-level indexing of grid points for some operations:
  ``big points`` and ``small points``. By default, big points are
  those for a mesh of half the density along every direction (that is,
  every other point of the grid). This means that there are eight
  ``small points`` for every big point. The arrays mentioned above are
  structured by keeping indexes relative to big points, and keeping
  for each a list for all the small-points. This was done originally
  to save space in indexing, but it complicates the code.  The number
  of big points is controlled by the parameter ``nsm`` (2 by
  default). Grid arrays might be in ``sequential`` mode, ordered in
  the natural Fortran sequence, considering all the small points, or
  in ``cluster`` mode, with lists of small-points ordered by their
  associated big-point. Conversion between the two is carried out by
  routine [[reord]].

## Parallelization issues

### Parallel distributions

When working in parallel, the distribution of mesh points among
processors is not unique during program execution, as different
operations have different load-balancing profiles:

* The fft involved in the Poisson equation spends the same time on
  every point. A ``uniform`` distribution is appropriate. Each process
  handles a parallelepipedic portion of the mesh box of the same size.

* The XC routine has no work to do on points where there is no charge
  density. A special algorithm determines the allocation of points to
  processors so that all have approximately the same workload. This is
  the ``linear`` distribution (this is a misnomer: it should be
  something like ``yes-no``).

* During computation of the charge density, and also during the
computation of the matrix elements of the 'grid' part of H, the work
on a given point depends on how many *pairs* or orbitals overlap over
it. The appropriate distribution is called ``quadratic``.

The work-load calculations, the nested-bisection algorithm used to
partition the mesh optimally, and all the communications involved in
redistributing mesh arrays, were implemented by Rogeli Grima of BSC
(Barcelona Supercomputing Center) in modules [[moremeshsubs]], [[meshcomm]], and [[schecomm]].

For more details see [this paper](http://doi.org/10.1007/s00214-010-0816-5).

### Re-distribution of orbital-based data to mesh form and back

In [[rhoofd]], one needs the density-matrix to compute the charge
density. Each non-zero entry in the DM corresponds to a pair of
orbitals that overlap. Those points in the overlap region will get a
contribution to the charge density. This is straightforward in serial
execution, but in parallel mode the distribution of the DM is
block-cyclic over orbitals, whereas the grid points are distributed
independently, and those points in the overlap region might not be handled by
the same processor that keeps the information on the DM entry.

Hence a re-distribution step is needed, coded in routine
[[matrixOtoM]] in module [[meshdscf]]. The new distribution of the DM
is such that all the entries relevant for a given mesh point are kept
in the same processor as the point itself.

The same problem appears, but in reverse, in [[vmat]], which computes
matrix elements of the Hamiltonian. The mesh arrays (typically potentials)
are to be sandwiched between two orbitals. Here an interim data structure
distributed in mesh form is created to hold the matrix elements, and then it
is converted to orbital form using routine [[matrixMtoO]].

In [[delk]], used in the Wannier90 interface, the interim structure is
converted and added to another array using routine [[matrixMtoOC]].

The communication operations needed by module [[meshdscf]] are
implemented in [[m_dscfcomm]], also coded by Rogeli Grima of BSC.  The
optimal sequence of point-to-point communications needed by these and
other mesh redistribution operations is determined by graph-coloring
techniques, implemented in [[scheComm]].

