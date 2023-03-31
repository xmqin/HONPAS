title: Reference-counted derived types

(Work in progress)

Derived types are a very useful addition to a programming
language. They allow to keep together in the same place a number or
related variables.

Sometimes, the same basic data in a derived-type variable needs to be
associated to larger data structures. The sparsity pattern of SIESTA
matrices is an example: every matrix in the program (for a given
geometry) uses the same sparse indexing arrays. Traditionally, there
was a set of global arrays (numh, listh, etc) that would take care of
the indexing throughout. With derived types, those arrays could be put in
a "sparsity_info" variable for easier handling.

Now consider what would happen if one wants to implement an
extrapolation scheme for density matrices in molecular-dynamics
runs. We must keep some past history of "geometry/density matrix"
pairs. Since different geometries need potentially different indexing
arrays, the old "array set" scheme would break down very soon (a
primitive extrapolation was implemented at some point, with only one
previous geometry: it used numh_old, listh_old, etc arrays...), and
the bookeeping was already a nightmare.

With the derived-type approach, one could keep an array of
sparsity_info elements, and associate them somehow to the relevant
density-matrix and geometry history (maybe implemented also as
arrays). This is do-able but also complicated.

Now imagine that the density matrix "object" contains within itself a
pointer to the relevant "sparsity_info" object, and that a
finite-stack of geometry/density-matrix objects is implemented. The
extrapolation algorithm can transparently use the stack metaphor
(which is the appropriate one for a finite-length re-cyclable
history), and forget about the array bookeeping. When the information
for a geometry step is no longer needed, the objects can be destroyed,
taking care that no other object in the program is using references to
the internal data. This can be easily achieved by implementing
[reference counting](https://en.wikipedia.org/wiki/Reference_counting).

In SIESTA, reference-counted derived types (and containers) are
implemented in modules whose name begins with `class` (this is an
unfortunate, probably temporary, choice), and use
[templates](https://en.wikipedia.org/wiki/Template_metaprogramming) to
support containers of various types.

For example, the sparsity pattern "bud" is defined in [[class_sparsity]], and
the containers involved in density-matrix extrapolation are

* A "sparse-2D matrix" bud, defined in [[class_dSpData2D]] (double precision version)

* A pair container, defined in [[class_Pair_Geometry_dSpData2D]] by
  including the ``Pair.T90`` template file.

* A finite-stack container, defined in [[class_Fstack_Pair_Geometry_dSpData2D]] by including
  the ``Fstack.T90`` template file.

The basic implementation of reference counting and the associated
functionality is in the [[basic_type.inc]] file.

