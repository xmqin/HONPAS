title: The core electronic-structure subsystem

Given a geometry, and information about the atomic
species involved (pseudopotentials, basis sets, etc, stored
in the appropriate
[tables](|page|/datastructures/10-species-tables.html)), SIESTA will:

* Determine the sparsity of the problem ([[hsparse]])
* Generate an overlap matrix S ([[overlap]])
* Initialize or read from file an starting density matrix (DM) ([[init_dm]])
* Generate a Hamiltonian (H) ([[setup_H0]] and [[setup_hamiltonian]])

and then, in what is called a self-consistent-field cycle:

* Solve the H/S generalized problem to generate a new density
  matrix ([[m_compute_dm:compute_dm]]).
* Maybe mix the input and output DMs to generate a new one ([[mixer]])
* Generate a new H and iterate until self-consistency ([[setup_hamiltonian]])

The last two points above might be replaced by the generation of a
new H from the output DM and the mixing of the input and output Hs.

