title: The outer loop

This functionality is not SIESTA-specific, but is implemented to
provide a more complete simulation package. The program has an outer
geometry loop: it computes the [electronic structure](|page|/core-es.html) (and
thus the forces and stresses) for a given geometry, updates the
atomic positions (and maybe the cell vectors) (see [[siesta_move]]) accordingly and moves on
to the next cycle.

The [[siesta_move]] routine is the main driver for all the changes in geometry needed to
implement various types of Molecular Dynamics (MD), relaxations, and dynamical-matrix calculations.

A new and more flexible method to implement extra functionality is based on Lua scripts: these
can read and modify the data structures in SIESTA (in particular the coordinates and forces). There
is a hook in [[siesta_move]] to handle control of the operations to the Lua interpreter. For more
information, see the manual.