title: TranSiesta details

This part of the documentation details the TranSiesta implementation.

TranSiesta is a solution method that calculates the (non) equilibrium density
by using open boundary conditions.

Currently TranSiesta does not implement non-collinear, nor spin-orbit calculations.

TranSiesta calculates the density matrix using the following procedure:

1. Setup all parameters that are to be used in the TranSiesta calculation.
   This is mainly the electrode Hamiltonian and self-energy calculations.
   The self-energies require the coupling along the semi-infinite direction
   to only couple to nearest neighbours, see [[m_ts_init:ts_init]].
2. Calculate an initial guess of the Hamiltonian and Fermi-level using
   any regular Siesta method, see [[compute_dm]].
3. When solving the Hartree potential it is vital that the boundary is
   "fixed" such that the potential is kept constant where the electrode
   are fixed. This is because the open boundary couples to a constant
   potential, the electrodes.
   This boundary is chosen as the electrode that has the largest cross
   section along the semi-infinite direction, see [[m_ts_hartree:ts_hartree_fix]].
4. Start the Green function method. In the
   TranSiesta SCF the electrode self-energies are taken from the TSGF*
   files and used to calculate the device Green function, see [[m_transiesta:transiesta]].
5. Mix the Hamiltonian or density matrix as in regular Siesta calculations.
6. Repeat step 3, 4 and 5.

Currently, the TranSiesta implemention can calculate the Green function
using 3 different methods:

1. The full Green function [[m_ts_fullg:ts_fullg]]/[[m_ts_fullk:ts_fullk]]
2. MUMPS library for sparse inversion [[m_ts_mumpsg:ts_mumpsg]]/[[m_ts_mumpsk:ts_mumpsk]]
3. Block tri diagonal inversion [[m_ts_trig:ts_trig]]/[[m_ts_trik:ts_trik]]

By far, the block tri diagonal inversion is the fastest and least memory
requiring method.
