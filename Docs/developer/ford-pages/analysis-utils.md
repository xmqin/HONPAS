title: Analysis and utilities

After the computation of the electronic structure, SIESTA can
optionally carry out several kinds of analysis, and
write the results either to the output file or to specific files ([[siesta_analysis]]):

* Output of the final (relaxed) coordinates in various formats
* Computation of wavefunctions at specific k-points ([[wwave]])
* COOP/COHP analysis ([[wwave]] with the scf loop k-point sampling).
* Band-structure calculation (with optional fat bands information) ([[bands]])
* Output of eigenvalues
* Partial density of states ([[pdos]] and subsidiary routines)
* Output of energy contributions and spin information
* Output of dipole moment (for molecules) or computation of bulk polarization via Berry's phase ([[ksv_pol]]).
* Calculation of optical conductivity ([[optical]]).
* Output of charge densities, potentials, etc. (a further call to [[dhscf]])
* Processing of information for Wannier90 ([[siesta2wannier90]]).
* Computation of local-density of states ([[local_DOS]]).
* Python utility tool [sisl](|page|/analysis-sisl.html)

In addition, several programs in the `Util` subdirectory of the SIESTA distribution are available
to carry out specific tasks:

* Bands: Tools for plotting band structures (including "fatbands"): [[gnubands]], [[eigfat2plot]]

* CMLComp: Tools to use the information contained in the CML file
         produced by the program (by Toby White, Andrew Walker,
         and others)

* Contour: grid2d: As Denchar but for any function defined in the 3D
         grid.  grid1d: Extracts a 1D line of data out of the 3D grid.
         Based on the 3D grid (by E. Artacho)

* Contrib: Code contributed by Siesta users (expanding collection). See
	 the individual documentation materials.

* COOP: Generation of COOP/COHP/PDOS curves for chemical analysis.
      Computation of data to generate "fatbands" plots: [[mprop]], [[fatband]].

* Denchar: Produces 2D and 3D plots of charge and spin density, and of
         wave-functions. Uses the density matrix and basis orbital
         information (by J. Junquera and P. Ordejon).

* DensityMatrix: Utilities to process and convert density-matrix files.

* Eig2DOS:   Estimation of the Density of States from the .EIG file.
           (by E. Artacho and A. Garcia)

* Gen-basis : Stand-alone program 'gen-basis' to generate basis sets and
            KB projectors, and 'ioncat' program for extraction of
            information from .ion files.

* Grid: Utilities for the manipulation of grid files.

* Grimme: Enable easy creation of the atomic potential block by reading
	the ChemicalSpecies block and printing out the relevant information.
	This makes it _very_ easy to create the correct units and quantities.
	(by N. Papior)

* Helpers: Some helper scripts for aiding script generation.

* HSX: Conversion tool between HS and HSX files (by A. Garcia)

* JobList:   A suite of programs to generate and dispatch multiple
           SIESTA jobs with varying options (by J. Soler)

* Macroave: Macroscopic averaging for interfaces and surfaces (by
	  J. Junquera)

* MD: Some sample scripts for the extraction of some MD information from
    the output file (by A. Garcia)

* MM_Examples: Force-field examples

* MPI_test:    Tests to help diagonose the interface to MPI.

* ON: Conversion of Order-N eigenstate to NetCDF file

* Optical: Calculation of optical properties (by D. Sanchez-Portal)

* Optimizer: General-purpose optimizer. Useful for basis-set and
	   pseudopotential optimization (by A. Garcia)

* pdosxml: Utility to process the PDOS file (in XML format). (by
	 A. Garcia)

* PEXSI: Utilites related to the output from PEXSI (by A. Garcia)

* Plrho: Plots in 3D the electron density and other functions
       calculated by siesta (by Jose M. Soler)

* Projections: Compute projections of electronic structure of a system
	     over the orbitals of a subsystem.

* PyAtom: Python scripts for plotting and data extraction (by A. Garcia)

* SCF: Python scripts for smaller stuff (by A. Garcia)

* Scripting: Experimental scripting modules in Python (by A. Garcia)

* sies2arc: Converts output coordinates to the arc movie format (by
	  J. Gale)

* SiestaSubroutine: Code and examples of the driving of Siesta by an
                  external agent (by J. Soler and A. Garcia)

* Sockets: Examples of the use of the f90 sockets interface (by M. Ceriotti)

* SpPivot: A utility to create a pivoting table for the sparse
	 patterns in Siesta. Can create GRAPHVIZ output for easy
	 display of the sparsity pattern (by N. Papior)

* STM/simple-stm:   Simple program for STM simulations (by P. Ordejon)

* STM/ol-stm:  Ordejon-Lorente STM image simulator

* TS: Contains several different utilities that are related to
    the NEGF code TranSiesta. (by N. Papior)
    See TS/README for details.

* VCA: Utilities to help in Virtual-Crystal calculations (by A. Garcia)

* Vibra: Package to compute phonon frequencies and modes (by P. Ordejon)

* WFS: Utilities for wavefunction-file manipulation
