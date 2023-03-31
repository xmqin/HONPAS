! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE atom_options

!
! Contains options for KB and basis set generation
! NOTE: This module is not MPI-safe. It should be used
!       only by the master node.
!
! The reason for this is to avoid a compilation cascade

  implicit none
  PUBLIC

  logical :: write_ion_plot_files    ! Write small auxiliary files?
  logical :: new_kb_reference_orbitals    ! New scheme for KB reference orbitals
  logical :: debug_kb_generation     ! Write auxiliary files for KB projectors

CONTAINS

  subroutine get_atom_options()
    use fdf

#ifdef MPI
    use mpi_siesta
    use sys, only: die

    integer MPIerror, Node
    call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
    if (Node /= 0) call die("Atom options can only be used by master node")
#endif

    write_ion_plot_files = fdf_boolean('WriteIonPlotFiles',.false.)
    debug_kb_generation  = fdf_boolean('Atom.Debug.KB.Generation',.false.)
    new_kb_reference_orbitals  = fdf_boolean('KB.New.Reference.Orbitals',.false.)
  end subroutine get_atom_options

END MODULE atom_options
