! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module meshphi

      use precision, only : grid_p, phi_grid_p

      implicit none

      public :: resetMeshPhi

! ----------------------------------------------------------------------
! Stores variables and arrays related to basis orbitals on the mesh
! ----------------------------------------------------------------------
! logical DirectPhi     : If true the phi is calculated on the fly
! integer endpht(0:nmpl): Last position occupied by points in lstpht
! integer lstpht(ntopl) : List of non-zero orbitals at point
! integer listp2(ntopl) : Maps orbital-mesh point to iop
! integer nphi          : Length of phi array second dimension
! real    phi(nsp,ntopl): Basis orbitals at mesh points (sparse)
! ----------------------------------------------------------------------
      logical,                   save :: DirectPhi
      integer,                   save :: nphi
      integer,          pointer, save :: endpht(:)
      integer,          pointer, save :: lstpht(:)
      integer,          pointer, save :: listp2(:)
      real(phi_grid_p), pointer, save :: phi(:,:)

      CONTAINS

      subroutine resetMeshPhi( )
!
! Deallocate arrays of module meshphi if necessary
!
      use alloc,         only : de_alloc

      implicit none

      if (associated(endpht))
     $     call de_alloc( endpht, 'endpht',  'PhiOnMesh' )
      if (associated(lstpht))
     $     call de_alloc( lstpht, 'lstpht',  'PhiOnMesh' )
      if (associated(listp2))
     $     call de_alloc( listp2, 'listp2',  'PhiOnMesh' )
      if (associated(phi))
     $     call de_alloc( phi, 'phi',  'PhiOnMesh' )

      end subroutine resetMeshPhi

      end module meshphi

!
