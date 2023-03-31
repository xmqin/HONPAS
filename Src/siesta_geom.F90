! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE siesta_geom
  use precision
  implicit none

  save
  
  ! Number of atoms in supercell, unit cell
  integer                         :: na_s, na_u

  !unit cell/supercell vectors by columns
  real(dp)                        :: ucell(3,3), ucell_last(3,3)
  real(dp)                        :: scell(3,3), scell_last(3,3)

  ! Shape of the system
  character(len=150)              :: shape*10

  ! Unit cell volume  (dangerous: the old code might have a BUG,
  ! as the volume is printed at the end 
  ! without being updated for the final cell change, and used in
  ! the calculation of Pmol in write_subs.

  real(dp)                        :: volume_of_some_cell

  !> Diagonal elements of supercell
  integer :: nsc(3) = 1
  !> The supercell offsets for the equivalent indices (size: `3, product(nsc)`)
  integer, pointer :: isc_off(:,:) => null()
   
  !> Previous geometry diagonal elements of supercell
  integer :: nsc_old(3) = 0

  ! Matrix of auxiliary supercell
  integer :: mscell(3,3)

  ! Unit cell "velocity" (time derivative)
  real(dp):: vcell(3,3)

  ! Atomic coordinates
  real(dp), pointer               :: xa(:,:)
  real(dp), pointer               :: xa_last(:,:)

  ! Atomic velocities
  real(dp), pointer               :: va(:,:)

  ! integer isa(na)           : Species index of each atom
  ! character cisa(na)        : Reference string for each atom
  ! NB cisa is this length in order to contain "siesta:e<isa>"
  ! where isa is the siesta element index, and we allow max 999
  ! such indices 
  integer,  pointer               :: isa(:)
  character(len=11), pointer      :: cisa(:) 

END MODULE siesta_geom
