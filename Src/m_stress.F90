! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE m_stress
  use precision
  implicit none

  ! Constrained stress tensor
  real(dp):: cstress(3,3)

  ! Kinetic contribution to stress tensor
  real(dp):: kin_stress(3,3) 

  ! Stress tensor without the intramolecular contribution
  real(dp):: mstress(3,3) 

  ! Total stress tensor, including kinetic components
  real(dp):: tstress(3,3) 

  ! Stress tensor = d_E/d_strain
  real(dp):: stress(3,3) 

END MODULE m_stress
