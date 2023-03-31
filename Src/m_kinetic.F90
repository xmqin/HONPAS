! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_kinetic

  use precision, only: dp
	implicit none

  public

  real(dp):: vn      = 0.0_dp   ! Velocity (time derivative) of the Nose thermostat
  real(dp):: vpr     = 0.0_dp   ! Velocity (time derivative) of the PR variables

  real(dp):: tempion = 0.0_dp   ! Ionic temperature

  real(dp):: kn       = 0.0_dp  ! Kinetic energy of the Nose' thermostat
  real(dp):: kpr      = 0.0_dp  ! Kinetic energy of the Parrinello-Rahman variables

end module m_kinetic



