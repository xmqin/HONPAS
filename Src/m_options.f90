! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_scf_options
!
use precision

real(dp), save  :: jinv0 = 0.2_dp
real(dp), save  :: kick_strength = 0.5_dp
real(dp), save  :: fluct_margin = 0.1_dp
integer, save   :: maxit = 5
integer, save   :: n_temp_jumps = 2
integer, save   :: n_temp_steps = 4
logical, save   :: initial_kick = .false.
logical, save   :: cycle_on_maxit = .true.
logical, save   :: variable_weight = .true.
integer, save   :: n_dnorm_samples = 4

end module m_scf_options

module m_options
!
!  Holds the Siesta run-time options
!
use m_scf_options
public

end module m_options
