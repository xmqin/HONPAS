! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
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
