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
