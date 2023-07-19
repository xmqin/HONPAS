module m_eo

  use precision, only: dp
	implicit none

  public

  real(dp), pointer, save :: eo(:,:,:)  ! Hamiltonian eigenvalues
  real(dp), pointer, save :: qo(:,:,:)  ! Occupations of eigenstates

end module m_eo





