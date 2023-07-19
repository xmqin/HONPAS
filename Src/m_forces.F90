module m_forces

  use precision, only: dp
	implicit none

  public

  real(dp), pointer, save :: fa(:,:)    ! Atomic forces
  integer:: ntcon   ! Total number of position constraints imposed
  real(dp), pointer, save :: cfa(:,:)   ! Atomic forces orthogonalized 
                                        ! to geometry constraints

end module m_forces





