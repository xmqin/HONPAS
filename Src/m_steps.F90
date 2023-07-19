module m_steps
  implicit none
  public
  integer:: inicoor           ! Initial step in geommetry iteration
  integer:: fincoor           ! Final step in geommetry iteration
  integer:: istp              ! Geommetry iteration step starting in istp=1
  logical:: final=.false.     ! Last geometry step?

end module m_steps
