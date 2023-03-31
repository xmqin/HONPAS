module m_matrix
  integer, parameter, private :: dp = selected_real_kind(10,100)

  type, public ::  matrix
     integer :: norbs
     integer :: no_l
     integer :: nnzl
     integer, pointer :: numcols(:) => null()
     integer, pointer :: cols(:)    => null()
     real, pointer    :: vals(:)    => null()
  end type matrix

end module m_matrix

