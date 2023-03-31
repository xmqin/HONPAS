module m_pivot_array

  implicit none

  public

  ! Current size of the pivoting arrays
  integer, save          :: Npiv = 0
  ! The pivoting array
  integer, save, allocatable :: ipiv(:)

  ! Counts of initializations
  integer, save :: N_init = 0

contains

  ! We initialize the pivoting array for rotating the inversion
  subroutine init_pivot(n)
    integer, intent(in) :: n

    N_init = N_init + 1

    if ( n > Npiv ) then
       Npiv = n
       if ( allocated(ipiv) ) deallocate(ipiv)
       allocate(ipiv(Npiv))
    end if

  end subroutine init_pivot

  subroutine clear_pivot()

    N_init = N_init - 1
    if ( N_init > 0 ) return

    ! Force it to be zero
    N_init = 0

    ! Deallocate the pivoting array
    deallocate(ipiv)
    Npiv = 0

  end subroutine clear_pivot

end module m_pivot_array
