! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module m_ts_global_vars
  
  implicit none
  
  save

  ! Whether transiesta is the solver
  logical :: TSmode = .false.

  ! Controls the change from diagon to transiesta solver
  logical :: TSinit = .false. , TSrun = .false.

  ! Whether this is an only overlap calculation
  logical :: onlyS = .false.

contains

  subroutine ts_method_init( start )

    use parallel, only : IONode
    use densematrix, only: resetDenseMatrix
    
    logical, intent(in) :: start

    ! If we are not in transiesta mode, return
    if ( .not. TSmode ) return

    if ( start ) then
       ! Quick return if we are already inside transiesta
       if ( TSrun ) return

       ! We will immediately start Transiesta
       TSinit = .false.
       TSrun  = .true.
       
       if ( IONode ) then
          write(*,'(/a)')'transiesta: Starting immediately'
       end if

       ! Be sure to not have *too* much memory allocated
       ! before entering transiesta
       ! NOTE, this should not be needed, but just in case
       ! future developments "forgets" about these arrays.
       call resetDenseMatrix()
       
       call ts_print_transiesta()

    else

       ! Tell transiesta to initialize the Hamiltonian with siesta

       TSinit = .true.
       TSrun  = .false.
       
       if ( IONode ) then
          write(*,'(/,a,/)') 'transiesta: Initialization run using siesta'
       end if

    end if
    
  end subroutine ts_method_init

  subroutine ts_print_transiesta()
    use parallel, only : IONode
    if ( IONode ) then
       write(*,'(/,t22,a)') repeat('*',27)
       write(*,  '(t22,a)') '*  WELCOME TO TRANSIESTA  *'
       write(*,'(t22,a,/)') repeat('*',27)
    end if
  end subroutine ts_print_transiesta

end module m_ts_global_vars
