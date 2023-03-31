! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2014, 2017.
!
! Please attribute the original author in case of dublication

! Enable the feature of manipulating with the
! Hamiltonian from outside.

! This is easily done by letting the user
! create a NetCDF file which contain dH elements
! which will enter the NEGF equations through:
!   H = H + dH
! We note that the user can input several different
!  dH 
! corresponding to different situations.
! If dH is complex, one can use dH as an onsite
! self-energy change, thus introducing spurious
! effects.
!
! The dH is effectively a change in the Hamiltonian and
! will thus enter in the equations for the bond-currents.
! If one has a dH term that should NOT contribute to the bond-currents
! please look at m_tbt_dSE.F90.

! There are 4 levels of usage:
!   1. constant dH
!      This will enter in a common equation
!   2. k-dependent dH
!      This will only enter for the corresponding
!      k-point.
!         H_k = H_k + dH_k' with dH_k' = 0 for k/=k'
!   3. Energy
!      This will as 2. only enter for specific
!      energy points.
!         H(E) = H(E) + dH(E') with dH(E') = 0 for E/=E'
!   4. k and E dependent dH
!         H_k(E) = H_k(E) + dH_k'(E') with dH_k'(E') = 0 for E/=E' and k/=k'

! Note that the highest level has precedence above the others,
! so specifying both 1. and 4. will only be the equivalent of
! using level 4.

! This is particularly handy for TB calculations
! but can prove just as useful for regular DFT
! calculations as one can use an already existing
! SIESTA.nc file, and then do several "case" studies
! on such a file, thus actually having FULL control
! over EVERYTHING.

module m_tbt_dH

  use m_tbt_delta, only: tDelta

  implicit none
 
  private

  ! The dH global variable
  type(tDelta), save :: dH

  ! Whether this is used or not
  logical :: use_dH = .false.

  public :: dH, use_dH

  public :: init_dH_options, print_dH_options
  public :: print_dH_warnings
  
contains

  subroutine init_dH_options( no_u )
    use m_tbt_delta, only: init_delta_options, delta_has_level

    integer, intent(in) :: no_u

    integer :: i

    use_dH = .false.
    call init_delta_options('dH', dH)

    do i = 1 , 4
       use_dH = use_dH .or. delta_has_level(dH, i)
    end do

#ifdef NCDF_4
    if ( use_dH ) call check_consecutive(no_u)
#endif

  end subroutine init_dH_options

#ifdef NCDF_4
  subroutine check_consecutive( no_u )

    use class_Sparsity

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif

    use m_tbt_delta, only: delta_has_level
    use m_tbt_delta, only: read_delta_Sp

    integer, intent(in) :: no_u

    type(Sparsity) :: sp
    integer :: il, lvls(4)

    integer :: io, ind
    integer, pointer :: ncol(:), col(:), ptr(:)

    ! Store the used levels
    lvls(:) = dH%lvls(:)

    do il = 1 , 4
       if ( .not. delta_has_level(dH, il) ) cycle

       ! Fake only reading one level
       ! This is to prevent a union of all the sparse patterns
       dH%lvls(:) = 0
       dH%lvls(il) = lvls(il)
       call read_delta_Sp(dH, no_u, sp)
       dH%lvls(:) = lvls(:)

       ! Get the sparse pattern
       call attach(sp, n_col=ncol,list_col=col,list_ptr=ptr)

       do io = 1, no_u
          do ind = ptr(io) + 2, ptr(io) + ncol(io)

             if ( col(ind-1) >= col(ind) ) then
                write(*,*)'tbt: delta Hamiltonian has a non-consecutive sparse pattern.'
                write(*,*)'     Sorted and consecutive columns for each row in the CSR matrix is'
                write(*,*)'     a requirement to calculate the correct bond currents from delta Hamiltonian'
                call die('tbt: delta Hamiltonian has a non-consecutive sparse pattern')
             end if

          end do
       end do

       call delete(sp)

    end do
    
  end subroutine check_consecutive
#endif
    

  subroutine print_dH_options( )

    use parallel, only : IONode

    character(len=*), parameter :: f10='(''tbt: '',a,t53,''='',tr4,a)'
    character(len=*), parameter :: f11='(''tbt: '',a)'
    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'

    if ( .not. IONode ) return

#ifdef NCDF_4
    if ( len_trim(dH%fname) == 0 ) then
       write(*,f11)'No delta Hamiltonian'
       return
    end if
    
    write(*,f10)'User selected dH file', trim(dH%fname)
#else
    write(*,f11)'delta Hamiltonian cannot be used (requires NetCDF4)'
#endif
    
  end subroutine print_dH_options

  subroutine print_dH_warnings( save_DATA )

    use dictionary
    use parallel, only : IONode

    type(dictionary_t), intent(in) :: save_DATA

    logical :: has

    if ( .not. IONode ) return

#ifdef NCDF_4
    if ( len_trim(dH%fname) == 0 ) return

    ! Check whether COHP is requested
    has = ('COHP-Gf' .in. save_DATA)
    has = has .or. ('COHP-A' .in. save_DATA)
    
    if ( has ) then
       write(*,'(a)') ' COHP curves are currently untested with dH terms!'
    end if

#endif
    
  end subroutine print_dH_warnings
  
end module m_tbt_dH
