! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2017.
!
! Please attribute the original author in case of dublication

! Enable the feature of manipulating with the
! Hamiltonian from outside.

! This is easily done by letting the user
! create a NetCDF file which contain dSE elements
! which will enter the NEGF equations through:
!   H = H + dSE
! We note that the user can input several different
!  dSE 
! corresponding to different situations.
! If dSE is complex, one can use dSE as an onsite
! self-energy change, thus introducing spurious
! effects.
!
! The dSE is effectively a change in the self-energies and
! will thus NOT enter in the equations for the bond-currents.
! If one has a dSE term that should contribute to the bond-currents
! please look at m_tbt_dH.F90.

! There are 4 levels of usage:
!   1. constant dSE
!      This will enter in a common equation
!   2. k-dependent dSE
!      This will only enter for the corresponding
!      k-point.
!         H_k = H_k + dSE_k' with dSE_k' = 0 for k/=k'
!   3. Energy
!      This will as 2. only enter for specific
!      energy points.
!         H(E) = H(E) + dSE(E') with dSE(E') = 0 for E/=E'
!   4. k and E dependent dSE
!         H_k(E) = H_k(E) + dSE_k'(E') with dSE_k'(E') = 0 for E/=E' and k/=k'

! Note that the highest level has precedence above the others,
! so specifying both 1. and 4. will only be the equivalent of
! using level 4.

! This is particularly handy for TB calculations
! but can prove just as useful for regular DFT
! calculations as one can use an already existing
! SIESTA.nc file, and then do several "case" studies
! on such a file, thus actually having FULL control
! over EVERYTHING.

module m_tbt_dSE

  use m_tbt_delta, only: tDelta

  implicit none
 
  private

  ! The dSE global variable
  type(tDelta), save :: dSE

  ! Whether this is used or not
  logical :: use_dSE = .false.

  public :: dSE, use_dSE

  public :: init_dSE_options, print_dSE_options

contains

  subroutine init_dSE_options( )
    use m_tbt_delta, only: init_delta_options, delta_has_level

    integer :: i

    use_dSE = .false.
    call init_delta_options('dSE', dSE)

    do i = 1 , 4
       use_dSE = use_dSE .or. delta_has_level(dSE, i)
    end do

  end subroutine init_dSE_options

  subroutine print_dSE_options( )

    use parallel, only : IONode

    character(len=*), parameter :: f10='(''tbt: '',a,t53,''='',tr4,a)'
    character(len=*), parameter :: f11='(''tbt: '',a)'
    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'

    if ( .not. IONode ) return

#ifdef NCDF_4
    if ( len_trim(dSE%fname) == 0 ) then
       write(*,f11)'No delta self-energy'
       return
    end if
    
    write(*,f10)'User selected dSE file', trim(dSE%fname)
#else
    write(*,f11)'delta self-energy cannot be used (requires NetCDF4)'
#endif
    
  end subroutine print_dSE_options
  
end module m_tbt_dSE
