!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_contour
!
! Routines that are used to setup the contour for integration of the GFs
! 
! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype
  use m_ts_chem_pot

  use m_ts_cctype
  use m_ts_io_contour

  use m_ts_contour_eq,  only : CONTOUR_EQ
  use m_ts_contour_neq, only : CONTOUR_NEQ

  use precision, only : dp

  implicit none

  private

  public :: read_contour_options
  public :: print_contour_options
  public :: print_contour_block
  public :: io_contour
  public :: sort_contour

contains

  subroutine io_contour(IsVolt, mus, slabel,suffix)
    use m_ts_contour_eq, only : io_contour_Eq
    use m_ts_contour_neq, only : io_contour_nEq
    logical, intent(in) :: IsVolt
    type(ts_mu), intent(in) :: mus(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

    call io_contour_Eq(mus,slabel,suffix)

    if ( IsVolt ) then
       call io_contour_nEq(slabel,suffix)
    end if

  end subroutine io_contour

  subroutine read_contour_options(N_Elec, Elecs, N_mu, mus, kT, IsVolt, Volt)
    use m_ts_contour_eq, only : read_contour_eq_options
    use m_ts_contour_neq, only : read_contour_neq_options

    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: N_mu
    type(ts_mu), intent(inout) :: mus(N_mu)
    ! SIESTA electronic temperature
    real(dp), intent(in) :: kT
    logical, intent(in) :: IsVolt
    real(dp), intent(in) :: Volt

    call read_contour_eq_options(N_mu,mus,Volt)

    if ( IsVolt ) then
       call read_contour_neq_options(N_Elec,Elecs,N_mu,mus,kT,Volt)
    end if

  end subroutine read_contour_options

  subroutine print_contour_options(prefix,IsVolt)
    use m_ts_contour_eq, only : print_contour_eq_options
    use m_ts_contour_neq, only : print_contour_neq_options

    character(len=*), intent(in) :: prefix
    logical, intent(in) :: IsVolt

    call print_contour_eq_options(prefix)

    if ( IsVolt ) then
       call print_contour_neq_options(prefix)
    end if

  end subroutine print_contour_options

  subroutine print_contour_block(prefix, IsVolt)
    use m_ts_contour_eq, only : print_contour_eq_block
    use m_ts_contour_neq, only : print_contour_neq_block

    character(len=*), intent(in) :: prefix
    logical, intent(in) :: IsVolt

    call print_contour_eq_block(prefix)

    if ( IsVolt ) then
       call print_contour_neq_block(prefix)
    end if

  end subroutine print_contour_block

  ! In order to retain the best numerical accuracy we introduce
  ! an ordering of the energy contour by weights.
  ! It is the only reasonable thing to access as the functional is
  ! non-deterministic.
  ! An example of the importance of this sorting:
  ! Take +200 voltage energy points. 
  ! The weights will be something like this:
  ! 1e-5,...,1e-4,...,1e-3,...,1e-4,...,1e-5
  ! Which means that the summation up till the half works great.
  ! But when we reach the last weight we could be in the situation where:
  !  1._dp + 1e-12_dp which will limit the accuracy obtained for that energy point.

  ! In order to circumvent this we simply sort by weights.
  subroutine sort_contour(NC,ce,cw)
    integer, intent(in) :: NC
    complex(dp), intent(inout) :: ce(NC), cw(NC)
    
    ! Local variables
    real(dp) :: tmp
    complex(dp) :: ctmp
    integer :: i,j
    ! As we will only do this once we dont need a fancy sorting
    ! algorithm...
    do i = 1 , NC - 1
       tmp = real(cw(i)*conjg(cw(i)),dp)
       do j = i+1, NC
          if ( tmp > real(cw(j)*conjg(cw(j))) ) then
             ctmp  = ce(i)
             ce(i) = ce(j)
             ce(j) = ctmp
             ctmp  = cw(i)
             cw(i) = cw(j)
             cw(j) = ctmp
          end if
       end do
    end do

  end subroutine sort_contour

end module m_ts_contour
