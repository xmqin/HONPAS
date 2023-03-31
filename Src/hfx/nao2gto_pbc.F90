! *** Module: nao2gto_pbc ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Periodic boundary conditions for Gaussian-type orbitals
!!
!! \author Honghui Shang
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_pbc

  use precision, only: dp

  implicit none

  private

  public :: cell_pbc, real_to_scaled, scaled_to_real, trans_pbc

  real(dp), parameter :: tiny = 1.0e-12_dp

  interface cell_pbc
    module procedure :: pbc1, pbc2, pbc3
  end interface

contains

! *****************************************************************************
!> \brief Computes a scaled position from a cartesian one
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[out] s: scaled position
!! \param[in] r: cartesian position
!! \param[in] cell: unit cell in reciprocal space
! *****************************************************************************
  subroutine real_to_scaled(s, r, rcell)

    implicit none

    ! Arguments
    real(dp), dimension(3), intent(out)  :: s
    real(dp), dimension(3), intent(in)   :: r
    real(dp), dimension(3,3), intent(in) :: rcell

    ! -------------------------------------------------------------------------

    s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
    s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
    s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)

  end subroutine real_to_scaled

! *****************************************************************************
!> \brief Retrieves a cartesian position from a scaled one
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[out] r: cartesian position
!! \param[in] s: scaled position
!! \param[in] cell: unit cell in real space
! *****************************************************************************
  subroutine scaled_to_real(r, s, cell)

    implicit none

    ! Arguments
    real(dp), dimension(3), intent(out)  :: r
    real(dp), dimension(3), intent(in)   :: s
    real(dp), dimension(3,3), intent(in) :: cell

    ! -------------------------------------------------------------------------

    r(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
    r(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
    r(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

  end subroutine scaled_to_real

! *****************************************************************************
!> \brief Applies periodic boundary conditions to a position vector r
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] r: position in supercell
!! \param[in] cell: unit cell in real space
!! \param[in] cell_inv: unit cell in reciprocal space
!! \param[out] r_pbc: position in unit cell
! *****************************************************************************
  subroutine trans_pbc(r, cell, cell_inv, r_pbc)

    implicit none

    ! Arguments
    real(dp), intent(in)  :: r(3), cell(3,3), cell_inv(3,3)
    real(dp), intent(out) :: r_pbc(3)

    ! Local variables
    real(dp) :: s(3)

    ! -------------------------------------------------------------------------

    s(1) = cell_inv(1,1)*r(1) + cell_inv(2,1)*r(2) + cell_inv(3,1)*r(3)
    s(2) = cell_inv(1,2)*r(1) + cell_inv(2,2)*r(2) + cell_inv(3,2)*r(3)
    s(3) = cell_inv(1,3)*r(1) + cell_inv(2,3)*r(2) + cell_inv(3,3)*r(3)

    if ( (dabs(s(1))-0.5d0) .gt. tiny ) s(1) = s(1) - anint(s(1))
    if ( (dabs(s(2))-0.5d0) .gt. tiny ) s(2) = s(2) - anint(s(2))
    if ( (dabs(s(3))-0.5d0) .gt. tiny ) s(3) = s(3) - anint(s(3))

    r_pbc(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
    r_pbc(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
    r_pbc(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

  end subroutine trans_pbc

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Computes a position in the unit cell from one in the supercell
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] r: position in supercell
!! \param[in] cell: unit cell in real space
!! \param[in] rcell: unit cell in reciprocal space
!! \return position in unit cell
! *****************************************************************************
  function pbc1(r, cell, rcell) result(r_pbc)

    implicit none

    ! Arguments
    real(dp), dimension(3), intent(in)   :: r
    real(dp), dimension(3,3), intent(in) :: cell, rcell

    ! Local variables
    real(dp), dimension(3) :: r_pbc
    real(dp), dimension(3) :: s

    ! -------------------------------------------------------------------------

    s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
    s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
    s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)
    s(1) = s(1) - anint(s(1))
    s(2) = s(2) - anint(s(2))
    s(3) = s(3) - anint(s(3))
    r_pbc(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
    r_pbc(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
    r_pbc(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

  end function pbc1

! *****************************************************************************
!> \brief Computes a position in the unit cell from one in the supercell,
!!        including a translation
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] r: position in supercell
!! \param[in] cell: unit cell in real space
!! \param[in] rcell: unit cell in reciprocal space
!! \param[in] nl: translations to apply in each direction
!! \return position in unit cell
! *****************************************************************************
  function pbc2(r, cell, rcell, nl) result(r_pbc)

    implicit none

    ! Arguments
    real(dp), dimension(3), intent(in)   :: r
    real(dp), dimension(3,3), intent(in) :: cell, rcell
    integer, dimension(3), intent(in)    :: nl

    ! Local variables
    real(dp), dimension(3) :: r_pbc, s

    ! -------------------------------------------------------------------------

    s(1) = rcell(1,1)*r(1) + rcell(1,2)*r(2) + rcell(1,3)*r(3)
    s(2) = rcell(2,1)*r(1) + rcell(2,2)*r(2) + rcell(2,3)*r(3)
    s(3) = rcell(3,1)*r(1) + rcell(3,2)*r(2) + rcell(3,3)*r(3)
    s(1) = s(1) - real(nint(s(1)) - nl(1), dp)
    s(2) = s(2) - real(nint(s(2)) - nl(2), dp)
    s(3) = s(3) - real(nint(s(3)) - nl(3), dp)
    r_pbc(1) = cell(1,1)*s(1) + cell(1,2)*s(2) + cell(1,3)*s(3)
    r_pbc(2) = cell(2,1)*s(1) + cell(2,2)*s(2) + cell(2,3)*s(3)
    r_pbc(3) = cell(3,1)*s(1) + cell(3,2)*s(2) + cell(3,3)*s(3)

  end function pbc2

! *****************************************************************************
!> \brief Computes a distance in the unit cell from one in the supercell
!!
!! \par History
!!      - 01.2010 Imported from CP2K [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] ra: first point in supercell
!! \param[in] rb: second point in supercell
!! \param[in] cell: unit cell in real space
!! \param[in] rcell: unit cell in reciprocal space
!! \return distance in unit cell
! *****************************************************************************
  function pbc3(ra, rb, cell, rcell) result(rab_pbc)

    implicit none

    ! Arguments
    real(dp), dimension(3), intent(in)   :: ra, rb
    real(dp), dimension(3,3), intent(in) :: cell, rcell

    ! Local variables
    integer                :: icell, jcell, kcell
    integer , dimension(3) :: periodic
    real(dp)               :: rab2, rab2_pbc
    real(dp), dimension(3) :: r, ra_pbc, rab, rb_image, rb_pbc, s2r, rab_pbc

    ! -------------------------------------------------------------------------

    ra_pbc(:) = pbc1(ra(:), cell, rcell)
    rb_pbc(:) = pbc1(rb(:), cell, rcell)

    rab2_pbc = huge(1.0_dp)

    do icell=-periodic(1),periodic(1)
      do jcell=-periodic(2),periodic(2)
        do kcell=-periodic(3),periodic(3)
          r = real((/icell,jcell,kcell/),dp)
          call scaled_to_real(s2r,r,cell)
          rb_image(:) = rb_pbc(:) + s2r
          rab(:) = rb_image(:) - ra_pbc(:)
          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
          if (rab2 < rab2_pbc) then
            rab2_pbc = rab2
            rab_pbc(:) = rab(:)
          end if
        end do
      end do
    end do

  end function pbc3

end module nao2gto_pbc
