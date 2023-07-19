! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module units
  ! Define various unit conversion factors from internal units.

  ! internally, siesta works with length: Bohr.
  !                               energy: Rydberg.
  !                                 time: femtosecond

  use precision, only : dp

  implicit none

!  The easy way to make sense of units conversion:

!  real(dp), parameter :: Bohr   = 1.0_dp
!  real(dp), parameter :: Rydberg = 1.0_dp
!  real(dp), parameter :: Femtosecond = 1.0_dp
!
!  Ang = Bohr / 0.529177
!   eV = Rydberg / 13.60580
!  Joule = eV / 1.6e-19_dp
!  Meter = Ang / 1.0e-10_dp
!  Pascal = Joule/Meter**2
!   kBar  = Pascal * 1.0e4
!   .... and so on.

  real(dp), parameter :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter :: eV     = 1._dp / 13.60580_dp
  real(dp), parameter :: kBar   = 1._dp / 1.47108e5_dp
  real(dp), parameter :: GPa    = kBar * 10
  real(dp), parameter :: Kelvin = eV / 11604.45_dp
  real(dp), parameter :: Debye  = 0.393430_dp
  real(dp), parameter :: amu    = 2.133107_dp

! pi to 50 digits
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter :: deg = pi / 180.0_dp

end module units
