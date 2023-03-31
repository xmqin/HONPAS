! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module trialorbitalclass

! Includes:
! Derived types:
!   trialorbital
! Variables:
!   cutoff: cutoff radii
!   T:      tolerances controling the cutoff radii
! Procedures:
!   function  gettrialwavefunction
!   real(dp)  function gettrialrcut( ofwhat )
!   integer   function gettriallmax( ofwhat )
!   function  gettrialcenter( ofwhat )
!   subroutine print_trialorb 

use precision,       only:dp

implicit none

type trialorbital
  real(dp),dimension(3) :: center ! Projection function center in 
                                  !   cristallographic coordinates relative 
                                  !   to the direct lattice vectors.
  real(dp),dimension(3) :: zaxis  ! Defines the axis from which the polar
                                  !   angle theta in spherical polar coordinates
                                  !   is measured.
                                  !   Default: (0.0 0.0 1.0)
  real(dp),dimension(3) :: xaxis  ! Defines the axis from which the azimuthal 
                                  !   angle phi in spherical coordinates is 
                                  !   measured.
                                  !   Must be orthogonal to z-axis.
  real(dp),dimension(3) :: yaxis  ! Angular momentum y-axis
  real(dp)              :: zovera ! z/a, diffusivity, spread. 
                                  !   Read from the nnkp file in Ang^-1
                                  !   Transformed later to Ang^-1
  integer               :: r      ! Radial quantum number
  integer               :: l      ! Angular momentum
  integer               :: mr     ! z-projection quantum number
  real(dp)              :: rcut   ! Siesta's cut-off radius: Bohr
  integer               :: lmax   ! Maximum total angular momentum
end type

! Cut-off radii in units of \alpha^-1
real(dp),parameter,dimension(3) ::                           &
                           cutoffs = (/6.934_dp,18.87_dp,35.44_dp/)
! Squared norm tolerance governing the cut-off radii
real(dp),parameter :: T = 0.0001_dp

CONTAINS

!
!<-----------------------WAVE FUNCTIONS----------------------------->
!
real(dp) function gettrialwavefunction( orbital, atpoint )
!
! Trial orbital as defined in Wannier90.
! See pages 35-38 of the Wannier90 Users Guide, Version 1.2
! Yields the function value at a given point relative to its center.
! It contains the library of Wannier90 trial orbitals in the form
! of statement functions.
!
! Argument: coordinates of the point in Bohrs
! Output unit: Bohr^{-3/2}
!

  use units,  only: pi    ! Value of pi

  implicit none
!
! Passed arguments
!

  real(dp), dimension(3), intent(in) :: atpoint
  type(trialorbital),     intent(in) :: orbital

  real(dp),dimension(3)              :: arg
! Angular-dependent factor of the wave-function
  real(dp)                           :: angular
!
! Statement functions: declaration
!
! x, y, and z coordinates
  real(dp)    :: x, y, z
! Distance with respect to the origin
  real(dp)    :: rr
! Spherical function
  real(dp)    :: sphere
! Angular functions associated with particular values of l (trialorbital%l),
! and m_r (trialorbital%mr) for l >= 0 
! s orbital
  real(dp)    :: s
! p orbitals
  real(dp)    :: px, py, pz
! d orbitals
  real(dp)    :: dz2, dxz, dyz, dxy, dx2y2
! f orbitals
  real(dp)    :: fz3, fxz2, fyz2, fzx2y2, fxyz, fxx23y2, fy3x2y2
! Angular functions associated with particular values of l (trialorbital%l),
! and m_r (trialorbital%mr) for l < 0 (hybrid functions)
! sp hybrids
  real(dp)    :: sp_1,    sp_2
! sp2 hybrids
  real(dp)    :: sp2_1,   sp2_2,   sp2_3
! sp3 hybrids
  real(dp)    :: sp3_1,   sp3_2,   sp3_3,   sp3_4
! sp3d hybrids
  real(dp)    :: sp3d_1,  sp3d_2,  sp3d_3,  sp3d_4,  sp3d_5
! sp3d2 hybrids
  real(dp)    :: sp3d2_1, sp3d2_2, sp3d2_3, sp3d2_4, sp3d2_5, sp3d2_6
  integer     :: rank
! Radial functions associated with different values of the r value
! defined in (trialorbital%r)
  real(dp)    :: R1, R2, R3

!
! Inverse square roots
!
  real(dp), parameter  :: rs2  = 1.0_dp/dsqrt(2.0_dp)
  real(dp), parameter  :: rs3  = 1.0_dp/dsqrt(3.0_dp)
  real(dp), parameter  :: rs6  = 1.0_dp/dsqrt(6.0_dp)
  real(dp), parameter  :: rs12 = 1.0_dp/dsqrt(12.0_dp)
!
! Constants depending on l and power of z
! See the prefactors before the angular functions in page 36
! of the Wannier90 Users Guide, Version 1.2
!
  real(dp), parameter  :: l0norm   = 1.0_dp/dsqrt(4.0_dp*pi)
  real(dp), parameter  :: l1norm   = dsqrt(3.0_dp)*l0norm
  real(dp), parameter  :: l2z2norm = dsqrt((5.0_dp/16.0_dp)/pi)
  real(dp), parameter  :: l2z1norm = dsqrt((15.0_dp/4.0_dp)/pi)
  real(dp), parameter  :: l2z0norm = l2z1norm*0.5_dp
  real(dp), parameter  :: l3z3norm = dsqrt(7.0_dp/pi)/4.0_dp
  real(dp), parameter  :: l3z2norm = dsqrt(21.0_dp/2.0_dp/pi)/4.0_dp
  real(dp), parameter  :: l3z1norm = dsqrt(105.0_dp/pi)/4.0_dp
  real(dp), parameter  :: l3z0norm = dsqrt(35.0_dp/2.0_dp/pi)/4.0_dp
!
! Real Spherical Harmonics
!
! To get a dimensionless spherical harmonics
!
  sphere(x,y,z,rank) = sqrt( x**2 + y**2 + z**2 )**rank
!
! s-orbital
!
       s(x,y,z)      = l0norm
!
! p-orbitals
!
      px(x,y,z)      = l1norm * x / sphere(x,y,z,1)
      py(x,y,z)      = l1norm * y / sphere(x,y,z,1)
      pz(x,y,z)      = l1norm * z / sphere(x,y,z,1)
!
! d-orbitals
!
     dz2(x,y,z)      = l2z2norm * (2.0_dp * z**2 - x**2 - y**2 )/sphere(x,y,z,2)
     dxz(x,y,z)      = l2z1norm * z * x / sphere(x,y,z,2)
     dyz(x,y,z)      = l2z1norm * z * y / sphere(x,y,z,2)
   dx2y2(x,y,z)      = l2z0norm * (x**2 - y**2)/ sphere(x,y,z,2)
     dxy(x,y,z)      = l2z1norm * x * y / sphere(x,y,z,2)
!
! f-orbitals
!
     fz3(x,y,z)      = l3z3norm *  &
 &                     (2.0_dp*z**2-3.0_dp*x**2-3.0_dp*y**2)*z/sphere(x,y,z,3)
    fxz2(x,y,z)      = l3z2norm *  &
 &                     (4.0_dp*z**2-x**2-y**2)*x/sphere(x,y,z,3)
    fyz2(x,y,z)      = l3z2norm *  &
 &                     (4.0_dp*z**2-x**2-y**2)*y/sphere(x,y,z,3)
  fzx2y2(x,y,z)      = l3z1norm *  &
 &                     z*(x**2-y**2)/sphere(x,y,z,3)
    fxyz(x,y,z)      = l3z1norm *  &
 &                     z*x*y*2.0_dp/sphere(x,y,z,3)
 fxx23y2(x,y,z)      = l3z0norm *  &
 &                     (x**2-3.0_dp*y**2)*x/sphere(x,y,z,3)
 fy3x2y2(x,y,z)      = l3z0norm *  &
 &                     (3.0_dp*x**2-y**2)*y/sphere(x,y,z,3)
!
! Hybrids. They follow the definitions given in page 37 
! of the Wannier90 Users Guide, Version 1.2
!
! sp hybrids
! 
  sp_1(x,y,z)   = rs2 * s(x,y,z) + rs2 * px(x,y,z)
  sp_2(x,y,z)   = rs2 * s(x,y,z) - rs2 * px(x,y,z)
!
! sp2 hybrids
! 
  sp2_1(x,y,z)  = rs3 * s(x,y,z) - rs6 * px(x,y,z) + rs2 * py(x,y,z)
  sp2_2(x,y,z)  = rs3 * s(x,y,z) - rs6 * px(x,y,z) - rs2 * py(x,y,z)
  sp2_3(x,y,z)  = rs3 * s(x,y,z) + rs6 * px(x,y,z) * 2.0_dp
!
! sp3 hybrids
! 
  sp3_1(x,y,z)  = 0.5_dp * ( s(x,y,z) + px(x,y,z) + py(x,y,z) + pz(x,y,z) )
  sp3_2(x,y,z)  = 0.5_dp * ( s(x,y,z) + px(x,y,z) - py(x,y,z) - pz(x,y,z) )
  sp3_3(x,y,z)  = 0.5_dp * ( s(x,y,z) - px(x,y,z) + py(x,y,z) - pz(x,y,z) )
  sp3_4(x,y,z)  = 0.5_dp * ( s(x,y,z) - px(x,y,z) - py(x,y,z) + pz(x,y,z) )
!
! sp3d hybrids
! 
  sp3d_1(x,y,z) =  rs3 * s(x,y,z)  - rs6 * px(x,y,z) + rs2 * py(x,y,z)
  sp3d_2(x,y,z) =  rs3 * s(x,y,z)  - rs6 * px(x,y,z) - rs2 * py(x,y,z)
  sp3d_3(x,y,z) =  rs3 * s(x,y,z)  + rs6 * px(x,y,z) * 2.0_dp
  sp3d_4(x,y,z) =  rs2 * pz(x,y,z) + rs2 * dz2(x,y,z)
  sp3d_5(x,y,z) = -rs2 * pz(x,y,z) + rs2 * dz2(x,y,z)
!
! sp3d2 hybrids
! (Bug corrected by J. Junquera in the definition of sp3d2_3 and sp3d2_4.
! In the original version by R. Korytar the px orbital was wrongly used 
! instead of the correct py orbital)
! 
  sp3d2_1(x,y,z) = rs6 * s(x,y,z) - rs2 * px(x,y,z) - rs12 * dz2(x,y,z) + &
 &                 dx2y2(x,y,z) / 2.0_dp
  sp3d2_2(x,y,z) = rs6 * s(x,y,z) + rs2 * px(x,y,z) - rs12 * dz2(x,y,z) + &
 &                 dx2y2(x,y,z) / 2.0_dp
  sp3d2_3(x,y,z) = rs6 * s(x,y,z) - rs2 * py(x,y,z) - rs12 * dz2(x,y,z) - & 
 &                 dx2y2(x,y,z) / 2.0_dp
  sp3d2_4(x,y,z) = rs6 * s(x,y,z) + rs2 * py(x,y,z) - rs12 * dz2(x,y,z) - &
 &                 dx2y2(x,y,z) / 2.0_dp
  sp3d2_5(x,y,z) = rs6 * s(x,y,z) - rs2 * pz(x,y,z) + rs3 * dz2(x,y,z)
  sp3d2_6(x,y,z) = rs6 * s(x,y,z) + rs2 * pz(x,y,z) + rs3 * dz2(x,y,z)
!
! Radial part 
! They follow the definitions given in page 38 
! of the Wannier90 Users Guide, Version 1.2
!
  R1(rr) = 2.0_dp * orbital%zovera**(3.0_dp/2.0_dp) *                     &
 &         exp( -orbital%zovera * rr )
  R2(rr) = 0.5_dp/sqrt(2.0_dp) * orbital%zovera**(3.0_dp/2.0_dp)  *       &
 &         (2.0_dp - orbital%zovera * rr ) * exp( -orbital%zovera * rr / 2.0_dp)
  R3(rr) = sqrt(4.0_dp/27.0_dp) * orbital%zovera**(3.0_dp/2.0_dp) *       &
 &         ( 1.0_dp - 2.0_dp * orbital%zovera * rr/3.0_dp   +             &
 &           2.0_dp * orbital%zovera**2 * rr**2 / 27.0_dp ) *             &
 &         exp( -orbital%zovera * rr / 3.0_dp )

!
! Executables                 !----------->
!
! argument is the relative position of a given point with respect the center
! of the trial function 
! arg = atpoint - orbital%center.
! Recall that this is done in phiatm()
  arg = atpoint

! Compute the x, y, and z components of arg
  x   = dot_product( orbital%xaxis, arg )
  y   = dot_product( orbital%yaxis, arg )
  z   = dot_product( orbital%zaxis, arg )
! Compute the distance between the point and the center of the trial function
  rr  = sphere(x,y,z,1)

!
! If out of the rcut sphere then vanish
!
  if ( rr .gt. orbital%rcut ) then
    gettrialwavefunction = 0.0_dp
    return
  endif

!
! Decipher arguments: 
! Compute the radial part of the trial function at the given point
!
  select case(orbital%r)
    case(1)
      gettrialwavefunction = R1(rr)
    case(2)
      gettrialwavefunction = R2(rr)
    case(3)
      gettrialwavefunction = R3(rr)
  end select
! Renormalize the wave function, since we cut it at R_c
  gettrialwavefunction = gettrialwavefunction / dsqrt(1.0_dp-T)

!
! Decipher arguments: 
! Compute the angular part of the trial function at the given point
!
  if ( rr .eq. 0.0_dp ) then
! Special treatment of the origin
    select case(orbital%l)
      case (0)
        angular = l0norm
! Hybrids are combinations of dz2 and s limits
      case (-1)
        angular = l0norm * rs2
      case (-2)
        angular = l0norm * rs3
      case (-3)
        angular = l0norm / 2.0_dp
      case (-4)
        if ( orbital%mr .lt. 4 ) then
!         angular = l0norm * rs3
        else
          angular = 0.0_dp
        endif
      case(-5)
          angular = l0norm * rs6
      case default
        angular = 0.0_dp
    end select

!! Special treatment of the origin
!    select case(orbital%l)
!      case (0)
!        angular = l0norm
!      case (2)
!        if(orbital%mr.eq.1) then
!! l=2,mr=1 doesn't vanish!
!!          angular = -l2z2norm  ... this is the limit in the z=0 plane
!! but other limiting value is +2*l2z2norm along the z axis
!! so there is a cusp and we take it's average: zero
!        else
!          angular = 0.0_dp
!        endif
!! Hybrids are combinations of dz2 and s limits
!      case (-1)
!        angular = l0norm/sqrt(2.0_dp)
!      case (-2)
!        angular = l0norm/sqrt(3.0_dp)
!      case (-3)
!        angular = l0norm/2.0_dp
!      case (-4)
!        if (orbital%mr.lt.4) then
!          angular = l0norm/sqrt(3.0_dp)
!        else
!          !angular = -l2z2norm/sqrt(2.0_dp)
!        endif
!      case(-5)
!        if (orbital%mr.lt.5) then
!          angular = l0norm/sqrt(6.0_dp)!+l2z2norm/sqrt(12.0_dp)
!        else
!          angular = l0norm/sqrt(6.0_dp)!-l2z2norm/sqrt(3.0_dp)
!        endif
!      case default
!        angular = 0.0_dp
!    end select

  else
!
! rr.not equal.0
!
    select case(orbital%l)
      case(0)
          angular = s(x,y,z)
      case(1)
        select case(orbital%mr)
          case(1)
            angular = pz(x,y,z)
          case(2)
            angular = px(x,y,z)
          case(3)
            angular = py(x,y,z)
        end select
      case(2)
        select case(orbital%mr)
          case(1)
            angular = dz2(x,y,z)
          case(2)
            angular = dxz(x,y,z)
          case(3)
            angular = dyz(x,y,z)
         case(4)
            angular = dx2y2(x,y,z)
          case(5)
            angular = dxy(x,y,z)
        end select
      case(3)
        select case(orbital%mr)
          case(1)
            angular = fz3(x,y,z)
          case(2)
            angular = fxz2(x,y,z)
          case(3)
            angular = fyz2(x,y,z)
          case(4)
            angular = fzx2y2(x,y,z)
          case(5)
            angular = fxyz(x,y,z)
          case(6)
            angular = fxx23y2(x,y,z)
          case(7)
            angular = fy3x2y2(x,y,z)
        end select
!
! Hybrids
!
      case(-1)
        select case(orbital%mr)
          case(1)
            angular = sp_1(x,y,z)
          case(2)
            angular = sp_2(x,y,z)
        end select
      case(-2)
        select case(orbital%mr)
          case(1)
            angular = sp2_1(x,y,z)
          case(2)
            angular = sp2_2(x,y,z)
          case(3)
            angular = sp2_3(x,y,z)
        end select
      case(-3)
        select case(orbital%mr)
          case(1)
            angular = sp3_1(x,y,z)
          case(2)
            angular = sp3_2(x,y,z)
          case(3)
            angular = sp3_3(x,y,z)
          case(4)
            angular = sp3_4(x,y,z)
        end select
      case(-4)
        select case(orbital%mr)
          case(1)
            angular = sp3d_1(x,y,z)
          case(2)
            angular = sp3d_2(x,y,z)
          case(3)
            angular = sp3d_3(x,y,z)
          case(4)
            angular = sp3d_4(x,y,z)
          case(5)
            angular = sp3d_5(x,y,z)
        end select
      case(-5)
        select case(orbital%mr)
          case(1)
            angular = sp3d2_1(x,y,z)
          case(2)
            angular = sp3d2_2(x,y,z)
          case(3)
            angular = sp3d2_3(x,y,z)
          case(4)
            angular = sp3d2_4(x,y,z)
          case(5)
            angular = sp3d2_5(x,y,z)
          case(6)
            angular = sp3d2_6(x,y,z)
        end select
    end select
  endif ! this ends the first if regarding rr = 0 or not

  !print *,3.0_dp/2_dp,1.5_dp,3.0_dp/2.0_dp,3.0_dp/2.0_dp
  gettrialwavefunction = gettrialwavefunction * angular
end function gettrialwavefunction

real(dp) function gettrialrcut( ofwhat )
  type(trialorbital),intent(in)    :: ofwhat
  gettrialrcut = ofwhat%rCut
end function

integer function gettriallmax( ofwhat )
  type(trialorbital),intent(in)    :: ofwhat
  gettriallmax = ofwhat%lMax
end function

function gettrialcenter( ofwhat )
  real(dp),dimension(3)            :: gettrialcenter
  type(trialorbital),intent(in)    :: ofwhat
  gettrialcenter = ofwhat%center
end function gettrialcenter

subroutine print_trialorb( what )
!
! This subroutine prints the information about the trial projection function
!
  type(trialorbital),intent(in)    :: what

  write(*,fmt='(a,3f8.3,a)') "print_trialorb: center = ",what%center," Bohr"
  write(*,fmt='(a,3f8.3)')   "print_trialorb: zaxis  = ",what%zaxis
  write(*,fmt='(a,3f8.3)')   "print_trialorb: xaxis  = ",what%xaxis
  write(*,fmt='(a,3f8.3)')   "print_trialorb: yaxis  = ",what%yaxis
  write(*,fmt='(a,1f8.3,a)') "print_trialorb: zovera = ",what%zovera," Bohr**-1"
  write(*,fmt='(a,i5)')      "print_trialorb: r      = ",what%r
  write(*,fmt='(a,i5)')      "print_trialorb: mr     = ",what%mr
  write(*,fmt='(a,i5)')      "print_trialorb: l      = ",what%l
  write(*,'(a,1f8.3,a)')     "print_trialorb: rcut   = ",what%rcut," Bohr"
  write(*,fmt='(a,i5)')      "print_trialorb: lmax   = ",what%lmax
end subroutine print_trialorb

endmodule trialorbitalclass
