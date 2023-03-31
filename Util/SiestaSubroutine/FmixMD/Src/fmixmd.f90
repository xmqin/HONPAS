! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine fmixmd( t0, tf, dtfast, dtconv, dtsamp, ffast, fconv, &
                   na, ma, xa, va, cell )

! Performs a mixed-force molecular dynamics simulation
! Ref: E.Anglada, J.Junquera, and J.M.Soler, PRE 68, 055701 (2003)
! Written by J.M.Soler. Nov.2003

  use fsiesta
  use sample_m

  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(dp),intent(in):: t0     ! Initial simulation time
  real(dp),intent(in):: tf     ! Final simulation time
  real(dp),intent(in):: dtfast ! Time step between fast-force evaluations
  real(dp),intent(in):: dtconv ! Step between converged-force evaluations
                               ! (must be a multiple of dtfast)
  real(dp),intent(in):: dtsamp ! Step between calls to sampling routine
                               ! (should be a multiple of dtconv)
  character(len=*),intent(in) :: ffast ! Fast force field label
  character(len=*),intent(in) :: fconv ! Coverged force field label
  integer, intent(in)   :: na        ! Number of atoms
  real(dp),intent(in)   :: ma(na)    ! Atomic masses
  real(dp),intent(inout):: xa(3,na)  ! Cartesian atomic coordinates
  real(dp),intent(inout):: va(3,na)  ! Atomic velocities
  real(dp),intent(inout):: cell(3,3) ! Simulation cell vectors

! Argument units are: fs, amu, Ang, Ang/fs
! Internal energy/force units are eV, eV/Ang
  integer, parameter :: utime = 10.1806_dp  ! (amu*Ang^2/eV) in fs

  integer  :: ia, it, nc, ns, nt
  real(dp) :: dt, e, fa(3,na), fc(3,na), fp(3,na), t

! Find some constants
  dt = dtfast/utime        ! Time step in units of amu*Ang^2/eV
  nt = nint(tf-t0)/dtfast  ! Total number of steps
  nc = nint(dtconv/dtfast) ! Number of steps between force corrections
  ns = nint(dtsamp/dtfast) ! Number of steps between sampling
  if (mod(ns,nc)/=0) print*, 'fmixmd: ERROR: ns not multiple of nc'

! Set physical units for communication with siesta
  call siesta_units( 'Ang', 'eV' )

! Launch siesta processes
  call siesta_launch( ffast )
  call siesta_launch( fconv )

! Find fast and converged forces for initial coordinates
  call siesta_forces( ffast, na, xa, cell, e, fa )
  call siesta_forces( fconv, na, xa, cell, e, fc )
  fa = fa + (fc-fa)*nc

! Time iteration
  do it = 1,nt
    t = t0 + it*dt

    ! New positions
    do ia = 1,na
      xa(:,ia) = xa(:,ia) + va(:,ia)*dt + fa(:,ia)/ma(ia) * dt**2/2
    end do

    ! New forces
    fp = fa       ! Save previous forces
    call siesta_forces( ffast, na, xa, cell, e, fa )
    if (mod(it,nc)==0) then
      call siesta_forces( fconv, na, xa, cell, e, fc )
      fa = fa + (fc-fa)*nc
    end if

    ! New velocities
    do ia = 1,na
      va(:,ia) = va(:,ia) + (fa(:,ia)+fp(:,ia))/(2*ma(ia)) * dt
    end do

    ! Sampling
    if (mod(it,ns)==0) &
      call sample( 'add', t*utime, cell, na, ma, xa, va, fa, e )

  end do ! it

! Stop siesta processes
  print *, "About to stop ffast..."
  call siesta_quit( ffast )
  print *, "About to stop fconv..."
  call siesta_quit( fconv )
  print *, "After call to siesta_quit"

end

