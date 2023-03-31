! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module sample_m

CONTAINS

subroutine sample( task, t, cell, na, ma, xa, va, fa, ep )
  implicit none
  integer, parameter :: dp = kind(1.d0)
  character(len=*),   intent(in) :: task      ! 'init'|'add'|'print'
  real(dp), optional, intent(in) :: t         ! time
  real(dp), optional, intent(in) :: cell(3,3) ! Cell vectors
  integer,  optional, intent(in) :: na        ! Number of atoms
  real(dp), optional, intent(in) :: ma(:)     ! Atomic masses
  real(dp), optional, intent(in) :: xa(:,:)   ! Atomic coordinates
  real(dp), optional, intent(in) :: va(:,:)   ! Atomic velocities
  real(dp), optional, intent(in) :: fa(:,:)   ! Atomic forces
  real(dp), optional, intent(in) :: ep        ! Potential energy

  logical, save :: first_time = .true.
  real(dp),save :: ecsum=0, ec2sum=0, epsum=0, ep2sum=0, etsum=0, et2sum=0
  integer, save :: nt=0

  integer  :: ia
  real(dp) :: ec, et

  select case(task)
  case('init')
    ecsum = 0._dp
    epsum = 0._dp
    etsum = 0._dp
    ec2sum = 0._dp
    ep2sum = 0._dp
    et2sum = 0._dp
    nt = 0
  case('add')
    ec = 0._dp
    do ia = 1,na
      ec = ec + ma(ia) * sum(va(:,ia)**2) / 2._dp
    end do
    et = ec + ep
    ecsum  = ecsum  + ec
    ec2sum = ec2sum + ec**2
    epsum  = epsum  + ep
    ep2sum = ep2sum + ep**2
    etsum  = etsum  + et
    et2sum = et2sum + et**2
    nt = nt + 1
  case('print')
    ecsum  = ecsum/nt
    ec2sum = ec2sum/nt
    ec2sum = dsqrt(ec2sum-ecsum**2)
    epsum  = epsum/nt
    ep2sum = ep2sum/nt
    ep2sum = dsqrt(ep2sum-epsum**2)
    etsum  = etsum/nt
    et2sum = et2sum/nt
    et2sum = dsqrt(et2sum-etsum**2)
    print*, 'sample: nt =', nt
    print*, 'sample: Ec =', ecsum, ' +/-', ec2sum
    print*, 'sample: Ep =', epsum, ' +/-', ep2sum
    print*, 'sample: Et =', etsum, ' +/-', et2sum
  case default
    print*, 'sample: ERROR: unknown task: ', trim(task)
    stop
  end select

end subroutine sample

end module sample_m
