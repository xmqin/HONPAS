! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_fire_para

use precision, only: dp
use parallel, only: Node
use m_mpi_utils, only: Globalize_sum

implicit none

  real(dp), parameter, public  :: FIRE_DEF_dt_inc = 1.1_dp
  real(dp), parameter, public  :: FIRE_DEF_dt_dec = 0.5_dp
  real(dp), parameter, public  :: FIRE_DEF_alphamax = 0.1_dp
  real(dp), parameter, public  :: FIRE_DEF_alpha_factor = 0.99_dp
  real(dp), parameter, public  :: FIRE_DEF_dtmax = 10.0_dp
  integer,  parameter, public  :: FIRE_DEF_nmin = 5

type, public :: fire_t
  integer   :: n
  integer   :: step = 0
  integer   :: npos = 0
  real(dp)  :: alpha 
  real(dp)  :: dt
  integer   :: nmin = FIRE_DEF_nmin
  real(dp)  :: dt_inc = FIRE_DEF_dt_inc
  real(dp)  :: dt_dec = FIRE_DEF_dt_dec
  real(dp)  :: alphamax = FIRE_DEF_alphamax
  real(dp)  :: alpha_factor = FIRE_DEF_alpha_factor
  real(dp)  :: dtmax = FIRE_DEF_dtmax
  real(dp), pointer :: v(:) => null()  ! velocities
  logical   :: debug = .false.
end type

private

public :: fire_setup, fire_step

CONTAINS

  subroutine fire_setup(b,n,dt,nmin,dt_inc,dt_dec,  &
                        alphamax,alpha_factor,debug)

    type(fire_t), intent(inout) :: b
    integer, intent(in)         :: n  ! size of problem
    real(dp), intent(in)        :: dt
    logical, intent(in), optional  :: debug
    integer, intent(in), optional  :: nmin
    real(dp), intent(in), optional :: dt_inc
    real(dp), intent(in), optional :: dt_dec
    real(dp), intent(in), optional :: alphamax
    real(dp), intent(in), optional :: alpha_factor

    if (present(debug)) then
       b%debug = debug
    endif
    if (present(nmin)) then
       b%nmin = nmin
    endif
    if (present(dt_inc)) then
       b%dt_inc = dt_inc
    endif
    if (present(dt_dec)) then
       b%dt_dec = dt_dec
    endif
    if (present(alphamax)) then
       b%alphamax = alphamax
    endif
    if (present(alpha_factor)) then
       b%alpha_factor = alpha_factor
    endif

    b%alpha = b%alphamax
    b%n  = n
    b%dt = dt
    b%npos = 0
    b%step = 0

    allocate(b%v(n))
    b%v = 0.0_dp

  end subroutine fire_setup

!----------------------------------------------------------

  subroutine fire_step(b,f,x,maxstep)
    type(fire_t), intent(inout)  :: b
    real(dp), intent(in)         :: f(:)
    real(dp), intent(inout)      :: x(:)
    real(dp), intent(in)         :: maxstep(:)

    integer :: i
    logical :: debug

! For a general problem, we should choose dt and the masses
! so that the first step does not travel too far compared to some
! "typical length". Say, delta_x = 0.1 bohr

      real(dp) dt, norm_v, norm_f, n_long, power
      real(dp) local_norm_v, local_norm_f, local_power, local_n_long

      real(dp), allocatable :: dx(:)
      integer, allocatable :: count(:)

      allocate(dx(size(x)))
      allocate(count(size(x)))
      count = 0

      debug = (b%debug .and. Node .eq. 0) 
      
      dt = b%dt
      b%step = b%step + 1

         Local_Power = sum(f*b%v)
         call Globalize_sum(Local_Power,Power)

         If (Power < 0) then
            if (debug) print *, "Power < 0"
            b%v = 0.0_dp
            b%dt = b%dt * b%dt_dec
            b%alpha = b%alphamax
            b%npos = 0
            norm_v = 0.0_dp
            if (b%debug) then
               local_norm_f = sqrt(sum(f*f))  ! For debugging -- fix later
               call Globalize_sum(Local_norm_f,norm_f)
            endif
         else
            if (debug) print *, "Power > 0"

            local_norm_v = sqrt(sum(b%v*b%v))
            call Globalize_sum(Local_norm_v,norm_v)
            local_norm_f = sqrt(sum(f*f))
            call Globalize_sum(Local_norm_f,norm_f)

            b%v = (1-b%alpha)*b%v + b%alpha * norm_v * f / norm_f

            if (b%npos > b%nmin) then
               b%dt = min(b%dt * b%dt_inc, b%dtmax)
               b%alpha = b%alpha * b%alpha_factor
            endif

            b%npos = b%npos + 1
         endif   ! Power

         ! Euler with midpoint correction
         dx =  b%v*dt/2     ! half of the mid-point rule (v(0))
         b%v = b%v + f*dt     ! v at dt
         dx = dx + b%v*dt/2  ! the other half  (v(t))

         where (abs(dx) <= maxstep)
            x = x + dx
         elsewhere
            x = x + sign(maxstep,dx)
            count = count + 1
         end where

         if (b%debug) then
            local_n_long = sum(count)
            call Globalize_sum(Local_n_long,n_long) ! Actually simple reduce
         endif

         if (debug) then
            write(6,"(a)") 'FIRE: step, dt, npos, alpha, magv, magf, long'
            write(6,"(a, i6, f10.4, i4, f10.4, 2f12.6, i5)")  &
                 "FIREdat:", b%step, b%dt, b%npos, b%alpha,  &
                 norm_v, norm_f, n_long
            if (b%n <=100) then
               write(6,"(a)") "FIRE: nvar, x, v, f"
               do i = 1, b%n
                  write(6,*) i, x(i), b%v(i), f(i)
               enddo
            endif
         endif

         deallocate(dx,count)

      end subroutine fire_step

      end module m_fire_para

