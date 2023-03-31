! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
! Swarm optimization
! Ref: Basic algorithm in Wikipedia article
! Needs lots of refinements, or it can be used
! to find good spots for further optimization
!
! Alberto Garcia, April 2007

!-----------------------------------------------------------------------

program swarm

use precision, only: dp
use vars_module, only: nvars, var, read_vars, constrained
use vars_module, only: generate_subs_file

use minimizer, only: energy=>objective_function

implicit none

real(dp), dimension(:,:), allocatable  :: x, best, v
real(dp), dimension(:), allocatable    :: y
real(dp), dimension(:), allocatable    :: r1, r2
real(dp), dimension(:), allocatable    :: gbest
real(dp), dimension(:), allocatable    :: best_value

real(dp) :: gbest_value

real(dp) :: dum, omega
real(dp) :: worst_value, best_of_this_iteration
logical  :: converged
integer  :: iter, iu, dum_array(1), min_part

! Parameters read from SWARM_PARAMS
integer       :: nparts      ! no of particles in swarm
real(dp)      :: ftol        ! Crude tolerance ...  0.001_dp
real(dp)      :: omega_max   !  0.90  original inertia
real(dp)      :: omega_min   !  0.40  final inertia
real(dp)      :: c1, c2      !  1.95,   2.10
integer       :: itmax       !  max no of iterations

integer :: i, j
external :: io_assign

!-------------------

   call read_vars("VARS")      ! Sets nvars

   call io_assign(iu)
   open(iu,file="SWARM_PARAMS",form="formatted",status="old")
   read(iu,*) nparts
   read(iu,*) omega_max
   read(iu,*) omega_min
   read(iu,*) c1
   read(iu,*) c2
   read(iu,*) ftol
   read(iu,*) itmax
   close(iu)

   if (nparts < nint(1.5 * nvars)) then
      print *, "Warning: No of particles rather small..."
   endif

     print *, "Number of variables: ", nvars
     print *, "Number of particles: ", nparts
     print *, "Omega, c1, c2: ", omega, c1, c2
     print *, "Max no of iterations, tolerance: ", itmax, ftol


   allocate(x(nvars,nparts),v(nvars,nparts))
   allocate(best(nvars,nparts))
   allocate(y(nparts),best_value(nparts))
   allocate(r1(nvars), r2(nvars))
   allocate(gbest(nvars))

   omega = omega_max
!
!    Initial random guesses
!
   do i = 1, nparts
     do j = 1, nvars
        call random_number(dum)
        x(j,i) = var(j)%min + dum * var(j)%range
     enddo
     v(:,i) = 0.0_dp
  enddo

!!!!!!!!!!!  open(99,file="INFO",form="formatted")

!$omp parallel do shared(x,y,nparts) private(i)
   do i=1, nparts
      y(i) = energy(x(:,i))
   enddo
!$omp end parallel do

!  Assign the best point to the global variable
   dum_array = minloc(y)
   min_part = dum_array(1)
   best_of_this_iteration = minval(y)
   gbest_value = best_of_this_iteration
   gbest(:) = x(:,min_part)
   call generate_subs_file(gbest(:),"best_so_far")

!  Individual bests now are just the actual values:
   best_value(1:nparts) = y(1:nparts)
   best(1:nvars,1:nparts) = x(1:nvars,1:nparts)

!  Crude convergence test for now
   worst_value = maxval(y)
   converged = ((worst_value - best_of_this_iteration) < ftol)

iter = 1   
do

!!$   do i = 1, nparts
!!$      write(99,"(a1,i2.2,10f12.6,:,20f12.6)") "x", i, x(:,i)
!!$      write(99,"(a1,i2.2,10f12.6,:,20f12.6)") "v", i, v(:,i)
!!$   enddo

   print *, "----------- iteration ", iter
   print *, "Best value, worst value: ", best_of_this_iteration, worst_value
   print *, "Best value so far: ", gbest_value
   do i=1, nvars
      var(i)%x = gbest(i)
      print *, "Best ", trim(var(i)%name), " :", var(i)%x
   enddo

   if (iter == itmax) then
      print *, "Max no of iterations reached"
      exit
   endif

   if (converged) then
      print *, "apparently converged"
      exit
   endif

   iter = iter + 1

   ! Update velocities
   do i = 1, nparts
      call random_number(r1(1:nvars))
      call random_number(r2(1:nvars))
      v(1:nvars,i) = omega * v(1:nvars,i) + &
                         c1*r1(1:nvars)*(best(:,i)-x(:,i)) + &
                         c2*r2(1:nvars)*(gbest(:)-x(:,i)) 
   enddo

   omega = omega - (omega_max-omega_min)/(itmax-1)
   print *, "New inertia: ", omega

   ! Update positions

   x(1:nvars,1:nparts) = x(1:nvars,1:nparts) + v(1:nvars,1:nparts)

   ! Check constraints

   do j = 1, nparts
      do i = 1, nvars

         ! Damping constraints: If the particle flies away,
         ! put particle in boundary, and reflect the
         ! velocity component with a damping factor
         ! See: Xu and Rahmat-Samii, 
         ! IEEE Trans. Antenna and Propag. Vol 55, no 3, March 2007

         if (x(i,j) < var(i)%min) then
            x(i,j) = var(i)%min
            call random_number(dum)
            v(i,j) = - dum * v(i,j)
         endif
         if (x(i,j) > var(i)%max) then
            x(i,j) = var(i)%max
            call random_number(dum)
            v(i,j) = - dum * v(i,j)
         endif
!!!!!!!         x(i,j) = constrained(x(i,j),i)
      enddo
   enddo

!$omp parallel do shared(x,y,nparts) private(i)
   do i=1, nparts
      y(i) = energy(x(:,i))
   enddo
!$omp end parallel do

   do i=1, nparts
      if (y(i) < best_value(i)) then
         best_value(i) = y(i)
         best(:,i) = x(:,i)
      endif
   enddo
   
   ! Update global best
   best_of_this_iteration = minval(y)
   if (best_of_this_iteration  < gbest_value) then
      gbest_value = best_of_this_iteration
      dum_array = minloc(y)
      min_part = dum_array(1)
      gbest(:) = x(:,min_part)
      call generate_subs_file(gbest(:),"best_so_far")
   endif

   worst_value = maxval(y)
   converged = ((worst_value - best_of_this_iteration) < ftol)

enddo
!

End Program swarm



