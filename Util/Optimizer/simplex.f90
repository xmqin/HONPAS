! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
! Simplex redesign
! Alberto Garcia, March 2007

!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------

program simplex
use precision, only: dp
use vars_module, only: nvars, var, read_vars, constrained
use vars_module, only: generate_subs_file

use minimizer, only: energy=>objective_function
use m_amoeba, only: amoeba

implicit none

real(dp), dimension(:,:), allocatable, save  :: p
real(dp), dimension(:), allocatable, save    :: y, x_final, pmin

real(dp)      :: ftol      ! = 0.0001_dp
real(dp)      :: lambda    ! = 0.4_dp
real(dp)      :: lambda_min, lambda_factor

integer :: mem_stat, iter, itmax, iu, imin, imin_dummy(1)

integer :: i

   call read_vars("VARS")      ! Sets nvars

   allocate(p(nvars+1,nvars),stat=mem_stat)
   if (mem_stat /=0 ) stop "alloc error"
   allocate(y(nvars+1),stat=mem_stat)
   if (mem_stat /=0 ) stop "alloc error"
   allocate(x_final(nvars),stat=mem_stat)
   if (mem_stat /=0 ) stop "alloc error"
   allocate(pmin(nvars),stat=mem_stat)
   if (mem_stat /=0 ) stop "alloc error"

!--------------------------------------------------------------------
   call io_assign(iu)
   open(iu,file="PARAMS",form="formatted",status="old")
   read(iu,*) lambda
   read(iu,*) lambda_factor
   read(iu,*) lambda_min
   read(iu,*) itmax
   read(iu,*) ftol
   close(iu)

     print *, lambda, lambda_factor, lambda_min
     print *, itmax, ftol
!--------------------------------------------------------------------

!
!    Initial guess at a random point in the corresponding
!    interval, or whatever  was specified by the user in the
!    VARS file (this should probably be made selectable in
!    the PARAMS file).
!
     do i = 1, nvars
        p(1,i) = var(i)%x
     enddo
     call generate_subs_file(p(1,:),"initial")

do
   call generate_subs_file(p(1,:),"best_so_far")

   if (lambda < lambda_min) then
      print *,  "Lambda < ", lambda_min,". Convergence presumably reached"
      exit
   endif
   print *, "Start of an amoeba cycle: [lambda= ", lambda, "]"

!
!  Simple-minded handling of constraints
! 
   do i=1, nvars
      p(1+i,:) = p(1,:)
      p(1+i,i) = constrained(p(1+i,i) + lambda*var(i)%range,i)
   enddo

!$omp parallel do shared(p,y,nvars) private(i)
   do i=1, nvars+1
      y(i) = energy(p(i,:))
   enddo
!$omp end parallel do

   print *, "Points in simplex and values:"
   do i=1, nvars+1
     print *, p(i,:), " --- ", y(i)
   enddo

   call amoeba(p,y,nvars,ftol,energy,iter,maxeval=itmax)

   print *, "Amoeba cycle finished in ", iter, " iterations."
   print *,  "Points in simplex and values:"
   do i=1, nvars+1
      print *, p(i,:), " --- ", y(i)
   enddo

   !
   ! Make sure that the minimum is at the first point, so
   ! that restarts build the new simplex around the minimum.
   !
   imin_dummy = minloc(y)
   imin = imin_dummy(1)
   pmin(:) = p(imin,:)     ! Swapping...
   p(imin,:) = p(1,:)
   p(1,:) = pmin(:)

   lambda = lambda * lambda_factor  ! Next simplex will be smaller
enddo
!
!     Take averages of positions in final simplex
!
      print *, "Final (average) point:"
      do i=1, nvars
         x_final(i) = sum(p(:,i)) / (nvars + 1)
         var(i)%x = x_final(i)
         print *, "Final (average)", trim(var(i)%name), " :", var(i)%x
      enddo
      call generate_subs_file(x_final,"final_average")
      print *, "Final (actual minimum):"
      do i=1, nvars
         var(i)%x = pmin(i)
         print *, "Final", trim(var(i)%name), " :", var(i)%x
      enddo
      call generate_subs_file(pmin(:),"final")

      deallocate(p,y,x_final,pmin)

End Program simplex



