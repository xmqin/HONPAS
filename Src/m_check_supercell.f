      module m_check_supercell

!     Fixes the supercell factors in those cases in which the
!     naive construction method does not work.

!     If the shortest distance between (super)lattice points is
!     smaller than the relevant cutoff, the factors connected
!     to the shortest supercell vector are considered for increase.
!     The procedure is repeated until the criterion is satisfied.
!
!     Two different algorithms are used. The one that gives the
!     smaller supercell is chosen.
!
      use precision, only : dp
      use m_minvec,  only : minvec
      use parallel, only  : IOnode

      implicit none

      public :: check_sc_factors
      private

      CONTAINS

      subroutine check_sc_factors(ucell, nsc, rmax)

      real(dp), dimension(3,3), intent(in) :: ucell
      integer, dimension(3), intent(inout) :: nsc
      real(dp),                 intent(in) :: rmax

      integer, dimension(3) :: original_nsc, nsc1, nsc2
      integer               :: original_vol, vol1, vol2, vol

      if (IOnode) write(6,"(a,3i5)")
     $                "Naive supercell factors: ", nsc

      original_nsc = nsc
      original_vol = nsc(1) * nsc(2) * nsc(3)      
!
!     Try both algorithms
!
      call check_sc_factors1(ucell, nsc, rmax)
      vol1 = nsc(1) * nsc(2) * nsc(3)      
      nsc1 = nsc
!
      nsc = original_nsc
      call check_sc_factors2(ucell, nsc, rmax)
      vol2 = nsc(1) * nsc(2) * nsc(3)      
      nsc2 = nsc

      if (vol1 > vol2) then
         nsc = nsc2
         vol = vol2
      else
         nsc = nsc1
         vol = vol1
      endif
!
!     Warn user 
!
      if (vol > original_vol)
     $     write(6,"(a)")
     $     "Auxiliary supercell needs to be larger than naive one..."
      if (vol > 2*original_vol)
     $     write(6,"(a)")
     $     "*** Auxiliary supercell is more than twice as big! "

      end subroutine check_sc_factors

      subroutine check_sc_factors1(ucell, nsc, rmax)

      real(dp), dimension(3,3), intent(in) :: ucell
      integer, dimension(3), intent(inout) :: nsc
      real(dp),                 intent(in) :: rmax

      real(dp), dimension(3,3) :: a, newcell, ctransf
      integer  :: i, j, dummy(1), mask(3)
      real(dp) :: minsize, rlen(3), coeffs(3)
!     real(dp) :: length(3)

!      print *, "rmax: ", rmax

      main_loop: do

         do i = 1, 3
            a(:,i) = ucell(:,i) * nsc(i)
         enddo

         call minvec(a, newcell, ctransf)
         minsize = minimum_vector_size(newcell)
!         print *, "minsize of rounded cell: ", minsize
         if (minsize > rmax) EXIT main_loop   ! we are done
!
!        Look for shortest vector in rounded cell, and increase
!        the diagonal factors of the vectors of the original cell
!        that enter into its linear combination.
!
         do i = 1, 3
            rlen(i) = sqrt(dot_product(newcell(:,i),newcell(:,i)))
         enddo
         dummy = minloc(rlen)
         j = dummy(1)            ! index of shortest vector in r cell
!         print *, "Shortest rcell vector: ", j
         coeffs(1:3) = ctransf(:,j)  ! linear combination
!         print *, "coeffs: ", coeffs
         mask(1:3) = 5000           ! large number (for minloc below)
         where ( coeffs /= 0)  mask = nsc
!         print *, "mask: ", mask
         dummy = minloc(mask)    ! index of smallest factor
         j = dummy(1)
!         print *, "increasing factor no. ", j
         nsc(j) = nsc(j) + 1
 
      enddo main_loop

      end subroutine check_sc_factors1
!
!     In the following version we increase the factor of
!     the smallest-length vector in the linear-combination.
!     
      subroutine check_sc_factors2(ucell, nsc, rmax)

      real(dp), dimension(3,3), intent(in) :: ucell
      integer, dimension(3), intent(inout) :: nsc
      real(dp),                 intent(in) :: rmax

      real(dp), dimension(3,3) :: a, newcell, ctransf
      integer  :: i, j, dummy(1)
      real(dp) :: minsize, rlen(3), coeffs(3)
      real(dp) :: length(3), mask(3)

!      print *, "rmax: ", rmax

      main_loop: do

!         print *, "Current nsc: ", nsc
         do i = 1, 3
            a(:,i) = ucell(:,i) * nsc(i)
            length(i) = sqrt(dot_product(a(:,i),a(:,i)))
         enddo
!         print "(a,3f10.4)", "Current lengths: ", length(:)

         call minvec(a, newcell, ctransf)
         minsize = minimum_vector_size(newcell)
!         print *, "minsize of rounded cell: ", minsize
         if (minsize > rmax ) EXIT main_loop
!
!        Look for shortest vector in rounded cell, and increase
!        the diagonal factors of the vectors of the original cell
!        that enter into its linear combination.
!
         do i = 1, 3
            rlen(i) = sqrt(dot_product(newcell(:,i),newcell(:,i)))
         enddo
         dummy = minloc(rlen)
         j = dummy(1)            ! index of shortest vector in r cell
!         print *, "Shortest rcell vector: ", j
         coeffs(1:3) = ctransf(:,j)  ! linear combination
!         print *, "coeffs: ", coeffs
         mask(1:3) = huge(1.0_dp)    ! large number (for minloc below)
         where ( coeffs /= 0)  mask = length*nsc
!         print *, "mask: ", mask
         dummy = minloc(mask)    ! index of smallest vector
         j = dummy(1)
!         print *, "increasing factor no. ", j
         nsc(j) = nsc(j) + 1
 
      enddo main_loop

      end subroutine check_sc_factors2


      function minimum_vector_size(cell) result(minsize)
      real(dp), intent(in)  :: cell(3,3)
      real(dp)              :: minsize

      real(dp) :: veclen
      integer :: j

      minsize = huge(1.0_dp)
      do j = 1,3
        veclen = sqrt(dot_product(cell(:,j),cell(:,j)))
        minsize = min(veclen, minsize)
      enddo

      end function minimum_vector_size

      end module m_check_supercell
