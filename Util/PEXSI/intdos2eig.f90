program intdos2eig
!
! This program processes the PEXSI_INTDOS file produced by
! the cumulative-DOS interface to PEXSI and produces a 
! "fake eigenvalues" file suitable to be processed by Util/Eig2DOS.
!
! The eigenvalues are placed at the jumps in the cumulative DOS, with
! as many copies as necessary to account for the jump (note that there
! is an implicit factor of two for the spin).
! 
!
implicit none
integer, parameter :: dp = selected_real_kind(10,100)

real(dp), allocatable :: e(:)
integer, allocatable :: incounts(:)
real(dp) :: delta, energy, ef, qtot
integer :: i, n, iostat, ic, njump, neigs, j

open(unit=1,file="PEXSI_INTDOS",form="formatted")

read(1,*) ef, qtot, n
print *, "Fermi level: ", ef

allocate(e(n), incounts(n))
do i = 1, n
   read(1,fmt=*) e(i), incounts(i)
enddo

close(1)

print *, "Read ", n, " values"

! Find number of fake eigenvalues
neigs = 0
do i = 2, n
   njump = (incounts(i)-incounts(i-1))
   if (mod(njump,2) /= 0) then
      write(0,*) " ** warning ** : odd jump in PEXSI_INTDOS"
   endif
   neigs = neigs + njump/2
enddo

! Write fake eigenvalues
! Simulate EIG format with 1 k-point, nspin=1
! Layout might be different, but equivalent if read list-directed,
! as in Eig2DOS.
!
open(unit=2,file="SIMEIG",form="formatted")
write(2,"(f10.4)") ef
write(2,"(3i5)") neigs, 1, 1
write(2,"(i5)") 1
!
delta = e(2) - e(1)
do i = 2, n
   njump = (incounts(i)-incounts(i-1))
   do j = 1, njump/2
      write(2,"(2f15.6)") e(i) - 0.5_dp*delta
   enddo
enddo
 
close(2)
deallocate(e, incounts)

end program intdos2eig

