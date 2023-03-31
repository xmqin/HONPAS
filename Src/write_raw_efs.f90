! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine write_raw_efs(s,n,f,e)
use precision, only: dp
use siesta_geom, only: isa
use atomlist, only: iza
use atmfuncs, only: labelfis

implicit none

integer, intent(in)  :: n     ! number of atoms
real(dp), intent(in) :: s(3,3)      ! stress
real(dp), intent(in) :: f(3,n)      ! forces
real(dp), intent(in) :: e           ! (Free-electronic) Energy

! Everything in internal Siesta units

integer i, j, iu

call io_assign(iu)
open(iu,file="FORCE_STRESS",form="formatted",status="unknown", &
     action="write",position="rewind")

write(iu,"(f30.10)") e
write(iu,'(3x,3f18.9)')  ((s(i,j), i=1,3),j=1,3)
write(iu,*) n
do i = 1,n
   write(iu,'(i3,i6,3f18.9,2x,a)') isa(i),iza(i),f(1:3,i), trim(labelfis(isa(i)))
enddo
call io_close(iu)

end subroutine write_raw_efs

   
