! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine ofc(fa, dx, na, has_constr, first)
! *******************************************************************
! Writes force constants matrix to file
! Input forces are in Ry/Bohr and input displacements are in Bohr.
! Force Constants written in file are in eV / Ang ** 2
! Written by P.Ordejon. August'98.
! Dynamic memory and save attribute for fres introduced by J.Gale
! Sept'99.
! Re-structured by Alberto Garcia, April 2007
! Added constrained by Nick Papior, June 2015
! Added displacement value by Nick Papior, June 2018
! ********* INPUT ***************************************************
! real*8 fa(3,na)             : atomic forces (in Ry / Bohr)
! real*8 dx                   : atomic displacements (in Bohr)
! integer na                  : number of atoms
! logical has_constr          : whether the forces are constrained or not
! logical first               : first call (initializes the force)
! ********** BEHAVIOUR **********************************************
! On the first call (undisplaced coordinates), the forces should be 
! zero (relaxed structure).
! However, since the relaxation is usually not perfect, some
! residual forces are obtained. These residual forces are 
! substracted from the forces on other steps, to calculate the
! force constants matrix
! *******************************************************************

  use precision
  use files, only : slabel, label_length
  use units, only : Ang, eV
  use alloc, only : re_alloc

  implicit none

  integer,  intent(in) :: na
  real(dp), intent(in) :: dx, fa(3,na)
  logical,  intent(in) :: has_constr, first

  external :: io_assign, io_close

  ! Local reference force for the un-displaced configuration.
  real(dp), dimension(:,:,:), pointer, save :: fres => null()

  ! Saved  variables and arrays
  character(len=label_length+4) :: fname
  real(dp), dimension(:,:), pointer :: fr
  real(dp) :: tmp
  
  integer :: i, ix, iu
  
  if ( .not. associated(fres) ) then
    call re_alloc( fres, 1, 3, 1, na, 1, 2, name='fres', routine='ofc' )
  end if
  
  if ( has_constr ) then
    fname = trim(slabel) // '.FCC'
    fr => fres(:,:,2)
  else
    fname = trim(slabel) // '.FC'
    fr => fres(:,:,1)
  end if

  if ( first ) then

    call io_assign(iu)
    open( iu, file=fname, status='unknown' )
    rewind(iu)

    tmp = dx / Ang
    if ( has_constr ) then
      write(iu,'(2a,i0,tr1,e22.16)') 'Force constants matrix (constrained). ', &
          'n_atoms, displacement [Ang]: ', na, tmp
    else
      write(iu,'(2a,i0,tr1,e22.16)') 'Force constants matrix. ', &
          'n_atoms, displacement [Ang]: ', na, tmp
    end if

    call io_close(iu)

    ! Copy over the residual (zero point) forces
    fr(:,:) = fa(:,:)

    ! We should not save anything now
    return

  end if

  call io_assign(iu)
  open( iu, file=fname, status='old',position="append", action="write")

  tmp = Ang ** 2 / eV / dx
  do i = 1 , na
    write(iu,'(3(tr1,e17.9))') ((-fa(ix,i)+fr(ix,i))*tmp,ix=1,3)
  enddo

  call io_close(iu)

end subroutine ofc
