! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine ofc(fa, dx, na)
C *******************************************************************
C Writes force constants matrix to file
C Input forces are in Ry/Bohr and input displacements are in Bohr.
C Force Constants written in file are in eV / Ang
C Written by P.Ordejon. August'98.
C Dynamic memory and save attribute for fres introduced by J.Gale
C Sept'99.
C Re-structured by Alberto Garcia, April 2007
C ********* INPUT ***************************************************
C real*8 fa(3,na)             : atomic forces (in Ry / Bohr)
C real*8 dx                   : atomic displacements (in Bohr)
C integer na                  : number of atoms
C ********** BEHAVIOUR **********************************************
C On the first call (undisplaced coordinates), the forces should be 
C zero (relaxed structure).
C However, since the relaxation is usually not perfect, some
C residual forces are obtained. These residual forces are 
C substracted from the forces on other steps, to calculate the
C force constants matrix
C *******************************************************************

      use precision
      use files, only : slabel, label_length
      use units, only : Ang, eV
      use alloc, only : re_alloc

      implicit          none

      integer, intent(in)  ::    na
      real(dp), intent(in) ::    dx, fa(3,na)

      external          io_assign, io_close

C Saved  variables and arrays
      character(len=label_length+3), save :: fname
      logical,                       save :: frstme = .true.
      real(dp), dimension(:,:), pointer, save :: fres

      integer    :: i, ix, unit1


      if (frstme) then

        fname = trim(slabel) // '.FC'
        nullify( fres )
        call re_alloc( fres, 1, 3, 1, na, name='fres', routine='ofc' )

        call io_assign(unit1)
        open( unit1, file=fname, status='unknown' )
        rewind(unit1)
        write(unit1,'(a)') 'Force constants matrix'
        call io_close(unit1)
        fres(:,:) = fa(:,:)
        frstme = .false.

      else

         call io_assign(unit1)
         open( unit1, file=fname, status='old',position="append",
     $         action="write")
         do i=1,na
            write(unit1,'(3f15.7)') ((-fa(ix,i)+fres(ix,i))*
     .           Ang**2/eV/dx, ix=1,3)
         enddo
         call io_close(unit1)

      endif

      end
