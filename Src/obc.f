! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine obc( polR, ucell, dx, nspin)
C *******************************************************************
C Writes the born effective charges to file 
C by Tom Archer
C Modified by Alberto Garcia, April 2007
C ********* INPUT ***************************************************
C real*8 dx             : atomic displacements (in Bohr)
C real*8 polR(3,nspin)  : polarization along lattice vectors(Bohr)
C real*8 ucell(3,3)     : cell vectors
C integer nspin         : spin polarized calculation flag
C ********** BEHAVIOUR **********************************************
C On the first call (undisplaced coordinates), initial polarization  
C is calculated.
C The born charges are calculated as the difference between the 
C polarization of the undisplaced and the displaced coordinates 
C divided by the displacement
C 
C a phase is introduced in the KSV calculation which is removed before
C the BC matrix is saved.  
C *******************************************************************

      use precision, only: dp
      use fdf,       only: fdf_string
      use alloc,     only: re_alloc

      implicit          none

      integer, intent(in)  :: nspin 
      real(dp), intent(in) :: polR(3,nspin),
     $                        ucell(3,3), dx
      external io_assign, io_close

c Internal

      integer    igrd, ispin, ix, unit1

      character(len=33), save :: fname
      logical, save  ::    first_time = .true.

      real(dp) :: dmod(3), phaseR(3), uR(3,3), pR(3), pxyz(3)

      double precision, dimension(:,:), pointer, save :: pres

!---------

      do ispin=1,nspin
         if (nspin.gt.1) then
            if (ispin.eq.1) write(6,'(/,a)')
     .           'obc: Macroscopic polarization for spin Up:'
            if (ispin.eq.2) write(6,'(/,a)')
     .           'obc: Macroscopic polarization for spin Down:'
         endif
         write(6,'(/,a,a)')
     .        'obc: Macroscopic polarization per unit cell',
     .        ' along the lattice vectors (Bohr): '
         write(6,'(a,3f12.6)') 'obc:',(polR(ix,ispin),ix=1,3)
      enddo 

      if (first_time) then

        fname = trim(fdf_string('SystemLabel','siesta')) // '.BC'

        nullify( pres )
        call re_alloc( pres, 1, 3, 1, nspin, name='pres',
     &                 routine='obc' )

        call io_assign(unit1)
        open( unit1, file=fname, status='unknown' )
        rewind(unit1)
        if (nspin.eq.1) then
           write(unit1,'(a)') 'BC matrix'
        else
           write(unit1,'(a)') 'BC matrix for spin-up + spin-down'
        endif  
        !  Set values of residual polarization
        do ispin=1,nspin
           do igrd=1,3
              pres(igrd,ispin) = polR(igrd,ispin)
           enddo
        enddo
        call io_close(unit1)

        first_time = .false.

      else

          !matrix to convert from lattice vectors to xyz
         do igrd=1,3
            dmod(igrd)=dot_product(ucell(:,igrd),ucell(:,igrd))
            dmod(igrd)=sqrt(dmod(igrd))
            do ix=1,3
               uR(ix,igrd)=ucell(ix,igrd)/dmod(igrd)
            enddo
         enddo

         !     write BC matrix in lattice vectors
         pR(:)=0.0_dp
         pxyz(:)=0.0_dp
         do ispin=1,nspin
            do igrd=1,3
               pR(igrd)=pR(igrd)+(polR(igrd,ispin)-pres(igrd,ispin))
            enddo
         enddo
         
         !remove extra phase
         do igrd=1,3
            phaseR(igrd) = nint((pR(igrd)/dmod(igrd))) * dmod(igrd)
            pR(igrd) = pR(igrd) - phaseR(igrd) 
         enddo

         !get polarization in cartesian coordinates
         do ix=1,3
            pxyz(ix) = uR(ix,1)*pR(1)+ uR(ix,2)*pR(2)+uR(ix,3)*pR(3)
         enddo

         call io_assign(unit1)
         open( unit1, file=fname, status='old',position="append",
     $         action="write")
         write(unit1,'(3f15.7)') (pxyz(ix)/dx,ix=1,3)
         call io_close(unit1)

      endif

      end
