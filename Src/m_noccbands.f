! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_noccbands
! This module finds the number of occupied bands


      implicit none

      integer, pointer, save, public, dimension(:) :: noccupied
! *********************************************************************
! integer noccupied: Number of occupied bands. 
!                    Its dimension will be determined by the number of
!                    different spin components.
! *********************************************************************

      contains 

      subroutine number_of_occupied_bands( qspin )

!
! Modules
!
      use precision
      use parallel,     only: IOnode
      use atmfuncs,     only: zvalfis
      use siesta_geom,  only: na_s, na_u, isa
      use alloc,        only: re_alloc
      use sys,          only: die
      use m_spin,       only: nspin

      real(dp), intent(in) :: qspin(2)
! *********************************************************************
! real*8  qspin  : Total population of spin up and down
! *********************************************************************

! Internal variables
      integer  :: ia, is, ispin, ntote
      real(dp) :: ntote_real, dq

! *********************************************************************
! real(dp) ntote_real  : Total number of electrons in the unit cell.
! *********************************************************************

! Reallocate arrays
      call re_alloc( noccupied, 1, nspin, name='noccupied', 
     .               routine='number_of_occupied_bands' )


! Compute the total number of electrons in the unit cell ---
!     If the calculation is spin unpolarized ---
      if (nspin .eq. 1) then
        ntote_real = 0.0_dp

        do ia = 1, na_u
          is = isa(ia) 
          ntote_real = ntote_real + zvalfis(is)
        enddo

        if ((nint(ntote_real) - ntote_real) .gt. 1.0e-6_dp) then
           if (IOnode)
     .       write(6,'(/,a,/,a,/a)')
     .      'noccbands: Non-integer number of electrons',
     .      'noccbands: This is hardly an insulator',
     .      'noccbands: No polarization calculation performed'
!           call die
        endif

        ntote = nint(ntote_real)

        if ( mod(ntote,2) .ne. 0) then
          if (IOnode)
     .       write(6,'(/,a,/,a,/a)')
     .      'noccbands: Odd total number of electrons',
     .      'noccbands: This is hardly an insulator'
!     .      'noccbands: No polarization calculation performed'
!          call die
        endif
        noccupied(1) = ntote / 2
        if (IOnode) then
          write(6,'(/a,i5)')
     .      'noccbands: Total number of electrons ',ntote
           write(6,'(a,i5)')
     .      'noccbands: Number of occupied bands ', noccupied(1)
        endif

!     If the calculation is spin polarized ---
      else

        dq = dabs( qspin(1) - qspin(2) )
        if ( dabs(dq - nint(dq)) .gt. 1.0d-2 ) then
          if (IOnode) 
     .     write(6,'(/,a,/,a,/a)')
     .  'noccbands: Spin polarization should have an integer value',
     .  'noccbands: This is not an insulator for both spin components'
!     .  'noccbands: No polarization calculation performed'
!          call die
        endif

        noccupied(1) = nint( qspin(1) )
        noccupied(2) = nint( qspin(2) )

        if (IOnode) then
          write(6,'(/a,i5)')
     .      'noccbands: Total number of electrons ', 
     -       noccupied(1) + noccupied(2)

          do ispin = 1, nspin
            write(6,'(a46,i5,a3,i5)')
     .        'noccbands: Number of occupied bands for spin ', ispin, 
     -        ' = ',  noccupied(ispin)
          enddo 
        endif 

      endif
      
      end subroutine number_of_occupied_bands

      end module m_noccbands
     


