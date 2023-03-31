! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module files
!
!     Contains the short system label, used to generate file names
!     slabel is currently set in reinit.
!
      integer, parameter, public                  :: label_length = 64
      character(len=label_length), save, public   :: slabel

      ! Standard files for output/input

      ! STDIN, the file that should be read fdf-options from
      character(len=label_length), save, public :: stdin_file
      ! STDOUT, the file that is printed to.
      !    If ' ', then regular STDOUT is used
      character(len=label_length), save, public :: stdout_file

      
      ! Derived type to hold some output file names
      type, public:: filesOut_t
        character(len=label_length+6)::
     &    rho   = ' ',  ! (pseudo)electron density
     &    drho  = ' ',  ! diff. between SCF and atomic electron densities
     &    rhoxc = ' ',  ! electron density including nonlinear core correction
     &    psch  = ' ',  ! soft diffuse ionic charge
     &    toch  = ' ',  ! total ionic+electronic charge
     &    vh    = ' ',  ! Hartree electrostatic potential
     &    vt    = ' ',  ! total effective potential
     &    vna   = ' '   ! neutral-atom potential
      end type filesOut_t

      private

      end module files
