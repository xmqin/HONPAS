! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
      module precision
!
!     Precision handling
!
      implicit none
      
      integer, parameter :: i8b = selected_int_kind(18)
      integer, parameter :: int_8 = selected_int_kind(10)

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
!
!     Select precision of certain parts of the program
!     These are set through preprocessor directives. The
!     default behavior is to use single-precision variables
!     for the values of the orbitals on the grid, the Broyden
!     mixing auxiliary arrays, and the O(N) arrays, and 
!     double precision for the grid function arrays.
!
#ifdef BROYDEN_DP
      integer, parameter :: broyden_p = dp
#else
      integer, parameter :: broyden_p = sp
#endif

#ifdef GRID_SP
      integer, parameter :: grid_p = sp
      integer, parameter :: phi_grid_p = sp
#elif defined(GRID_DP)
      integer, parameter :: grid_p = dp
      integer, parameter :: phi_grid_p = dp
#else
      integer, parameter :: grid_p = dp
      integer, parameter :: phi_grid_p = sp
#endif

#ifdef ON_DP
      integer, parameter :: on_p = dp
#else
      integer, parameter :: on_p = sp
#endif

! For future use in trying to limit the memory usage by
! reducing precision in Transiesta
! Preprocessor flag could be: 
! dp = TRANSIESTA_DP
! sp = TRANSIESTA_SP
      integer, parameter :: ts_p = dp

      public

      end module precision

