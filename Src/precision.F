!!@LICENSE
!
      module precision
!
!     Precision handling
!
      implicit none

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
!
!     Select precision of certain parts of the program
!     These are set through preprocessor directives. The
!     default behavior is to use single-precision variables
!     
!     If you use GRID_DP, please note that it is advantageous
!     to enable also PHI_GRID_SP, since the phi array on the
!     grid can safely be kept in single precision, with significant
!     savings in memory and negligible numerical changes.

#ifdef GRID_DP
      integer, parameter :: grid_p = dp

#ifdef PHI_GRID_SP
      integer, parameter :: phi_grid_p = sp
#else
      integer, parameter :: phi_grid_p = dp
#endif

#else
      integer, parameter :: grid_p = sp
      integer, parameter :: phi_grid_p = sp
#endif

#ifdef ON_DP
      integer, parameter :: on_p = dp
#else
      integer, parameter :: on_p = sp
#endif

#ifdef BROYDEN_DP
      integer, parameter :: broyden_p = dp
#else
      integer, parameter :: broyden_p = sp
#endif


      public

      end module precision

