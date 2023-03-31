! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      function length(a)
*
*     This function calculates the length of a vector defined in
*    cartesian coordinates.
*
      use precision
      real(dp), intent(in) :: a(3)
      real(dp) :: length
*
*     a(3) : cartesian coordinates of the vector.  
*
      length = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

      return
      end

