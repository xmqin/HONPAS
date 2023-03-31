! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program hs2hsx

! Generates an HSX file  (more compact) from an HS file

use hsx_m, only: write_hsx_file, hsx_t, read_hs_file

type(hsx_t)  :: h

call read_hs_file(h,"HS") 

call write_hsx_file(h,"HSX_OUT")

end program hs2hsx

