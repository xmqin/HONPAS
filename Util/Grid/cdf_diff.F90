! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------
program cdf_diff
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

use m_grid

type(grid_t)   :: gp_total, gp_molecule, gp_slab
type(grid_t)   :: gp

integer :: n(3)

call get_cdf_grid("TotalSystem.Rho.grid.nc",gp_total)
call get_cdf_grid("SingleMolecule.Rho.grid.nc",gp_molecule)
call get_cdf_grid("Slab.Rho.grid.nc",gp_slab)

n(:) = gp_total%n(:)
if (any(gp_molecule%n /= n)) STOP "molecule"
if (any(gp_slab%n /= n)) STOP "slab"

allocate(gp%grid(n(1),n(2),n(3)))

gp%n = n
gp%cell = gp_total%cell
gp%grid = gp_total%grid - gp_slab%grid - gp_molecule%grid

print *, gp%cell
print *, gp%n
print *, gp%grid(1,1,1)

call put_cdf_grid(gp,"Diff.Rho.grid.nc")

#endif /* CDF */
end program cdf_diff

