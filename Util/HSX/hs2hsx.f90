program hs2hsx

! Generates an HSX file  (more compact) from an HS file

use hsx_m, only: write_hsx_file, hsx_t, read_hs_file

type(hsx_t)  :: h

call read_hs_file(h,"HS") 

call write_hsx_file(h,"HSX_OUT")

end program hs2hsx

