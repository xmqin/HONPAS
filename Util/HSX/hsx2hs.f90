program hsx2hs

! Generates an HS file from an HSX file (more compact)

use hsx_m, only: read_hsx_file, hsx_t, write_hs_file

type(hsx_t)  :: h

call read_hsx_file(h,"HSX") 

call write_hs_file(h, "HS")

end program hsx2hs
