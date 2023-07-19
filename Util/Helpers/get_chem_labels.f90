program get_chem_labels
!
! Prints to standard output a list of chemical species labels in an fdf file.
! Alberto Garcia, Sep 2009
!
  use f2kcli     
  use chemical, only: read_chemical_types, number_of_species, species_label
  use fdf, only : fdf_init

  integer            :: i, n
  logical            :: silent

  integer             :: narg
  character(len=256)  :: filein

  narg = command_argument_count()
  if (narg == 1) then
     call get_command_argument(1,filein)
  else
     stop "Usage: get_chem_labels FDF_FILE"
  endif

  call fdf_init(filein,".tmp_fdflog")

  call read_chemical_types(silent=.true.)
  n = number_of_species()
  do i = 1, n
     write(*,"(a)", advance="no") trim(species_label(i)) // " "
  enddo
  write(*,*)

end program get_chem_labels
