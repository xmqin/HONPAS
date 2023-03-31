program tst_ncdf

  use dictionary
  use netcdf_ncdf

  use tst_ncdf_utils

  implicit none

  character(len=30) :: fname

  integer :: Node, Nodes
  type(hNCDF) :: ncdf
  type(dictionary_t) :: dic, lv, atts

  call tst_mpi_init(Node,Nodes)

  fname = 'NCDF_CRT_3.nc'
  ! Create the netcdf-file
  call ncdf_create(ncdf,fname,&
       mode=NF90_64BIT_OFFSET, overwrite=.true.)

  ! Lets try and create a netcdf file, purely from a 
  ! dictionary.

  ! Lets first create the top level
  dic = ('DIMno'.kv. 10) // ('DIMunlim'.kv.NF90_UNLIMITED) // &
       ('DIMna'.kv. 100) // ('DIMxyz'.kv.3)
  ! Create the variable xa(xyz,na,unlim) (NF90_DOUBLE), atts
  atts = ('info'.kv.'Coordinates')//('unit'.kv.'Ang')
  lv = ('atts'.kvp.atts) // ('type'.kv.NF90_DOUBLE)
  lv = lv // ('dims'.kv.'xyz,na,unlim')
  dic = dic // ('VARv'.kvp.lv)
  call print(dic)

  ! Create another variable
  call nullify(lv)   ! we will retrieve them later on..
  call nullify(atts)
  ! Note that the atts dictionary now "floats" in dic
  atts = ('info'.kv.'Hello')//('unit'.kv.'ASNTHeosN')
  lv = ('atts'.kvp.atts) // ('type'.kv.NF90_DOUBLE_COMPLEX)
  lv = lv // ('dims'.kv.'xyz,na')
  dic = dic // ('VARcheck'.kvp.lv)
  call print(dic)

  call ncdf_crt(ncdf,dic)

  call ncdf_close(ncdf)

  call check_nc(fname)

#ifdef NCDF_PARALLEL
  call MPI_Finalize(Nodes)
#endif
  
end program tst_ncdf

