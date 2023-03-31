program tst_ncdf_3
#ifdef NCDF_PARALLEL
  use mpi
#endif
  use dictionary
  use netcdf_ncdf

  use tst_ncdf_utils

  implicit none

  character(len=30) :: fname

  type(hNCDF) :: ncdf
  integer :: Node, Nodes, i
  character(len=1) :: ci
  type(dictionary_t) :: dic
  logical :: assert
  integer :: ilist1(1,10), ilist2(2,10)

  call tst_mpi_init(Node,Nodes)

  ! Create the netcdf-file
  if ( Nodes > 1 ) then
     fname = 'NCDF_par_3.nc'
     call ncdf_create(ncdf,fname, &
          mode=NF90_MPIIO, &
          overwrite=.true., comm=MPI_Comm_World)
  else
     fname = 'NCDF_seq_3.nc'
     call ncdf_create(ncdf,fname,&
          mode=NF90_64BIT_OFFSET, overwrite=.true.)
  end if

  ilist1 = 0
  ilist2 = 0
  
  ! Define dimensions
  call ncdf_def_dim(ncdf,'x',1)
  call ncdf_def_dim(ncdf,'y',NF90_UNLIMITED)
  
  dic = ('unit'.kv.'m')//('Name'.kv.'What is this')
  call ncdf_def_var(ncdf,'v',NF90_INT,(/'x','y'/), atts=dic)
  call delete(dic) ! ensure deletion

  call ncdf_print(ncdf)
  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'v',i,start=(/1,i/))
     end if
  end do

  ! redefining a NetCDF file after already ending the definition step is ONLY 
  ! allowed in NetCDF 3 formats...
  call ncdf_def_dim(ncdf,'z',2)
  dic = ('unit'.kv.'m')//('Name'.kv.'Height')
  dic = dic//('ATT_DELETE'.kv.1)
  call ncdf_def_var(ncdf,'h',NF90_INT,(/'z','y'/), atts=dic)
  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'h',(/i,i*2/),start=(/1,i/))
     end if
  end do

  ! We assert that the dimensions exist
  dic = ('x'.kv.1)//('y'.kv.10)
  call ncdf_assert(ncdf,assert,dims=dic)
  call delete(dic)
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion...',Node
  end if

  call ncdf_get_var(ncdf,'v',ilist1)
  assert = .true.
  do i = 1 , 10
     if ( ilist1(1,i) /= i ) then
        assert = .false.
     end if
  end do
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion...',Node
  end if


  call ncdf_get_var(ncdf,'h',ilist2)
  assert = .true.
  do i = 1 , 10
     if ( any(ilist2(:,i) /= (/i,i*2/)) ) then
        assert = .false.
     end if
  end do
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion...',Node
  end if

  call ncdf_close(ncdf)

  call check_nc(fname)

#ifdef NCDF_PARALLEL
  call MPI_Finalize(Nodes)
#endif

end program tst_ncdf_3

