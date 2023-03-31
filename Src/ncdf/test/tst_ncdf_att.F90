program tst_ncdf
  use netcdf_ncdf

  use tst_ncdf_utils

  implicit none

  character(len=30) :: fname

  type(hNCDF) :: ncdf
  integer :: Node, Nodes, i
  type(dictionary_t) :: dic, ld
  character(len=255) :: a1, a2
  integer, pointer :: i0, i1(:)
  real(dp), pointer :: d1(:)
  type(variable_t) :: v

  call tst_mpi_init(Node,Nodes)

  ! Create the netcdf-file
  if ( Nodes > 1 ) then
     fname = 'NCDF_par_3.nc'
     call ncdf_create(ncdf,'NCDF_par_3.nc',&
          mode=NF90_MPIIO, &
          overwrite=.true., comm=MPI_Comm_World)
  else
     fname = 'NCDF_seq_3.nc'
     call ncdf_create(ncdf,'NCDF_seq_3.nc',&
          mode=NF90_64BIT_OFFSET, overwrite=.true.)
  end if


  ! Define dimensions
  call ncdf_def_dim(ncdf,'x',1)
  call ncdf_def_dim(ncdf,'y',NF90_UNLIMITED)
  
  dic = ('unit'.kv.'m')//('Name'.kv.'What is this')
  call ncdf_def_var(ncdf,'v',NF90_DOUBLE,(/'x','y'/), atts=dic)
  call delete(dic) ! ensure deletion


  dic = ('g1'.kv.1)//('g2'.kv.(/1,2/))//('g3'.kv.(/1._dp,2._dp/))
  dic = dic // ('g4'.kv.'hello')
  call ncdf_put_gatt(ncdf,atts=dic)
  call delete(dic)

  ! We must only redefine access patterns after
  ! enddef...
  call ncdf_enddef(ncdf)

  call ncdf_par_access(ncdf,access=NF90_COLLECTIVE)

  call ncdf_print(ncdf)
  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'v',real(i,8),start=(/1,i/))
     end if
  end do

  ! redefining a NetCDF file after already ending the definition step is ONLY 
  ! allowed in NetCDF 3 formats...
  call ncdf_def_dim(ncdf,'z',2)
  dic = ('unit'.kv.'m')//('Name'.kv.'Height')
  dic = dic//('ATT_DELETE'.kv.1)
  call ncdf_def_var(ncdf,'h',NF90_DOUBLE,(/'z','y'/), atts=dic)
  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'h',(/real(i,8),real(i*2,8)/),start=(/1,i/))
     end if
  end do

  call ncdf_close(ncdf)
  call check_nc(fname)

  if ( Nodes > 1 ) then
     call ncdf_open(ncdf,'NCDF_par_3.nc')
  else
     call ncdf_open(ncdf,'NCDF_seq_3.nc')
  end if

  call ncdf_get_gatt(ncdf,atts=dic)
  
  ! Loop attributes
  ld = .first. dic
  do 
     if ( .empty. ld ) exit
     
     ! check contents
     a1 = .key. ld
     v  = .val. ld
     select case ( a1 )
     case ( 'g1' )
        call associate(i0,v)
        if ( i0 /= 1 ) then
           print *,'ERROR on g1',i0
           stop 9
        end if
     case ( 'g2' )
        call associate(i1,v)
        if ( any(i1 /= (/1,2/)) ) then
           print *,'ERROR on g2'
           stop 9
        end if
     case ( 'g3' )
        call associate(d1,v)
        if ( maxval(abs(d1 - (/1._dp,2._dp/))) > 1.e-10_dp ) then
           print *,'ERROR on g3'
           stop 9
        end if
     case ( 'g4' )
        call assign(a2,v)
        if ( a2 /= 'hello' ) then
           print *,'ERROR on g4'
           stop 9
        end if
     end select
     ld = .next. ld
  end do
  
  call ncdf_close(ncdf)

  call check_nc(fname)

#ifdef NCDF_PARALLEL
  call MPI_Finalize(Nodes)
#endif
  
end program tst_ncdf

