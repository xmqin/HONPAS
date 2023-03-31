program tst_ncdf_4
#ifdef NCDF_PARALLEL
  use mpi
#endif
  use dictionary
  use netcdf_ncdf

  use tst_ncdf_utils

  implicit none

  character(len=30) :: fname

  type(hNCDF) :: ncdf, grp1, grp2
  integer :: Node, Nodes, i, comp_lvl
  type(dictionary_t) :: dic
  logical :: assert

  call tst_mpi_init(Node,Nodes)

  ! Create the netcdf-file
  if ( Nodes > 1 ) then
     fname = 'NCDF_par_4.nc'
     call ncdf_create(ncdf,fname, &
          mode=ior(NF90_NETCDF4, NF90_MPIIO), overwrite=.true., &
          parallel=.true., &
          comm=MPI_Comm_World)
     ! parallel writes are not allowed with compression
     ! Offset positions are not well defined.
     comp_lvl = 0
  else
     fname = 'NCDF_seq_4.nc'
     call ncdf_create(ncdf,fname, &
          mode=NF90_NETCDF4, overwrite=.true.)
     comp_lvl = 3
  end if

  call ncdf_def_dim(ncdf,'x',1)
  call ncdf_def_dim(ncdf,'y',NF90_UNLIMITED)
  call ncdf_def_dim(ncdf,'z',2)

  dic = ('unit'.kv.'m')//('Name'.kv.'What is this')
  call ncdf_def_var(ncdf,'v',NF90_DOUBLE,(/'x','y'/), &
       atts=dic,compress_lvl=comp_lvl)
  dic = dic // ('Name'.kv.'Height') // ('unit'.kv. 'm')
  call ncdf_def_var(ncdf,'h',NF90_DOUBLE,(/'z','y'/), &
       atts=dic,compress_lvl=comp_lvl)
  call delete(dic)

  call ncdf_default(ncdf,access=NF90_COLLECTIVE)

  call ncdf_def_grp(ncdf,'info',grp1)
  call ncdf_def_dim(grp1,'N',Nodes*2)
  call ncdf_def_dim(grp1,'Np1',Nodes*2+1)
  call ncdf_def_var(grp1,'Nodes',NF90_INT,(/'N'/))
  call ncdf_def_var(grp1,'Inde',NF90_INT,(/'Np1'/))
  call ncdf_default(grp1,access=NF90_COLLECTIVE)
  call ncdf_par_access(grp1,name='Inde',access=NF90_INDEPENDENT)
  call ncdf_def_grp(grp1,'scndlevel',grp2)
  call ncdf_def_dim(grp2,'N',Nodes*2)
  call ncdf_def_var(grp2,'Nodes',NF90_INT,(/'N'/))
  call ncdf_default(grp2,access=NF90_COLLECTIVE)

  ! print out the leveled netcdf
  call ncdf_print(ncdf)
  call ncdf_print(grp1)
  call ncdf_print(grp2)

  do i = 1 , Nodes * 2
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(grp1,'Nodes',i,start=(/i/))
        call ncdf_put_var(grp2,'Nodes',i,start=(/i/))
     end if
  end do

  do i = 1 , Nodes * 2 + 1
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(grp1,'Inde',i,start=(/i/))
     end if
  end do

  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'v',real(i,8),start=(/1,i/))
     end if
  end do
  do i = 1 , 10
     if ( mod(i,Nodes) == Node ) then
        call ncdf_put_var(ncdf,'h',(/real(i,8),real(i*2,8)/),start=(/1,i/))
     end if
  end do

  ! We assert that the dimensions exist
  dic = ('x'.kv.1)//('y'.kv.10)//('z'.kv.2)
  call ncdf_assert(ncdf,assert,dims=dic)
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion...',Node
  end if

  ! Groups inherit "unknown" dimensions from their
  ! parent. Hence we just check them also
  dic = dic // ('N'.kv.Nodes * 2)
  call ncdf_assert(grp1,assert,dims=dic)
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED, grp1',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion..., grp1',Node
  end if
  call ncdf_assert(grp2,assert,dims=dic)
  call delete(dic)
  if ( .not. assert ) then
     write(*,*) 'ASSERTION NOT FULFILLED, grp2',Node
     stop 9
  else
     write(*,*) 'Fulfilled assertion..., grp2',Node
  end if
  
  call ncdf_close(ncdf)

  call check_nc(fname)

#ifdef NCDF_PARALLEL
  call MPI_Finalize(Nodes)
#endif

end program tst_ncdf_4
