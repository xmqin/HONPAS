module tst_ncdf_utils

  use iso_c_binding
#ifdef NCDF_PARALLEL
  use mpi
#endif
  use variable
  use dictionary

  implicit none

  integer, parameter :: sp = selected_real_kind(p=6)
  integer, parameter :: dp = selected_real_kind(p=14)
#ifndef NCDF_PARALLEL
  integer, parameter :: MPI_Comm_World = -1
#endif

  integer, private, save :: Node, Nodes

  interface c_interface
     function nf90_set_log_level(level) bind (C, name = "nc_set_log_level")
       use iso_c_binding
       implicit none
       integer(c_int) :: nf90_set_log_level
       integer(c_int), intent (in) :: level
     end function nf90_set_log_level
  end interface c_interface

contains

  subroutine tst_mpi_init(oNode,oNodes)
#ifdef NCDF_PARALLEL
    use mpi
#endif
    use netcdf_ncdf

    integer, intent(out) :: oNode, oNodes
    
#ifdef NCDF_PARALLEL
    integer :: MPIerror
    
    call MPI_Init(MPIerror)
    
    call MPI_Comm_rank(MPI_COMM_WORLD,Node,MPIerror)
    call MPI_Comm_size(MPI_COMM_WORLD,Nodes,MPIerror)
#else
    Node = 0
    Nodes = 1
#endif

    oNode = Node
    oNodes = Nodes

    call ncdf_IOnode(Node == 0)
    
  end subroutine tst_mpi_init
  
  subroutine msg(c)
    character(len=*), intent(in) :: c
    if ( Node == 0 ) then
       write(*,*) c
    end if
  end subroutine msg
  
  subroutine goto_dir(a)
    character(len=*), intent(in) :: a
    if ( Node == 0 ) then
       call system('mkdir -p '//a)
    end if
    call chdir(a) 
  end subroutine goto_dir
  
  subroutine check_nc(file)
#ifdef NCDF_PARALLEL
    use mpi
#endif
    character(len=*), intent(in) :: file
    logical :: exist
#ifdef NCDF_PARALLEL
    integer :: MPIerror
#endif
    inquire(file=file,exist=exist)
#ifdef NCDF_PARALLEL
    call MPI_barrier(MPI_COMM_world,MPIerror)
#endif
    if ( exist ) then
       if ( Node == 0 ) then
          call system('ncdump -v v '//file)
          call system('ncdump -k '//file)
       end if
    else if ( Node == 0 ) then
       write(*,'(a)')'File: '//trim(file)//' does not exist!'
       stop 9
    end if
  end subroutine check_nc

end module tst_ncdf_utils
