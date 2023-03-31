subroutine memory_all(str,comm)
  use m_rusage, only :  rss_max
#ifdef MPI  
  use mpi_siesta
#endif
  character(len=*), intent(in) :: str
  integer, intent(in)          :: comm

  integer :: mpierror
  real    :: mem(2), tmem(2)
  integer :: nprocs, myrank

#ifdef MPI
  call MPI_Comm_Size( Comm, nprocs, MPIerror )
  call MPI_Comm_Rank( Comm, myrank, MPIerror )
#else
  nprocs = 1
  myrank = 0
#endif

  mem(2) = rss_max()
  mem(1) = - mem(2)

#ifdef MPI
  call MPI_Reduce(mem,tmem,2,MPI_Real,MPI_max,0,comm,MPIerror)
  mem = tmem
#endif

  if (myrank == 0) then
     write(6,"(a,2f12.2)") " &m -- Peak memory (Mb) " // trim(str) // " (max,min): ",  mem(2), -mem(1)
  endif

end subroutine memory_all
