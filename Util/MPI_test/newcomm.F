
      program main

      use mpi_siesta

      implicit none

      integer, parameter :: dp = selected_real_kind(14,100)

      real(dp)  mypi, pi, pisub
      integer n, myid, numprocs, i, rc, ntotal, ierr

      integer :: comm_sub, group_world, group_sub, myid_sub
      logical :: worker


      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (myid  == 0) then
         print *, "Using ", numprocs, " procs in World."
      endif

      worker = (myid <= 1)
!
      call MPI_COMM_GROUP(MPI_COMM_WORLD, group_world, ierr) 
      call MPI_Group_incl(group_world, 2, (/0,1/), group_sub, ierr)
      call MPI_Comm_create(MPI_COMM_WORLD, group_sub, comm_sub, ierr)
!
      MPI_COMM_WORLD = comm_sub

      if (worker) then
         call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid_sub, ierr )
         mypi = (1+myid_sub) * 2.0
         call MPI_REDUCE(mypi,pi,1,mpi_double_precision,MPI_SUM,0,           
     &     MPI_COMM_WORLD,ierr)
         if (myid_sub == 0) then
            print *, "Using ", numprocs, " procs. pi: ", pi
         endif
      endif

      call sub(pi)

 30   call MPI_FINALIZE(rc)

      CONTAINS

      subroutine sub(x)
      real(dp), intent(in) :: x

      integer :: myid
      real(dp) :: y, ytot

      y = x
      call MPI_COMM_RANK( true_MPI_COMM_WORLD, myid, ierr )
      call MPI_Bcast(y,1,mpi_double_precision,0,
     $               true_MPI_COMM_WORLD,ierr)     
      call MPI_REDUCE(y,ytot,1,mpi_double_precision,MPI_SUM,0,           
     &     true_MPI_COMM_WORLD,ierr)
      if (myid == 0) then
         print *, "Total value using all procs: ", ytot
      endif
      
      end subroutine sub

      end




