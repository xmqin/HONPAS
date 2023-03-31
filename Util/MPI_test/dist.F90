
      program disttest

        ! Redistribution of orbital data.
        ! Two different distributions: 
        !   dist1: block-cyclic (as in Siesta)
        !   dist2: one block per processor, with fat last block (as in PEXSI)

      use mpi
      use class_Dist
      use m_matrix, only: matrix
      use m_redist, only: redistribute

      implicit none

      integer nprocs, i, j, ierr

      integer :: group_world, group1, group2
      integer :: nprocs1, nprocs2, bs, pbs, norbs, io, ncomms
      integer :: ig, ibeg, iend

      real    :: x


      integer ::  myrank1, myrank2, myid
      logical :: proc_in_set2, proc_in_set1
      integer, allocatable :: ranks(:)

      type(Dist) :: dist1, dist2

      type(matrix) :: m1, m2
         
!--------------------------------
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      if (myid  == 0) then
         print *, "Using ", nprocs, " procs in World."
         write(*,fmt="(a)",advance="no") "Enter number of procs in first set: "
         read *, nprocs1
         write(*,fmt="(a)",advance="no") "Enter number of procs in 2nd set: "
         read *, nprocs2
         write(*,fmt="(a)",advance="no") "Enter number of orbitals: "
         read *, norbs
         write(*,fmt="(a)",advance="no") "Enter blocksize for default dist: "
         read *, bs
      endif
      call MPI_Bcast(nprocs1,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(nprocs2,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(norbs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     
      call MPI_Bcast(bs,1,mpi_integer,0,MPI_COMM_WORLD,ierr)     

      call MPI_COMM_GROUP(MPI_COMM_WORLD, group_world, ierr) 

      ! New group, just for cleanliness
      ! Count from the end in this case
      allocate(ranks(nprocs1))
      do i = 1, nprocs1
         ranks(i) = nprocs - i
      end do
      call MPI_Group_incl(group_world, nprocs1, ranks, group1, ierr)
      call newDistribution(bs,group1,dist1,TYPE_BLOCK_CYCLIC,"bc dist")
      deallocate(ranks)

      ! New group, just for cleanliness
      allocate(ranks(nprocs2))
      do i = 1, nprocs2
         ranks(i) = i-1
      end do
      call MPI_Group_incl(group_world, nprocs2, ranks, group2, ierr)
      pbs = norbs/nprocs2
      call newDistribution(pbs,group2,dist2,TYPE_PEXSI,"px dist")
      deallocate(ranks)

      print *, "Group 1, Group 2: ", group(dist1), group(dist2)

      call mpi_group_rank(group(dist1),myrank1,ierr)
      call mpi_group_rank(group(dist2),myrank2,ierr)
      proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
      proc_in_set2 = (myrank2 /= MPI_UNDEFINED)

      ! Create source matrix
      if (proc_in_set1) then
         m1%norbs = norbs
         m1%no_l = num_local_elements(dist1,norbs,myrank1)
         allocate(m1%numcols(m1%no_l))
         do io = 1, m1%no_l
            call random_number(x)
            m1%numcols(io) = x*(0.4*norbs) + 1
         enddo
         m1%nnzl = sum(m1%numcols(1:m1%no_l))
         allocate(m1%cols(m1%nnzl))
         do j = 1, m1%nnzl
            call random_number(x)
            m1%cols(j) = x*norbs + 1
         enddo
      endif

      call redistribute(norbs,m1,dist1,m2,dist2,mpi_comm_world)


      if (proc_in_set1) then
         ibeg = 1
         do io = 1, m1%no_l
            iend = ibeg + m1%numcols(io) - 1
            ig = index_local_to_global(dist1,io,myrank1)
            print "(a,i4,a,2i5,2x,i5,2x,10i3)", "Src: ", myrank1, &
                 " il, ig, ncols, cols: ", io, ig, &
                 m1%numcols(io), m1%cols(ibeg:iend)
            ibeg = iend + 1
         enddo
      endif

      call MPI_Barrier(mpi_comm_world,ierr)

      if (proc_in_set2) then
         ibeg = 1
         do io = 1, m2%no_l
            iend = ibeg + m2%numcols(io) - 1
            ig = index_local_to_global(dist2,io,myrank2)
            print "(a,i4,a,2i5,2x,i5,2x,10i3)", "Dst: ", myrank2, &
              " il, ig, ncols, cols: ", io, ig, &
              m2%numcols(io), m2%cols(ibeg:iend)
            ibeg = iend + 1
         enddo
      endif

      call MPI_FINALIZE(ierr)

    end program disttest




