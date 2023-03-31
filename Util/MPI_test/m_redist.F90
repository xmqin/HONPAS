module m_redist
public :: redistribute
CONTAINS
  subroutine redistribute(norbs,m1,dist1,m2,dist2,mpi_comm)

      use mpi
      use class_Dist
      use m_matrix, only: matrix
      use m_transfers, only: comm, do_transfers

      implicit none

    integer, intent(in)       :: norbs
    type(matrix) :: m1
    type(matrix) :: m2
    type(Dist) :: dist1, dist2
    integer, intent(in)   :: mpi_comm

    type(comm), dimension(:), allocatable, target :: comms
    type(comm), dimension(:), allocatable, target :: commsnnz
    type(comm), pointer :: c, cnnz

    integer ::  myrank1, myrank2, myid
    logical ::  proc_in_set1, proc_in_set2
    integer ::  ierr

    integer ::  i, io, g1, g2

      g1 = group(dist1)
      g2 = group(dist2)
      call mpi_group_rank(g1,myrank1,ierr)
      call mpi_group_rank(g2,myrank2,ierr)
      proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
      proc_in_set2 = (myrank2 /= MPI_UNDEFINED)

      call mpi_comm_rank(mpi_comm,myid,ierr)
      
      ! Figure out the communication needs
      call analyze_comms()

      ! In preparation for the transfer, we allocate
      ! storage for the second group of processors
      ! Note that m2%numcols (and, in general, any of the 2nd set 
      ! of arrays), will not be allocated by those processors
      ! not in the second set.

      if (proc_in_set2) then
         m2%no_l = num_local_elements(dist2,norbs,myrank2)
         allocate(m2%numcols(m2%no_l))
      endif

      call do_transfers(comms,m1%numcols,m2%numcols, &
                        g1,g2,mpi_comm)

      ! Now we can figure out how many non-zeros there are
      if (proc_in_set2) then
         m2%nnzl = sum(m2%numcols(1:m2%no_l))
         allocate(m2%cols(m2%nnzl))
         allocate(m2%vals(m2%nnzl))
      endif

      ! Generate a new comms-structure with new start/count indexes

      allocate(commsnnz(size(comms)))
      do i = 1, size(comms)
         c => comms(i)
         cnnz => commsnnz(i)

         cnnz%src = c%src
         cnnz%dst = c%dst
         if (myrank1 == c%src) then
            ! Starting position at source: previous cols plus 1
            cnnz%i1 = sum(m1%numcols(1:(c%i1-1))) + 1
            ! Number of items transmitted: total number of cols
            cnnz%nitems = sum(m1%numcols(c%i1 : c%i1 + c%nitems -1))
         endif
         if (myrank2 == c%dst) then
            ! Starting position at destination: previous cols plus 1
            cnnz%i2 = sum(m2%numcols(1 : (c%i2-1))) + 1
            ! Number of items transmitted: total number of cols
            cnnz%nitems = sum(m2%numcols(c%i2 : c%i2 + c%nitems -1))
         endif
      end do

      ! Transfer the cols arrays
      call do_transfers(commsnnz,m1%cols,m2%cols, &
                        g1, g2, mpi_comm)

      ! Transfer the values arrays
!      call do_transfers(commsnnz,m1%vals,m2%vals, &
!                        g1,g2,mpi_comm)

      CONTAINS

!-----------------------------------------------------
   subroutine analyze_comms()

      integer, allocatable, dimension(:) :: p1, p2, isrc, idst
      integer :: ncomms

      ! Find the communication needs for each orbital
      ! This information is replicated in every processor
      ! (Note that the indexing functions are able to find
      !  out the information for any processor. For the
      ! block-cyclic and "pexsi" distributions, this is quite
      ! easy. For others, the underlying indexing arrays might
      ! be large...)

      ! It might not be necessary to have this in memory. It 
      ! can be done on the fly
      allocate(p1(norbs),p2(norbs),isrc(norbs),idst(norbs))

      if (myid == 0) then
         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
      endif
      do io = 1, norbs
         p1(io) = node_handling_element(dist1,io)
         p2(io) = node_handling_element(dist2,io)
         isrc(io) = index_global_to_local(dist1,io,p1(io))
         idst(io) = index_global_to_local(dist2,io,p2(io))
         if (myid == 0) then
            if ((norbs < 1000) .or. (mod(io,12) == 0)) then
               write(6,"(5i10)") io, p1(io), isrc(io), p2(io), idst(io)
            endif
         endif
      enddo

      ! Aggregate communications
      ! First pass: find out how many there are, on the basis
      ! of groups of orbitals that share the same source and
      ! destination. Due to the form of the distributions, the
      ! local indexes are also correlative in that case, so we
      ! only need to check for p1 and p2. (Check whether this
      ! applies to every possible distribution...)

      ncomms = 1
      do io = 2, norbs
         if ((p1(io) /= p1(io-1)) .or. (p2(io) /= p2(io-1))) then
            ncomms = ncomms + 1
         else
            !
         endif
      enddo

      allocate(comms(ncomms))

      ! Second pass: Fill in the data structures
      ncomms = 1
      c => comms(ncomms)
      io = 1
      c%src = p1(io)
      c%dst = p2(io)
      c%i1  = isrc(io)
      c%i2  = idst(io)
      c%nitems = 1
      do io = 2, norbs
         if ((p1(io) /= p1(io-1)) .or. (p2(io) /= p2(io-1))) then
            ! end of group -- new communication
            ncomms = ncomms + 1
            c => comms(ncomms)
            c%src = p1(io)
            c%dst = p2(io)
            c%i1  = isrc(io)
            c%i2  = idst(io)
            c%nitems = 1
         else
            ! we stay in the same communication
            c%nitems = c%nitems + 1
         endif
      enddo

      if (myid == 0) then
         do i = 1, ncomms
            c => comms(i)
            print "(a,i5,a,2i5,2i7,i5)", &
                 "comm: ", i, " src, dst, i1, i2, n:", &
                 c%src, c%dst, c%i1, c%i2, c%nitems
         enddo
      endif

      deallocate(p1,p2,isrc,idst)

    end subroutine analyze_comms

  end subroutine redistribute

end module m_redist
