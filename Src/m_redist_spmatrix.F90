
! --- Tangled code
module m_redist_spmatrix
#ifdef SIESTA__PEXSI
 implicit none
 type, public :: comm_t
    integer :: src, dst, i1, i2, nitems
 end type comm_t
 
 integer, parameter, private :: dp = selected_real_kind(10,100)
 
 type, public :: dp_pointer
    ! boxed array pointer type
    real(dp), pointer :: data(:) => null()
 end type dp_pointer
 
 type, public ::  aux_matrix
    integer :: norbs = -1
    integer :: no_l  = -1
    integer :: nnzl  = -1
    integer, pointer :: numcols(:) => null()
    integer, pointer :: cols(:)    => null()
    ! array of 1D pointers
    type(dp_pointer), dimension(:), pointer :: vals(:) => null()
 end type aux_matrix
 public :: redistribute_spmatrix
CONTAINS
 subroutine redistribute_spmatrix(norbs,m1,dist1,m2,dist2,bridge_comm)
 
   use mpi
   use class_Distribution
   use alloc,       only: re_alloc, de_alloc
 
   implicit none
 
   integer, intent(in)       :: norbs   ! Overall number of rows
   type(aux_matrix) :: m1               ! Source matrix
   type(aux_matrix) :: m2               ! Destination matrix -- it is allocated
   type(Distribution) :: dist1, dist2           ! Distributions
   integer, intent(in)   :: bridge_comm    ! Umbrella Communicator
 
   type(comm_t), dimension(:), allocatable, target :: comms
   type(comm_t), dimension(:), allocatable, target :: commsnnz
   type(comm_t), pointer :: c, cnnz
 
   integer ::  myrank1, myrank2, myid, gg
   logical ::  proc_in_set1, proc_in_set2
   integer ::  ierr
 
   integer ::  i, io, g1, g2, j, nvals
   integer ::  comparison, n1, n2, c1, c2
   integer, parameter :: dp = selected_real_kind(10,100)
   real(dp), dimension(:), pointer  :: data1 => null(), data2 => null()
 
   integer, allocatable :: ranks1(:), ranks2(:)
 
   ! The communicators are a sanity check on the ranks
 
   call mpi_comm_rank(bridge_comm,myid,ierr)
   c1 = ref_comm(dist1)
   c2 = ref_comm(dist2)
   call MPI_Comm_Compare(c1,c2,comparison,ierr)
 
   select case (comparison)
   case (MPI_IDENT)
      ! Communicators are identical
   case (MPI_CONGRUENT)
      ! Communicators have the same group and rank order, but they differ in context
   case (MPI_SIMILAR)
      ! Rank order is different
      call die("Different rank order in communicators")
   case (MPI_UNEQUAL)
      ! Big mess
      call die("Incompatible distribution reference communicators")
   end select
 
   ! Now check congruence with the provided bridge_comm
   
   call MPI_Comm_Compare(c1,bridge_comm,comparison,ierr)
   select case (comparison)
   case (MPI_IDENT)
      ! Communicators are identical
   case (MPI_CONGRUENT)
      ! Communicators have the same group and rank order, but they differ in context
      ! We will use bridge_comm
   case (MPI_SIMILAR)
      ! Rank order is different
      call die("Different rank order in dist communicators and bridge comm")
   case (MPI_UNEQUAL)
      ! Big mess
      call die("Incompatible bridge and dist communicators")
   end select
 
   ! Now create groups g1 and g2.
   ! (DO NOT trust the internal handles)
   call MPI_Comm_Group(bridge_comm,gg,ierr)
   call get_ranks_in_ref_comm(dist1, ranks1)
   call get_ranks_in_ref_comm(dist2, ranks2)
   n1 = size(ranks1)
   n2 = size(ranks2)
   call MPI_Group_Incl(gg,n1,ranks1,g1,ierr)
   call MPI_Group_Incl(gg,n2,ranks2,g2,ierr)
 
   ! The rest is the same as before
 
   call mpi_group_rank(g1,myrank1,ierr)
   call mpi_group_rank(g2,myrank2,ierr)
 
   proc_in_set1 = (myrank1 /= MPI_UNDEFINED)
   proc_in_set2 = (myrank2 /= MPI_UNDEFINED)
 
   if (proc_in_set1 .or. proc_in_set2) then
     print "(a,3i6,2l2)", "world_rank, rank1, rank2, ing1?, ing2?", myid,  &
        myrank1, myrank2, proc_in_set1, proc_in_set2
   endif
 
   ! Figure out the communication needs
   call analyze_comms()
 
   ! In preparation for the transfer, we allocate
   ! storage for the second group of processors
   ! Note that m2%numcols (and, in general, any of the 2nd set 
   ! of arrays), will not be allocated by those processors
   ! not in the second set.
 
 
   if (proc_in_set2) then
      m2%norbs = norbs
      m2%no_l = num_local_elements(dist2,norbs,myrank2)
      call re_alloc(m2%numcols,1,m2%no_l,"m2%numcols","redistribute_spmatrix")
   endif
 
   if (myid == 0) print *, "About to transfer numcols..."
   call do_transfers_int(comms,m1%numcols,m2%numcols, &
        g1,g2,bridge_comm)
   if (myid == 0) print *, "Transferred numcols."
 
   ! We need to tell the processes in set 2 how many
   ! "vals" to expect.
   if (proc_in_set1) then
      if (associated(m1%vals)) then
         nvals = size(m1%vals)
      else
         nvals = 0
      endif
   endif
   ! Now do a broadcast within bridge_comm, using as root one
   ! process in the first set. Let's say the one with rank 0
   ! in g1, the first in the set, which will have rank=ranks1(1)
   ! in bridge_comm
   call MPI_Bcast(nvals,1,MPI_Integer,ranks1(1),bridge_comm,ierr)
 
   ! Now we can figure out how many non-zeros there are
   if (proc_in_set2) then
      m2%nnzl = sum(m2%numcols(1:m2%no_l))
      call re_alloc(m2%cols,1,m2%nnzl,"m2%cols","redistribute_spmatrix")
 
      if (nvals > 0) then
         allocate(m2%vals(nvals))
         do j=1,nvals
            call re_alloc(m2%vals(j)%data,1,m2%nnzl, &
                 "m2%vals(j)%data","redistribute_spmatrix")
         enddo
      endif
 
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
 
 !!$         do i = 1, size(comms)
 !!$            c => commsnnz(i)
 !!$            if (myrank1 == c%src) then
 !!$               print "(a,i5,a,2i5,2i7,i5)", &
 !!$                 "commnnz(src): ", i, " src, dst, i1, (), n:", &
 !!$                 c%src, c%dst, c%i1, -1, c%nitems
 !!$            endif
 !!$            if (myrank2 == c%dst) then
 !!$               print "(a,i5,a,2i5,2i7,i5)", &
 !!$                 "commnnz(dst): ", i, " src, dst, (), i2, n:", &
 !!$                 c%src, c%dst, -1, c%i2, c%nitems
 !!$            endif
 !!$         enddo
 
   if (myid == 0) print *, "About to transfer cols..."
   ! Transfer the cols arrays
   call do_transfers_int(commsnnz,m1%cols,m2%cols, &
        g1, g2, bridge_comm)
 
   if (myid == 0) print *, "About to transfer values..."
   ! Transfer the values arrays
   do j=1, nvals
      if (proc_in_set1) data1 => m1%vals(j)%data
      if (proc_in_set2) data2 => m2%vals(j)%data
      call do_transfers_dp(commsnnz,data1,data2, &
           g1,g2,bridge_comm)
   enddo
   nullify(data1,data2)
   if (myid == 0) print *, "Done transfers."
 
   deallocate(commsnnz)
   deallocate(comms)
   deallocate(ranks1, ranks2)
 
   call MPI_group_free(gg,ierr)
   call MPI_group_free(g1,ierr)
   call MPI_group_free(g2,ierr)
 
 CONTAINS
 
   
   !-----------------------------------------------------
      subroutine analyze_comms()
   
         integer, allocatable, dimension(:) :: p1, p2, isrc, idst
         integer :: ncomms
   
         ! To turn on debug printing, set this to .true.
         logical, save :: comms_not_printed = .false. 
   
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
   
   !      if (myid == 0) then
   !         write(6,"(5a10)") "Orb", "p1", "i1", "p2", "i2"
   !      endif
         do io = 1, norbs
            p1(io) = node_handling_element(dist1,io)
            p2(io) = node_handling_element(dist2,io)
            isrc(io) = index_global_to_local(dist1,io,p1(io))
            idst(io) = index_global_to_local(dist2,io,p2(io))
   !         if (myid == 0) then
   !            if ((norbs < 1000) .or. (mod(io,12) == 0)) then
   !               write(6,"(5i10)") io, p1(io), isrc(io), p2(io), idst(io)
   !            endif
   !        endif
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
   
         if (myid == 0 .and. comms_not_printed) then
            do i = 1, ncomms
               c => comms(i)
               write(6,"(a,i5,a,2i5,2i7,i5)") &
                    "comm: ", i, " src, dst, i1, i2, n:", &
                    c%src, c%dst, c%i1, c%i2, c%nitems
            enddo
            comms_not_printed = .false.
         endif
   
         deallocate(p1,p2,isrc,idst)
   
       end subroutine analyze_comms
 
 end subroutine redistribute_spmatrix
 !--------------------------------------------------
    subroutine do_transfers_int(comms,data1,data2,g1,g2,bridge_comm)
 
      use mpi
      type(comm_t), intent(in), target     :: comms(:)
      integer, dimension(:), pointer  :: data1
      integer, dimension(:), pointer  :: data2
      integer, intent(in)                :: g1
      integer, intent(in)                :: g2
      integer, intent(in)                :: bridge_comm
 
      integer                 :: basegroup, nsize1, nsize2, ierr
      integer, allocatable    :: comm_rank1(:), comm_rank2(:)
 
 
      integer :: ncomms
      integer :: i
      integer :: nrecvs_local, nsends_local
      integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)
      integer :: src_in_comm, dst_in_comm
      integer :: myrank1, myrank2, myrank
      type(comm_t), pointer :: c
 
 
       ! Find the rank correspondences, in case
       ! there is implicit renumbering at the time of group creation
 
       call  MPI_Comm_group( bridge_comm, basegroup, ierr )
       call  MPI_Comm_Rank( bridge_comm, myrank, ierr )
 
       call  MPI_Group_Size( g1, nsize1, ierr )
       call  MPI_Group_Size( g2, nsize2, ierr )
 
       allocate(comm_rank1(0:nsize1-1))
       call MPI_Group_translate_ranks( g1, nsize1, (/ (i,i=0,nsize1-1) /), &
                                       basegroup, comm_rank1, ierr )
 !      print "(i4,a,10i3)", myrank, ":Ranks of g1 in base group:", comm_rank1
 
       allocate(comm_rank2(0:nsize2-1))
       call MPI_Group_translate_ranks( g2, nsize2, (/ (i,i=0,nsize2-1) /), &
                                       basegroup, comm_rank2, ierr )
 !      print "(i4,a,10i3)", myrank,":Ranks of g2 in base group:", comm_rank2
 
       call mpi_group_rank(g1,myrank1,ierr)
       call mpi_group_rank(g2,myrank2,ierr)
       
 !      print "(i4,a,2i6)", myrank,": Ranks in g1 and g2: ", myrank1, myrank2
 !      print "(i4,a,2i3)", myrank,": g1 and g2: ", g1, g2
 
 
       ! Do the actual transfers. 
       ! This version with non-blocking communications
 
      ncomms = size(comms)
 
       ! Some bookkeeping for the requests
       nrecvs_local = 0
       nsends_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank2 == c%dst) then
             nrecvs_local = nrecvs_local + 1
          endif
          if (myrank1 == c%src) then
             nsends_local = nsends_local + 1
          endif
       enddo
       allocate(local_reqR(nrecvs_local))
       allocate(local_reqS(nsends_local))
       allocate(statuses(mpi_status_size,nrecvs_local))
 
       ! First, post the receives
       nrecvs_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank2 == c%dst) then
             nrecvs_local = nrecvs_local + 1
             src_in_comm = comm_rank1(c%src)
             call MPI_irecv(data2(c%i2),c%nitems,MPI_integer,src_in_comm, &
                            i,bridge_comm,local_reqR(nrecvs_local),ierr)
          endif
       enddo
 
       ! Post the sends
       nsends_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank1 == c%src) then
             nsends_local = nsends_local + 1
             dst_in_comm = comm_rank2(c%dst)
             call MPI_isend(data1(c%i1),c%nitems,MPI_integer,dst_in_comm, &
                         i,bridge_comm,local_reqS(nsends_local),ierr)
          endif
       enddo
 
       ! A former loop of waits can be substituted by a "waitall",
       ! with every processor keeping track of the actual number of 
       ! requests in which it is involved.
 
       ! Should we wait also on the sends?
 
       call MPI_waitall(nrecvs_local, local_reqR, statuses, ierr)
 
 
       ! This barrier is needed, I think
       call MPI_Barrier(bridge_comm,ierr)
 
       deallocate(local_reqR, local_reqS, statuses)
 
     end subroutine do_transfers_int
 
 !--------------------------------------------------
    subroutine do_transfers_dp(comms,data1,data2,g1,g2,bridge_comm)
 
      use mpi
      integer, parameter :: dp = selected_real_kind(10,100)
 
      type(comm_t), intent(in), target     :: comms(:)
      real(dp), dimension(:), pointer :: data1
      real(dp), dimension(:), pointer :: data2
      integer, intent(in)                :: g1
      integer, intent(in)                :: g2
      integer, intent(in)                :: bridge_comm
 
      integer                 :: basegroup, nsize1, nsize2, ierr
      integer, allocatable    :: comm_rank1(:), comm_rank2(:)
 
 
      integer :: ncomms
      integer :: i
      integer :: nrecvs_local, nsends_local
      integer, allocatable :: statuses(:,:), local_reqR(:), local_reqS(:)
      integer :: src_in_comm, dst_in_comm
      integer :: myrank1, myrank2, myid
      type(comm_t), pointer :: c
 
      call  MPI_Comm_Rank( bridge_comm, myid, ierr )
 !     print *, "Entering transfer_dp"
 !     print *, "rank, Associated data1: ", myid, associated(data1)
 !     print *, "rank, Associated data2: ", myid, associated(data2)
 
       ! Find the rank correspondences, in case
       ! there is implicit renumbering at the time of group creation
 
       call  MPI_Comm_group( bridge_comm, basegroup, ierr )
       call  MPI_Group_Size( g1, nsize1, ierr )
       call  MPI_Group_Size( g2, nsize2, ierr )
       allocate(comm_rank1(0:nsize1-1))
       call MPI_Group_translate_ranks( g1, nsize1, (/ (i,i=0,nsize1-1) /), &
                                       basegroup, comm_rank1, ierr )
 !      print "(a,10i3)", "Ranks of g1 in base group:", comm_rank1
       allocate(comm_rank2(0:nsize2-1))
       call MPI_Group_translate_ranks( g2, nsize2, (/ (i,i=0,nsize2-1) /), &
                                       basegroup, comm_rank2, ierr )
 !      print "(a,10i3)", "Ranks of g2 in base group:", comm_rank2
 
       call mpi_group_rank(g1,myrank1,ierr)
       call mpi_group_rank(g2,myrank2,ierr)
 
       ! Do the actual transfers. 
       ! This version with non-blocking communications
 
      ncomms = size(comms)
 
       ! Some bookkeeping for the requests
       nrecvs_local = 0
       nsends_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank2 == c%dst) then
             nrecvs_local = nrecvs_local + 1
          endif
          if (myrank1 == c%src) then
             nsends_local = nsends_local + 1
          endif
       enddo
       allocate(local_reqR(nrecvs_local))
       allocate(local_reqS(nsends_local))
       allocate(statuses(mpi_status_size,nrecvs_local))
 
       ! First, post the receives
       nrecvs_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank2 == c%dst) then
             nrecvs_local = nrecvs_local + 1
             src_in_comm = comm_rank1(c%src)
             call MPI_irecv(data2(c%i2),c%nitems,MPI_Double_Precision,src_in_comm, &
                            i,bridge_comm,local_reqR(nrecvs_local),ierr)
          endif
       enddo
 
       ! Post the sends
       nsends_local = 0
       do i=1,ncomms
          c => comms(i)
          if (myrank1 == c%src) then
             nsends_local = nsends_local + 1
             dst_in_comm = comm_rank2(c%dst)
             call MPI_isend(data1(c%i1),c%nitems,MPI_Double_Precision,dst_in_comm, &
                         i,bridge_comm,local_reqS(nsends_local),ierr)
          endif
       enddo
 
       ! A former loop of waits can be substituted by a "waitall",
       ! with every processor keeping track of the actual number of 
       ! requests in which it is involved.
 
       ! Should we wait also on the sends?
 
       call MPI_waitall(nrecvs_local, local_reqR, statuses, ierr)
 
 
       ! This barrier is needed, I think
       call MPI_Barrier(bridge_comm,ierr)
 
       deallocate(local_reqR, local_reqS, statuses)
 
     end subroutine do_transfers_dp
#endif
end module m_redist_spmatrix
! --- End of tangled code
